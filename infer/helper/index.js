

var gammaCof = [
  76.18009172947146,
  -86.50532032941677,
  24.01409824083091,
  -1.231739572450155,
  0.1208650973866179e-2,
  -0.5395239384953e-5];

var logGamma = function(xx) {
  var x = xx - 1.0;
  var tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  var ser = 1.000000000190015;
  for (var j = 0; j <= 5; j++) {
    x += 1;
    ser += gammaCof[j] / x;
  }
  return -tmp + Math.log(2.5066282746310005 * ser);
}

// HT https://en.wikipedia.org/wiki/Digamma_function#Computation_and_approximation
var digamma = function(x) {
  if (x < 6) {
    return digamma(x + 1) - 1 / x;
  }
  return Math.log(x) -
    1 / (2 * x) -
    1 / (12 * Math.pow(x, 2)) +
    1 / (120 * Math.pow(x, 4)) -
    1 / (252 * Math.pow(x, 6)) +
    1 / (240 * Math.pow(x, 8)) -
    5 / (660 * Math.pow(x, 10)) +
    691 / (32760 * Math.pow(x, 12)) -
    1 / (12 * Math.pow(x, 14));
}

var product = function(xs) {
  var ret = 1
  for (var i = 0; i < xs.length; i++) {
    ret *= xs[i]
  }
  return ret
}

var choose = function(n,k) {
  // optimizations for common cases
  if (k == 0) { return 1 }
  if (k == 1) { return n }
  if (n < k) { return 0 }
  var botMax = ((k > n-k) ? k : n-k);
  var botMin = ((k > n-k) ? n-k : k);
  var top = _.range(botMax+1, n+1)
  var bot = _.range(1, botMin + 1);
  return product(top)/product(bot)
}

var add = function(xs, ys) {
  return xs.map(function(x,i) {
    return x + ys[i]
  })
}


var mul = function(xs, ys) {
  return xs.map(function(x,i) {
    return x * ys[i]
  })
}

function randomExponential(rate) {
  // Allow to pass a random uniform value or function
  // Default to Math.random()
  var U = Math.random();

  return -Math.log(U)/rate;
}

function randomDiscrete(ps) {
  var bins = [];
  var runningTotal = 0;
  var U = Math.random();
  for(var i = 0; i < ps.length; i++) {
    var binLeft = runningTotal, binRight = runningTotal + ps[i];
    if (U > binLeft && U < binRight) {
      return i
    }
    runningTotal += ps[i]
  }

  return bins
}


var NS_PER_SEC = 1e9;
var ns_total = 0;

var gillespie = function(rates, state, tStart, tEnd, lawInputs, laws) {
  var jumps = [];
  var tAcc = tStart;

  var accState = state.slice()

  var numSpecies = laws[0].length
  var numLaws = laws.length

  while (tAcc < tEnd) {
    // calculate hazards and sum hazards in an inlined loop rather than function call
    var hazards = [];
    var totalHazard = 0;
    for(var i = 0; i < numLaws; i++) {
      var input = lawInputs[i]
      var haz = rates[i];
      for(var j = 0; j < numSpecies; j++) {
        haz *= choose(accState[j], input[j])
      }
      hazards.push(haz)
      totalHazard += haz
    }

    if (totalHazard == 0) {
      return jumps
    }
    var dt = randomExponential(totalHazard)
    var hazardsNormalized = hazards.map(function(h) { return h / totalHazard })
    var lawNum = randomDiscrete(hazardsNormalized)

    var stateUpdate = laws[lawNum];

    tAcc = tAcc + dt

    accState = add(accState, stateUpdate);

    var jump = {dt: dt,
                t: tAcc,
                lawNum: lawNum,
                state: global.T.fromScalars(accState)}

    jumps.push(jump)
  }

  return jumps
}

var nGillespie = function(n, rates, state, tStart, tEnd, lawInputs, laws) {
  var ret = []
  //var t = Date.now()
  ns_total = 0
  for(var i = 0; i < n; i++ ) {
    ret.push(gillespie(rates, state, tStart, tEnd, lawInputs, laws))
  }
  // console.log('total ' + (Date.now() - t))
  // console.log(ns_total / NS_PER_SEC)

  return ret;
}

module.exports = {
  logGamma: logGamma,
  digamma: digamma,
  gillespie: gillespie,
  nGillespie: nGillespie
}
