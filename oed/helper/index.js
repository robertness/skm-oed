

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

var gillespie = function(rates, vState, tStart, tEnd,
                         lawInputs, speciesNums, lawNums, vLaws
                        ) {
  var jumps = [];
  var tAcc = tStart;

  var accState = vState;

  var getGValue = function(i) {
    var input = lawInputs[i]
    return product(speciesNums.map(function(j) {
      var state_j = global.T.get(accState, j)
      var input_j = input[j]
      return choose(state_j, input_j)
    }))
  }

  while (tAcc < tEnd) {

    var gValues = lawNums.map(getGValue)

    var hazards = mul(rates, gValues),
        totalHazard = global._.sum(hazards);
    if (totalHazard == 0) {
      return jumps
    }
    var dtDist = new global.dists.Exponential({a: totalHazard})
    var dt = dtDist.sample()
    var lawDist = new global.dists.Discrete({ps: hazards})
    var lawNum = lawDist.sample()
    var stateUpdate = vLaws[lawNum];
    accState = global.T.add(accState, stateUpdate);

    tAcc = tAcc + dt

    var jump = {dt: dt,
                t: tAcc,
                lawNum: lawNum,
                prevState: vState,
                state: accState,
                gValues: gValues
               }
    jumps.push(jump)
  }

  return jumps
}

module.exports = {
  logGamma: logGamma,
  digamma: digamma,
  gillespie: gillespie
}
