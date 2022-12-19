const main = require('./main');
const generic = require('./generic');
const process = require('./process');

// Pair Production ee->WW
$('.ex-eeww').click(() => {
  main.cleanAll();
  process.addProcess('"e+", "e-"', '"W+", "W-"');
  process.ProcessList[0].setSqrts(500);
  generic.messageGUI('Example loaded.', 'alert-success');
});
