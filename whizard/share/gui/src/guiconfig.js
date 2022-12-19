const log = require('js-logger');


function setDefaults() {
  const context = {};
  // Directory to store sindarin file and run whizard
  context.outputDir = 'output';

  // Default port
  context.port = '3000';

  if (typeof(process.env.WHZ_GUI_PORT) !== 'undefined') {
    context.port = process.env.WHZ_GUI_PORT;
  }
  context.port = Number(context.port);

  // Time in ms to update histograms/log file while whizard is running.
  context.updateTime = 10000;

  // Use Google Chart API for Latex images
  // Internet connection is required.
  context.useGoogleLatex = false;

  // Set to log.WARN for less verbose
  context.logLevel = log.DEBUG;
  return context;
}


export const context = setDefaults();
