const simulation = require('./simulation');
const generic = require('./generic');
const scan = require('./scan');

export const ProcessList = [];


function SimulationData() {
  this.nEvents = 0;
}


function IntegrationData() {
  this.nlo = false;
  // Default values
  this.sqrts = 0;
  this.nCalls = 10000;
  this.nIter = 5;
}


function ScanData() {
  this.Sets = [];
  this.type = '';
  this.title = '';
  this.xlabel = '';
  this.ylabel = '';
  this.xmin = '';
  this.xmax = '';
}


function SindarinProcessIsNlo() {
  return this.integrationData.nlo;
}


function SindarinProcessSetNlo(value) {
  this.integrationData.nlo = value;
}


function SindarinProcessGetSqrts() {
  return this.integrationData.sqrts;
}


function SindarinProcessSetSqrts(sqrts) {
  this.integrationData.sqrts = sqrts;
}


function SindarinProcessGetNIter() {
  return this.integrationData.nIter;
}


function SindarinProcessSetNIter(nIt) {
  this.integrationData.nIter = nIt;
}


function SindarinProcessGetNCalls() {
  return this.integrationData.nCalls;
}


function SindarinProcessSetNCalls(nCalls) {
  this.integrationData.nCalls = nCalls;
}


function SindarinProcessSetNEvents(nEv) {
  this.simulationData.nEvents = nEv;
}


function SindarinProcessGetNEvents() {
  return this.simulationData.nEvents;
}


// TODO: (bcn 2016-07-09) have to test this together with scan.ScansList
export function extAssignScans() {
  for (let i = 0; i < ProcessList.length; i++) {
    ProcessList[i].importScanData(i);
  }
}


function SindarinProcessToString() {
  try {
    if (this.counter <= 0) {
      throw new Error('Invalid process counter');
    }
    let str = 'process proc_' + this.counter + ' = '
      + this.incoming + ' => ' + this.outgoing;
    if (this.isNlo()) {
      str += ' {nlo_calculation = "Full"}';
    }
    return str;
  } catch (err) {
    console.error(err);
  }
}


function ProcessDisplayName() {
  return 'proc_' + this.counter + ' = ' + this.incoming + ' => ' + this.outgoing;
}


function importScanData(processID) {
  this.ScanData.Sets = scan.ScansList[processID].ScansContainer;
  this.ScanData.type = scan.ScansList[processID].type;
  this.ScanData.title = scan.ScansList[processID].title;
  this.ScanData.xlabel = scan.ScansList[processID].xlabel;
  this.ScanData.ylabel = scan.ScansList[processID].ylabel;
  this.ScanData.xmin = scan.ScansList[processID].xmin;
  this.ScanData.xmax = scan.ScansList[processID].xmax;
}


function processWriteToSindarin(i) {
  let str = this.toString() + '\n';
  // If scans defined overwrite
  if (this.ScanData.Sets.length > 0) {
    str += '#Plot data' + '\n';
    str += '$x_label = "' + this.ScanData.xlabel + '"' + '\n';
    str += '$y_label = "' + this.ScanData.ylabel + '"' + '\n';
    str += '$title = "' + this.ScanData.title + '"' + '\n';
    str += 'plot lineshape_' + i + ' { x_min = ' +
      this.ScanData.xmin + ' x_max = ' + this.ScanData.xmax + ' }' + '\n';
    str += 'scan sqrts = (';
    for (let j = 0; j < this.ScanData.Sets.length; j++) {
      const e = this.ScanData.Sets[j];
      str += '(' + e.min + ' => ' + e.max + ' /+ ' + e.inc + '),';
    }
    str = str.substring(0, str.length - 1);
    str += ') {' + '\n';
    str += '\tintegrate (proc_' + this.counter +
        ') { iterations = 2:1000:"gw", 1:2000 }' + '\n';
    str += 'record lineshape_' + i +
      '(sqrts, integral (proc_' + i + ') / 1000)' + '\n';
    str += '}\ncompile_analysis ';
  } else {
    str += 'sqrts = ' + this.getSqrts() + '\n';
    if (this.getNCalls() > 0 && this.getNIter() > 0) {
      str += 'integrate(proc_' + this.counter + ') {iterations=' +
        this.getNIter() + ':' + this.getNCalls() + ':"gw"}\n';
    }
  }
  if (this.getNEvents() > 0) {
    str += 'simulate(proc_' + this.counter +
        ') {n_events=' + this.getNEvents() + '}\n';
  }
  return str;
}


export function SindarinProcess(incoming, outgoing) {
  this.counter = 0;
  this.incoming = incoming;
  this.outgoing = outgoing;
  this.integrationData = new IntegrationData();
  this.simulationData = new SimulationData();
  this.ScanData = new ScanData();
  this.isNlo = SindarinProcessIsNlo;
  this.setNlo = SindarinProcessSetNlo;
  this.setSqrts = SindarinProcessSetSqrts;
  this.getSqrts = SindarinProcessGetSqrts;
  this.getNIter = SindarinProcessGetNIter;
  this.setNIter = SindarinProcessSetNIter;
  this.getNCalls = SindarinProcessGetNCalls;
  this.setNCalls = SindarinProcessSetNCalls;
  this.getNEvents = SindarinProcessGetNEvents;
  this.setNEvents = SindarinProcessSetNEvents;
  this.toString = SindarinProcessToString;
  this.writeToSindarin = processWriteToSindarin;
  this.name = ProcessDisplayName;
  this.importScanData = importScanData;
}


export function rebuildProcessList() {
  $('#pop_process').empty();
  $('#pop_process').append('<div class="row">');
  let procIndex = 1;
  for (let i = 0; i < ProcessList.length; i++) {
    if (ProcessList[i] instanceof SindarinProcess) {
      ProcessList[i].counter = procIndex;
      procIndex++;
      $('#pop_process').append(
          '<div class="col-md-10"><a href="javascript:;" class="process">' +
          ProcessList[i].name() + '</a></div><div class="col-md-2">' +
          '<a href="javascript:;" class="process-remove" processid=' + i +
          '><span class="glyphicon glyphicon-remove-sign" aria-hidden="true">' +
          '</span></a></div>');
    }
    $('#pop_process').append('</div>');
  }
}


// Add a new process
// TODO: (bcn 2016-06-26) Adding a process should not change simulation and scan
export function addProcess(incoming, outgoing) {
  ProcessList.push(new SindarinProcess(incoming, outgoing));
  simulation.addSimulation();
  scan.Scan.newProcess();
  rebuildProcessList();
}


function suggestAddingProcessIfNoneAdded(list) {
  if (list.filter((value) => value !== null).length === 0) {
    $('#simulate-process-list').html('Please add a process.');
    $('#integrate-process-list').html('Please add a process.');
    $('#scan-process-list').html('Please add a process.');
    $('.simulate-right, .integrate-right, .scan-right').hide();
  }
}

function constructProcessList(processList, list, jQuerySelector, type, id) {
  $(jQuerySelector).empty();
  for (let i = 0; i < processList.length; i++) {
    const CSSClass = list[i].status ? 'label-success' : 'label-default';
    const Text = list[i].status ? 'On' : 'Off';
    const name = generic.texImageOrPlain(processList[i].name());
    $(jQuerySelector).append(
        '<a href="#" class="list-group-item process-entry-' + type + '" process-id="' +
        i + '">' + name + '<br><span id="' + id + '_' + i +
        '" class="label ' + CSSClass + '">' + Text + '</span></a>');
  }
}

function constructIntegrationList(processList) {
  $('#integrate-process-list').empty();
  for (let i = 0; i < processList.length; i++) {
    const name = generic.texImageOrPlain(ProcessList[i].name());
    $('#integrate-process-list').append(
        '<a href="#" class="list-group-item process-entry" process-id="' +
        i + '">' + name + '</a>');
  }
}

// Generate process list to choose setups from
export function displayProcessList() {
  constructIntegrationList(ProcessList);
  constructProcessList(ProcessList, simulation.SimulateList,
      '#simulate-process-list', 'sim', 'proc_indicator');
  constructProcessList(ProcessList, scan.ScansList,
      '#scan-process-list', 'scan', 'proc_indicator_scan');
  suggestAddingProcessIfNoneAdded(ProcessList);
}


export function setupJquery() {
  $('#button-add-process').click(() => {
    const incoming = $('#conf-process-in').val();
    const outgoing = $('#conf-process-out').val();
    if (incoming && outgoing) {
      addProcess(generic.parseParticleNameString(incoming),
        generic.parseParticleNameString(outgoing));
      generic.messageGUI('New process is added.', 'alert-success');
    } else {
      let incomingMissing = '';
      let outgoingMissing = '';
      if (!incoming) incomingMissing = 'No incoming particles';
      if (!outgoing) outgoingMissing = 'No outgoing particles';
      generic.messageGUI('Adding process failed! ' + incomingMissing + ' and ' +
        outgoingMissing, 'alert-danger');
    }
  });

  // Mini-button: Remove process
  $(document).on('click', '.process-remove', function processRemoveHandler() {
    const id = $(this).attr('process-id');
    ProcessList.splice(id, 1);
    rebuildProcessList();
    simulation.removeSimulateElement(id);
  });
}
