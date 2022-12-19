import 'babel-polyfill';
import '../public/bootstrap.min';
const generic = require('./generic');
const cuts = require('./cuts');
// const guiconfig = require('./guiconfig');
const constructSindarin = require('./constructSindarin');
const models = require('./models');
import './examples';
const alias = require('./alias');
const process = require('./process');
const scan = require('./scan');
const simulation = require('./simulation');
const context = require('./guiconfig').context;

const ToolbarColumns = 4;

export function rebuildSindarin() {
  let SindarinList = [];
  const modelString = $('#conf-model').text();
  SindarinList.push(new models.SindarinModelData(modelString));

  if ($('#conf-additional').val()) {
    const AdditionalCode = new
      constructSindarin.SindarinAdditionalCode($('#conf-additional').val());
    SindarinList.push(AdditionalCode);
  }

  if ($('#conf-beams').val()) {
    SindarinList.push(new constructSindarin.SindarinAssignment('beams',
          $('#conf-beams').val() + ' => ' + $('#conf-pdf').text()));
  }

  const NewLineStarter = '\n\t and ';
  const cutsList = cuts.cutsClosure.getCutsArray();
  if (cutsList.length > 0) {
    let CutsRHS = '';
    for (let i = 0; i < cuts.Instance.length; i++) {
      CutsRHS += cutsList[i] + NewLineStarter;
    }
    CutsRHS = CutsRHS.substring(0, CutsRHS.length - NewLineStarter.length);
    SindarinList.push(new cuts.SindarinCuts(CutsRHS));
  }

  // Access Scans data
  process.extAssignScans();

  SindarinList = SindarinList.concat(alias.aliasList);
  SindarinList = SindarinList.concat(simulation.SimulateList);

  const a = new constructSindarin.SindarinGenerator(SindarinList, process.ProcessList);
  return a.construct();
}


export function cleanAll() {
  $('input[type="text"]').val('');
  $('#conf-additional').val('');
  alias.cleanAlias();
  cuts.cutsClosure.clean();
  scan.Scan.clean();
  process.ProcessList = [];
  simulation.SimulateList = [];
  scan.ScansList = [];
}


export function rebuildPreviewTab() {
  const SindarinScript = rebuildSindarin();
  $('#preview').html('<pre>' + SindarinScript + '</pre>');
}


function hideOptionalFields() {
  $('.outputcontainer').hide();
  $('#pbar').hide();
  $('#form-events, #form-calls, #form-iterations').hide();
  $('#gui-box').hide();
  $('.simulate-right, .integrate-right').hide();
  $('#struct-scan').hide();
}


$(document).ready(() => {
  alias.setupJquery();
  models.setupJquery(ToolbarColumns);
  cuts.setupJquery();
  simulation.setupJquery();
  scan.setupJquery();
  process.setupJquery();
  hideOptionalFields();

  // Selecting Tabs:Integration > Process
  $(document).on('click', '.process-entry', function processEntryHandler() {
    $('.process-entry').removeClass('active');
    $(this).addClass('active');
  });

  // Integrate checked, show #form-iterations and #for-calls
  $('#conf-integrate').change(function showIterationsCalls() {
    if ($(this).prop('checked')) {
      $('#form-iterations').fadeIn('fast');
      $('#form-calls').fadeIn('fast');
    } else {
      $('#form-iterations').fadeOut('fast');
      $('#form-calls').fadeOut('fast');
    }
  });

  $('#conf-int-nlo').change(function nloSetter() {
    try {
      if (simulation.getActiveProcessId() < 0) throw new {error: 'Please select a process'};
      if ($(this).prop('checked')) {
        process.ProcessList[simulation.getActiveProcessId()].setNlo(true);
      } else {
        process.ProcessList[simulation.getActiveProcessId()].setNlo(false);
      }
    } catch (err) {
      generic.messageGUI(err, 'alert-danger');
    }
  });

  $('#conf-int-sqrts').change(function sqrtsSetter() {
    try {
      if (simulation.getActiveProcessId() < 0) throw new {error: 'Please select a process'};
      process.ProcessList[simulation.getActiveProcessId()].setSqrts($(this).val());
    } catch (err) {
      generic.messageGUI(err, 'alert-danger');
    }
  });

  $('#conf-int-itt').change(function iterSetter() {
    try {
      if (simulation.getActiveProcessId() < 0) throw new {error: 'Please select a process'};
      process.ProcessList[simulation.getActiveProcessId()].setNIter($(this).val());
    } catch (err) {
      generic.messageGUI(err, 'alert-danger');
    }
  });

  $('#conf-int-cpi').change(function callsSetter() {
    try {
      if (simulation.getActiveProcessId() < 0) throw new {error: 'Please select a process'};
      process.ProcessList[simulation.getActiveProcessId()].setNCalls($(this).val());
    } catch (err) {
      generic.messageGUI(err, 'alert-danger');
    }
  });

  $(document).on('click', '.process-entry', function processEntryHandler() {
    simulation.setActiveProcessId($(this).attr('process-id'));
    const p = process.ProcessList[simulation.getActiveProcessId()];
    if (p.isNlo()) {
      $('#conf-int-nlo').prop('checked', true);
    } else {
      $('#conf-int-nlo').prop('checked', false);
    }
    $('#conf-int-itt').val(p.getNIter());
    $('#conf-int-cpi').val(p.getNCalls());
    $('#conf-int-sqrts').val(p.getSqrts());
    $('.integrate-right').fadeIn('fast');
  });

  $(document).on('click', '.process-entry-sim', function processEntrySimHandler() {
    simulation.setActiveProcessId($(this).attr('process-id'));
    const p = simulation.SimulateList[simulation.getActiveProcessId()];
    // Fill simulate fields
    $('#conf-sim-sim').prop('checked', p.getStatus());
    $('#conf-sim-events').val(p.getEvents());
    // Fill histogram fields
    simulation.Simulate.fillHistogramFieldsHTML();
    // Process selected show right column
    $('.simulate-right').fadeIn('fast');
  });

  // Simulate checked, show form-events
  $('#conf-simulate').change(function showFormEvents() {
    if ($(this).prop('checked')) {
      $('#form-events').fadeIn('fast');
    } else {
      $('#form-events').fadeOut('fast');
    }
  });

  // Tab preview clicked, generate script
  $('#tab_button_preview').click(() => {
    rebuildPreviewTab();
  });

  //  Changing tab, rebuild process list
  $('#tab_button_integrate, #tab_button_simulate, #tab_button_scan').click(() => {
    process.displayProcessList();
  });

  // Tab Cuts clicked, generate particles list
  $('#tab_button_cuts').click(() => {
    cuts.cutsClosure.rebuildParticlesHTML();
  });

  // Tab Simulate clicked, generate particles popup list
  $('#tab_button_simulate').click(() => {
    simulation.Simulate.rebuildParticlesHTML();
  });

  // Clicking on the model
  $(document).on('click', '.model', function showModelOnClick() {
    $('#conf-model').html($(this).text() + ' <span class="caret"></span>');
  });

  //  Remove Alias
  $(document).on('click', '.alias-remove', function removeAlias() {
    const id = $(this).attr('alias-id');
    alias.removeAlias(id);
  });

  // Button: Save Sindarin
  $(document).on('click', '.savesin', () => {
    const SindarinScript = rebuildSindarin();
    $.post('/savesin',
        {src: SindarinScript}, (data) => {
          generic.messageGUI(data, 'alert-success');
        });
  });

  // Button: Run Whizard
  $('.runwhiz').click(function whizardRunner() {
    // Run option: [--rebuild-events, --rebuild-grids, --rebuild]
    const option = (typeof $(this).attr('opt') !== 'undefined') ?
      $(this).attr('opt') : '';

    // Animation
    $('#pbar').show();
    $('#whizoutput').fadeOut('fast');
    $('.outputcontainer').fadeOut('fast');

    // Functionality
    $('.runwhiz, .runarrow').attr('disabled', 'disabled');
    const SindarinScript = rebuildSindarin();
    let whizRunning = true;
    generic.monitorLogChanges(whizRunning);

    // eslint-disable-next-line no-unused-vars
    $.post('/runwhiz', {src: SindarinScript, option}, (data) => {
      // Animation
      $('.outputcontainer').fadeIn('fast');
      $('#pbar').fadeOut('fast');

      // Functionality
      whizRunning = false;
      $('.runwhiz, .runarrow').removeAttr('disabled');
      // Display output whizard file
      $('#whizoutput').load('/' + context.outputDir + '/whizard.log').fadeIn('fast');
      // Display pdf (assuming there exists one for now)
      $('#out_hist').html(
          '<embed src="whizard_analysis.pdf" width="100%" height="700px">');

      // Whiz->GUI error parser
      // AM: check for other keywords, change method
      const CritKeyword = 'FATAL ERROR:';
      const str = $.ajax({url: '/' + context.outputDir + '/whizard.log',
        async: false}).responseText;

      if (str.indexOf(CritKeyword) > -1) {
        const s = str.substring(str.lastIndexOf(CritKeyword),
            str.lastIndexOf('*')).replace(/\*/g, '');
        if (s) generic.messageGUI(s, 'alert-danger');
      }
    });
  });

  // Design stuff
  $('[rel=popover]').popover({
    html: true,
    content: () => $('#pop_models').html(),
  });

  $('[rel=popover_aliases]').popover({
    html: true,
    content: () => $('#pop_aliases').html(),
  });

  $('[rel=popover_process]').popover({
    html: true,
    content: () => $('#pop_process').html(),
  });

  // Popover for histogram
  $('[data-toggle=popover_simulation_hist]').popover({
    html: true,
    content: () => $('#pop_sim_subevent').html(),
  });
});
