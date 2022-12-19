const cuts = require('./cuts');

let activeProcessId = -1;
export const SimulateList = [];


export function setActiveProcessId(id) {
  activeProcessId = id;
}


export function getActiveProcessId() {
  return activeProcessId;
}

function SimulateSetEvents(events) {
  this.events = events;
}


function SimulateGetEvents() {
  return this.events;
}


function SimulateSetStatus(status) {
  this.status = status;
}


function SimulateGetStatus() {
  return this.status;
}


function SindarinSimulateToString() {
  let s = '';
  if (this.Histogram.doHistogram) {
    s += '# New histogram \n';
    s += 'histogram distribution_' + this.procid + ' (' + this.Histogram.xmin +
        ', ' + this.Histogram.xmax + ', ' +
        (this.Histogram.xmax - this.Histogram.xmin) + ' / ' + this.Histogram.ticks
        + '.) {' + '\n';
    s += '$title = " ' + this.Histogram.title + '"' + '\n';
    s += '$x_label = "' + this.Histogram.xlabel + '"' + '\n';
    s += '}' + '\n';
    s += 'analysis = record distribution_' + this.procid + ' (eval ' +
        this.Histogram.analysis + ' ["' + this.Histogram.subevent + '"])' + '\n';
  }
  s += 'simulate(proc_' + this.procid + ') { n_events = ' +
    this.events + ' }' + '\n';
  if (this.Histogram.doHistogram) {
    s += 'compile_analysis' + '\n';
  }
  return s;
}


function setHistTitle(title) {
  this.title = title;
}


function setHistXmin(x) {
  this.xmin = x;
}


function setHistXmax(x) {
  this.xmax = x;
}


function setHistTicks(ticks) {
  this.ticks = ticks;
}


function setHistXLabel(l) {
  this.xlabel = l;
}


function setHistSubevent(s) {
  this.subevent = s;
}


function setHistStatus(s) {
  this.doHistogram = s;
}


function setHistAnalysis(s) {
  this.analysis = s;
}


function HistogramData() {
  this.doHistogram = false;
  this.xmin = -1;
  this.xmax = 1;
  this.ticks = 30;
  this.title = '';
  this.xlabel = '';
  this.subevent = '';
  this.analysis = 'cos (Theta)';

  this.setHistStatus = setHistStatus;
  this.setHistTitle = setHistTitle;
  this.setHistXLabel = setHistXLabel;
  this.setHistSubevent = setHistSubevent;
  this.setHistXmin = setHistXmin;
  this.setHistXmax = setHistXmax;
  this.setHistTicks = setHistTicks;
  this.setHistAnalysis = setHistAnalysis;
}


export function SindarinSimulate() {
  this.procid = SimulateList.length + 1;
  this.status = 0;
  this.events = 0;

  this.setEvents = SimulateSetEvents;
  this.getEvents = SimulateGetEvents;
  this.getStatus = SimulateGetStatus;
  this.setStatus = SimulateSetStatus;
  this.toString = SindarinSimulateToString;

  this.Histogram = new HistogramData();
  this.writeToSindarin = () => {
    if (this.getStatus() && this.getEvents() > 0) {
      return this.toString();
    }
    return '';
  };
}


function rebuildSimulateList() {
  let NewIndex = 1;
  for (let i = 0; i < SimulateList.length; i++) {
    if (SimulateList[i] instanceof SindarinSimulate) {
      if (SimulateList[i] === null) continue;
      SimulateList[i].procid = NewIndex++;
    }
  }
}


// Removes element from the list and rebuilds index list
export function removeSimulateElement(id) {
  SimulateList[id] = null;
  SimulateList.splice(id, 1);
  rebuildSimulateList();
}


export function addSimulation() {
  SimulateList.push(new SindarinSimulate());
  rebuildSimulateList();
}


export const Simulate = {
  rebuildParticlesHTML: () => {
    $('#pop_sim_subevent_list').html('');
    const particles = cuts.cutsClosure.getActiveParticles();
    for (let i = 0; i < particles.length; i++) {
      $('#pop_sim_subevent_list').append(
        '<li role="presentation"><a href="#" class="simulation-particles-click">'
        + particles[i] + '</a></li>');
    }
  },
  fillHistogramFieldsHTML: () => {
    $('#conf-sim-hist').prop('checked',
        SimulateList[activeProcessId].Histogram.doHistogram);
    if ($('#conf-sim-hist').prop('checked')) $('#struct-sim-hist').fadeIn('fast');
    else $('#struct-sim-hist').fadeOut('fast');
    $('#conf-sim-hist-title').val(SimulateList[activeProcessId].Histogram.title);
    $('#conf-sim-hist-x').val(SimulateList[activeProcessId].Histogram.xlabel);
    $('#conf-sim-hist-minx').val(SimulateList[activeProcessId].Histogram.xmin);
    $('#conf-sim-hist-maxx').val(SimulateList[activeProcessId].Histogram.xmax);
    $('#conf-sim-hist-ticks').val(SimulateList[activeProcessId].Histogram.ticks);
    $('#conf-sim-hist-subevent').val(
        SimulateList[activeProcessId].Histogram.subevent);
    $('#sim-hist-analysis').html(
        SimulateList[activeProcessId].Histogram.analysis +
        ' <span class="caret"></span>');
  },
};


export function setupJquery() {
  // [Appearance]
  // Selecting Tabs:Integration > Process
  $(document).on('click', '.process-entry-sim', function selectProcessEntry() {
    $('.process-entry-sim').removeClass('active');
    $(this).addClass('active');
  });

  //  Checkbox Simulation->Histograms
  $('#conf-sim-hist').change(function transferHistStatus() {
    if ($(this).prop('checked')) $('#struct-sim-hist').fadeIn('fast');
    else $('#struct-sim-hist').fadeOut('fast');
    SimulateList[activeProcessId].Histogram.setHistStatus($(this).prop('checked'));
  });

  // Clicking on Simulation->Histograms->Subevent->Particle
  $(document).on('click', '.simulation-particles-click', function transferSubevent() {
    const Input = $('#conf-sim-hist-subevent').val();
    $('#conf-sim-hist-subevent').val(Input + '' + $(this).text());
    SimulateList[activeProcessId].Histogram.setHistSubevent(
        $('#conf-sim-hist-subevent').val());
  });

  // Histogram Data from HTML forms
  $('#conf-sim-hist-title').change(function transferTitle() {
    SimulateList[activeProcessId].Histogram.setHistTitle($(this).val());
  });

  $('#conf-sim-hist-x').change(function transferXLabel() {
    SimulateList[activeProcessId].Histogram.setHistXLabel($(this).val());
  });

  $('#conf-sim-hist-subevent').change(() => {
    SimulateList[activeProcessId].Histogram.setHistSubevent(
        $('#conf-sim-hist-subevent').val());
  });

  $('#conf-sim-hist-minx').change(function transferXMin() {
    SimulateList[activeProcessId].Histogram.setHistXmin($(this).val());
  });

  $('#conf-sim-hist-maxx').change(function transferXMax() {
    SimulateList[activeProcessId].Histogram.setHistXmax($(this).val());
  });

  $('#conf-sim-hist-ticks').change(function transferHistTicks() {
    SimulateList[activeProcessId].Histogram.setHistTicks($(this).val());
  });

  $('.sim-hist-analysis').click(function transferAnalysis() {
    $('#sim-hist-analysis').html($(this).text() + ' <span class="caret"></span>');
    SimulateList[activeProcessId].Histogram.setHistAnalysis($(this).attr('key'));
  });

  // Simulate Data from HTML forms
  $('#conf-sim-sim').change(function transferIndicator() {
    SimulateList[activeProcessId].setStatus($(this).prop('checked'));
    // Changing On/Off indicator
    if (SimulateList[activeProcessId].status) {
      $('#proc_indicator_' + activeProcessId).removeClass(
          'label-default label-success').addClass('label-success').text('On');
    } else {
      $('#proc_indicator_' + activeProcessId).removeClass(
          'label-default label-success').addClass('label-default').text('Off');
    }
  });

  $('#conf-sim-events').change(function transferEvents() {
    SimulateList[activeProcessId].setEvents($(this).val());
  });

  // Hiding optional fields
  $('#struct-sim-hist').hide();
}
