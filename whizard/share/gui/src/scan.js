const simulation = require('./simulation');

// Object structure:
// ScansList[i] : contain all Scan information for a i'th process
// ScansList[i].ScansContainer : contains union of subintervals for a i'th process
// ScansList[i].ScansContainer[j] : contains j'th individual subinterval for
//                          i'th process ex: (88.0 GeV => 90.0 GeV /+ 0.5 GeV)
export const ScansList = [];


function ScanElement(min, max, inc) {
  this.min = min;
  this.max = max;
  this.inc = inc;
}


function addScanElement(min, max, inc) {
  this.ScansContainer.push(new ScanElement(min, max, inc));
}


function ScanProcess() {
  this.ScansContainer = [];
  this.status = false;
  this.type = '';
  this.title = '';
  this.xlabel = '';
  this.ylabel = '';
  this.xmin = '';
  this.xmax = '';
  this.addScanElement = addScanElement;
}


export const Scan = {
  newProcess: () => {
    ScansList.push(new ScanProcess());
  },

  // Rebuilding Scan configuration for a selected process
  fillHTMLFields: () => {
    Scan.clean();
    for (let i = 0; i <
        ScansList[simulation.getActiveProcessId()].ScansContainer.length; i++) {
      const elem = ScansList[simulation.getActiveProcessId()].ScansContainer[i];
      Scan.addNew(elem.min, elem.max, elem.inc, i);
    }
    $('#conf-scan-check').prop('checked', ScansList[simulation.getActiveProcessId()].status);
    $('#scan-plot-title').val(ScansList[simulation.getActiveProcessId()].title);
    $('#scan-plot-xlabel').val(ScansList[simulation.getActiveProcessId()].xlabel);
    $('#scan-plot-ylabel').val(ScansList[simulation.getActiveProcessId()].ylabel);
    $('#scan-plot-xmin').val(ScansList[simulation.getActiveProcessId()].xmin);
    $('#scan-plot-xmax').val(ScansList[simulation.getActiveProcessId()].xmax);
  },

  getScanCodeHTML: (min, max, inc, sid) =>
    '<div class="row" rel="scan-elem" scanid="' + sid + '"> ' +
    '  <div class="col-sm-4"> ' +
    '    <div class="form-group"> ' +
    '      <input type="text" class="form-control conf-scan-min" ' +
    '      id="" placeholder="Min value" value="' + min + '"> ' +
    '    </div> ' +
    '  </div>  ' +
    '  <div class="col-sm-4"> ' +
    '    <div class="form-group"> ' +
    '      <input type="text" class="form-control conf-scan-max" ' +
    '      id="" placeholder="Max value" value="' + max + '"> ' +
    '    </div> ' +
    '  </div>         ' +
    '  <div class="col-sm-4"> ' +
    '    <div class="form-group"> ' +
    '      <input type="text" class="form-control conf-scan-inc" ' +
    '      id="" placeholder="Increment" value="' + inc + '"> ' +
    '    </div> ' +
    '  </div>  ' +
    '</div>',

  addNew: (min, max, inc, sid) => {
    $('#scansContainer').append(Scan.getScanCodeHTML(min, max, inc, sid));
  },

  clean: () => {
    $('#scansContainer').html('');
  },
};


export function setupJquery() {
  // Selecting: Scan > Process
  $(document).on('click', '.process-entry-scan', function selectScan() {
    $('.scan-right').fadeIn('fast');
    $('.process-entry-scan').removeClass('active');
    $(this).addClass('active');
    simulation.setActiveProcessId($(this).attr('process-id'));
    Scan.fillHTMLFields();
    if (ScansList[simulation.getActiveProcessId()].status) {
      $('#struct-scan').fadeIn('fast');
    } else {
      $('#struct-scan').fadeOut('fast');
    }
  });

  // Button: New Scan subinterval
  $('.scan-newscan').click(() => {
    Scan.addNew('', '', '',
        ScansList[simulation.getActiveProcessId()].ScansContainer.length);
    ScansList[simulation.getActiveProcessId()].addScanElement(0, 0, 0);
  });

  // Button: clean Scans
  $('.scan-clean').click(() => {
    ScansList[simulation.getActiveProcessId()].ScansContainer = [];
    Scan.clean();
  });

  // Checkbox button: Enable scan for this process
  $('#conf-scan-check').click(function enableScan() {
    ScansList[simulation.getActiveProcessId()].status = $(this).prop('checked');
    if ($(this).prop('checked')) {
      $('#struct-scan').fadeIn('fast');
    } else {
      $('#struct-scan').fadeOut('fast');
    }
    // Changing On/Off indicator
    if (ScansList[simulation.getActiveProcessId()].status) {
      $('#proc_indicator_scan_' + simulation.getActiveProcessId()).removeClass(
          'label-default label-success').addClass('label-success').text('On');
    } else {
      $('#proc_indicator_scan_' + simulation.getActiveProcessId()).removeClass(
          'label-default label-success').addClass('label-default').text('Off');
    }
  });

  // Modify field: Scan-minimum value
  $(document).on('change', '.conf-scan-min', function transferScanMin() {
    const min = $(this).val();
    const intervalID = $(this).parent().parent().parent().attr('scanid');
    ScansList[simulation.getActiveProcessId()].ScansContainer[intervalID].min = min;
  });

  // Modify field: Scan-maximum value
  $(document).on('change', '.conf-scan-max', function transferScanMax() {
    const max = $(this).val();
    const intervalID = $(this).parent().parent().parent().attr('scanid');
    ScansList[simulation.getActiveProcessId()].ScansContainer[intervalID].max = max;
  });

  // Modify field: Scan-inc value
  $(document).on('change', '.conf-scan-inc', function transferScanInc() {
    const inc = $(this).val();
    const intervalID = $(this).parent().parent().parent().attr('scanid');
    ScansList[simulation.getActiveProcessId()].ScansContainer[intervalID].inc = inc;
  });

  // Plot fields modification
  $('#scan-plot-title').change(() => {
    ScansList[simulation.getActiveProcessId()].title = $('#scan-plot-title').val();
  });

  $('#scan-plot-xlabel').change(() => {
    ScansList[simulation.getActiveProcessId()].xlabel = $('#scan-plot-xlabel').val();
  });

  $('#scan-plot-ylabel').change(() => {
    ScansList[simulation.getActiveProcessId()].ylabel = $('#scan-plot-ylabel').val();
  });

  $('#scan-plot-xmin').change(() => {
    ScansList[simulation.getActiveProcessId()].xmin = $('#scan-plot-xmin').val();
  });

  $('#scan-plot-xmax').change(() => {
    ScansList[simulation.getActiveProcessId()].xmax = $('#scan-plot-xmax').val();
  });
}
