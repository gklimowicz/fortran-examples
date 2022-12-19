const alias = require('./alias');
const process = require('./process');
const generic = require('./generic');


export function SindarinCuts(cuts) {
  this.CutsData = cuts;
  this.toString = () => 'cuts = ' + this.CutsData;
  this.writeToSindarin = () => {
    if (this.CutsData.length > 0) {
      return this.toString();
    }
    return '# No cuts defined';
  };
}


// TODO: (bcn 2016-06-14) reduce complexity
export const cutsClosure = (() => {
  let LastClickedCutName;
  let LastClickedCutEq;
  let ActiveInputParticleElement = null;
  const CutNames = ['Pt', 'E', 'M', 'M2', 'abs(Eta)', 'abs(cos(Theta))'];
  const CutEq = ['>', '<'];

  const Public = {};
  // LastActiveInputElement (particles in cuts)
  Public.setLastActiveInputElement = (element) => {
    ActiveInputParticleElement = element;
  };
  Public.getLastActiveInputElement = () => ActiveInputParticleElement;

  // LastClickedCutName (name of cut)
  Public.setLastClickedCutName = (element) => {
    LastClickedCutName = element;
  };
  Public.getLastClickedCutName = () => LastClickedCutName;

  // LastClickedCutEq (inequality in cut)
  Public.setLastClickedCutEq = (element) => {
    LastClickedCutEq = element;
  };
  Public.getLastClickedCutEq = () => LastClickedCutEq;

  // Returns array of particles used in GUI
  Public.getActiveParticles = () => {
    let ParticlesList = [];

    // Get particles from alias list
    for (let i = 0; i < alias.aliasList.length; i++) {
      if (alias.aliasList[i] instanceof alias.SindarinAlias) {
        const thisAlias = alias.aliasList[i].alias;
        ParticlesList = generic.arrayUnique(ParticlesList.concat(thisAlias.split(':')));
      }
    }

    // Construct particle list from process definitions
    for (let i = 0; i < process.ProcessList.length; i++) {
      if (process.ProcessList[i] === null) continue; // Check if process was removed
      let Process = process.ProcessList[i].incoming + ' ' +
        process.ProcessList[i].outgoing;
      Process = Process.replace(/"/g, '').replace(/'/g, '')
        .replace(/\(|\)/g, '').replace(/,/g, ' ');
      Process = Process.split(' ').filter((n) => n !== '');
      ParticlesList = generic.arrayUnique(ParticlesList.concat(Process));
    }

    return ParticlesList;
  };

  // getCutsArray: Returns array of used cuts in sindarin format
  Public.getCutsArray = () => {
    const CutsList = [];
    $('[rel=cut-elem]').each(function getCuts(i, obj) { // eslint-disable-line no-unused-vars
      const cutName = $(this).find('.cut-name').text().replace(/ /g, '');
      const cutEq = $(this).find('.cut-eq').text().replace(/ /g, '');
      const cutVal = $(this).find('.cut-val').val();

      let cutAssignment = $(this).find('.cut-assignment').val();
      // Contains an array of individual particles used in cut
      if (cutAssignment === undefined) return;

      cutAssignment = cutAssignment.split(' ').remove(' ').remove('');
      // Constructing string of form: p1,p2,p3...
      if (cutAssignment.length && cutVal.length &&
          cutEq.length && cutName.length) {
        let ParticlesString = '';
        for (let j = 0; j < cutAssignment.length; j++) {
          ParticlesString += generic.parseParticleName(cutAssignment[j]) + ':';
        }
        ParticlesString = ParticlesString.slice(0, -1);

        // Cut format going into Whizard
        CutsList.push('all ' + cutName + ' ' + cutEq + ' ' +
            cutVal + ' [collect[' + ParticlesString + ']]');
      }
    });

    return CutsList;
  };

  // Rebuilds particles list html code in cuts tab
  Public.rebuildParticlesHTML = () => {
    // cleaning current list and building new
    $('#cuts-html-particles-list').html('');
    const particles = Public.getActiveParticles();
    for (let i = 0; i < particles.length; i++) {
      $('#cuts-html-particles-list').append('<li role="presentation">' +
          '<a href="#" class="cuts-particles-click">' +
          particles[i] + '</a></li>');
    }
  };

  // cleans all cuts data in html
  Public.clean = () => {
    $('#cutsContainer').html('');
  };

  Public.addNewCut = (name, eq, cutValue, cutAssignment) => {
    $('#cutsContainer').append(
        Public.getCutHTML(name, eq, cutValue, cutAssignment));
  };

  Public.getCutHTML = (name, eq, cutValue, cutAssignment) => {
    let CutCode = '<div class="row" rel="cut-elem"> ' +
                  '<div class="col-md-12">' +
                  '<div class="btn-group">' +
                  '<button type="button" class=' +
                  '"btn btn-default dropdown-toggle cut-name" ' +
                  'data-toggle="dropdown"> ' + name +
                  ' <span class="caret"></span> </button> ' +
                  '<ul class="dropdown-menu scrollable-menu" role="menu">';

    for (let i = 0; i < CutNames.length; i++) {
      CutCode += '<li><a href="#" class="cuts-select-name">' +
        CutNames[i] + '</a></li>';
    }

    CutCode += '<li><a href="#" class="cuts-select-delete">' +
      '<span class="glyphicon glyphicon-remove-sign" aria-hidden="true">' +
      '</span> Delete</a></li></ul></div><div class="btn-group">' +
      '<button type="button" class="btn btn-default dropdown-toggle cut-eq"' +
      'data-toggle="dropdown">' + eq + ' <span class="caret"></span>' +
      '</button><ul class="dropdown-menu" role="menu">';

    for (let i = 0; i < CutEq.length; i++) {
      CutCode += '<li><a href="#" class="cuts-select-eq">' +
        CutEq[i] + '</a></li>';
    }

    CutCode += '</ul> </div> <div class="btn-group">' +
      '<input type="text" class="form-control cut-val"' +
      'placeholder="500 GeV" value="' + cutValue + '"> </div>' +
      '<div class="btn-group pull-right"> <form class="form-inline">' +
      '<div class="form-group"><label for="exampleInputName2">Union:</label>' +
      '<input type="text" rel="" class="form-control cuts-particles-active ' +
      'cut-assignment" value="' + cutAssignment + '">' +
      '</div> </form> </div> </div> </div>';

    return CutCode;
  };

  return Public;
})();


export function setupJquery() {
  const cuts = cutsClosure;
  // Clicking on Cuts->Particles
  // Input is added for the last active input field
  $(document).on('click', '.cuts-particles-click', () => {
    const oldList = cuts.getLastActiveInputElement().val();
    cuts.getLastActiveInputElement().val(oldList + ' ' + $(this).text());
  });

  // Clicking on the input field making it active
  $(document).on('click', '.cuts-particles-active', () => {
    cuts.setLastActiveInputElement($(this));
  });

  // Checking focus to show/hide particle Cuts -> Particles List
  $(document).on('click', 'body', () => {
    if ($('.cuts-particles-active').is(':focus') ||
        $('.cuts-particles-click').is(':focus')) {
      $('#cuts-html-particles').fadeIn('fast');
    } else {
      $('#cuts-html-particles').fadeOut('fast');
    }
  });

  // Selecting cut name (Pt, M)
  $(document).on('click', '.cut-name', () => {
    cuts.setLastClickedCutName($(this));
  });

  $(document).on('click', '.cuts-select-name', () => {
    cuts.getLastClickedCutName().html($(this).text() + ' <span class="caret"></span>');
  });

  // Selecting cut inequality sign (>, <)
  $(document).on('click', '.cut-eq', () => {
    cuts.setLastClickedCutEq($(this));
  });

  $(document).on('click', '.cuts-select-eq', () => {
    cuts.getLastClickedCutEq().html($(this).text() + ' <span class="caret"></span>');
  });

  // Selecting to remove cut
  $(document).on('click', '.cuts-select-delete', () => {
    cuts.getLastClickedCutName().parent().parent().unbind().remove();
  });

  // Button: Cuts > New Cut
  $('.cuts-newcut').click(() => {
    cuts.addNewCut('Pt', '>', '', '');
  });
}
