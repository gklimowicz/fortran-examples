const generic = require('./generic');


export let aliasList = [];


export function SindarinAlias(str, alias) {
  this.name = str;
  this.alias = alias;
  this.writeToSindarin = () => 'alias ' + this.name + ' = ' + this.alias;
}


function rebuildAliasList() {
  $('#pop_aliases').empty();
  $('#pop_aliases').append('<div class="row">');
  for (let i = 0; i < aliasList.length; i++) {
    const alias = aliasList[i].writeToSindarin();
    $('#pop_aliases').append('<div class="col-md-10">' +
        '<a class="alias">' + alias + '</a></div>' +
        '<div class="col-md-2"><a href="javascript:;" class="alias-remove" alias-id='
        + i + '><span class="glyphicon glyphicon-remove-sign" ' +
        'aria-hidden="true"></span></a></div>');
    $('#pop_aliases').append('</div>');
  }
}


export function addAlias(name, str) {
  aliasList.push(new SindarinAlias(name, str));
  rebuildAliasList();
}


export function removeAlias(id) {
  if (parseInt(Number(id), 10) == id) {  // eslint-disable-line eqeqeq
    aliasList.splice(id, 1);
    rebuildAliasList();
  } else {
    console.error('Did not get an integer id but ', id);
  }
}


export function cleanAlias() {
  aliasList = [];
  rebuildAliasList();
}


export function setupJquery() {
  $('#button_alias').click(() => {
    // Checking if both fields are non-empty
    if ($('#conf-alias-lhs').val() && $('#conf-alias-rhs').val()) {
      addAlias($('#conf-alias-lhs').val(), $('#conf-alias-rhs').val());
      generic.messageGUI('New alias is added.', 'alert-success');
      // $('#conf-alias-lhs').val('');
      // $('#conf-alias-rhs').val('');
    }
  });
}
