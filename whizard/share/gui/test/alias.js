const alias = require('../lib/alias.js');
const expect = require('chai').expect;
const chai = require('chai');
const chaiAsPromised = require('chai-as-promised');
chai.use(chaiAsPromised);
GLOBAL.$ = require('jquery');
// const chaiJquery = require('chai-jquery');

describe('SindarinAlias', () => {
  const testAliasName = 'quark';
  const testAlias = 'u:U';
  it('should be constructable', () => {
    const sindarinAlias = new alias.SindarinAlias(testAliasName, testAlias);
    expect(sindarinAlias.alias).to.equal(testAlias);
    expect(sindarinAlias.alias).to.equal(testAlias);
  });
  it('should write to the expected string', () => {
    const sindarinAlias = new alias.SindarinAlias(testAliasName, testAlias);
    expect(sindarinAlias.writeToSindarin()).to.equal('alias quark = u:U');
  });
  it('should use jquery to build an empty alias list', () => {
    alias.cleanAlias();
    // expect($('#pop_aliases').get).to.have.length.above(0);
    // expect($('#pop_aliases').val()).to.contain('text');
  });
  it('should be able to add an alias', () => {
    alias.addAlias('foo', 'bar');
    expect(alias.aliasList).to.have.length(1);
    alias.cleanAlias();
    expect(alias.aliasList).to.have.length(0);
  });
  it('should be able to remove an alias', () => {
    alias.addAlias('foo', 'bar');
    alias.removeAlias(0);
    expect(alias.aliasList).to.have.length(0);
  });
});
