const utils = require('../lib/utils');
const expect = require('chai').expect;
const chai = require('chai');
const chaiAsPromised = require('chai-as-promised');
chai.use(chaiAsPromised);

describe('utils', () => {
  const testDir = 'test-dir';
  before(() => utils.rmdir(testDir));
  it('should create folders', () => {
    const promise = utils.mkdir(testDir);
    const result = expect(promise).to
      .eventually.equal('Successfully created folder');
    return result;
  });
  it('should fail to create them twice', () => {
    const promise = utils.mkdir(testDir);
    const result = expect(promise).to.eventually.equal('Folder ' + testDir +
        ' already exists');
    return result;
  });
  it('should remove folders', () => {
    const promise = utils.rmdir(testDir);
    const result = expect(promise).to.eventually.equal('Successfully removed folder');
    return result;
  });
  it('fail to remove them twice', () => {
    const promise = utils.rmdir(testDir);
    const result = expect(promise).to.eventually.equal('Folder ' + testDir +
        ' does not exist');
    return result;
  });
});
