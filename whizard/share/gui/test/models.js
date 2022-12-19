const models = require('../lib/models.js');
const chai = require('chai');
const expect = require('chai').expect;
const chaiAsPromised = require('chai-as-promised');
chai.use(chaiAsPromised);

describe('sindarinModel', () => {
  const modelName = 'THDM';
  const description = 'Two-Higgs Doublet Model';
  it('should be constructable', () => {
    const test = new models.SindarinModel(modelName, description);
    expect(test.toString()).to.equal('model = THDM');
  });
});
