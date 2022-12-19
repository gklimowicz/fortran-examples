const cuts = require('../lib/cuts');
const chai = require('chai');
const chaiAsPromised = require('chai-as-promised');
chai.use(chaiAsPromised);

describe('cuts.cutsClosure', () => {
  it('should be empty', () => {
    const test = cuts.cutsClosure.getActiveParticles();
    return test === [];
  });
});

