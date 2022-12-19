const process = require('../lib/process.js');
const expect = require('chai').expect;
const chai = require('chai');
const chaiAsPromised = require('chai-as-promised');
chai.use(chaiAsPromised);

describe('SindarinProcess', () => {
  it('should be constructable', () => {
    const incoming = 'e1, E1';
    const outgoing = 'e2, E2';
    const sindarinProcess = new process.SindarinProcess(incoming, outgoing);
    expect(sindarinProcess.counter).to.equal(0);
    expect(sindarinProcess.incoming).to.equal(incoming);
    expect(sindarinProcess.outgoing).to.equal(outgoing);
    expect(sindarinProcess.isNlo()).to.equal(false);
    sindarinProcess.setNlo(true);
    expect(sindarinProcess.isNlo()).to.equal(true);
    sindarinProcess.setNEvents(17);
    expect(sindarinProcess.getNEvents()).to.equal(17);
  });
  it('should be able to write itself as sindarin', () => {
    const incoming = 'e1, E1';
    const outgoing = 'e2, E2';
    const sindarinProcess = new process.SindarinProcess(incoming, outgoing);
    sindarinProcess.counter = 2;
    sindarinProcess.setNEvents(3);
    const expectName = 'proc_2 = e1, E1 => e2, E2';
    expect(sindarinProcess.name()).to.equal(expectName);
    sindarinProcess.setSqrts('500 GeV');
    sindarinProcess.setNIter(5);
    sindarinProcess.setNCalls('100');
    expect(sindarinProcess.writeToSindarin()).to.equal('process ' + expectName +
        '\nsqrts = 500 GeV' +
        '\nintegrate(proc_2) {iterations=5:100:"gw"}' +
        '\nsimulate(proc_2) {n_events=3}\n');
  });
});

describe('ProcessList', () => {
  it('should allow to add a process', () => {
    process.addProcess('foo', 'bar');
  });
  it('should be rebuildable', () => {
    process.rebuildProcessList();
  });
  it('should be displayable', () => {
    process.displayProcessList();
  });
});
