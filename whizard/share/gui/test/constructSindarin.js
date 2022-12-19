const construct = require('../lib/constructSindarin.js');
const expect = require('chai').expect;

describe('SindarinAssignment', () => {
  it('should be constructable', () => {
    const testSindarinAssignment = new construct.SindarinAssignment('foo', 17);
  });
  it('should write the expected string', () => {
    const testSindarinAssignment = new construct.SindarinAssignment('foo', 17);
    expect(testSindarinAssignment.writeToSindarin()).to.equal('foo = 17');
  });
});

describe('SindarinCommand', () => {
  it('should just write the input string (no logic here)', () => {
    const testSindarinCommand = new construct.SindarinCommand('foo');
    expect(testSindarinCommand.writeToSindarin()).to.equal('foo');
  });
});

describe('SindarinAdditionalCode', () => {
  it('should just write the input string (no logic here)', () => {
    const testSindarinAdditionalCode = new construct.SindarinAdditionalCode('foo');
    expect(testSindarinAdditionalCode.writeToSindarin()).to.equal('foo');
  });
});

describe('SindarinWriteHeader', () => {
  it('should write a header', () => {
    const testSindarinWriteHeader = construct.sindarinWriteHeader();
    expect(testSindarinWriteHeader).to.contain('automatically generated');
  });
});

describe('SindarinGenerator', () => {
  const testSindarinList = [new construct.SindarinAdditionalCode('foo')];
  const testProcessList = [new construct.SindarinAdditionalCode('bar')];
  it('should be constructable', () => {
    const testSindarinGenerator = new
      construct.SindarinGenerator(testSindarinList, testProcessList);
  });
  it('should be construct a sindarin', () => {
    const testSindarinGenerator = new
      construct.SindarinGenerator(testSindarinList, testProcessList);
    const result = testSindarinGenerator.construct();
    expect(result).to.contain('Date:');
    expect(result).to.contain('bar');
    expect(result).to.contain('foo');
  });
});
