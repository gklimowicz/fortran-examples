const promisify = require('es6-promisify');
const fs = require('fs');

export const mkdir = promisify(fs.mkdir, function cb(err) {
  if (err) {
    if (err.code === 'EEXIST') {
      return this.resolve('Folder ' + err.path + ' already exists');
    }
    return this.reject(err);
  }
  return this.resolve('Successfully created folder');
});

export const rmdir = promisify(fs.rmdir, function cb(err) {
  if (err) {
    if (err.code === 'ENOENT') {
      return this.resolve('Folder ' + err.path + ' does not exist');
    }
    return this.reject(err);
  }
  return this.resolve('Successfully removed folder');
});
