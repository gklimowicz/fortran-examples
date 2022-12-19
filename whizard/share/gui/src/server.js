const express = require('express');
const bodyParser = require('body-parser');
const log = require('js-logger');
const utils = require('./utils');
const context = require('./guiconfig').context; // This should parse a settings JSON/YAML
const app = express();
const startServer = () => {
  app.set('view engine', 'ejs');
  app.use('/public', express.static(process.cwd() + '/public'));
  app.use('/browser', express.static(process.cwd() + '/browser'));
  app.use('/output', express.static(process.cwd() + '/output'));
  // app.use(express.logger('dev'));
  app.use(bodyParser.json());       // to support JSON-encoded bodies
  app.use(bodyParser.urlencoded({     // to support URL-encoded bodies
    extended: true,
  }));
  app.get('/', (req, res) => {
    res.render('index');
  });
  app.listen(context.port, () => {
    console.warn('Listening on port: ' + context.port);
  });
};
log.useDefaults();
log.setLevel(context.logLevel);

utils.mkdir(context.outputDir)
  .then(log.info).catch(log.error);


startServer();
// let mongo = require('mongodb').MongoClient;
// const promisify = require('es6-promisify');
// const mongoPort = 27017;
  // console.log('MongoDB successfully connected on port: ' + mongoPort);
// It is probably better to save the state in a mongo db but it also adds a
// further dependency. For now disabled
// promisify(mongo.connect)('mongodb://localhost:' + mongoPort + '/whizard')
  // .then(function(mongodb) {
    // startServer(mongodb);
  // })
  // .catch(console.log.bind(console));

// Button: Run Sindarin
app.post('/runwhiz', (req, res) => {
  const fs = require('fs');
  console.log('Run Whizard clicked.');
  fs.writeFile(context.outputDir + '/gui-generated.sin', req.body.src, (err) => {
    if (err) {
      console.log(err);
    } else {
      console.log('The file was saved!');
    }
    return err;
  });
  const exec = require('child_process').exec;
  const cmd = ' cd ' + context.outputDir +
        ' && whizard ' + req.body.option + ' gui-generated.sin';
  console.log('Executing:' + cmd);
  exec(cmd, (error, stdout, stderr) => {
    if (stdout !== '') {
      console.log('---------stdout: ---------\n' + stdout);
    }
    if (stderr !== '') {
      console.log('---------stderr: ---------\n' + stderr);
    }
    if (error !== null) {
      console.log('---------exec error: ---------\n[' + error + ']');
    }
    res.end('finished');
  });
});


// Button: Save Sindarin
app.post('/savesin', (req, res) => {
  const fs = require('fs');
  console.log('Save Sindarin clicked.');
  fs.writeFile(context.outputDir + '/gui-generated.sin', req.body.src,
      (err) => {
        if (err) {
          res.end('Error saving.');
          console.log(err);
        } else {
          console.log('The file was saved!');
          res.end('Saved Succesfully.');
        }
        return err;
      });
});


// Action: Check file timestamp
app.post('/checktimestamp', (req, res) => {
  const fs = require('fs');
  fs.stat(req.body.filename, (err, stats) => {
    if (err) {
      res.end('Could not read filestamp.');
      console.log(err);
    } else {
      res.end(stats.mtime + '');
    }
    return err;
  });
});
