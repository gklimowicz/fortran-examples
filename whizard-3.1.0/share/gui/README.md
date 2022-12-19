# WHIZARD-GUI

## Usage

Run `node lib/server.js` (or `nodejs lib/server.js` under Ubuntu or `npm start`)
in your terminal. Then go to your favorite browser and open `localhost:PORT`,
where PORT can be specified as environment variable `WHZ_GUI_PORT`, e.g.
`WHZ_GUI_PORT=2000 npm start`. If no port is given, the default value is 3000.

### Usage with ssh
In .ssh/config, add the line "LocalForward PORT localhost:PORT" to the entry
for your machine.

## Dependencies

### node

#### Ubuntu packages
Just `apt-get install npm nodejs`.

#### Locally from source
See `https://nodejs.org`.

### Packages from npm
`npm install` (in package directory, no arguments) will install dependencies in
local `node_modules` folder.

## Developer Guidelines

### Used technologies
- EJS - EmbeddedJavaScript is the view engine
- Bootstrap - ?
- Babel - allows to write in ES6, which is transpiled to ES5
- ESLint - configurable linter that helps to write correct code

### Directory structure
```
├── .babelrc            # configuration file for babel
├── .eslintrc           # configuration file for eslint
├── .istanbul.yml       # configuration file for coverage report
├── browser             # resulting files from webpack, used for client
├── lib                 # resulting files from babel, used for node and tests
├── src                 # we can use ES6 here, will be transpiled to lib
│   ├── guiconfig.js    # user configuration
│   ├── index.js        # express routes to serve sites
│   ├── server.js       # main that starts the server (use the one in lib though)
│   └── utils.js        # utility functions
├── node_modules        # this is created by npm for other modules
├── output              # output folder for sindarins and runs
├── package.json        # package dependencies for npm and scripts
├── public              # static information
│   ├── css             # stylesheet
│   ├── fonts           # extra fonts
│   └── images          # used images
├── test                # unit and integration tests
└── views               # pages that are interpreted by EJS
    ├── footer.ejs
    ├── header.ejs      # header loads css, fonts and sets meta data
    ├── index.ejs       # main page
    └── navbar.ejs
```

### Final tips
Be careful when changing version numbers in `package.json`. Things might break and
our tests are still insufficient to detect this fast.
