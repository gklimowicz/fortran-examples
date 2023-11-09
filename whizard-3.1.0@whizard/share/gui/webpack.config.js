var path = require('path');
var webpack = require('webpack');
var build_dir = path.resolve(__dirname, 'browser');
var src_dir = path.join(__dirname, 'src');

module.exports = {
  // Makes sure errors in console map to the correct file and line number
  name: 'browser',
  entry: [
    'babel-polyfill',
    './src/main',
  ],
  output: {
    path: build_dir,
    publicPath: '/',
    filename: '[name].js',
  },
  devServer: {
    contentBase: build_dir,
  },
  debug: true,
  devtool: 'source-map',
  module: {
    loaders: [
      {
        test: src_dir,
        include: src_dir,
        loader: 'babel-loader',
        query: {
          plugins: ['transform-runtime'],
          presets: ['es2015'],
        },
      },
    ],
  },
  stats: {
    colors: true,
  },
};
