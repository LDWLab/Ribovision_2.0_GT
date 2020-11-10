var path = require("path");
var webpack = require('webpack');
var BundleTracker = require('webpack-bundle-tracker');
const VueLoaderPlugin = require('vue-loader/lib/plugin')
var CleanObsoleteChunks = require('webpack-clean-obsolete-chunks');

module.exports = {
  context: __dirname,
  devtool: "source-map",

  entry: './static/js/index',

  output: {
      path: path.resolve('./static/bundles/'),
      filename: "[name]-[hash].js",
  },

  plugins: [
    new BundleTracker({filename: './webpack-stats.json'}),
    new VueLoaderPlugin(),
    new CleanObsoleteChunks({verbose: false,})
  ],
  module: {
    rules: [
      {
        test: /\.js$/,
        exclude: /node_modules/,
        use: ['babel-loader']
      },{
        test: /\.vue$/,
        loader: 'vue-loader'
      },
    ]
  },
  resolve: {
    extensions: ['*', '.js', '.jsx'],
    alias: {
      'vue$': 'vue/dist/vue.esm.js' // 'vue/dist/vue.common.js' for webpack 1
    }
  }

};