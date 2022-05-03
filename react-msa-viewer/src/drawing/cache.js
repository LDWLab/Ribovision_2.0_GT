/**
* Copyright 2018, Plotly, Inc.
* All rights reserved.
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

class CanvasCharCache {
  constructor(g) {
    this.cache = {};
    this.cacheHeight = 0;
    this.cacheWidth = 0;
  }

  // returns a cached canvas
  getFontTile(letter, width, height, font) {
    // validate cache
    if (width !== this.cacheWidth || height !== this.cacheHeight || font !== this.font) {
      this.updateDimensions(width, height);
      this.font = font;
    }

    if (this.cache[letter] === undefined) {
      this.createTile(letter, width, height);
    }

    return this.cache[letter];
  }

  // creates a canvas with a single letter
  // (for the fast font cache)
  createTile(letter, width, height, font) {
    const canvas = this.cache[letter] = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    this.ctx = canvas.getContext('2d');
    this.ctx.font = this.font + "px mono";

    this.ctx.textBaseline = 'middle';
    this.ctx.textAlign = "center";

    return this.ctx.fillText(letter, width / 2, height / 2, width, font);
  }
  updateDimensions(width, height) {
    this.invalidate();
    this.cacheWidth = width;
    this.cacheHeight = height;
  }

  invalidate() {
    // TODO: destroy the old canvas elements
    this.cache = {};
  }
}

export default CanvasCharCache;
