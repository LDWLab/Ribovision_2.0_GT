/**
* Copyright 2018, Plotly, Inc.
* All rights reserved.
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

// TODO: use more abstract, non-canvas-like API
class DrawingBase {
  constructor(el) {
    this.el = el;
    this.state = {};
  }

  updateEl(el) {
    this.el = el;
  }

  // props
  font(fontName) {
    this.state.font = fontName;
  }

  globalAlpha(globalAlpha) {
    this.state.globalAlpha = globalAlpha;
  }

  startDrawingFrame() {
    this.clear();
  }
  endDrawingFrame() {}
  save(){}
  restore(){}
}

export default DrawingBase;
