import {Tooltip} from './Tooltip.js'
import {XYDispatch} from './PositionDispatch.js'
import { MSAViewer, SequenceViewer, Labels, } from '@plotly/react-msa-viewer';
import React, { Component } from "react";

var AlnViewer = class RV3AlnViewer extends Component {
  state = { 
      tileWidth: 18,
      tileHeight: 18,
      aaPos: 0,
      seqPos: 0,
      width: (window.innerWidth - 300) * 0.7,
      height: ((window.innerHeight - 171)/2) * 0.9,
      highlight: null 
  };
  handleResize = () => {
      this.setState({
          width: (window.innerWidth - 300) * 0.7,
          height: ((window.innerHeight - 171)/2) * 0.9
      });
      var style = document.querySelector('[data="rv3_style"]');
      style.innerHTML += ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}"
  };
  componentDidMount() {
      window.addEventListener("resize", this.handleResize);
      var style = document.querySelector('[data="rv3_style"]');
      style.innerHTML += ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}"
  };
  componentWillUnmount() {
      window.removeEventListener("resize", this.handleResize);
  };
  onResidueMouseEnter = e => {
      if (vm.topology_loaded == 'True'){
          let resiPos = vm.structure_mapping[e.position];
          if (resiPos !== undefined){
              viewerInstanceTop.pluginInstance.highlight(resiPos, resiPos);
              viewerInstance.visual.highlight({
                  data:[{
                          entity_id:vm.entityId,
                          start_residue_number:resiPos,
                          end_residue_number:resiPos,
                      },],
              });
          }
      }
      if (!window.ajaxRun){
          
          if (e.position !== undefined){
              window.ajaxRun = true;
              registerHoverResiData(e, this, vm);
          }
      } else {
          this.setState({ fold: undefined, phase: undefined });
      }
   };
  onResidueMouseLeave = e => {
      if (vm.topology_loaded == 'True'){
          viewerInstanceTop.pluginInstance.clearHighlight();
          viewerInstance.visual.clearHighlight();
      }
      this.setState({ fold: undefined, phase: undefined });
  };
  highlightRegion = () => {
      const highlight = {
          sequences: {
            from: 0,
            to: 2
          },
          residues: {
            from: 2,
            to: 13
          }
        };
      this.setState({ highlight });
  };
  removeHighlightRegion = () => {
      this.setState({ highlight: null });
  };
  render() {
      const xPos = this.state.tileWidth * (this.state.aaPos);
      const yPos = this.state.tileHeight * (this.state.seqPos);
      const maxXpos = window.aaFreqs.length - Math.round((((window.innerWidth - 300) * 0.7)/this.state.tileWidth))+2;
      const maxYpos = vm.fastaSeqNames.length - Math.round(((((window.innerHeight - 171)/2) * 0.9)/this.state.tileHeight))+2;
      return (
      <div style={{ display: "flex" }}>
          <div>
            <input
              style = {{ 
                  width: (window.innerWidth - 300) * 0.7+"px",
                  position: "relative",
                  left: (window.innerWidth - 300) * 0.2+"px"
              }}
              type="range"
              min="0"
              max={maxXpos}
              value={this.state.aaPos}
              onChange={(evt) => this.setState({ aaPos: evt.target.value })}
              className="slider"
              id="xPosSlider"
              />
          <MSAViewer 
          {...msaOptions}
          ref={(ref) => (this.el = ref)}
          highlight={this.state.highlight}
          width={this.state.width}
          height={this.state.height}
          tileWidth={this.state.tileWidth}
          tileHeight={this.state.tileHeight}
          position={{ xPos, yPos }}
          >
          <div style={{ position: "relative", display: "flex"}}>
          <Labels style={{
              width: (window.innerWidth - 300) * 0.2
              }}/>
          <div>
              <SequenceViewer
                onResidueMouseEnter={this.onResidueMouseEnter}
                onResidueMouseLeave={this.onResidueMouseLeave}
              />
              {this.state.fold && (
                <div
                  style={{
                    position: "absolute",
                    opacity: 0.8,
                    ...this.state.tooltipPosition,
                  }}
                >
                  <Tooltip>
                    Fold: {this.state.fold} <br></br>
                    Phase: {this.state.phase}
                  </Tooltip>
                </div>
              )}
              </div>
          </div>
          <XYDispatch parent_state={this.state} />
          </MSAViewer>
      </div>
      <input
          style={{ 
              width: ((window.innerHeight - 171)/2)*0.9+"px",
          }}
          type="range"
          min="0"
          max={maxYpos}
          value={this.state.seqPos}
          onChange={(evt) => this.setState({ seqPos: evt.target.value })}
          className="slider"
          id="yPosSlider"
          />
      </div>
      );
  }
}

export {AlnViewer}