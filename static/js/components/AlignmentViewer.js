import {Tooltip} from './Tooltip.js'
import {XYDispatch} from './PositionDispatch.js'
//import { MSAViewer, SequenceViewer, Labels, } from '@plotly/react-msa-viewer';
import { MSAViewer, SequenceViewer, Labels, } from './MSAV.umd.js';
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
        var style = document.querySelector('[data="rv3_style"]');
        style.innerHTML += ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}";

        var handleMolStarTopViewHovers = function (alnViewerClass, residueNumber){
            var alignmentNumber = Number(_.invert(vm.structure_mapping)[residueNumber]);
            var numVisibleTiles = Math.round(alnViewerClass.state.width/alnViewerClass.state.tileWidth);
            if (!isNaN(alignmentNumber)){
                if (alnViewerClass.state.aaPos > alignmentNumber || alignmentNumber > alnViewerClass.state.aaPos+numVisibleTiles){
                    let visiblePos = alignmentNumber-Math.round(numVisibleTiles/2);
                    if (visiblePos < 0) {visiblePos = 0};
                    alnViewerClass.setState({ aaPos: visiblePos })
                }
                alnViewerClass.highlightRegion({
                    sequences: {from: 0, to: vm.fastaSeqNames.length},
                    residues: {from: alignmentNumber, to: alignmentNumber}
                });
            }
        }
        window.addEventListener("resize", this.handleResize);
        document.addEventListener('PDB.topologyViewer.mouseover', (e) => {
            handleMolStarTopViewHovers(this, e.eventData.residueNumber);
        });
        document.addEventListener('PDB.topologyViewer.mouseout', () => {
            this.removeHighlightRegion();
        });
        document.addEventListener('PDB.molstar.mouseover', (e) => {
            handleMolStarTopViewHovers(this, e.eventData.seq_id);
        });
        document.addEventListener('PDB.molstar.mouseout', () => {
            this.removeHighlightRegion();
        });
    };
    componentWillUnmount() {
        window.removeEventListener("resize", this.handleResize);
    };
    onResidueMouseEnter = e => {
        this.highlightRegion({
            sequences: {from: 0, to: vm.fastaSeqNames.length},
            residues: {from: e.position+1, to: e.position+1}
        })
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
    highlightRegion = (highlight) => {
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
                    <Labels 
                      style = {{
                        width: (window.innerWidth - 300) * 0.2
                        }}
                    />
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
            <div style={{width: "30px"}}>
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
        </div>
        );
    }
}

export {AlnViewer}