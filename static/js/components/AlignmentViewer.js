import {Tooltip} from './Tooltip.js'
import {XSlider, YSlider} from './Sliders.js'
import {XYDispatch} from './PositionDispatch.js'
//import { MSAViewer, SequenceViewer, Labels, } from '@plotly/react-msa-viewer';
import { MSAViewer, 
        SequenceViewer, 
        Labels, 
        PositionBar, } from './MSAV.umd.js';
import React, { Component } from "react";

var AlnViewer = class RV3AlnViewer extends Component {
    state = { 
        tileWidth: 17,
        tileHeight: 17,
        aaPos: 0,
        seqPos: 0,
        width: (window.innerWidth - 300) * 0.7,
        height: ((window.innerHeight - 171)/2) * 0.8,
        highlight: null,
        colorScheme: vm.colorScheme,
    };
    handleResize = () => {
        this.setState({
            width: (window.innerWidth - 300) * 0.7,
            height: ((window.innerHeight - 171)/2) * 0.8
        });
        var style = document.querySelector('[data="rv3_style"]');
        style.innerHTML += ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}"
    };
    componentDidMount() {
        vm.colorScheme = 'clustal2';
        var style = document.querySelector('[data="rv3_style"]');
        style.innerHTML += ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}";
        window.ajaxRun = false;
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
        $('#alnSequenceViewer').mouseleave(function () {
            window.ajaxRun = false;
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
        if (vm.topology_loaded){
            let resiPos = vm.structure_mapping[e.position+1];
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
        //if (prop_loaded){
        //    let seqIx = vm.fastaSeqNames.indexOf(e.sequence.name.split(' ').slice(1,).join(' '))
        //    Plotly.Fx.hover('total',[
        //        {curveNumber:5, pointNumber:5},
        //        {curveNumber:6, pointNumber:6}
        //    ]);
        //}
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
        if (vm.topology_loaded){
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
        var maxXpos = vm.fasta_data.split('\n')[1].length - Math.round((((window.innerWidth - 300) * 0.7)/this.state.tileWidth))+2;
        var maxYpos = vm.fastaSeqNames.length - Math.round(((((window.innerHeight - 171)/2) * 0.8)/this.state.tileHeight))+2;
        var alnViewerAdjHeight = ((window.innerHeight - 171)/2) * 0.8;
        var alnViewerAdjWidth = (window.innerWidth - 300) * 0.7;
        if (maxYpos < 0){ maxYpos = 0 };
        if (maxXpos < 0){ maxXpos = 0 };
        return (
        <div style={{ display: "flex" }}>
            <div>
                <XSlider 
                  alnViewerAdjWidth={alnViewerAdjWidth}
                  maxXpos={maxXpos}
                  MSAVObject={this}
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
                  colorScheme={this.state.colorScheme}
                >
                <div style={{ position: "relative", display: "flex", height:this.state.height+this.state.tileHeight}}>
                        <Labels 
                          id="alnViewerLabels"
                          style = {{
                            width: (window.innerWidth - 300) * 0.2,
                            paddingTop: 13.6,
                            marginRight: 3,
                            }}
                        />
                    <div>
                        <PositionBar 
                          markerSteps={5} 
                          startIndex={0} 
                        />
                        <SequenceViewer
                          id="alnSequenceViewer"
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
            <YSlider 
              alnViewerAdjHeight={alnViewerAdjHeight}
              maxYpos={maxYpos}
              MSAVObject={this}
            />
        </div>
        );
    }
}

export {AlnViewer}