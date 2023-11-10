import {Tooltip} from './Tooltip.js'
import {XSlider, YSlider} from './Sliders.js'
import {XYDispatch} from './PositionDispatch.js'
import ScrollBooster from 'scrollbooster';
//import { MSAViewer, SequenceViewer, Labels, } from '@plotly/react-msa-viewer';
import { MSAViewer, 
        SequenceViewer, 
        Labels, 
        PositionBar,
        OverviewBar } from './MSAV.umd.js';
import React, { Component } from "react";

var AlnViewer = class RV3AlnViewer extends Component {
    unSelectNucleotide = window.unSelectNucleotide;
    state = { 
        tileWidth: 17,
        tileHeight: 17,
        aaPos: vm.aaPos,
        seqPos: vm.seqPos,
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
        style.innerHTML = ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}"
    };
    componentDidMount() {
        var style = document.querySelector('[data="rv3_style"]');
        style.innerHTML = ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}";
        window.ajaxRun = false;
        var handleMolStarTopViewHovers = function (alnViewerClass, residueNumber, eventEntityID){
            var alignmentNumber = Number(_.invert(vm.structure_mapping)[residueNumber]);
            var numVisibleTiles = Math.round(alnViewerClass.state.width/alnViewerClass.state.tileWidth);
            if (!isNaN(alignmentNumber)&&eventEntityID==vm.entityID){
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
        document.addEventListener('PDB.RNA.viewer.mouseover', (e) => {
            handleMolStarTopViewHovers(this, e.eventData.label_seq_id, e.eventData.entityId);
        });
        document.addEventListener('PDB.RNA.viewer.mouseout', () => {
            this.removeHighlightRegion();
        });
        document.addEventListener('PDB.molstar.mouseover', (e) => {
            handleMolStarTopViewHovers(this, e.eventData.seq_id, e.eventData.entity_id);
        });
        document.addEventListener('PDB.molstar.mouseout', () => {
            this.removeHighlightRegion();
        });
        $('#alnSequenceViewer').mouseleave(function () {
            window.ajaxRun = false;
        });
        new ScrollBooster({
            viewport: document.querySelector("#alnViewerLabels").firstElementChild,
            scrollMode: 'native',
            direction: 'horizontal',
            bounce: false,
        });
        vm.msavWillMount = true;
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
                viewerInstanceTop.viewInstance.selectResidue(resiPos);
                viewerInstance.visual.highlighting({
                    data:[{
                            entity_id:`${vm.entityID}`,
                            //start_residue_number:resiPos,
                            //end_residue_number:resiPos,
                            residue_number:resiPos,
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
            //window.clearHighlight());
            let resiPos = vm.structure_mapping[e.position+1];
            if (resiPos !== undefined){
                viewerInstanceTop.viewInstance.clearSelection(resiPos);
            }
            viewerInstance.visual.clearHighlighting();
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
        var maxXpos = window.aaFreqs.length - Math.round((((window.innerWidth - 300) * 0.7)/this.state.tileWidth))+2;
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
                  {...window.msaOptions}
                  id = "MSAViewer"
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
                    <div>
                        <Labels 
                          id="alnViewerLabels"
                          style = {{
                            width: (window.innerWidth - 300) * 0.2,
                            paddingTop: 13.6,
                            marginRight: 3,
                            }}
                        />
                        
                    </div>
                    <div>
                        <PositionBar 
                          markerSteps={5} 
                          startIndex={1} 
                        />
                        <SequenceViewer
                          id="alnSequenceViewer"
                          onResidueMouseEnter={this.onResidueMouseEnter}
                          onResidueMouseLeave={this.onResidueMouseLeave}
                        />
                        <OverviewBar id="conservationBar" method='proteovision'/>
                        {this.state.fold && (
                          <div
                            style={{
                              position: "absolute",
                              opacity: 0.8,
                              ...this.state.tooltipPosition,
                            }}
                          >
                            {<Tooltip>
                              Fold: {this.state.fold} <br></br>
                              Phase: {this.state.phase}
                        </Tooltip>}
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