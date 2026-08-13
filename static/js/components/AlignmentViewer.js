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
        maskFeatures: [],
        hideSequences: null,
        hideColumnMap: null,
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
        this._isMounted = true;
        var style = document.querySelector('[data="rv3_style"]');
        style.innerHTML = ".slider::-webkit-slider-thumb { width: "+(window.innerWidth - 300)*0.05+"px}";
        window.ajaxRun = false;
        var handleMolStarTopViewHovers = function (alnViewerClass, residueNumber, eventEntityID){
            var alignmentNumber = Number(_.invert(vm.structure_mapping)[residueNumber]);
            var renderedPosition = alnViewerClass.resolveRenderedPosition(alignmentNumber - 1);
            var renderedAlignmentNumber = renderedPosition === undefined ? undefined : renderedPosition + 1;
            var numVisibleTiles = Math.round(alnViewerClass.state.width/alnViewerClass.state.tileWidth);
            if (renderedAlignmentNumber !== undefined && !isNaN(renderedAlignmentNumber) && eventEntityID==vm.entityID){
                if (alnViewerClass.state.aaPos > renderedAlignmentNumber || renderedAlignmentNumber > alnViewerClass.state.aaPos+numVisibleTiles){
                    let visiblePos = renderedAlignmentNumber-Math.round(numVisibleTiles/2);
                    if (visiblePos < 0) {visiblePos = 0};
                    alnViewerClass.setState({ aaPos: visiblePos })
                }
                alnViewerClass.highlightRegion({
                    sequences: {from: 0, to: vm.fastaSeqNames.length},
                    residues: {from: renderedAlignmentNumber, to: renderedAlignmentNumber}
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
        this._isMounted = false;
        window.removeEventListener("resize", this.handleResize);
    };
    // When hideSequences/hideColumnMap are set (masking Focus/Hide mode), the
    // rendered sequence has had columns removed, so a rendered column index
    // no longer equals the real alignment position. Translate it back here
    // so tooltips/3D highlighting stay correctly aligned.
    resolveAlnPosition = (renderedPosition) => {
        if (this.state.hideColumnMap) {
            return this.state.hideColumnMap[renderedPosition];
        }
        return renderedPosition;
    };
    resolveRenderedPosition = (alignmentPosition) => {
        if (this.state.hideColumnMap) {
            var renderedPosition = this.state.hideColumnMap.indexOf(alignmentPosition);
            return renderedPosition === -1 ? undefined : renderedPosition;
        }
        return alignmentPosition;
    };
    renderPositionMarker = ({index}) => {
        var originalPosition = this.state.hideColumnMap[index];
        if (originalPosition === undefined) {
            return <div style={{width: this.state.tileWidth, display: 'inline-block', textAlign: 'center'}}>.</div>;
        }
        var alignmentNumber = originalPosition + 1;
        var previousAlignmentNumber = index > 0 ? this.state.hideColumnMap[index - 1] + 1 : undefined;
        var startsFragment = previousAlignmentNumber === undefined || alignmentNumber !== previousAlignmentNumber + 1;
        var label = startsFragment || alignmentNumber % 5 === 0 ? alignmentNumber : '.';
        return <div style={{width: this.state.tileWidth, display: 'inline-block', textAlign: 'center'}}>{label}</div>;
    };
    onResidueMouseEnter = e => {
        this.highlightRegion({
            sequences: {from: 0, to: vm.fastaSeqNames.length},
            residues: {from: e.position+1, to: e.position+1}
        })
        if (vm.topology_loaded && !window.mask3DUpdateInProgress && window.viewerInstance && viewerInstance.plugin){
            let alnPosition = this.resolveAlnPosition(e.position);
            if (alnPosition === undefined) { return; }
            let resiPos = vm.structure_mapping[alnPosition+1];
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
    };
    onResidueMouseLeave = e => {
        if (vm.topology_loaded && !window.mask3DUpdateInProgress && window.viewerInstance && viewerInstance.plugin){
            //window.clearHighlight());
            let alnPosition = this.resolveAlnPosition(e.position);
            let resiPos = alnPosition !== undefined ? vm.structure_mapping[alnPosition+1] : undefined;
            if (resiPos !== undefined){
                viewerInstanceTop.viewInstance.clearSelection(resiPos);
            }
            viewerInstance.visual.clearHighlighting();
        }
        if (!this._isMounted) { return; }
        this.setState({ fold: undefined, phase: undefined });
    };
    highlightRegion = (highlight) => {
        if (!this._isMounted) { return; }
        this.setState({ highlight });
    };
    removeHighlightRegion = () => {
        if (!this._isMounted) { return; }
        this.setState({ highlight: null });
    };
    render() {
        const xPos = this.state.tileWidth * (this.state.aaPos);
        const yPos = this.state.tileHeight * (this.state.seqPos);
        var renderedSequences = this.state.hideSequences || window.msaOptions.sequences;
        var renderedLength = renderedSequences.length ? renderedSequences[0].sequence.length : 0;
        var maxXpos = renderedLength - Math.round((((window.innerWidth - 300) * 0.7)/this.state.tileWidth))+2;
        if (vm.fastaSeqNames) {
            var maxYpos = vm.fastaSeqNames.length - Math.round(((((window.innerHeight - 171)/2) * 0.8)/this.state.tileHeight));
        } else {
            maxYpos = 0
        }
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
                  sequences={renderedSequences}
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
                          key={this.state.hideColumnMap ? this.state.hideColumnMap.join(',') : 'full'}
                          markerSteps={5}
                          startIndex={1}
                          markerComponent={this.state.hideColumnMap ? this.renderPositionMarker : undefined}
                        />
                        <SequenceViewer
                          id="alnSequenceViewer"
                          onResidueMouseEnter={this.onResidueMouseEnter}
                          onResidueMouseLeave={this.onResidueMouseLeave}
                          features={this.state.maskFeatures}
                        />
                        <OverviewBar
                          key={this.state.hideColumnMap ? this.state.hideColumnMap.join(',') + '-' + vm.selected_property : 'full-' + vm.selected_property}
                          id="conservationBar"
                          method='proteovision'
                          columnMap={this.state.hideColumnMap}
                        />
                        {this.state.fold && (
                          <div
                            style={{
                              position: "absolute",
                              opacity: 0.8,
                              ...this.state.tooltipPosition,
                            }}
                          >
                            {/*<Tooltip>
                              Fold: {this.state.fold} <br></br>
                              Phase: {this.state.phase}
                        </Tooltip>*/}
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