import React from "react";
export function YSlider (props){
  const {alnViewerAdjHeight, maxYpos, MSAVObject} = props;
  if (alnViewerAdjHeight < MSAVObject.state.tileHeight*vm.fastaSeqNames.length) {
    return(
      <div style={{width: "30px"}}>
      <input
        style={{ 
            width: alnViewerAdjHeight+"px",
        }}
        type="range"
        min="1"
        max={maxYpos}
        value={MSAVObject.state.seqPos}
        onChange={(evt) => MSAVObject.setState({ seqPos: evt.target.value })}
        className="slider"
        id="yPosSlider"
      />
      </div>
    );
  }else{
    return (null);
  }
}

export function XSlider (props){
  const {alnViewerAdjWidth, maxXpos, MSAVObject} = props;
  if (alnViewerAdjWidth < MSAVObject.state.tileWidth*window.aaFreqs.length) {
    return(
      <div>
        <input
          style = {{ 
              width: (window.innerWidth - 300) * 0.7+"px",
              position: "relative",
              left: ((window.innerWidth - 300) * 0.2)+3+"px"
              }}
          type="range"
          min="1"
          max={maxXpos}
          value={MSAVObject.state.aaPos}
          onChange={(evt) => MSAVObject.setState({ aaPos: evt.target.value })}
          className="slider"
          id="xPosSlider"
        />
      </div>
    );
  }else{
    return (null);
  }
}