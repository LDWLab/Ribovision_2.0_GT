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