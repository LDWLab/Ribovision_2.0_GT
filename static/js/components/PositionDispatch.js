import { msaConnect, withPositionStore, } from './MSAV.umd.js';
import { Component } from "react";

class MSAPluginComponent extends Component {
  // called on every position update (e.g. mouse movement or scrolling)
  shouldRerender(newPosition) {
    this.props.parent_state.seqPos = newPosition.currentViewSequence;
    this.props.parent_state.aaPos = newPosition.currentViewSequencePosition;
    return true;
  }
  render() {
    return null;
  }
}

const MSAPluginConnected = withPositionStore(
  MSAPluginComponent
);

// select attributes from the main redux store
var mapStateToProps = (state) => {
  return {
    height: state.props.height,
    sequences: state.sequences,
  };
};

// subscribe to the main redux store
const XYDispatch = msaConnect(mapStateToProps)(
  MSAPluginConnected
);

export {XYDispatch}