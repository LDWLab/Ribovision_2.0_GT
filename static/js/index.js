import { actions, MSAViewer, SequenceViewer } from '@plotly/react-msa-viewer';
import React, { Component } from "react";
import ReactDOM, { render } from 'react-dom';
//import PdbTopologyViewer from 'pdb-topology-viewer/src/web-component/index'
import PdbTopologyViewerPlugin from 'pdb-topology-viewer/build/pdb-topology-viewer-plugin'

function MyTopView() {
  return (
    class PdbTopologyViewer extends HTMLElement {

      static get observedAttributes() {
        return ['entry-id', 'entity-id', 'chain-id', 'display-style', 'error-style', 'menu-style', 'subscribe-events'];
      }

      constructor() {
        super();
      }

      validateParams() {
        if(typeof this.entryId == 'undefined' || typeof this.entityId == 'undefined' || this.entryId == null || this.entityId == null) return false;
        return true
      }

      invokePlugin(){
        let paramValidatity = this.validateParams();
        if(!paramValidatity) return

        // create an instance of the plugin
        if(typeof this.pluginInstance == 'undefined') this.pluginInstance = new PdbTopologyViewerPlugin();

        let options = {
          entryId: this.entryId,
          entityId: this.entityId,
        }

        if(typeof this.chainId !== 'undefined' && this.chainId !== null) options['chainId'] = this.chainId;
        if(typeof this.displayStyle !== 'undefined' && this.displayStyle !== null) options['displayStyle'] = this.displayStyle;
        if(typeof this.errorStyle !== 'undefined' && this.errorStyle !== null) options['errorStyle'] = this.errorStyle;
        if(typeof this.menuStyle !== 'undefined' && this.menuStyle !== null) options['menuStyle'] = this.menuStyle;
        if(typeof this.subscribeEvents !== 'undefined' && this.subscribeEvents !== null) options['subscribeEvents'] = this.subscribeEvents;

        this.pluginInstance.render(this, options);

      }

      attributeChangedCallback() {
        this.entryId = this.getAttribute("entry-id");
        this.entityId = this.getAttribute("entity-id");
        this.chainId = this.getAttribute("chain-id");
        this.displayStyle = this.getAttribute("display-style");
        this.errorStyle = this.getAttribute("error-style");
        this.menuStyle = this.getAttribute("menu-style");
        this.subscribeEvents = this.getAttribute("subscribe-events");
        this.invokePlugin();
      }

    }
  )
}


const options = {
  sequences: [
    {
      name: "seq.1",
      sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.7",
      sequence: "MEEPQSDPSIEP-PLSQ------WKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.2",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.3",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.4",
      sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.5",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.6",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    }
  ],
  height: 60,
};

function Tooltip(props) {
  const { direction, style, children, ...otherProps } = props;
  const containerStyle = {
    display: "inline-block"
  };
  const tooltipStyle = {
    position: "relative",
    width: "160px"
  };
  const textStyle = {
    color: "#fff",
    fontSize: "14px",
    lineHeight: 1.2,
    textAlign: "center",
    backgroundColor: "#000",
    borderRadius: "3px",
    padding: "7px"
  };
  const triangleStyle = {
    position: "absolute",
    width: 0,
    fontSize: 0,
    lineHeight: 0,
    visibility: "visible",
    opacity: 1
  };
  switch (direction) {
    case "up":
    case "down":
      triangleStyle.borderLeft = "5px solid transparent";
      triangleStyle.borderRight = "5px solid transparent";
      triangleStyle.left = "50%";
      triangleStyle.marginLeft = "-5px";
      break;
    case "left":
    case "right":
      triangleStyle.borderTop = "5px solid transparent";
      triangleStyle.borderBottom = "5px solid transparent";
      triangleStyle.top = "50%";
      triangleStyle.marginTop = "-5px";
      break;
    default:
  }
  switch (direction) {
    case "down":
      triangleStyle.borderTop = "5px solid #000";
      triangleStyle.top = "100%";
      containerStyle.paddingBottom = "5px";
      break;
    case "up":
      triangleStyle.borderBottom = "5px solid #000";
      triangleStyle.top = "0%";
      triangleStyle.marginTop = "-5px";
      containerStyle.paddingTop = "5px";
      break;
    case "left":
      triangleStyle.borderRight = "5px solid #000";
      triangleStyle.marginLeft = "-5px";
      containerStyle.paddingLeft = "5px";
      break;
    case "right":
      triangleStyle.left = "100%";
      triangleStyle.borderLeft = "5px solid #000";
      containerStyle.paddingRight = "5px";
      break;
    default:
  }
  return (
    <div style={{ ...containerStyle, ...style }} {...otherProps}>
      <div style={tooltipStyle}>
        <div style={textStyle}>{children}</div>
        <div style={triangleStyle} />
      </div>
    </div>
  );
}
Tooltip.defaultProps = {
  style: {},
  direction: "down"
};

function MyMSA() {
  class SimpleTooltip extends Component {
      state = {};
      onResidueMouseEnter = e => {
        let direction, tooltipPosition;
        if (e.position < 10) {
          direction = "left";
          tooltipPosition = {
            top: (e.i - 0.3) * 20 + "px",
            left: (e.position + 1) * 20 + "px"
          };
        } else {
          direction = "right";
          tooltipPosition = {
            top: (e.i - 0.3) * 20 + "px",
            left: e.position * 20 - 165 + "px"
          };
        }
        this.setState({
          lastEvent: e,
          tooltipPosition,
          direction
        });
      };
      onResidueMouseLeave = e => {
        this.setState({ lastEvent: undefined });
      };
      render() {
        return (
          <div>
            <MSAViewer {...options}>
              <div style={{ position: "relative" }}>
                <SequenceViewer
                  onResidueMouseEnter={this.onResidueMouseEnter}
                  onResidueMouseLeave={this.onResidueMouseLeave}
                />
                {this.state.lastEvent && (
                  <div
                    style={{
                      position: "absolute",
                      opacity: 0.8,
                      ...this.state.tooltipPosition
                    }}
                  >
                    <Tooltip direction={this.state.direction}>
                      {this.state.lastEvent.residue}
                    </Tooltip>
                  </div>
                )}
              </div>
            </MSAViewer>
          </div>
        );
      }
  }
  return <SimpleTooltip />;
};



ReactDOM.render(<MyMSA />, document.getElementById("target"));
ReactDOM.render(<MyTopView />, document.getElementById("target"));