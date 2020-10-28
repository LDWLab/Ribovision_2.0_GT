<template>
    <div>
      <p><div id="v-step-0"></div></p>
      <treeselect id="treeselect" v-model="value" :multiple="true" :options="options" v-on:input="loadAlignment(value)"/>
    </div>
</template>


<script>
import { actions, MSAViewer, SequenceViewer } from '@plotly/react-msa-viewer';
import React, { Component } from "react";
import ReactDOM, { render } from 'react-dom';
(function() {
  var mousePos;
  document.onmousemove = handleMouseMove;
  function handleMouseMove(event) {
      var eventDoc, doc, body;
      event = event || window.event; // IE-ism
      // If pageX/Y aren't available and clientX/Y are,
      // calculate pageX/Y - logic taken from jQuery.
      // (This is to support old IE)
      if (event.pageX == null && event.clientX != null) {
          eventDoc = (event.target && event.target.ownerDocument) || document;
          doc = eventDoc.documentElement;
          body = eventDoc.body;
          event.pageX = event.clientX +
            (doc && doc.scrollLeft || body && body.scrollLeft || 0) -
            (doc && doc.clientLeft || body && body.clientLeft || 0);
          event.pageY = event.clientY +
            (doc && doc.scrollTop  || body && body.scrollTop  || 0) -
            (doc && doc.clientTop  || body && body.clientTop  || 0 );
      }
      mousePos = {
        x: event.pageX,
        y: event.pageY
    };
    window.mousePos = mousePos;
  }
})();

const options = {
  sequences: [
    {
      name: "seq.1",
      sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.2",
      sequence: "MEEPQSDPSIEP-PLSQ------WKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.3",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.4",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.5",
      sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
    },
    {
      name: "seq.6",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    },
    {
      name: "seq.7",
      sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
    }
  ],
  height: 120,
};

function Tooltip(props) {
  const { style, children, ...otherProps } = props;
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
  return (
    <div style={{ ...containerStyle, ...style }} {...otherProps}>
      <div style={tooltipStyle}>
        <div style={textStyle}>{children}</div>
      </div>
    </div>
  );
}
Tooltip.defaultProps = {
  style: {},
};
function ajax(url) {
  return new Promise((resolve, reject) => {
      $.ajax({
          url: url,
          type: 'GET',
          dataType: "json",
          success: function(data) {
              resolve(data)
          },
          error: function(error) {
              console.log(`Error ${error}`);
              reject(error)
          }
      })
  })
}


window.ajaxRun = false;
function MyMSA() {
  class SimpleTooltip extends Component {
      state = {};
      onResidueMouseEnter = e => {
        if (!window.ajaxRun){
        window.ajaxRun = true;
        const strainQuery = '&res__poldata__strain__strain=';
        if (e.position !== undefined){
          var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(e.position) + 1)}&aln=10${strainQuery}Escherichia coli str. K-12 substr. MG1655`
            ajax(url).then(alnpos_data => {
                if (alnpos_data.count != 0){
                  ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
                    var alnViewEle = document.querySelector("canvas");
                    let boundingBox = alnViewEle.getBoundingClientRect();
                    if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
                      let tooltipPosition;
                      tooltipPosition = {
                        top: mousePos.y +"px",
                        left: mousePos.x +"px",
                      };
                      if (resiData["Structural fold"][0] !== undefined && resiData["Associated data"][0] !== undefined){
                        this.setState({
                          fold: resiData["Structural fold"][0][1],
                          phase: resiData["Associated data"][0][1],
                          tooltipPosition,
                        });
                      }else{
                        this.setState({
                          fold: 'NA',
                          phase: 'NA',
                          tooltipPosition,
                        });
                      }
                    }
                    window.ajaxRun = false;
                  });
                }else{
                  window.ajaxRun = false;
                }
            }).catch(error => {
                console.log(error);
            })
          }
        }else{
          this.setState({ fold: undefined, phase: undefined });
        }
      };
      onResidueMouseLeave = e => {
        this.setState({ fold: undefined, phase: undefined });
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
                {this.state.fold && (
                  <div
                    style={{
                      position: "absolute",
                      opacity: 0.8,
                      ...this.state.tooltipPosition
                    }}
                  >
                    <Tooltip>
                      Fold: {this.state.fold} <br></br>
                      Phase: {this.state.phase}
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
  import Treeselect from '@riophae/vue-treeselect'
  export default {
      // register the component
      components: { Treeselect },
      data() {
        return {
          // define the default value
          value: null,
          // define options
          options: [ {
            id: 'a',
            label: 'a',
            children: [ {
              id: 'aa',
              label: 'aa',
            }, {
              id: 'ab',
              label: 'ab',
            } ],
          }, {
            id: 'b',
            label: 'b',
          }, {
            id: 'c',
            label: 'c',
          } ],
        }
      },
      methods:{
        loadAlignment(value){
          if (value == 'a'){
          ReactDOM.render(<MyMSA />, document.getElementById("v-step-0"));
          }
        }
      }
    }
</script>