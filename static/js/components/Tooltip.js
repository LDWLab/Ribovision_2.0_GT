import React from "react";
export function Tooltip (props) {
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
    <div id="tooltip" style={{ ...containerStyle, ...style }} {...otherProps}>
      <div style={tooltipStyle}>
        <div style={textStyle}>{children}</div>
      </div>
    </div>
  );
}
