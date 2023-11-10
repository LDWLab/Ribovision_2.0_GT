var RE = /\.md$/;
export default (function (input) {
  var _input$split = input.split('#'),
      path = _input$split[0],
      hash = _input$split[1];

  if (RE.test(path)) {
    path = path.replace(RE, '');
  }

  return "" + path + (hash ? "#" + hash : '');
});