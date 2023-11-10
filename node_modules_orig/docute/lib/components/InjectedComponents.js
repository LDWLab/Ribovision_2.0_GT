export default {
  name: 'InjectedComponents',
  functional: true,
  props: {
    position: {
      type: String,
      required: true
    }
  },
  render: function render(h, _ref) {
    var props = _ref.props,
        parent = _ref.parent;
    var components = parent.$pluginApi.getComponents(props.position);
    if (components.length === 0) return;
    return h('div', {
      class: 'InjectedComponents',
      attrs: {
        'data-position': props.position
      }
    }, components.map(function (_ref2) {
      var component = _ref2.component,
          props = _ref2.props;
      return h(component, {
        props: props
      });
    }));
  }
};