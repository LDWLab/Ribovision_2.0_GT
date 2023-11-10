var getComponent = function getComponent(tag) {
  return {
    functional: true,
    render: function render(h, ctx) {
      return h(tag, ctx.data, ctx.children);
    }
  };
};

export default (function (Vue) {
  Vue.component('v-style', getComponent('style'));
  Vue.component('v-script', getComponent('script'));
});