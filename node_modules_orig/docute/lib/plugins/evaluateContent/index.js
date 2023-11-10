export default {
  name: 'hoistTags',
  extend: function extend(api) {
    api.extendMarkedRenderer(function (renderer) {
      var hoistedTagsRe = /^<(script|style)(?=(\s|>|$))/i;

      renderer.html = function (html) {
        html = html.trim();

        if (hoistedTagsRe.test(html)) {
          return html.replace(/^<(script|style)/, '<v-$1').replace(/<\/(script|style)>$/, '</v-$1>');
        }

        return html;
      };
    });
  }
};