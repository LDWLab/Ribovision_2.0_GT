var getComponent = function getComponent(str, className) {
  return typeof str === 'string' ? {
    template: "<div class=\"" + className + "\">" + str + "</div>"
  } : str;
};

export default {
  name: 'banner-footer',
  extend: function extend(api) {
    var _api$store$getters$co = api.store.getters.config,
        banner = _api$store$getters$co.banner,
        footer = _api$store$getters$co.footer;

    if (banner) {
      api.registerComponent('content:start', getComponent(banner, 'docute-banner'));
    }

    if (footer) {
      api.registerComponent('content:end', getComponent(footer, 'docute-footer'));
    }
  }
};