<script>
import { isExternalLink } from '../utils';
export default {
  functional: true,
  props: ['openInNewTab', 'externalLinkIcon'],
  render: function render(h, _ref) {
    var data = _ref.data,
        children = _ref.children,
        _ref$props = _ref.props,
        openInNewTab = _ref$props.openInNewTab,
        externalLinkIcon = _ref$props.externalLinkIcon;
    var attrs = Object.assign({}, data.attrs);
    var to = attrs.to;

    if (isExternalLink(to)) {
      delete attrs.to;
      delete attrs.prefetchFiles;
      return h('a', Object.assign({}, data, {
        class: [data.class, 'is-external-link'],
        attrs: Object.assign({}, attrs, {
          href: to,
          target: openInNewTab === false ? '_self' : '_blank'
        })
      }), [].concat(children, [openInNewTab === false || externalLinkIcon === false ? null : h('external-link-icon', {
        class: 'external-link-icon'
      })]));
    }

    return h('router-link', data, children);
  }
};
</script>

<style>
.external-link-icon {
  margin-left: 3px;
}
</style>