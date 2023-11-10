<script>
export default {
  name: 'DocuteRoot',
  render: function render(h) {
    return h('div', {
      attrs: {
        id: this.$store.getters.target,
        class: 'Root'
      }
    }, [h('router-view')]);
  },
  created: function created() {
    this.insertStyle();
  },
  computed: {
    css: function css() {
      var cssVariables = this.$store.getters.cssVariables;
      var content = Object.keys(cssVariables).reduce(function (res, key) {
        res += "--" + key.replace(/[A-Z]/g, function (m) {
          return "-" + m.toLowerCase();
        }) + ":" + cssVariables[key] + ";";
        return res;
      }, '');
      return ":root{" + content + "}";
    }
  },
  watch: {
    css: function css() {
      this.insertStyle();
    }
  },
  methods: {
    insertStyle: function insertStyle() {
      if (this.$ssrContext) {
        this.$ssrContext.insertedStyle = this.css;
        return;
      }

      var ID = 'docute-inserted-style';
      var style = document.getElementById(ID);

      if (style) {
        style.innerHTML = this.css;
      } else {
        style = document.createElement('style');
        style.id = ID;
        style.innerHTML = this.css;
        document.head.insertBefore(style, document.head.firstChild);
      }
    }
  }
};
</script>

<style src="../css/global.css"></style>