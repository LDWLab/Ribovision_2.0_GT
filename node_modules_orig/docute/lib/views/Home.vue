<template>  
  <div class="Page" :class="{[`layout-${$store.getters.config.layout}`]: true}">
    <SiteHeader />
    <div class="Wrap">
      <Sidebar />
      <SidebarMask />
      <div class="Main">
        <div class="Content" v-if="$store.state.fetchingFile">
          <content-loader
            :height="160"
            :width="400"
            :speed="2"
            :primaryColor="$store.getters.cssVariables.loaderPrimaryColor"
            :secondaryColor="$store.getters.cssVariables.loaderSecondaryColor"
          >
            <rect x="0" y="5" rx="4" ry="4" width="117" height="6.4" />
            <rect x="0" y="25" rx="3" ry="3" width="85" height="6.4" />
            <rect x="0" y="60" rx="3" ry="3" width="350" height="6.4" />
            <rect x="0" y="80" rx="3" ry="3" width="380" height="6.4" />
            <rect x="0" y="100" rx="3" ry="3" width="201" height="6.4" />
          </content-loader>
        </div>
        <div class="Content" v-else>
          <InjectedComponents position="content:start" />
          <component v-if="pageTitle" :is="MarkdownTitle" class="page-title" />
          <component :class="{'has-page-title': pageTitle}" :is="PageContent" />
          <EditLink />
          <PrevNextLinks />
          <InjectedComponents position="content:end" />
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import Vue from 'vue';
import jump from 'jump.js';
import { ContentLoader } from 'vue-content-loader';
import Sidebar from '../components/Sidebar.vue';
import SidebarMask from '../components/SidebarMask.vue';
import SiteHeader from '../components/Header.vue';
import PrevNextLinks from '../components/PrevNextLinks.vue';
import EditLink from '../components/EditLink.vue';
import { INITIAL_STATE_NAME } from '../utils/constants';
import hooks from '../hooks';

function _await(value, then, direct) {
  if (direct) {
    return then ? then(value) : value;
  }

  if (!value || !value.then) {
    value = Promise.resolve(value);
  }

  return then ? value.then(then) : value;
}

export default {
  name: 'PageHome',
  components: {
    ContentLoader: ContentLoader,
    Sidebar: Sidebar,
    SidebarMask: SidebarMask,
    SiteHeader: SiteHeader,
    PrevNextLinks: PrevNextLinks,
    EditLink: EditLink
  },
  serverPrefetch: function serverPrefetch() {
    try {
      var _this2 = this;

      return _await(_this2.fetchFile(_this2.$route.path), function () {
        _this2.setTitle();
      });
    } catch (e) {
      return Promise.reject(e);
    }
  },
  mounted: function mounted() {
    if (!window[INITIAL_STATE_NAME]) {
      this.fetchFile(this.$route.path).then(this.setInitialState);
    }
  },
  beforeRouteUpdate: function beforeRouteUpdate(to, from, next) {
    next();

    if (to.path !== from.path) {
      this.fetchFile(to.path);
    }
  },
  watch: {
    '$route.hash': function $routeHash() {
      var _this3 = this;

      this.$nextTick(function () {
        _this3.jumpToHash();
      });
    },
    pageTitle: function pageTitle() {
      this.setTitle();
    }
  },
  computed: {
    pageTitle: function pageTitle() {
      return this.$store.state.page.title;
    },
    MarkdownTitle: function MarkdownTitle() {
      return {
        name: 'MarkdownTitle',
        template: "<h1>" + this.pageTitle + "</h1>"
      };
    },
    PageContent: function PageContent() {
      var env = this.$store.state.env;
      var _this$$store$getters$ = this.$store.getters.config.componentMixins,
          componentMixins = _this$$store$getters$ === void 0 ? [] : _this$$store$getters$;
      var component = {
        mixins: [].concat(componentMixins, env.mixins.map(function (mixin) {
          // eslint-disable-next-line no-new-func
          var fn = new Function('Vue', "return " + mixin.trim());
          return fn(Vue);
        })),
        name: 'PageContent',
        template: "<div class=\"page-content\">" + this.$store.state.page.content + "</div>"
      };
      hooks.process('extendMarkdownComponent', component);
      return component;
    }
  },
  methods: {
    fetchFile: function fetchFile(path) {
      try {
        var _this5 = this;

        return _await(_this5.$store.dispatch('fetchFile', path), function () {
          hooks.invoke('onContentWillUpdate', _this5);
          return _await(_this5.$nextTick(), function () {
            hooks.invoke('onContentUpdated', _this5);

            _this5.jumpToHash();
          });
        });
      } catch (e) {
        return Promise.reject(e);
      }
    },
    jumpToHash: function jumpToHash() {
      var hash = decodeURI(this.$route.hash);

      if (hash) {
        var el = document.querySelector(hash);

        if (el) {
          var header = document.querySelector('.Header');
          jump(el, {
            a11y: true,
            duration: 0,
            offset: -(header.clientHeight + 30)
          });
        }
      }
    },
    setInitialState: function setInitialState() {
      if (/(Prerender|jsdom|PhantomJS)/i.test(navigator.userAgent)) {
        var script = document.createElement('script');
        script.textContent = "window." + INITIAL_STATE_NAME + " = " + JSON.stringify({
          page: this.$store.state.page,
          env: this.$store.state.env,
          fetchingFile: false
        });
        document.head.appendChild(script);
      }
    },
    setTitle: function setTitle() {
      var path = this.$route.path;
      var _this$$store$getters = this.$store.getters,
          config = _this$$store$getters.config,
          homePaths = _this$$store$getters.homePaths;
      var title = homePaths.indexOf(path) > -1 ? config.title : this.pageTitle + " - " + config.title; // Strip HTML tags

      title = title.replace(/<(?:.|\n)*?>/gm, '');

      if (this.$ssrContext) {
        this.$ssrContext.title = title;
      } else {
        document.title = title;
      }
    }
  }
};
</script>

<style src="../css/prism.css"></style>

<style src="../css/page-content.css"></style>

<style scoped>
.Main {
  padding-left: var(--sidebar-width);
  padding-top: calc(var(--header-height) + 40px);
  padding-bottom: 2rem;
  background: var(--main-background)
}

@media screen and (max-width: 768px) {

.Main {
    padding-left: 0
}
  }

@media print {

.Main {
    padding-left: 0;
    padding-top: 30px
}
  }

.Content {
  padding: 0 20px 0 80px
}

@media screen and (max-width: 768px) {

.Content {
    padding: 0 20px
}
  }

.layout-wide .Content {
  max-width: 770px;
  margin: 0 auto;
  padding: 0 2.5rem
}

@media screen and (max-width: 768px) {

.layout-wide .Content {
    max-width: 100%;
    padding: 0 20px
}
  }

@media print {

.layout-wide .Content {
    padding: 0
}
  }

.page-title {
  font-size: 3rem;
  margin: 0;
  margin-bottom: 1.4rem;
  font-weight: 300;
  line-height: 1.1;
}
</style>