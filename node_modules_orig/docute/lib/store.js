/* globals "1.17.1" */
import Vue from 'vue';
import Vuex from 'vuex';
import marked from './utils/marked';
import highlight from './utils/highlight';
import { getFilenameByPath, getFileUrl, isExternalLink, inBrowser } from './utils';
import markedRenderer from './utils/markedRenderer';
import hooks from './hooks';
import load from './utils/load';
import prismLanguages from './utils/prismLanguages';
import { defaultCssVariables, darkCssVariables } from './utils/cssVariables';
import { INITIAL_STATE_NAME } from './utils/constants';

function _await(value, then, direct) {
  if (direct) {
    return then ? then(value) : value;
  }

  if (!value || !value.then) {
    value = Promise.resolve(value);
  }

  return then ? value.then(then) : value;
}

Vue.use(Vuex);
var initialState = inBrowser && window[INITIAL_STATE_NAME];

var getDefaultTheme = function getDefaultTheme(store, _ref) {
  var theme = _ref.theme,
      detectSystemDarkTheme = _ref.detectSystemDarkTheme;

  if (!inBrowser || !detectSystemDarkTheme) {
    return theme || 'default';
  }

  var mq = window.matchMedia('(prefers-color-scheme: dark)');
  mq.addListener(function () {
    store.commit('SET_THEME', mq.matches ? 'dark' : 'default');
  });
  return theme || (mq.matches ? 'dark' : 'default');
};

var store = new Vuex.Store({
  state: Object.assign({
    originalConfig: {},
    page: {
      title: null,
      headings: null,
      html: ''
    },
    env: {},
    showSidebar: false,
    fetchingFile: true
  }, initialState),
  mutations: {
    SET_CONFIG: function SET_CONFIG(state, config) {
      if (config === void 0) {
        config = {};
      }

      config.layout = config.layout || 'narrow'; // TODO: remove `centerContent` in next major version

      if (config.centerContent) {
        config.layout = 'narrow';
      }

      config.theme = getDefaultTheme(store, config);
      state.originalConfig = config;
    },
    SET_PAGE: function SET_PAGE(state, page) {
      state.page = page;
    },
    TOGGLE_SIDEBAR: function TOGGLE_SIDEBAR(state, show) {
      state.showSidebar = typeof show === 'boolean' ? show : !state.showSidebar;
    },
    SET_FETCHING: function SET_FETCHING(state, fetching) {
      state.fetchingFile = fetching;
    },
    SET_ENV: function SET_ENV(state, env) {
      state.env = env;
    },
    SET_THEME: function SET_THEME(state, theme) {
      state.originalConfig.theme = theme;
    }
  },
  actions: {
    fetchFile: function fetchFile(_ref2, path) {
      var commit = _ref2.commit,
          getters = _ref2.getters,
          dispatch = _ref2.dispatch;

      try {
        commit('TOGGLE_SIDEBAR', false);
        commit('SET_FETCHING', true);
        var page = Object.assign({
          markdown: true
        }, getters.config.routes && getters.config.routes[path]);

        if (!page.content && !page.file) {
          var filename = getFilenameByPath(path);
          page.file = getFileUrl(getters.config.sourcePath, filename);
          page.editLink = getters.config.editLinkBase && getFileUrl(getters.config.editLinkBase, filename);
        }

        return _await(Promise.all([!page.content && fetch(page.file, getters.config.fetchOptions).then(function (res) {
          return res.text();
        }).then(function (res) {
          page.content = res;
        }), dispatch('fetchPrismLanguages')]), function () {
          return _await(hooks.processPromise('processMarkdown', page.content), function (_hooks$processPromise) {
            page.content = _hooks$processPromise;
            return _await(hooks.processPromise('processPage', page), function (_hooks$processPromise2) {
              page = _hooks$processPromise2;
              var env = {
                headings: [],
                mixins: [],
                config: getters.config
              };

              if (page.markdown) {
                page.content = marked(page.content, {
                  renderer: markedRenderer(hooks),
                  highlight: highlight,
                  env: env
                });
              }

              return _await(hooks.processPromise('processHTML', page.content), function (_hooks$processPromise3) {
                page.content = _hooks$processPromise3;
                page.headings = env.headings;

                if (!page.title) {
                  page.title = env.title;
                }

                commit('SET_PAGE', page);
                commit('SET_ENV', env);
                commit('SET_FETCHING', false);
              });
            });
          });
        });
      } catch (e) {
        return Promise.reject(e);
      }
    },
    fetchPrismLanguages: function fetchPrismLanguages(_ref3) {
      var getters = _ref3.getters;
      var langs = getters.config.highlight;

      if (!langs || langs.length === 0) {
        return Promise.resolve();
      }

      return load(langs.reduce(function (res, lang) {
        if (prismLanguages[lang]) {
          res = res.concat(prismLanguages[lang]);
        }

        res.push(lang);
        return res;
      }, []).filter(function (lang, i, arr) {
        // Dedupe
        return arr.indexOf(lang) === i && prismLanguages.builtin.indexOf(lang) === -1;
      }).map(function (lang) {
        return "https://unpkg.com/prismjs@" + "1.17.1" + "/components/prism-" + lang + ".js";
      }), 'prism-languages');
    }
  },
  getters: {
    target: function target(_ref4) {
      var _target = _ref4.originalConfig.target;
      if (!_target) return 'docute';
      if (_target[0] === '#') return _target.slice(1);
      return _target;
    },
    languageOverrides: function languageOverrides(_ref5) {
      var originalConfig = _ref5.originalConfig;
      // `locales` is for legacy support
      var overrides = originalConfig.overrides || originalConfig.locales;
      return overrides && Object.keys(overrides).reduce(function (res, path) {
        if (overrides[path].language) {
          res[path] = overrides[path];
        }

        return res;
      }, {});
    },
    currentLocalePath: function currentLocalePath(_ref6, _ref7) {
      var route = _ref6.route;
      var languageOverrides = _ref7.languageOverrides;

      if (languageOverrides) {
        // Is it a locale?
        for (var _i = 0, _Object$keys = Object.keys(languageOverrides); _i < _Object$keys.length; _i++) {
          var localePath = _Object$keys[_i];

          if (localePath !== '/') {
            var RE = new RegExp("^" + localePath);

            if (RE.test(route.path)) {
              return localePath;
            }
          }
        }
      }

      return '/';
    },
    config: function config(_ref8, _ref9) {
      var originalConfig = _ref8.originalConfig;
      var currentLocalePath = _ref9.currentLocalePath,
          languageOverrides = _ref9.languageOverrides;
      return languageOverrides ? Object.assign({}, originalConfig, {}, languageOverrides[currentLocalePath]) : originalConfig;
    },
    homePaths: function homePaths(_, _ref10) {
      var languageOverrides = _ref10.languageOverrides;
      var localePaths = languageOverrides ? Object.keys(languageOverrides) : [];
      return [].concat(localePaths, ['/']);
    },
    sidebarLinks: function sidebarLinks(_, _ref11) {
      var sidebar = _ref11.sidebar;
      return sidebar ? sidebar.reduce(function (res, next) {
        // backward compabillity
        var children = next.children || next.links || [];
        return [].concat(res, children);
      }, []).filter(function (item) {
        return !isExternalLink(item.link);
      }) : [];
    },
    sidebar: function sidebar(_, _ref12) {
      var config = _ref12.config;
      var sidebar = config.sidebar || [];
      return typeof sidebar === 'function' ? sidebar(store) : sidebar;
    },
    cssVariables: function cssVariables(_, _ref13) {
      var config = _ref13.config;
      return Object.assign({}, config.theme === 'dark' ? darkCssVariables : defaultCssVariables, {}, typeof config.cssVariables === 'function' ? config.cssVariables(config.theme) : config.cssVariables);
    }
  }
});

if (process.env.NODE_ENV === 'development' && inBrowser) {
  window.store = store;
}

export default store;