<template>  
  <header class="Header">
    <div class="Wrap">
      <div class="header-inner">
        <div class="header-left">
          <h1 class="site-title">
            <SidebarToggle />
            <router-link to="/"> <component :is="Logo" /> </router-link>
          </h1>
          <HeaderNav v-if="leftNav" :nav="leftNav" />
        </div>
        <div class="header-right">
          <InjectedComponents position="header-right:start" />
          <HeaderNav v-if="rightNav" :nav="rightNav" />
          <InjectedComponents position="header-right:end" />
        </div>
      </div>
    </div>
  </header>
</template>

<script>
import HeaderNav from './HeaderNav.vue';
import SidebarToggle from './SidebarToggle.vue';
export default {
  components: {
    HeaderNav: HeaderNav,
    SidebarToggle: SidebarToggle
  },
  computed: {
    leftNav: function leftNav() {
      var nav = this.$store.getters.config.nav;
      return nav && nav.filter(function (item) {
        return item.position === 'left';
      });
    },
    rightNav: function rightNav() {
      var nav = this.$store.getters.config.nav;
      return nav && nav.filter(function (item) {
        return item.position === 'right' || !item.position;
      });
    },
    Logo: function Logo() {
      var _this$$store$getters$ = this.$store.getters.config,
          title = _this$$store$getters$.title,
          logo = _this$$store$getters$.logo;

      if (typeof logo === 'string') {
        return {
          template: logo
        };
      }

      if (typeof logo === 'object') {
        return logo;
      }

      return {
        render: function render(h) {
          return h('span', null, [title]);
        }
      };
    }
  }
};
</script>

<style scoped>
.Header {
  height: var(--header-height);
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  z-index: 33;
  border-bottom: 1px solid var(--border-color);
  background: var(--header-background);
  color: var(--header-text-color)
}

@media print {

.Header {
    position: static
}
  }

.Wrap {
  height: 100%;
}

.site-title {
  font-weight: normal;
  margin: 0 25px 0 0;
  font-size: 1.2rem;
  display: flex;
  align-items: center;
  white-space: nowrap
}

.site-title a {
    color: inherit;
    text-decoration: none;
  }

.header-inner {
  height: 100%;
  padding: 0 20px;
  position: relative;
  display: flex;
  align-items: center
}

@media print {

.header-inner {
    padding: 0
}
  }

.header-left {
  display: flex;
}

.header-right {
  display: flex;
  position: absolute;
  right: 20px;
  top: 0;
  height: 100%;
  background: var(--header-background)
}

@media print {

.header-right {
    display: flex;
    right: 0;
    padding-right: 0
}

    .header-right /deep/ .header-nav {
      display: flex;
    }
  }
</style>