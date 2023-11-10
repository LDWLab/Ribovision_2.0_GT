<template>  
  <div class="prev-next-links" v-if="prevLinkItem || nextLinkItem">
    <div class="prev-link" v-if="prevLinkItem">
      <router-link :to="prevLinkItem.link">
        ← {{ prevLinkItem.title }}
      </router-link>
    </div>
  
    <div class="next-link" v-if="nextLinkItem">
      <router-link :to="nextLinkItem.link">
        {{ nextLinkItem.title }} →
      </router-link>
    </div>
  </div>
</template>

<script>
import { mapGetters } from 'vuex';
export default {
  computed: Object.assign({}, mapGetters(['sidebarLinks']), {
    currentLink: function currentLink() {
      return this.$route.path;
    },
    currentLinkIndex: function currentLinkIndex() {
      // Related:
      // - https://github.com/vuejs/vue/issues/8728
      // - https://github.com/egoist/docute/pull/171
      var sidebarLinks = this.sidebarLinks;

      for (var i = 0; i < sidebarLinks.length; i++) {
        var item = sidebarLinks[i];

        if (item.link === this.currentLink) {
          return i;
        }
      }

      return false;
    },
    prevLinkItem: function prevLinkItem() {
      return typeof this.currentLinkIndex === 'number' ? this.sidebarLinks[this.currentLinkIndex - 1] : null;
    },
    nextLinkItem: function nextLinkItem() {
      return typeof this.currentLinkIndex === 'number' ? this.sidebarLinks[this.currentLinkIndex + 1] : null;
    }
  })
};
</script>

<style scoped>
.prev-next-links {
  overflow: auto;
  margin-top: 40px;
  padding-top: 30px;
  border-top: 1px solid var(--border-color)
}

@media print {

.prev-next-links {
    display: none
}
  }

.prev-link {
  float: left;
}

.next-link {
  float: right;
}
</style>