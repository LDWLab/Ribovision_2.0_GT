<template>  
  <div class="VersionsSelector">
    <DocuteSelect
      @change="handleChange"
      :value="currentVersionLink"
      v-slot="{value}"
    >
      <option disabled>Choose Version</option>
      <option
        v-for="(version, key) in $store.getters.config.versions"
        :key="key"
        :value="version.link"
        :selected="value === version.link"
      >
        {{ key }}
      </option>
    </DocuteSelect>
  </div>
</template>

<script>
import { isExternalLink } from '../../utils';
export default {
  methods: {
    handleChange: function handleChange(link) {
      if (isExternalLink(link)) {
        location.href = link;
      } else {
        this.$router.push(link);
      }
    }
  },
  computed: {
    currentVersionLink: function currentVersionLink() {
      var versions = this.$store.getters.config.versions;

      for (var _i = 0, _Object$keys = Object.keys(versions); _i < _Object$keys.length; _i++) {
        var version = _Object$keys[_i];
        var link = versions[version].link;

        if (link !== '/' && this.$route.path.startsWith(link)) {
          return link;
        }
      }

      return '/';
    }
  }
};
</script>

<style scoped>
.VersionsSelector {
  padding: 0 20px;
  margin-top: 10px;
}
</style>