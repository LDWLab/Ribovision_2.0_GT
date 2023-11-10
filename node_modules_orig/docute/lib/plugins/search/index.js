import SearchBar from './SearchBar.vue';
export default {
  name: 'search',
  extend: function extend(api) {
    api.registerComponent('header-right:start', SearchBar);
  }
};