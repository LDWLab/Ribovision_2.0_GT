<template>
  <div id="app">
    <h1>Test</h1>
    <Treeselect/>
    <RVGuide/>
  </div>
</template>

<script>
  import RVGuide from './components/RVGuide.vue'
  import Treeselect from './components/TreeSelect.vue'
  export default {
    components: {
      RVGuide,
      Treeselect
    }
  }
</script>

