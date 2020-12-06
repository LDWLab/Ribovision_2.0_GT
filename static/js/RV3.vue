<template>
  <div>
    <RVGuide/>
    <DropDownTree/>
  </div>
</template>

<script>
  import RVGuide from './components/RVGuide.vue'
  import DropDownTree from './components/DropDownTree.vue'
  require('vue-tour/dist/vue-tour.css')
  export default {
    components: {
      RVGuide,
      DropDownTree
    }
  }
</script>

