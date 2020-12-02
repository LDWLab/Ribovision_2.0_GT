import Vue from 'vue';
import RV3 from './RV3.vue';
import VueTour from 'vue-tour'



//From here https://github.com/pulsardev/vue-tour
Vue.use(VueTour)

var vm_both = new Vue({
  el: '#app',
  template: '<RV3/>',
  components: { RV3 }
})

window.vm = vm_both.$children[0].$children[1];