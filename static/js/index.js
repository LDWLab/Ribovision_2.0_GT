import Vue from 'vue';
import App from './App.vue';
import VueTour from 'vue-tour'



//From here https://github.com/pulsardev/vue-tour
Vue.use(VueTour)

var vm = new Vue({
  el: '#app',
  template: '<App/>',
  components: { App }
})


