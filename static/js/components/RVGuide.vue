<template>
  <div>
    <header class="pink section">DESIRE: DatabasE for Study and Imaging of Ribosomal Evolution
    <button v-on:click="startTour();" style="float: right;">Help</button></header>
    <v-tour name="myTour" :steps="steps" :options="{ highlight: true }"></v-tour>
  </div>
</template>

<script>
  export default {
    name: 'my-tour',
    data () {
      return {
        steps: [
          {
            target: '#v-step-0',  // We're using document.querySelector() under the hood
            header: {
              title: 'Alignment viewer',
            },
            content: `Alignment viewer!`
          },
          {
            target: '#treeselect',
            content: 'Selection!'
          },
          {
            target: '[data-v-step="2"]',
            content: 'Try it, you\'ll love it!<br>You can put HTML in the steps and completely customize the DOM to suit your needs.',
            params: {
              placement: 'top' // Any valid Popper.js placement. See https://popper.js.org/popper-documentation.html#Popper.placements
            }
          }
        ]
      }
    },
    methods: {
      startTour(){
        this.$tours['myTour'].start()
      }
    }
  }
</script>