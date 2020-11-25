<template>
    <div>
        <header class="pink section">DESIRE: DatabasE for Study and Imaging of Ribosomal Evolution
        <button v-on:click="startTour();" style="float: right;">Help</button></header>
        <v-tour name="myTour" :steps="steps" :options="{ highlight: true }"></v-tour>
    </div>
</template>

<script>
    const tourSteps = [
        {
            target: 'header',
            header: {
                title: 'Welcome to RiboVision3!',
            },
            content: ``
        },{
            target: '#tree_type',
            header: {
                title: 'Mode of operation',
            },
            content: `Select on three possible modes of operation.</br>
            <b>Orthologs</b> retrieves orthologous alignments from our database.</br>
            <b>Paralogs</b> retrieves paralogous alignments from our database.</br>
            <b>Upload</b> allows you to upload your own fasta formatted alignment.`,
        },{
            target: '#treeselect',
            header: {
                title: 'Phylogenetic selection',
            },
            content: `Select a phylogenetic group. Supports searching and multiple groups.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var treeselectEl = vm.$refs["treeselect"]
                resolve (
                    treeselectEl.$emit('input', [2]),
                    treeselectEl.openMenu()
                )
            })
        },{
            target: '#selectaln',
            header: {
                title: 'Alignment selection',
            },
            content: `Select an alignment from our database.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var selectAlnEl = document.querySelector("#selectaln")
                resolve (
                    selectAlnEl.click()
                )
            })
        }

    ]

    export default {
        name: 'my-tour',
        data () {
            return {
                steps: tourSteps
                //[
                //    {
                //        target: '#v-step-0',    // We're using document.querySelector() under the hood
                //        header: {
                //            title: 'Alignment viewer',
                //        },
                //        content: `Alignment viewer!`
                //    },
                //    {
                //        target: '#treeselect',
                //        content: 'Selection!'
                //    },
                //    {
                //        target: '[data-v-step="2"]',
                //        content: 'Try it, you\'ll love it!<br>You can put HTML in the steps and completely customize the DOM to suit your needs.',
                //        params: {
                //            placement: 'top' // Any valid Popper.js placement. See https://popper.js.org/popper-documentation.html#Popper.placements
                //        }
                //    }
                //]
            }
        },
        methods: {
            startTour(){
                this.$tours['myTour'].start()
            }
        },
        mounted: function () {
            this.$tours['myTour'].start()
        }
    }
</script>