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
            content: `Select on three possible modes of operation.<br/>
            <b>Orthologs</b> retrieves orthologous alignments.<br/>
            <b>Paralogs</b> retrieves paralogous alignments.<br/>
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
                resolve (
                    vm.alnobj = {id: 1, text: "uL02"},
                )
            })
        },{
            target: '.alignment_section',
            header: {
                title: 'Alignment viewer',
            },
            content: `This is the alignment viewer. 
            Hover over residue to reveal additional data for it.<br/>
            The viewer window can be moved by dragging or by using the scrollbars.`,
        },
    ]

    export default {
        name: 'my-tour',
        data () {
            return {
                steps: tourSteps
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