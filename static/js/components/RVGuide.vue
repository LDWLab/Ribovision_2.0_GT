<template>
    <div>
        <header class="pink section">
            <span class="title">RiboVision3 Webserver for Study and Imaging of Ribosomal Protein Evolution </span>
            <button class="btn btn-outline-dark" v-on:click="startTour();" style="float: right;">Help</button>
            <p style="padding:5px;float: right;"></p>
            <button class="btn btn-outline-dark" id="resetButton" v-on:click="resetRV3State();" style="float: right;">Reset</button>
            <p style="padding:5px;float: right;"></p>
            <button class="btn btn-outline-dark" id="saveButton" v-on:click="saveRV3State();" style="float: right;">Save session</button>
            <p style="padding:5px;float: right;"></p>
            <label for="inputRV3State" id="rv3-state-upload" class="btn btn-outline-dark">
                Upload session
            </label>
            <input id="inputRV3State" type="file" accept=".json" ref="rv3_state_file" v-on:change="loadRV3State()"/>
        </header>
        <v-tour 
          name="myTour"
          :steps="steps" 
          :callbacks="myCallbacks"
          :options="{ highlight: true }"
          >
        </v-tour>
    </div>
</template>

<script>
    import {initialState} from './DropDownTreeVars.js'
    import {readLoadRV3State} from './loadRV3State.js'
    export default {
        name: 'my-tour',
        data () {
            return {
                steps: tourSteps,
                myOptions: {
                    useKeyboardNavigation: false,
                },
                myCallbacks: {
                    onSkip: this.skipTour,
                    onFinish: this.stopTour
                }
            }
        },
        methods: {
            startTour(){
                this.$tours['myTour'].start()
            },
            stopTour(){
                Object.assign(vm.$data, initialState());
            },
            skipTour(){
                Object.assign(vm.$data, initialState());
            },
            resetRV3State(){
                Object.assign(vm.$data, initialState());
            },
            saveRV3State(){
                let anchor = document.createElement('a');
                vm.uploadSession=true;
                var saveData = vm.$data;
                saveData["window.selectSections_RV1"]=window.selectSections_RV1;
                anchor.href = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(saveData, replacer));
                anchor.target = '_blank';
                anchor.download = "rv3State.json";
                anchor.click();
                vm.uploadSession=false;
            },
            loadRV3State(){
                if (this.$refs.rv3_state_file.files.length == 0){return;}
                Object.assign(vm.$data, initialState());
                readLoadRV3State(this.$refs.rv3_state_file.files[0]);
                this.$refs.rv3_state_file.files = null;
            }
        },
        mounted: function () {
            if (localStorage.getItem("hasCodeRunBefore") === null) {
                this.$tours['myTour'].start();
                localStorage.setItem("hasCodeRunBefore", true);
            }
        },
    }

    var getExampleFasta = function(){
        $.ajax({
            url: `static/alignments/EFTU_example.fas`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                vm.file = new File([data], "EFTU_example.fas", {});
            },
        })
    };
    const tourSteps = [
        {
            target: 'header',
            header: {
                title: 'Welcome to RiboVision3!',
            },
            content: `RiboVision3 is a visualization tool for ribosomal proteins 
            designed to visualize phylogenetic, structural, and physicochemical 
            properties in primary, secondary, and tertiary representations.`
        },{
            target: '#tree_type',
            header: {
                title: 'Mode of operation',
            },
            content: `Select either of two possible modes of operation.<br/>
            <b>DESIRE</b> retrieves alignments from the DatabasE for Study and Imaging of Ribosomal Evolution.<br/>
            <b>User upload</b> allows you to upload your own fasta formatted alignment.`,
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.type_tree="orth",
                )
            })
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
                vm.type_tree="orth";
                var treeselectEl = vm.$refs["treeselect"];
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
        },{
            target: '#downloadFastaBtn',
            header: {
                title: 'Download alignment',
            },
            content: `Download the generated alignment in a fasta format.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#downloadAlnImageBtn',
            header: {
                title: 'Download alignment image',
            },
            content: `Download the visible part of the alignment as a png image.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#selectAlnColorScheme',
            header: {
                title: 'Select alignment color scheme.',
            },
            content: `Select an alignment color scheme from the available ones.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#pdb_input',
            header: {
                title: 'Select PDB id for structure display',
            },
            content: `Select a PDB from the available ones in the dropdown menu.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.pdbid = "4v9d",
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select polymer for structure display',
            },
            content: `Select one of the filtered polymers to see 2D and 3D structures.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var polSele = document.querySelector("#polymerSelect")
                resolve (
                    vm.chainid = ["CC"],
                    vm.$nextTick(function(){
                        polSele.lastElementChild.click();
                    }),
                )
            })
        },{
            target: '.topology_section',
            header: {
                title: 'Topology viewer',
            },
            content: `This is the topology viewer that shows secondary protein structure.`,
        },{
            target: '.molstar_section',
            header: {
                title: '3D viewer',
            },
            content: `This is the 3D viewer that shows tertiary protein structure.<br/>
            The alignment, topology, and 3D viewers have integrated hover effects.`,
        },{
            target: '.menuSelectbox',
            header: {
                title: 'Annotation data',
            },
            content: `Calculated annotation data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".menuSelectbox");
                var exampleData = topviewer.pluginInstance.domainTypes[4];
                resolve (
                    topviewer.pluginInstance.updateTheme(exampleData.data),
                    window.viewerInstance.visual.select({data: selectSections_RV1.get(exampleData.label), nonSelectedColor: {r:255,g:255,b:255}}),
                    annotationSelect.selectedIndex=4,
                )
            })
        },{
            target: '.saveSVG',
            header: {
                title: 'Save an SVG image.',
            },
            content: `Saves the current view of the topology viewer`,
            params: {
              placement: 'left'
            },
        },{
            target: '.resetIcon',
            header: {
                title: 'Reset the view.',
            },
            content: `Resets the current view of the topology viewer`,
            params: {
              placement: 'right'
            },
        },{
            target: '#downloadDataBtn',
            header: {
                title: 'Download annotation data',
            },
            content: `Downloads annotation data calculated from the 
            alignment and mapped on the structure residues in csv format.`,
        },{
            target: '#maskingSection',
            header: {
                title: 'Masking ranges',
            },
            content: `Here you can specify ranges that mask mapped data on the topology and 3D viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_filter=true,
                )
            })
        },{
            target: '#filterSection',
            header: {
                title: '3D viewer',
            },
            content: `Here you can specify ranges that truncate the structure shown on the topology and 3D viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_filter=false,
                    vm.checked_selection=true,
                )
            })
        },{
            target: '#customDataSection',
            header: {
                title: 'Upload custom data',
            },
            content: `Upload custom data in csv format to be maped on the topology and 3D viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_selection=false,
                    vm.checked_customMap=true,
                )
            })
        },{
            target: '#tree_type',
            header: {
                title: 'Upload custom alignment',
            },
            content: `Using a custom alignment is the other mode of operation.<br/>
            Changing between modes clears the viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_customMap=false,
                    vm.type_tree="upload",
                    vm.cleanTreeOpts(),
                    document.getElementById("pdbeMolstarView").textContent = null,
                    document.getElementById("topview").textContent = null,
                )
            })
        },{
            target: '#inputUploadFasta',
            header: {
                title: 'Input custom alignment',
            },
            content: `Select a fasta format alignment from your computer.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    getExampleFasta(),
                )
            })
        },{
            target: '#uploadShowFasta',
            header: {
                title: 'Upload the chosen alignment.',
            },
            content: `The alignment will be sent to our server, but it won't be stored there. <br/>
            Our server will calculate amino-acid propensities and check the format.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#uploadShowFasta")
                resolve (
                    uploadButton.click(),
                )
            })
        },{
            target: '#pdb_input_custom',
            header: {
                title: 'Write a PDB ID for structure display',
            },
            content: `In the case of uploaded alignment we let you write in any PDB ID of length 4.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.pdbid = "1efu",
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select polymer for structure display',
            },
            content: `In the case of uploaded alignment we do not filter the avaialable PDB chains. <br/>
            You can select any polymer from the PDB structure.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var polSele = document.querySelector("#polymerSelect")
                resolve (
                    vm.chainid = ["B"],
                    vm.$nextTick(function(){
                        polSele.lastElementChild.click();
                    }),
                )
            })
        },{
            target: '#warningPoorStructureAln',
            header: {
                title: 'Warning for poor alignment.',
            },
            content: `A warning will be displayed here when the selected structure and alignment sequences 
            have poor alignment.<br/>
            The number of misaligned positions will be indicated.<br/>
            The user can input a different pdb or select a new chain or restart with a new alignment.`,
        },{
            target: '#pdb_input_custom',
            header: {
                title: 'Write a different PDB ID for structure display',
            },
            content: `Writing a new PDB id will clear all data related to the old PDB.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.pdbid = "1eft",
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select polymer for structure display',
            },
            content: `Selecting a polymer that produces good alignment with the sequence alignment does not raise a warning.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var polSele = document.querySelector("#polymerSelect")
                resolve (
                    vm.chainid = ["A"],
                    vm.$nextTick(function(){
                        polSele.lastElementChild.click();
                    }),
                )
            })
        },{
            target: '#resetButton',
            header: {
                title: 'Reset session',
            },
            content: `Reset the current RiboVision3 session.<br/>
            All loaded data will be removed.`,
        },{
            target: '#saveButton',
            header: {
                title: 'Save the session',
            },
            content: `Downloads a RiboVision3 session file.<br>
            The state of current alignment, structure and frequency viewers will be saved.<br>
            Masking ranges, truncation ranges, or custom mapping data will not be saved.`,
        },{
            target: '#rv3-state-upload',
            header: {
                title: 'Load a session',
            },
            content: `Upload a RiboVision3 session file.<br>
            The file will load a previously saved RiboVision3 session.`,
        },{
            target: 'footer',
            header: {
                title: 'Thank you',
            },
            content: `Thank you for reading our guide!
            Ending this guide will reset the session.`,
        },
    ]

</script>