<template>
    <div>
        <header class="pink section" style="display:flex;">
            <div style="padding-top:10px">
                <span class="title" >ProteoVision: Advanced Visualization of Proteins </span>
            </div>
            <div class="headerOptions" style="margin-left: auto;padding-top:10px;">
                <span title="Start an interactive guide">
                    <button class="btn btn-outline-dark" v-on:click="startTour();" style="float: right;">Help</button>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Go to ProteoVision documentation">
                    <button class="btn btn-outline-dark" id="aboutButton" v-on:click="goToAboutPage();" style="float: right;">About</button>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Go to DESIRE api">
                    <a href="/desire-api/" target="_blank" id="desireAPIButton" class="btn btn-outline-dark" style="float: right;">API</a>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Reset the current session">
                    <button class="btn btn-outline-dark" id="resetButton" v-on:click="resetRV3State();" style="float: right;">Reset</button>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Save a ProteoVision session file">
                    <button class="btn btn-outline-dark" id="saveButton" v-on:click="saveRV3State();" style="float: right;">Save session</button>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Load a ProteoVision session file">
                    <label for="inputRV3State" id="rv3-state-upload" class="btn btn-outline-dark">Load session</label>
                    <input id="inputRV3State" type="file" accept=".json" ref="rv3_state_file" v-on:change="loadRV3State()"/>
                </span>
            </div>
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
    import {replacer} from './jsonMapReplacer.js'
    import {readLoadRV3State} from './loadRV3State.js'
    import {cloneDeep} from 'lodash'
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
                },
            }
        },
        methods: {
            startTour(){
                this.resetRV3State();
                vm.guideOff = false;
                this.$tours['myTour'].start()
            },
            stopTour(){
                vm.guideOff = true;
                this.resetRV3State();
            },
            skipTour(){
                vm.guideOff = true;
                this.resetRV3State();
            },
            resetRV3State(){
                Object.assign(vm.$data, initialState());
                window.tempCSVdata = null;
                clearInputFile(document.getElementById('inputRV3State'));
                vm.flushDjangoSession();
            },
            saveRV3State(){
                let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
                let anchor = document.createElement('a');
                vm.uploadSession=true;
                var saveData = _.cloneDeep(vm.$data);
                saveData.topology_loaded = false;
                saveData.postedPDBEntities = false;
                saveData.domain_or_selection = null;
                saveData.checked_domain = false;
                saveData.checked_selection = false;
                saveData.filter_range = null;
                saveData.selected_domain = [];
                saveData["window.selectSections_RV1"]=window.selectSections_RV1;
                saveData["window.aaFreqs"]=window.aaFreqs;
                saveData["window.barColors"]=window.barColors;
                anchor.href = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(saveData, replacer));
                anchor.target = '_blank';
                anchor.download = `PVState-${month}-${date}-${year}.json`;
                anchor.click();
                vm.uploadSession=false;
            },
            loadRV3State(){
                if (this.$refs.rv3_state_file.files.length == 0){return;}
                Object.assign(vm.$data, initialState());
                window.tempCSVdata = null;
                window.selectSections_RV1 = null;
                window.barColors = null;
                window.aaFreqs = null;
                readLoadRV3State(this.$refs.rv3_state_file.files[0]);
                clearInputFile(document.getElementById('inputRV3State'));
            }, downloadAboutDoc(){
                var req = new XMLHttpRequest();
                req.open("GET", `static/alignments/About.pdf`, true);
                req.responseType = "blob";
                req.onload = function (event) {
                  var blob = req.response;
                  console.log(blob.size);
                  var link=document.createElement('a');
                  link.href=window.URL.createObjectURL(blob);
                  link.download="AboutProteoVision.pdf";
                  link.click();
                };
                req.send();
            }, goToAboutPage(){
                window.open("https://apollo2.chemistry.gatech.edu/AboutProteoVision/", "_blank"); 
            },
        },
        mounted: function () {
            if (localStorage.getItem("hasCodeRunBefore") === null) {
                tourSteps[0].content += '<br><b>First time users are advised to complete this guide by only clicking the Next button ▼</b>';
                tourSteps.unshift(cookieNotice);
                tourSteps[1].before = function before(type) {
                    return new Promise((resolve, reject) => {
                        resolve (
                            localStorage.setItem("hasCodeRunBefore", true),
                            vm.guideOff = false,
                        )
                    })
                }
                this.$tours['myTour'].start();
            }
        },
    }

    const cookieNotice = {
        target: 'header',
            header: {
                title: 'Privacy Notice',
            },
        content: `ProteoVision uses two essential cookies.
            One ensures you do not see this message every time you visit the website;<br>
            the other ensures our server can validate a secure connection to your browser.<br>
            We do not store any other data from you. Uploaded CSV files are kept in your browser memory and are not stored between sessions. 
            Uploaded alignments are sent to our server for processing, however they are deleted immediately on completion of the job.<br>
            Continuing to use our website indicates you consent to this data processing.`
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
    var tourSteps = [
        {
            target: 'header',
            header: {
                title: 'Welcome to ProteoVision!',
            },
            content: `ProteoVision is a visualization tool for ribosomal proteins 
            designed to display phylogenetic, structural, and physicochemical 
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
                title: 'Phylogenetic browser',
            },
            content: `Select a phylogenetic group. This menu supports searching and selection of multiple groups.<br>
            Initially only the three major phylogenetic branches are shown, deeper branches take a few seconds to load. 
            Deeper branches can be opened by clicking on the triangle ▸ next to each parent branch.`,
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
            content: `Select an alignment from the DESIRE database.`,
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
            The viewer window can be moved by dragging or by using the scrollbars.<br/>
            Hover over residue to reveal associated data for it.`,
        },{
            target: '#conservationBar',
            header: {
                title: 'Sequence conservation',
            },
            content: `Bar representation of conservation for each alignment position.
            The conservation is calculated as Shannon entropy and represented by the column heights.
            The columns can be colored with custom data or with ProteoVision calculated data.`,
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
            content: `Download the visible part of the alignment or the entire alignment as a png image. 
            Large alignments might take some time so please be patient.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#selectColorMappingProps',
            header: {
                title: 'Select calculated data.',
            },
            content: `Select a property to color the conservation bar. Same coloring will be applied to the 
            2D and 3D viewers if they are enabled.`,
            params: {
              placement: 'right'
            }
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
            target: '#showFrequencies',
            header: {
                title: 'Show amino-acid frequencies',
            },
            content: `Select this checkbox to show amino-acid frequencies, calculated from the loaded alignment.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_propensities = true,
                    vm.handlePropensities(true)
                )
            })
        },{
            target: '#downloadFreqsBtn',
            header: {
                title: 'Download amino-acid frequencies',
            },
            content: `Download amino-acid frequencies showed in the graph as a csv file.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#total',
            header: {
                title: 'Frequency graph',
            },
            content: `Frequencies for the 20 amino-acids are shown as a boxplot figure. Each species is represented with a datapoint. <br>
            Hovering on a datapoint highlights the corresponding species in the alignment viewer.`,
        },{
            target: '#pdb_input',
            header: {
                title: 'Select a structure for 2D and 3D display',
            },
            content: `Type any PDB ID or select a one from the dropdown menu. PDBs can be searched by the species name or PDB ID.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_propensities = false,
                    vm.pdbid = "4v9d",
                    vm.$children[1].search = "4v9d",
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select a polymer for structure display',
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
                title: '3D MolStar viewer',
            },
            content: `This is the 3D viewer that shows tertiary protein structure.<br/>
            The alignment, topology, and 3D viewers have integrated hover effects.`,
        },{
            target: '.menuSelectbox',
            header: {
                title: 'Mapping data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".menuSelectbox");
                var exampleData = topviewer.pluginInstance.domainTypes[4];
                resolve (
                    vm.selected_property = "Polarity",
                    topviewer.pluginInstance.updateTheme(exampleData.data),
                    window.viewerInstance.visual.select({data: selectSections_RV1.get(exampleData.label), nonSelectedColor: {r:255,g:255,b:255}}),
                    annotationSelect.selectedIndex=4,
                )
            })
        },{
            target: '.gradient_section',
            header: {
                title: 'Data range',
            },
            content: `Gradient bar showing the range of values for the currently selected mapping data.`,
            params: {
              placement: 'left'
            },
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
                title: 'Download calculated data',
            },
            content: `Downloads data calculated from the alignment and 
            mapped on the structure residues as a PyMOL script or in csv format.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#showRNAcontext',
            header: {
                title: 'Show ribosomal context in 3D Viewer',
            },
            content: `When selected, the 3D viewer shows the entire ribosomal structure in gray
            and the protein of interest is highlighted. Color mapping is possible but can 
            take a while, since the entire ribosome is being recolored.
            The user is advised to exercise patience as the web server might be unresponsive when colormapping with this option.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#domainSelectionSection',
            header: {
                title: 'ECOD domains',
            },
            content: `ECOD annotations for the selected pdb and chain are retrieved and displayed here.
            You can select a domain by which the 2D and 3D representations will be truncated.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.domain_or_selection="domain"
                )
            })
        },{
            target: '#filterSection',
            header: {
                title: 'Custom truncation range',
            },
            content: `Here you can specify a range to truncate the structure shown on the topology and 3D viewers.
            <br>You can either truncate the structure by range or by ECOD domain.
            Amino acid frequencies will be recalculated based on the active selection between these two.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.domain_or_selection="selection"
                )
            })
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
                    vm.domain_or_selection=null,
                    vm.checked_filter=true,
                )
            })
        },{
            target: '#customDataSection',
            header: {
                title: 'Upload custom data',
            },
            content: `Upload custom data in csv format to be mapped on the alignment conservation bar 
            as well as the topology and 3D viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_filter=false,
                    vm.checked_customMap=true,
                )
            })
        },{
            target: '#downloadExampleCSV',
            header: {
                title: 'Download a mapping data example',
            },
            content: `CSV format example that is supported by ProteoVision. The file must have a header row with column labeled
            <b>Index</b> indicating the structure residues. At least one more header is necessary which 
            labels the data column. The viridis colormap is used to map the datapoints with colors.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#propensitiesSubstructure',
            header: {
                title: 'Recalculate AA frequencies',
            },
            content: `Once a structure has been defined, AA frequencies can be recalculated by a given secondary structure.<br>
            The options are Coil, Strand, or Helix residues.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_customMap=false,
                    vm.checked_propensities = true
                )
            })
        },{
            target: '#tree_type',
            header: {
                title: 'Upload a custom alignment',
            },
            content: `Using a custom alignment is the other mode of operation.<br/>
            Changing between modes clears the viewers.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                vm.topology_loaded = false;
                resolve (
                    vm.checked_propensities = false,
                    vm.type_tree="upload",
                    document.getElementById('tree_type').children[1].click(),
                    vm.cleanTreeOpts(),
                    document.getElementById("pdbeMolstarView").textContent = null,
                    document.getElementById("topview").textContent = null,
                    vm.topology_loaded = false,
                )
            })
        },{
            target: '#downloadExampleFasta',
            header: {
                title: 'Download example alignment.',
            },
            content: `Download an example of a fasta format alignment.`,
            params: {
              placement: 'right'
            },
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
            Our server will calculate amino-acid frequencies and check the format.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#uploadShowFasta")
                resolve (
                    uploadButton.click(),
                    vm.fetchingPDBwithCustomAln=true,
                )
            })
        },{
            target: '#warningCDHITtruncation',
            header: {
                title: 'Warning for clustering of alignment.',
            },
            content: `A warning will be displayed here when the uploaded alignment sequences 
            have been clustered by CD-HIT. This is done to ensure there is no overrepresentation of certain sequences.<br/>
            The number of clustered sequences at 90% identity threshold will be indicated.<br/>
            The user can input a different PDB or select a new polymer or restart with a new alignment.`,
        },{
            target: '#cdHITResults',
            header: {
                title: 'CD-HIT options.',
            },
            content: `The user can select to use their original unclustered alignment from this dropdown menu.
            The user can also download the CD-HIT report from their alignment.`,
        },{
            target: '#pdb-upload',
            header: {
                title: 'Upload custom PDB file.',
            },
            content: `When using custom alignment ProteoVision supports a custom PDB structure file.
            The structure file must be in PDB format and must contain only a single chain.`,
            params: {
              placement: 'right'
            },
        },{
            target: '.autocomplete',
            header: {
                title: 'Write a PDB ID for structure display',
            },
            content: `For uploaded alignment you can write in any PDB ID of length 4.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.pdbid = "1efu",
                    vm.$children[0].search = "1efu",
                )
            })
        },{
            target: '#blastingPDBsMSG',
            header: {
                title: 'Background BLAST search',
            },
            content: `The first sequence of the uploaded alignment will be BLASTed against the PDB database 
            to find structures with similar chains. This process takes some time so you are free to input any 4 letter PDB while waiting.
            <br>When the BLAST finishes a searchable dropdown menu will be available with the PDB results.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select polymer for structure display',
            },
            content: `For uploaded alignment while BLAST is running we do not filter the available PDB chains and 
            you can select any polymer from the PDB structure.<br>
            When the BLAST finishes, only chains with high similarity will be shown here.`,
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
            The user can input a different PDB or select a new polymer or restart with a new alignment.`,
        },{
            target: '#completeBLASTsMSG',
            header: {
                title: 'Background BLAST complete',
            },
            content: `A message will be shown here when the BLAST search is complete.
            PDB IDs with chains that are highly similar to the first sequence of the alignment will be populated in the
            <b>Input PDB</b> box as searchable dropdown menu.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var tempMap = new Map();
                tempMap.set("1EFT", ["A"]);
                tempMap.set("1EFU", ["A", "C"]);
                vm.fetchingPDBwithCustomAln = 'complete';
                resolve (
                    vm.blastPDBresult = [{id:"1EFT", name:"1EFT"},
                                        {id:"1EFU", name:"1EFU"}],
                    vm.blastMAPresult = tempMap,
                );
            })
        },{
            target: '.autocomplete',
            header: {
                title: 'Select a different PDB ID for structure display',
            },
            content: `Writing a new PDB ID will clear all data related to the old PDB.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                vm.pdbid = "1eft";
                vm.$children[0].search = "1eft";
                var pdbinput = document.querySelector('.input-group-text');
                var autoresult = document.querySelector('#autocomplete-results');
                autoresult.firstElementChild.click();
                resolve (
                    vm.$nextTick(function(){
                        vm.$children[0].isOpen=true;
                        pdbinput.click();
                    }),
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select a polymer for structure display',
            },
            content: `If you select a polymer that produces good alignment with 
            the sequence alignment, it won't raise a warning.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var polSele = document.querySelector("#polymerSelect")
                resolve (
                    vm.$children[0].isOpen=false,
                    vm.chainid = ["A"],
                    vm.$nextTick(function(){
                        polSele.lastElementChild.click();
                    }),
                )
            })
        },{
            target: '#aboutButton',
            header: {
                title: 'About ProteoVision',
            },
            content: `Redirects to a comprehensive online documentation that describes
            all functions and features of ProteoVision.`,
        },{
            target: '#desireAPIButton',
            header: {
                title: 'DESIRE API',
            },
            content: `This is a link to the DESIRE REST API that supports this webserver.`,
        },{
            target: '#resetButton',
            header: {
                title: 'Reset the session',
            },
            content: `Reset the current ProteoVision session.<br/>
            All loaded data will be removed.`,
        },{
            target: '#saveButton',
            header: {
                title: 'Save the session',
            },
            content: `Downloads a ProteoVision session file.<br>
            The state of current alignment, structure, and frequency viewers will be saved.<br>
            Masking ranges and truncation ranges will not be saved.`,
        },{
            target: '#rv3-state-upload',
            header: {
                title: 'Load a session',
            },
            content: `Upload a ProteoVision session file.<br>
            The file will load a previously saved ProteoVision session.`,
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