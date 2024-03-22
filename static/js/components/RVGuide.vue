<template>
    <div>
        <header class="pink section" style="display:flex;">
            <div style="padding-top:10px">
                <span class="title" >RiboVision 2.0: Advanced Visualization of RNA molecules </span>
            </div>
            
            <div class="headerOptions" style="margin-left: auto;padding-top:10px;">
                <p style="padding:2px;float: right;"></p>
                <span title="Go to RiboVision 2.0 documentation">
                    <button class="btn btn-outline-dark" id="aboutButton" v-on:click="goToAboutPage();" style="float: right;">About</button>
                </span>
                <span title="Start an interactive guide">
                    <button class="btn btn-outline-dark" v-on:click="startUploadTour();" style="float: right;">User Upload Tour</button>
                    <button class="btn btn-outline-dark" v-on:click="startTour();" style="float: right;">RiboVision Tour</button>
                </span>
                
                
                <!--
                <p style="padding:2px;float: right;"></p>
                <span title="Go to DESIRE api">
                    <a href="/desire-api/" target="_blank" id="desireAPIButton" class="btn btn-outline-dark" style="float: right;">API</a>
                </span>
               
                
                <p style="padding:2px;float: right;"></p>
                <span title="Reset the current session">
                    <button class="btn btn-outline-dark" id="resetButton" v-on:click="resetRV3State();" style="float: right;">Reset</button>
                </span>
                -->
                <!--<p style="padding:2px;float: right;"></p>
                <span title="Save a RiboVision 2.0 session file">
                    <button class="btn btn-outline-dark" id="saveButton" v-on:click="saveRV3State();" style="float: right;">Save session</button>
                </span>
                <p style="padding:2px;float: right;"></p>
                <span title="Load a RiboVision 2.0 session file">
                    <label for="inputRV3State" id="rv3-state-upload" class="btn btn-outline-dark">Load session</label>
                    <input id="inputRV3State" type="file" accept=".json" ref="rv3_state_file" v-on:change="loadRV3State()"/>
                </span>-->
            </div>
            
        </header>
        <v-tour 
          name="myTour"
          :steps="steps" 
          :callbacks="myCallbacks"
          :options="{ highlight: true }"
          >
        </v-tour>
        <v-tour 
          name="userUploadTour"
          :steps="userUploadSteps" 
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
    import {uploadCustomFullSequence} from './handleUploadCustomFullSequence.js'
    import {cloneDeep} from 'lodash'
    var riboTourFinished = false
    export default {
        name: 'my-tour',
        data () {
            return {
                userUploadSteps: uploadTourSteps,
                steps: tourSteps,
                myOptions: {
                    useKeyboardNavigation: false,
                },
                myCallbacks: {
                    onSkip: this.skipTour,
                    onFinish: this.stopTour,
                },
                uploadCallbacks: {
                    onSkip: this.skipTour,
                    onFinish: this.stopTour,
                }
            }
        },
        methods: {
            startTour(){
                this.resetRV3State();
                vm.guideOff = false;
                vm.tourType = "Ribovision"
                this.$tours['myTour'].start()
            },
            startUploadTour(){
                this.resetRV3State();
                vm.guideOff = false;
                vm.tourType = "Upload"
                this.$tours['userUploadTour'].start()
            },
            stopTour(){
                if (localStorage.getItem("uploadTourCompleted") != null || vm.tourType != null) {
                    vm.guideOff = true;
                    vm.tourType = null;
                    this.resetRV3State();
                } else {
                    this.$tours['userUploadTour'].start();
                    localStorage.setItem("uploadTourCompleted", true)
                }
            },
            skipTour(){
                vm.guideOff = true;
                this.resetRV3State();
            },
            resetRV3State(){
                Object.assign(vm.$data, initialState());
                window.tempCSVdata = null;
                //clearInputFile(document.getElementById('inputRV3State'));
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
                window.open("https://apollo2.chemistry.gatech.edu/AboutRiboVision2/about/", "_blank"); 
            },
        },
        mounted: function () {
            //vm.guideOff = false
            //localStorage.setItem("hasCodeRunBefore", false);
            //if (localStorage.getItem("hasCodeRunBefore") !== 'true') {
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
                this.$tours['myTour'].start()

                /*this.$tours['myTour'].start(() => {
                    this.$tours['userUploadTour'].start();
                });*/
            }
            //localStorage.setItem("hasCodeRunBefore", true)
            
        },
    }

    const cookieNotice = {
        target: 'header',
            header: {
                title: 'Privacy Notice',
            },
        content: `RiboVision 2.0 uses two essential cookies.
            One ensures you do not see this message every time you visit the website;<br>
            the other ensures our server can validate a secure connection to your browser.<br>
            We do not store any other data from you. Uploaded CSV files are kept in your browser memory and are not stored between sessions. 
            Uploaded alignments are sent to our server for processing; however, they are deleted immediately upon completion of the job.<br>
            Continuing to use our website indicates you consent to this data processing.`
    }

    var getExampleFasta = function(){
        $.ajax({
            url: `static/alignments/RNAseP_example.fas`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                vm.file = new File([data], "RNAseP_example.fas", {});
            },
        })
    };
    var getExampleCIF = function(){
        $.ajax({
            url: `static/alignments/3Q1Q.cif`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                vm.$refs.customCIFfile.files = [new File([data], "3Q1Q.cif", {})];
            },
        })
    };
    var getExamplePDB = function(){
        $.ajax({
            url: `static/alignments/3Q1Q_B.pdb`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                let list = new DataTransfer();
                let file =new File([data], "3Q1Q_B.pdb", {});
                list.items.add(file);

                let myFileList = list.files;
                vm.$refs.customPDBfile.files = myFileList;
            },
        })
    };
    var getExampleFullSeqPDB = function(){
        $.ajax({
            url: `static/alignments/3Q1Q_B.fas`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                let list = new DataTransfer();
                let file =new File([data], "3Q1Q_B.fas", {});
                list.items.add(file);

                let myFileList = list.files;
                vm.$refs.customFullSequence.files = myFileList;
                uploadCustomFullSequence();
            },
        })
    };
    var getExampleDataCSV = function(){
        $.ajax({
            url: `static/alignments/custom_data.csv`,
            type: 'GET',
            dataType: "text",
            success: function(data) {
                let list = new DataTransfer();
                let file =new File([data], "custom_data.csv", {});
                list.items.add(file);

                let myFileList = list.files;
                vm.$refs.custom_csv_file.files = myFileList;
                vm.csv_data=new File([data], "custom_data.csv", {});
                handleCustomMappingData();
            },
        })
    };
    var uploadTourSteps = [
        /*{
            target: 'header',
            header: {
                title: 'User-Upload Mode',
            },
            content: `In user-upload mode, you can visualize custom structures and alignments.`
        },*/
        {
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
                    vm.type_tree="upload",
                    document.getElementById('tree_type').children[1].click(),
                    vm.cleanTreeOpts(),
                    //document.getElementById("pdbeMolstarView").textContent = null,
                    //document.getElementById("topview").textContent = null,
                    vm.topology_loaded = false,
                )
            })
        },{
            target: '#downloadExampleFasta',
            header: {
                title: 'Download example alignment.',
            },
            content: `Download an example of a fasta-format alignment.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#inputUploadFasta',
            header: {
                title: 'Input custom alignment',
            },
            content: `Select a fasta-format alignment from your computer.`,
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
            Our server will calculate amino-acid frequencies and check the format.
            WAIT UNTIL MSA IS LOADED BEFORE CLICKING "Next".`,
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
        },
        {
            target: '#warningCDHITtruncation',
            header: {
                title: 'Warning for clustering of alignment.',
            },
            content: `A warning will be displayed here when the uploaded alignment sequences 
            have been clustered by CD-HIT. This is done to ensure there is no overrepresentation of certain sequences.<br/>
            The number of clustered sequences at 90% identity threshold will be indicated.<br/>
            The user can input a different PDB or select a new polymer or restart with a new alignment.`,
        },
        {
            target: '#cdHITResults',
            header: {
                title: 'CD-HIT options.',
            },
            content: `The user can select to use their original unclustered alignment from this dropdown menu.
            The user can also download the CD-HIT report from their alignment.`,
        },
        {
            target: '#radioCIF',
            header: {
                title: 'Select CIF format for the custom 3D structure. ',
            },
            content: `Users have an option to select the format of the 3D custom structure (CIF or PDB).`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#radioCIF")
                resolve (
                    uploadButton.click(),
                    
                )
            })
        },{
            target: '#cif-upload',
            header: {
                title: 'Upload custom CIF file.',
            },
            content: `If the CIF option is selected, the users will need to provide a 3D structure in the CIF format and the entity ID of the desired RNA chain. The additional requirements for the CIF file format are detailed in the documentation.`,
            params: {
              placement: 'right'
            },
        
        },{
            target: '#radioPDB',
            header: {
                title: 'Select PDB format for the custom 3D structure.',
            },
            content: `RiboVision 2.0 also supports a custom structure in the PDB format.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    uploadButton.click(),
                    
                )
            })
        },{
            target: '#pdb-upload',
            header: {
                title: 'Upload custom PDB file.',
            },
            content: `The PDB structure file must contain only a single RNA chain.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    getExamplePDB(),
                    vm.pdbFileUploadedFlag=true,
                    
                )
            })
        },

        {
            target: '#full-sequence-upload',
            header: {
                title: 'Upload a complete RNA sequence',
            },
            content: `In addition to the PDB structure, a complete RNA sequence for the structure in the PDB file must also be uploaded separately. This sequence is required to generate the complete 2D diagram without omission of unresolved regions.
            WAIT until 2D and 3D structures are generated BEFORE CLICKING "Next"`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    getExampleFullSeqPDB(),
                    
                    
                )
            })
        },{
            target: '#downloadDataBtn',
            header: {
                title: 'Download Custom Data.',
            },
            content: `Custom mode also allows the users to map and  visualize custom data. The data must be uploaded as a CSV file.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    //vm.checked_customMap = true,
                    //vm.cleanCustomMap(true),
                    
                    
                )
            })
        },
        {
            target: '#uploadCustomData',
            header: {
                title: 'Choose to upload custom data from CSV file.',
            },
            content: `The data in the CSV file are organized in a specific way. Please download an example of the custom CSV file here.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                
                resolve (
                    vm.checked_customMap = true,
                    vm.cleanCustomMap(true),
                    
                    
                )
            })
        },
        {
            target: '#inputUploadCSV',
            header: {
                title: 'Select and upload  CSV data file.',
            },
            content: `CSV format example that is supported by RiboVision 2.0. The file must have a header row with column labeled Index indicating the structure residues. At least one more header is necessary which labels the data column. The viridis colormap is used to map the datapoints with colors.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                
                resolve (

                    getExampleDataCSV()
                    
                    
                )
            })
        },
        {
            target: '.menuSelectbox',
            header: {
                title: 'Representation of Custom Data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".menuSelectbox");
                var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');

                
                resolve (
                    vm.selected_property = "circle",
                    
                    selectBoxEle.value="2",
                    selectBoxEle.dispatchEvent(new Event('change')),

                )
            })
        },
        {
            target: '.mappingSelectbox',
            header: {
                title: 'Custom Mapping data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".mappingSelectbox");
                var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
                var exampleData = topviewer.viewInstance.uiTemplateService.domainTypes[2];
                
                
                resolve (
                    vm.selected_property = "Custom Data",
                    topviewer.viewInstance.uiTemplateService.domainTypes[2],
                    
                    annotationSelect.selectedIndex=1
                )
            })
        },
       
        {
            target: '#aboutButton',
            header: {
                title: 'About RiboVision 2.0',
            },
            content: `Redirects to a comprehensive online documentation that describes
            all functions and features of RiboVision 2.0.`,
        },
        //{
            //target: '#resetButton',
            //header: { 
            //    title: 'Reset the session',
            //},
            //content: `Reset the current ProteoVision session.<br/>
            //All loaded data will be removed.`,
        //},
        /*{
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
        },*/{
            target: 'footer',
            header: {
                title: 'Thank you',
            },
            content: `Thank you for reading our guide!
            Ending this guide will reset the session.`,
        },
    ]
    var tourSteps = [
        {
            target: 'header',
            header: {
                title: 'Welcome to RiboVision 2.0!',
            },
            content: `RiboVision 2.0 is a web server for visualization of (ribosomal) RNAs 
            designed to display phylogenetic, structural, and evolutionary 
            properties in primary, secondary, and tertiary representations.`
        }, 
        
        {
            target: '#tree_type',
            header: {
                title: 'Mode of operation',
            },
            content: `Select either of two possible modes of operation.<br/>
            <b>RiboVision 2.0 </b> retrieves alignments from the DatabasE for Study and Imaging of Ribosomal Evolution (DESIRE).<br/>
            <b>User upload</b> allows you to upload your own fasta-formatted alignment and 3D structure.`,
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
            target: '#select_protein_type',
            header: {
                title: 'RNA family',
            },
            content: `Select an RNA family.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var selectEl = document.querySelector('#select_protein_type');
                resolve (
        
                    
                    selectEl.value = "LSU-rRNA",
                    selectEl.dispatchEvent(new Event('change'))
                    
    
                )
            })
        },{
            target: '#selectaln',
            header: {
                title: 'Alignment selection',
            },
            content: `Select an alignment from the RiboVision database.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
               
                resolve (
                    vm.alnobj = {id: 256, text: "5S"},
                    
                    
                )
            })
        },{
            target: '.alignment_section',
            header: {
                title: 'Alignment viewer',
            },
            content: `This is the MSA alignment viewer. 
            The viewer window can be moved by dragging or by using the scrollbars.<br/>
            Hover over a residue to reveal its associated data.`,
        },{
            target: '#conservationBar',
            header: {
                title: 'Sequence conservation',
            },
            content: `Bar representation of conservation for each alignment position.
            The conservation is calculated as Shannon entropy and is represented by the column heights.
            These columns can be colored with custom data or with RiboVision-2.0-calculated data.`,
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
            Large alignments might take some time, so please be patient.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#selectColorMappingProps',
            header: {
                title: 'Select calculated data.',
            },
            content: `Select a property to color the conservation bar. The same coloring will be applied to the 
            2D and 3D viewers if they are enabled.`,
            params: {
              placement: 'right'
            }
        },{
            target: '#selectAlnColorScheme',
            header: {
                title: 'Select alignment color scheme.',
            },
            content: `Select an alignment color scheme from those available.`,
            params: {
              placement: 'right'
            },
        },
        
        {
            target: '#pdb_input',
            header: {
                title: 'Select a structure for 2D and 3D display',
            },
            content: `Type any PDB ID or select one from the dropdown menu. PDBs can be searched by the species name or PDB ID.`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                resolve (
                    vm.checked_propensities = false,
                    vm.pdbid = "7k00",
                    vm.$children[1].search = "7k00",
                )
            })
        },{
            target: '#polymerSelect',
            header: {
                title: 'Select a polymer for structure display',
            },
            content: `Select one of the filtered polymers to see 2D and 3D structures. WAIT UNTIL 2D and 3D STRUCTURES LOADED BEFORE CLICKING "Next".`,
            params: { 
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                var polSele = document.querySelector("#polymerSelect")
                vm.RVGuideEntityId = 23;
                resolve (
                    vm.chainid = ["b"],
                    vm.$nextTick(function(){
                        polSele.lastElementChild.click();
                    }),
                )
            })
        },{
            target: '.topology_section',
            header: {
                title: 'RNA topology viewer',
            },
            content: `This is the RNA topology viewer that depicts 2D RNA layout and base pairing.`,
        },{
            target: '.molstar_section',
            header: {
                title: '3D MolStar viewer',
            },
            content: `This is the 3D viewer that shows tertiary protein structure.<br/>
            The alignment, topology, and 3D viewers feature integrated hover effects.`,
        },{
            target: '.mappingSelectbox',
            header: {
                title: 'Mapping data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".mappingSelectbox");
                
               
                var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
                var exampleData = topviewer.viewInstance.uiTemplateService.domainTypes[1];
                
                resolve (
                    vm.selected_property = "Shannon entropy",
                    topviewer.viewInstance.uiTemplateService.domainTypes[1],
                    //window.viewerInstance.visual.select({data: selectSections_RV1.get(exampleData.label), nonSelectedColor: {r:255,g:255,b:255}}),
                    annotationSelect.selectedIndex=1
                )
            })
        },{
            target: '#basePairingSelectElement',
            header: {
                title: 'Base-pairing data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var menuToClick = document.querySelector("#basePairingSelectElement");
                menuToClick.click();
                var checkBoxAll = document.querySelector("#Checkbox_cWS");
                checkBoxAll.checked = true;
                let x = document.getElementById("PdbeTopViewer");
                x.viewInstance.uiTemplateService.changeBP("cWS", false);
                resolve (
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
            before: type => new Promise((resolve, reject) => {
                var menuToClick = document.querySelector("#basePairingSelectElement");
                menuToClick.click();
                resolve();
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
            target: '#rnaTopologyReset-7k00',
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
            The user is advised to exercise patience, as the web server might be unresponsive when colormapping with this option.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#customDataSection',
            header: {
                title: 'Upload custom data',
            },
            content: `Upload custom data in CSV format to be mapped on the alignment conservation bar 
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
            content: `CSV format example that is supported by RiboVision 2.0. The file must have a header row with column labeled
            <b>Index</b> indicating the structure residues. At least one more header is necessary which 
            labels the data column. The viridis colormap is used to map the datapoints with colors.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#proteinSelect',
            header : {
                title: 'Select RNA-protein contacts to visualize'
            },
            content: 'The RNA contacts with selected proteins are represented by coloring respective RNA nucleotides within the 2D and 3D structure viewers. The 3D viewer also shows the structure of the selected protein.',
            params: {
                placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                const newValue = "50S ribosomal protein L25";
                //var polymerSelect2 = document.querySelector("#polymerSelect2");

                var topviewer = document.getElementById("PdbeTopViewer");
                
                
                
                //polymerSelect2.value = newValue;
                //polymerSelect2.dispatchEvent(new Event('change'))
                vm.selectedProteins = ['f'];
                vm.pchainid = ['f'];
                const newEntityID = 23;
                
                vm.entityID = newEntityID;
                topviewer.viewInstance.uiTemplateService.colorMapContacts(); 
                topviewer.viewInstance.uiTemplateService.colorMapContacts(); 
                showModificationsAndContactsHelper("" + newEntityID);
                //setTimeout(function() {
                    //secondFunctionCall();
                    //showModificationsAndContactsHelper("" + newEntityID);
                    //}, 5100);
                resolve();
            }),/*
            after: type => new Promise((resolve, reject) => {
                document.getElementsByClassName("v-step__button v-step__button-stop")[0].addEventListener('click', function() {
                    riboTourFinished = true;
                });
                resolve();
            })*/
           
        },
        /*
        {
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
                    vm.type_tree="upload",
                    document.getElementById('tree_type').children[1].click(),
                    vm.cleanTreeOpts(),
                    //document.getElementById("pdbeMolstarView").textContent = null,
                    //document.getElementById("topview").textContent = null,
                    vm.topology_loaded = false,
                )
            })
        },{
            target: '#downloadExampleFasta',
            header: {
                title: 'Download example alignment.',
            },
            content: `Download an example of a fasta-format alignment.`,
            params: {
              placement: 'right'
            },
        },{
            target: '#inputUploadFasta',
            header: {
                title: 'Input custom alignment',
            },
            content: `Select a fasta-format alignment from your computer.`,
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
            Our server will calculate amino-acid frequencies and check the format.
            WAIT UNTIL MSA IS LOADED BEFORE CLICKING "Next".`,
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
        },
        {
            target: '#warningCDHITtruncation',
            header: {
                title: 'Warning for clustering of alignment.',
            },
            content: `A warning will be displayed here when the uploaded alignment sequences 
            have been clustered by CD-HIT. This is done to ensure there is no overrepresentation of certain sequences.<br/>
            The number of clustered sequences at 90% identity threshold will be indicated.<br/>
            The user can input a different PDB or select a new polymer or restart with a new alignment.`,
        },
        {
            target: '#cdHITResults',
            header: {
                title: 'CD-HIT options.',
            },
            content: `The user can select to use their original unclustered alignment from this dropdown menu.
            The user can also download the CD-HIT report from their alignment.`,
        },
        {
            target: '#radioCIF',
            header: {
                title: 'Select CIF format for the custom 3D structure. ',
            },
            content: `Users have an option to select the format of the 3D custom structure (CIF or PDB).`,
            params: {
              placement: 'right'
            },
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#radioCIF")
                resolve (
                    uploadButton.click(),
                    
                )
            })
        },{
            target: '#cif-upload',
            header: {
                title: 'Upload custom CIF file.',
            },
            content: `If the CIF option is selected, the users will need to provide a 3D structure in the CIF format and the entity ID of the desired RNA chain. The additional requirements for the CIF file format are detailed in the documentation.`,
            params: {
              placement: 'right'
            },
        
        },{
            target: '#radioPDB',
            header: {
                title: 'Select PDB format for the custom 3D structure.',
            },
            content: `RiboVision 2.0 also supports a custom structure in the PDB format.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    uploadButton.click(),
                    
                )
            })
        },{
            target: '#pdb-upload',
            header: {
                title: 'Upload custom PDB file.',
            },
            content: `The PDB structure file must contain only a single RNA chain.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    getExamplePDB(),
                    vm.pdbFileUploadedFlag=true,
                    
                )
            })
        },

        {
            target: '#full-sequence-upload',
            header: {
                title: 'Upload a complete RNA sequence',
            },
            content: `In addition to the PDB structure, a complete RNA sequence for the structure in the PDB file must also be uploaded separately. This sequence is required to generate the complete 2D diagram without omission of unresolved regions.
            WAIT until 2D and 3D structures are generated BEFORE CLICKING "Next"`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    getExampleFullSeqPDB(),
                    
                    
                )
            })
        },{
            target: '#downloadDataBtn',
            header: {
                title: 'Download Custom Data.',
            },
            content: `Custom mode also allows the users to map and  visualize custom data. The data must be uploaded as a CSV file.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                //let uploadButton = document.querySelector("#radioPDB")
                resolve (
                    //vm.checked_customMap = true,
                    //vm.cleanCustomMap(true),
                    
                    
                )
            })
        },
        {
            target: '#uploadCustomData',
            header: {
                title: 'Choose to upload custom data from CSV file.',
            },
            content: `The data in the CSV file are orgainized in a specific way. Please download an example of the custom CSV file here.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                
                resolve (
                    vm.checked_customMap = true,
                    vm.cleanCustomMap(true),
                    
                    
                )
            })
        },
        {
            target: '#inputUploadCSV',
            header: {
                title: 'Select and upload  CSV data file.',
            },
            content: `CSV format example that is supported by RiboVision 2.0. The file must have a header row with column labeled Index indicating the structure residues. At least one more header is necessary which labels the data column. The viridis colormap is used to map the datapoints with colors.`,
            params: {
              placement: 'right'
            },
        
            before: type => new Promise((resolve, reject) => {
                
                resolve (

                    getExampleDataCSV()
                    
                    
                )
            })
        },
        {
            target: '.menuSelectbox',
            header: {
                title: 'Representation of Custom Data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".menuSelectbox");
                var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');

                
                resolve (
                    vm.selected_property = "circle",
                    
                    selectBoxEle.value="2",
                    selectBoxEle.dispatchEvent(new Event('change')),

                )
            })
        },
        {
            target: '.mappingSelectbox',
            header: {
                title: 'Custom Mapping data',
            },
            content: `Calculated mapping data from the alignment can be selected from this dropdown menu.<br/>
            The data gets mapped on the alignment conservation bar, as well as the topology and 3D viewers.`,
            params: {
              placement: 'left'
            },
            before: type => new Promise((resolve, reject) => {
                var topviewer = document.getElementById("PdbeTopViewer");
                var annotationSelect = document.querySelector(".mappingSelectbox");
                var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
                var exampleData = topviewer.viewInstance.uiTemplateService.domainTypes[2];
                
                
                resolve (
                    vm.selected_property = "Custom Data",
                    topviewer.viewInstance.uiTemplateService.domainTypes[2],
                    
                    annotationSelect.selectedIndex=1
                )
            })
        },
       
        {
            target: '#aboutButton',
            header: {
                title: 'About RiboVision 2.0',
            },
            content: `Redirects to a comprehensive online documentation that describes
            all functions and features of RiboVision 2.0.`,
        },
        //{
            //target: '#resetButton',
            //header: { 
            //    title: 'Reset the session',
            //},
            //content: `Reset the current ProteoVision session.<br/>
            //All loaded data will be removed.`,
        //},
        /*{
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
        },*//*{
            target: 'footer',
            header: {
                title: 'Thank you',
            },
            content: `Thank you for reading our guide!
            Ending this guide will reset the session.`,
        },*/
    ]

</script>