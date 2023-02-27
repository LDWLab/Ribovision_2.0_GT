<template>
    <div id="phylo_tree_dropdown">
        <div class="left-sidebar">
            <div id="tree_type" class="btn-group btn-group-toggle" data-toggle="buttons">
                
                <label class="btn btn-outline-dark" style="margin: 0 1% 0 0;width:100%;" for="orthologs" >
                    <input type="radio" id="orthologs" value="orth" v-model="type_tree" v-on:input="cleanTreeOpts()" checked>
                    RiboVision
                </label>
                <!--<label class="btn btn-outline-dark" for="paralogs">
                    <input type="radio" id="paralogs" value="para" v-model="type_tree" v-on:input="cleanTreeOpts()">
                    Paralogs
                </label>-->
                
                    <label class="btn btn-outline-dark" style="margin: 0 0 0 1%;width:50%;" for="upload">
                    <input type="radio" id="upload" value="upload" v-model="type_tree" v-on:input="cleanTreeOpts()">
                    User upload
                </label>
                <!--
                -->
            </div>
            <div id="treeselect" v-if="type_tree=='para'|type_tree=='orth'">
            <treeselect ref="treeselect"
              :load-options="loadOptions"
              v-model="tax_id" 
              v-on:input="loadProteinTypes(tax_id, type_tree)"
              placeholder="Select a phylogenetic group"
              no-children-text="Loading... or no children"
              :multiple="true" 
              :options="options" 
              :flat="true"
              :limit="3"
              :default-expand-level="1"
              >Loading phylogenetic tree <img src="static/img/loading.gif" alt="Loading phylogenetic tree" style="height:25px;">
              </treeselect>
            </div>
            <div v-else>
                Select an alignment file:
                <p><input id="inputUploadFasta" class="btn btn-outline-dark" type = "file" accept=".fasta,.fas,.fa" ref="custom_aln_file" v-on:change="handleFileUpload()"/></p>
                <p><button id="uploadShowFasta" class="btn btn-outline-dark" v-on:click="submitCustomAlignment()" v-if="file&&type_tree=='upload'">Upload the alignment</button></p>
                <p><button id="downloadExampleFasta" class="btn btn-outline-dark" v-on:click="getExampleFile(`static/alignments/EFTU_example.fas`, `PVExampleAlignment.fas`)" v-if="!file&&type_tree=='upload'">Download example alignment</button></p>
            </div>
            <p>
                <select class="btn btn-outline-dark dropdown-toggle" id="select_protein_type" v-if="tax_id" v-model="protein_type_obj" v-on:change="loadData(tax_id, type_tree)">
                    <option :value="null" selected disabled hidden>Select RNA type</option>
                    <option v-for="proteinType in proteinTypes" v-bind:key="proteinType">{{ proteinType }}</option>
                </select>
            </p>
            <p>
                <select class="btn btn-outline-dark dropdown-toggle" id="selectaln" v-if="protein_type_obj && tax_id" v-model="alnobj">
                    <option :value="null" selected disabled hidden>Select an alignment</option>
                    <option v-for="aln in alignments" v-bind:key = aln.txt v-bind:value="{ id: aln.value, text: aln.text }">{{ aln.text }}</option>
                </select>
            </p>
                <!--<span v-if="alnobj&&alnobj!='custom'">Select structure for mapping:</span>-->
                <div v-if="alnobj&&alnobj=='custom'&&file&&type_tree=='upload'">
                    <label for="uploadCustomPDB" id="pdb-upload" class="btn btn-outline-dark">Upload a custom PDB</label>
                    <input id="uploadCustomPDB" class="btn btn-outline-dark" type="file" accept=".pdb" ref="customPDBfile" v-on:change="uploadCustomPDB()"/>
                    <!--OR<br>-->
                </div>
                <span v-if="alnobj">Select/type PDB entry:</span>
                <select class="btn btn-outline-dark dropdown-toggle" id="pdb_input" v-if="alnobj&&alnobj!='custom'" v-model="pdbid">
                    <option :value="null" selected disabled hidden>Select PDB entry</option>
                    <option v-for="pdb in pdbs" v-bind:value="pdb.id">{{pdb.name}}</option>
                </select>
                <p>
                <autocomplete id="pdb_input" isAsync:true :items="pdbs" v-if="alnobj&&alnobj!='custom'" v-model="pdbid"></autocomplete>
                <autocomplete isAsync:true :items="blastPDBresult" v-if="alnobj&&alnobj=='custom'" v-model="pdbid"></autocomplete>
                <!--
                <div id="blastingPDBsMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln&&fetchingPDBwithCustomAln==true">
                    <b>BLASTing first alignment sequence against PDB sequences</b>
                    <img src="static/img/loading.gif" alt="BLASTing available PDBs" style="height:25px;">
                </div>
                <div id="blastedPDBsNoneMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln&&fetchingPDBwithCustomAln=='none'">
                    <b>BLAST didn't match any PDBs!<br>You can still type in any PDB.</b>
                </div>
                <div id="blastedPDBsNoneMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln&&fetchingPDBwithCustomAln=='error'">
                    <b>BLAST error! Please try different alignment or refreshing the page!</b>
                </div>
                <span id="completeBLASTsMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln=='complete'"><b>Completed BLAST for similar PDBs.</b></span>
                <div v-if="hide_chains" id="onFailedChains">Looking for available polymers <img src='static/img/loading.gif' alt='Searching available polymers' style='height:25px;'></div>
                </p>-->
            <p><select multiple class="form-control btn-outline-dark" id="polymerSelect" v-bind:style="{ resize: 'both'}"  v-if="chains&&fasta_data&&pdbid||uploadSession" v-model="chainid" >
                <option :value ="null" selected disabled>Select an RNA chain to visualize</option>
                <option v-for="chain in chains" v-bind:key="chain.value" v-bind:value="chain.value" @click="postStructureData(pdbid, chainid); calculateProteinContacts(pdbid, chainid); populateECODranges(pdbid, chainid); showPDBViewer(pdbid, chainid, chain.entityID);">{{ chain.text }}</option>
            </select></p>
            <!-- 
            -->
            <!-- 
            -->
            <div v-if="structure_mapping">
                <select id="downloadDataBtn" class="btn btn-outline-dark dropdown-toggle" v-model="downloadMapDataOpt" v-if="topology_loaded">
                    <option :value="null" selected disabled>Download mapped data</option>
                    <option value='csv'>As CSV file</option>
                    <option value='pymol'>As PyMOL script</option>
                </select>
            </div>
            <p><div v-if="topology_loaded&&type_tree=='orth'" class="checkbox" id="showRNAcontext">
                <label><input type="checkbox" v-model="checkedRNA" v-on:change="updateMolStarWithRibosome(checkedRNA)">
                    Show ribosomal context in 3D</label>
            </p></div>
            <div v-if="topology_loaded&&!checkedRNA">  
                <div id="customDataSection">
                <p><div class="checkbox">
                        <label><input type="checkbox" v-model="checked_customMap" v-on:change="cleanCustomMap(checked_customMap)">
                        Upload custom mapping data</label>
                        <p><input class="btn btn-outline-dark" id="inputUploadCSV" v-if="checked_customMap" type="file" accept=".csv" ref="custom_csv_file" v-on:change="handleCustomMappingData()"/></p>
                        <p v-if="raiseCustomCSVWarn" v-html="raiseCustomCSVWarn"></p>
                        <p><button class="btn btn-outline-dark" id="downloadExampleCSV" v-if="checked_customMap&&csv_data==null" type="button" v-on:click="getExampleFile(`static/alignments/rv3_example_cusom_mapping.csv`, `PVExampleCustomMapping.csv`)">
                        Download example mapping data
                        </button></p>
                    </div>
                </p></div>
                <div v-if="topology_loaded&&protein_contacts">
                    <p><select multiple class="form-control btn-outline-dark" id="polymerSelect2" v-bind:style="{ resize: 'both'}" v-model="pchainid">
                    <label>Select RNA-protein contacts to view in 3D</label>
                    <option :value ="null" selected disabled></option>
                    <option v-for="chain in protein_chains" v-bind:value="chain.value" v-bind:key="chain.key" @click="showContacts();">{{ chain.text }}</option>
                    </select></p>
                </div>   
                <p><select multiple class="form-control btn-outline-dark" id="polymerSelect3" v-bind:style="{ resize: 'both'}" v-model="modifications" v-if="modified">
                <label>Select modified residues to highlight</label>
                <option :value ="null" selected disabled></option>
                <option v-for="[text, k] of modified_residues.entries()" v-bind:value="text" v-bind:key="k" @click="showModifications();">{{ text }}</option>
                </select></p>
            </div>
        </div>
        <div class="alignment_section">
            <div id="alnif" v-if="alnobj">
                <div id="alnMenu" style="display: flex;">
                    <button id="downloadFastaBtn" class="btn btn-outline-dark" style="margin: 0 1%;" v-if="msavWillMount" type="button" v-on:click="downloadAlignmentData()">
                        Download alignment
                    </button>
                    <select id="cdHITResults" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="cdhitSelectedOpt" v-if="cdHITReport">
                        <option :value="null" selected disabled>See cdhit options</option>
                        <option v-for="prop in cdhitOpts" :value="prop.value" :key="prop.Name">{{ prop.Name }}</option>
                    </select>
                    <select id="downloadAlnImageBtn" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="downloadAlignmentOpt" v-if="msavWillMount">
                        <option :value="null" selected disabled>Download alignment image</option>
                        <option value='full'>Full alignment</option>
                        <option value='visible'>Visible alignment</option>
                    </select>
                    <select id="selectColorMappingProps" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="selected_property" v-if="msavWillMount">
                        <option :value="null" selected disabled>Select data</option>
                        <option value="Select data">Clear data</option>
                        <option v-for="prop in available_properties" :key="prop.Name">{{ prop.Name }}</option>
                    </select>
                    <select id="selectAlnColorScheme" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="colorScheme" v-if="msavWillMount">
                        <option :value="null" selected disabled>Select a colorscheme</option>
                        <option v-for="colorscheme in availColorschemes" :key="colorscheme">{{ colorscheme }}</option>
                    </select>
                </div>
                <div id="alnDiv">Loading alignment <img src="static/img/loading.gif" alt="Loading alignment" style="height:25px;"></div>
            </div>
        </div>
        <div class="warningSection">
            <div id="warningCDHITtruncation" v-if="cdHITReport&&didCDHit_truncate" >
                <b>Warning, your alignment sequences were clustered by cdhit! See dropdown menu above the alignment for options.<br/>
                Original alignment had {{this.cdHITnums[0]}} sequences, which were clustered in {{this.cdHITnums[1]}} groups using threshold of 90% identity.</b>
            </div>
            <div id="warningPoorStructureAln" v-if="poor_structure_map" >
                <b>Warning, poor alignment between the structure and sequences!!!<br/>
                Found {{poor_structure_map}} poorly aligned residues.
                Proceed with caution or try a different structure.</b>
            </div>
        </div>
        <div class="topology_section">
            <span id="topif" v-if="chainid.length>0||customPDBsuccess">
                <div v-if="!topology_loaded">
                    Wait for alignment-structure mapping <img src="static/img/loading.gif" alt="Loading topology viewer" style="height:25px;">
                </div>
                <div id="topview"></div>
            </span>
        </div>
        <div class = "gradient_section" v-if = "topology_loaded">
            <img id = 'gradientSVG' 
                v-for="prop in available_properties" 
                :key="prop.Name"
                v-if = "selected_property == prop.Name"
                :src="prop.url"
        </div>
        <div class="molstar_section">
            <div v-if="PDBparsing==true">
                Parsing PDB structure <img src="static/img/loading.gif" alt="Parsing PDB structure" style="height:25px;">
            </div>
            <div v-if="PDBparsing=='error'">
                Failed to parse the PDB structure!!! Try a different structure.
            </div>
            <span id="molif" v-if="chainid.length>0||customPDBsuccess">
                <div id ="pdbeMolstarView">
                    Loading Molstar Component <img src="static/img/loading.gif" alt="Loading MolStar" style="height:25px;">
                </div>
            </span>
        </div>
        <footer>
            <div id="footerDiv" style="display: flex;"></div>
        </footer>
    </div>
</template>


<script>
  import schemes from './msaColorSchemes/index.js'
  import {ajaxProper} from './ajaxProper.js'
  import {addFooterImages} from './Footer.js'
  import {initialState} from './DropDownTreeVars.js'
  import {filterAvailablePolymers} from './filterRiboChains.js'
  import {parseRNAchains} from './handleRNAchains.js'
  import {generateChainsFromLiteMol} from './handleChainData.js'
  import {colorByMSAColorScheme} from './handleMSAbasedColoring.js'
  import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'
  import {loadAlignmentViewer} from './loadAlignmentViewer.js'
  import {customCSVhandler} from './handleCSVdata.js'
  import {AlnViewer} from './AlignmentViewer.js'
  import {updateProperty} from './handleCSVdata.js'
  import {populatePDBsFromCustomAln} from './populatePDBsFromCustomAln.js'
  import {populateECODranges} from './populateECODranges.js'
  import {postCIFdata} from './postCustomStruct.js'
  import {uploadCustomPDB} from './handleUploadPDB_RNA.js'
  import {loadViewersWithCustomUploadStructure} from './handleViewersWithUploadPDB.js'
  import ReactDOM, { render } from 'react-dom';
  import React, { Component } from "react";
  import Treeselect from '@riophae/vue-treeselect'
  import Autocomplete from './Autocomplete.vue'
  import { intersection } from 'lodash';
  import {downloadPyMOLscript} from './handlePyMOLrequest.js'
  //import {parseRNAchains} from './handleRNAchains.js'
  export default {
      // register the component
      components: { Treeselect, Autocomplete },
      data: function () {
        return initialState();
      },
      
      watch: {
        type_tree: function (type_tree){
            if (this.type_tree == "orth"){
                document.getElementById('tree_type').children[0].click();
            }else if (this.type_tree == "upload"){
                document.getElementById('tree_type').children[1].click();
            }
        },csv_data: function(csv_data){
            customCSVhandler(csv_data);
        },alnobj: function (data){
            if (!data){return;}
            if (data == "custom"){
                this.showAlignment(null, null, vm.type_tree);
            }else{
                this.populatePDBs(data);
                this.showAlignment(data.id, vm.tax_id, vm.type_tree);
            }
        },pdbid: function (pdbid){
            if (!pdbid){return;}
            if (vm.type_tree == "upload"){
                this.getPDBchains(pdbid, null);
            }else{
                this.getPDBchains(pdbid, vm.alnobj.id);
            }
        },unfilteredChains: function(chain_list){
            if (!chain_list){return;}
            if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
            if (this.type_tree != "upload") {
                //parseRNAchains(chain_list);
                filterAvailablePolymers(chain_list, this.alnobj.id, this);
            } else if (vm.blastMAPresult == null){
                let chain_options = []
                for (let i = 0; i < chain_list.length; i++) {
                    let chain_listI = chain_list[i]
                    if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
                    if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
                    if (typeof(chain_listI.source[0]) === "undefined") {continue;}
                    chain_options = pushChainData(chain_options, chain_listI);
                }
                if (chain_options.length === 0) {
                    chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
                }
                vm.chains = chain_options;
                this.hide_chains = null;
            } else {
                let chain_options = [];
                var chainsFromBlast = vm.blastMAPresult.get(this.pdbid);
                for (let i = 0; i < chain_list.length; i++) {
                    let chain_listI = chain_list[i]
                    if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
                    if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
                    if (typeof(chain_listI.source[0]) === "undefined") {continue;}
                    if (!chainsFromBlast){
                        chain_options = pushChainData(chain_options, chain_listI);
                    } else {
                        let intersectedChains = _.intersection(chainsFromBlast, chain_listI["in_chains"]);
                        intersectedChains.forEach(function(chainVal){
                            chain_options.push({
                                text: `${chainVal} ${chain_listI["molecule_name"][0]}`,
                                value: chainVal,
                                sequence: chain_listI["sequence"],
                                entityID: chain_listI["entity_id"],
                                startIndex: chain_listI.source[0].mappings[0].start.residue_number,
                                endIndex: chain_listI.source[0].mappings[0].end.residue_number
                            });
                        });
                    }
                }
                if (chain_options.length === 0) {
                    chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
                }
                vm.chains = chain_options;
                this.hide_chains = null;
            }
        },colorScheme: function (scheme){
            if (window.PVAlnViewer){
                window.PVAlnViewer.setState({colorScheme:scheme});
            }
            if (this.topology_loaded == true){
                colorByMSAColorScheme(scheme, this);
            }
        },postedPDBEntities: function (successPost){
            if (successPost){
                this.showTopologyViewer(this.pdbid, this.chainid, this.fasta_data);
            } else {
                const topview_item = document.getElementById("topview");
                if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading Structure Data ", true)}
            }
        },customPDBsuccess:function(successPost){
            if (successPost){
                this.$nextTick(function(){
                    loadViewersWithCustomUploadStructure();
                })
            }
        },topology_loaded: function(topology_loaded){
            if (window.tempCSVdata!= null && this.topology_loaded){
                vm.csv_data = window.tempCSVdata;
                window.tempCSVdata = null;
                customCSVhandler(vm.csv_data);
            }
            if (this.selected_property){
                this.$nextTick(function(){
                    recolorTopStar(this.selected_property);
                });
            }
        },downloadAlignmentOpt: function(opt){
            if (this.uploadSession){return;}
            else if (opt == 'visible'){
                downloadAlignmentImage();
            } else if (opt == 'full'){
                downloadFullAlignmentImage();
            }
            this.downloadAlignmentOpt = null;
        },downloadMapDataOpt: function(opt){
            if (this.uploadSession){return;}
            else if (opt == 'csv'){
                downloadCSVData();
                this.downloadMapDataOpt = null;
            } else if (opt == 'pymol'){
                downloadPyMOLscript();
                this.downloadMapDataOpt = null;
            }
        },selected_property: function(name){
            if (this.uploadSession){return;}
            if (!name){return;}
            if(this.colorSchemeData){this.colorSchemeData = null;}
            var updatedBarColors = [];
            if(name == "Select data") {
                window.aaFreqs.forEach(function(){
                        updatedBarColors.push("#808080")
                    });
            } else {
                let min = Math.min(...aaPropertyConstants.get(name));
                let max = Math.max(...aaPropertyConstants.get(name));
                let colormapArray = aaColorData.get(name);
                let propData = this.aa_properties.get(name);
                var separatedData = [];
                if (this.aa_properties.has(name)){
                    propData.forEach(function(data, index){
                        separatedData.push([index+1, Number(math.sum(data).toFixed(2))]);
                    })
                } else if (name == 'TwinCons'){
                    separatedData = this.unmappedTWCdata;
                } else {
                    //assume custom data
                    if (this.structure_mapping && window.custom_prop){
                        var customProp = window.custom_prop.get(name);
                        window.aaFreqs.forEach(function(aaFr, alnIx){
                            var strucIx = vm.structure_mapping[alnIx+1];
                            if (strucIx && customProp[strucIx-1]){
                                customProp.forEach(function(customData){
                                    if (customData[0] == strucIx){
                                        separatedData.push([alnIx+1, Number(customData[1])]);
                                    }
                                });
                            } else {
                                separatedData.push([alnIx+1, NaN]);
                            }
                        });
                    }
                }
                const [rgbMap, MappingData] = parsePVData(separatedData, min, max, colormapArray);
                rgbMap.forEach(function (data){
                    updatedBarColors.push(rgbToHex(...data[0]))
                })
            }
            window.barColors = updatedBarColors;
            var alnDiv = document.querySelector('#alnDiv');
            window.msaOptions.colorScheme = this.colorScheme;
            this.aaPos = window.PVAlnViewer.state.aaPos;
            this.seqPos = window.PVAlnViewer.state.seqPos;
            ReactDOM.unmountComponentAtNode(alnDiv);
            this.msavWillMount = null;
            this.$nextTick(function(){
                //const root = ReactDOM.createRoot(alnDiv)
                //root.render(<AlnViewer ref={(PVAlnViewer) => {window.PVAlnViewer = PVAlnViewer}}/>)
                ReactDOM.render(
                    <AlnViewer ref={(PVAlnViewer) => {window.PVAlnViewer = PVAlnViewer}}/>, alnDiv
                );
            });
            if (this.topology_loaded){
                recolorTopStar(name);
            }
        },cdhitSelectedOpt: function(opt){
            if (opt=="untrunc"){
                this.cdhitOpts = this.cdhitOpts.filter(function( obj ) {
                    return obj.value !== 'untrunc';
                });
                if (this.cdhitOpts.filter(e => e.value === 'trunc').length === 0) {
                    this.cdhitOpts.push({Name:'Reload truncated alignment', value:'trunc'});
                }
                cleanupOnNewAlignment(vm, "Loading alignment...");
                vm.showAlignment(null, null, "upload");
                vm.didCDHit_truncate = false;
            }
            if (opt=="trunc"){
                this.cdhitOpts = this.cdhitOpts.filter(function( obj ) {
                    return obj.value !== 'trunc';
                });
                if (this.cdhitOpts.filter(e => e.value === 'untrunc').length === 0) {
                    this.cdhitOpts.push({Name:'Reload original alignment', value:'untrunc'})
                }
                cleanupOnNewAlignment(vm, "Loading alignment...");
                vm.showAlignment(null, null, "upload");
                vm.didCDHit_truncate = true;
            }
            if (opt=="download"){
                let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
                let anchor = document.createElement('a');
                anchor.href = 'data:text;charset=utf-8,' + encodeURIComponent(vm.cdHITReport);
                anchor.target = '_blank';
                anchor.download = `PVcdhitReport-${month}-${date}-${year}.txt`;
                anchor.click();
            }
            if (!this.opt){
                this.cdhitSelectedOpt = null;
            }
        }
    },methods: {
        handleFileUpload(){
            this.file = this.$refs.custom_aln_file.files[0];
            if (this.tax_id != null){this.tax_id = null;}
        },
        submitCustomAlignment(){
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(this);}
            if (vm.fasta_data){
                let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
                vm.fasta_data = cleanFasta;
            }
            let formData = new FormData();
            var fr = new FileReader();
            var uploadedFile = this.file;
            fr.onload = function(){
                if (validateFasta(fr.result)){
                    vm.blastPDBresult = [];
                    vm.blastMAPresult = null;
                    let firstSeq = parseFastaString(fr.result)[1].replace(/-/g,'');
                    vm.populatePDBsFromCustomAln(firstSeq);
                    formData.append('custom_aln_file', uploadedFile)
                    cleanupOnNewAlignment(vm, "Loading alignment...");
                    $.ajax({
                        url: '/custom-aln-data',
                        data: formData,
                        cache: false,
                        contentType: false,
                        processData: false,
                        method: 'POST',
                        type: 'POST', // For jQuery < 1.9
                        success: function(data){
                            if (data == "Success!"){
                                vm.didCDHit_truncate = true;
                            } else {
                                vm.didCDHit_truncate = false;
                            }
                            vm.alnobj = "custom";
                            vm.showAlignment(null, null, "upload");
                        },
                        error: function(error) {
                            alert(`${error.responseText}`);
                        }
                    });
                }
            };
            fr.readAsText(this.file)
        },
        cleanTreeOpts() {
            if (this.uploadSession){return;}
            Object.assign(vm.$data, initialState());
            this.type_tree="upload";
            this.topology_loaded=false;
            this.schemesMgr = new schemes();
            //cleanupOnNewAlignment(this, "Select new alignment!");
            //[this.options, this.tax_id, this.alnobj, this.chainid] = [null, null, null, null];
            //var topview = document.getElementById("PdbeTopViewer");
            //var molstar = document.querySelector(".msp-plugin");
            //if (molstar) {molstar.textContent = null}
            //if (topview) {topview.textContent = null}
        }, loadOptions({ action, callback }) {
            if (this.type_tree == "orth"){
                if (action === "LOAD_CHILDREN_OPTIONS") {
                    action = "";
                    callback();
                //     When they figure out LOAD_CHILDREN_OPTIONS with async search
                //     ajax(`/alignments/showTaxonomy-api/${parentNode.id}`).then(data => {
                //         let fetched_data = [data]
                //         parentNode.children = fetched_data[0].children
                //         callback()
                //     }).catch(error => {
                //         parentNode.children = []
                //         console.log(error)
                //         callback(new Error(`Failed to load options: network error: ${error}`))
                //     })
                };
                if (action === "LOAD_ROOT_OPTIONS") {
                    ajax(`/alignments/showTaxonomy-api/0`).then(data => {
                        data.isDisabled = true;
                        this.options = [data];
                        callback();
                    }).catch(error => {
                        console.log(error)
                    })
                };
                if (action === "LOAD_ROOT_OPTIONS") {
                    ajax('/alignments/showTaxonomy').then(data => {
                        if (this.type_tree == "orth"){
                            this.options = null;
                            data.isDisabled = true;
                            this.options = [data];
                            callback();
                        }
                    }).catch(error => {
                        console.log(error)
                    })
                };
            }
            if (this.type_tree == "para"){
                loadParaOptions(action, callback, this);
            }
        }, 
        loadProteinTypes (tax_id, type_tree) {
            if (type_tree == "orth") {
                this.alignments = null;
                this.proteinTypes = null;
                var url = '/proteinTypes';
                // let taxIDs = tax_id;
                let taxIDs = '';
                vm.tax_id.forEach(element => {
                    taxIDs += element + ",";
                });
                ajax(url, {taxIDs}).then(data => {
                    let results = data["results"];
                    let numSublists = results.length;
                    let intersection = results[0];
                    for (let i = 1; i < numSublists; i++) {
                        intersection = intersection.filter(value => results[i].includes(value));
                    }
                    vm.proteinTypes = intersection;
                    // let proteinTypes = data["proteinTypesList"];
                    // vm.proteinTypes = proteinTypes;
                });
            }
        },
        showModifications() {
            viewerInstanceTop.viewInstance.uiTemplateService.colorMapModifications();  
            showModificationsAndContactsHelper("" + this.entityID);
        },
        showContacts() {
            viewerInstanceTop.viewInstance.uiTemplateService.colorMapContacts();  
            showModificationsAndContactsHelper("" + this.entityID);
        },
        getR2DT(sequence) {
            //console.log("getR2DT.vue", sequence);
            //var url = `r2dt/${sequence}`
            
            //ajax(url).then((returnedObject) => {
                    //let returnedData = returnedObject["RNA_2D_json"]; 
                    //console.log("r2dt_json", returnedData);   
                    //viewerInstanceTop.viewInstance.uiTemplateService.render(returnedData);  
        //})
        
        
        },
        loadData (value, type_tree) {
            if (this.uploadSession){return;}
            if (type_tree == "upload"){this.tax_id = null; return;}
            if (value.length == 0){this.tax_id = null; return;}
            cleanupOnNewAlignment(this, "Select new alignment!");
            vm.chainid = [];
            if (this.alnobj != null) {this.alnobj = null;}
            if (type_tree == "orth"){
                this.alignments = null;
                /*
                var url = '/desire-api/taxonomic-groups/?format=json&taxgroup_id__in=' + value
                ajax(url).then(data => {
                    loadOrthAlns(data, this);
                });
                */
                var url = '/getAlignmentsFilterByProteinTypeAndTaxIds';
                let selectedProteinType = vm.protein_type_obj;
                let taxIDs = '';
                vm.tax_id.forEach(element => {
                    taxIDs += element + ",";
                });
                ajax(url, {selectedProteinType, taxIDs}).then(data => {
                // ajax(url).then(data => {
                    let results = data["results"];
                    let intersection = results[0];
                    let numSublists = results.length;
                    for (let sublistIndex = 1; sublistIndex < numSublists; sublistIndex++) {
                        intersection = intersection.filter(intersectionSublistEntry => {
                            let
                                sublist = results[sublistIndex],
                                sublistLength = sublist.length,
                                inclusionFlag = false;
                            for (let sublistEntryIndex = 0; sublistEntryIndex < sublistLength; sublistEntryIndex++) {
                                if (sublist[sublistEntryIndex][0] == intersectionSublistEntry[0]) {
                                    inclusionFlag = true;
                                    break;
                                }
                            }
                            return inclusionFlag;
                        });
                    }
                    vm.alignments = []
                    intersection.forEach(intersectionEntry => {
                        let alignment = new Object();
                        alignment.text = intersectionEntry[0];
                        alignment.value = intersectionEntry[1];
                        vm.alignments.push(alignment);
                    });
                    // let alignmentNamesAndPrimaryKeys = data["alignmentNamesAndPrimaryKeys"];
                    // alignmentNamesAndPrimaryKeys.forEach(alignmentNamesAndPrimaryKey => {
                    //     let alignment = new Object();
                    //     alignment.text = alignmentNamesAndPrimaryKey[0];
                    //     alignment.value = alignmentNamesAndPrimaryKey[1];
                    //     vm.alignments.push(alignment);
                    // });
                    
                    // loadOrthAlns(data, this);
                });
                // var url = '/desire-api/taxonomic-groups/?format=json&taxgroup_id__in=' + value;
                // ajax(url).then(data => {
                //     loadOrthAlns(data, this);
                // });

            }
            if (type_tree == "para"){
                loadParaAlns (value, this)
            }
        }, getPDBchains(pdbid, aln_id) {
            if (this.uploadSession){return;}
            if (pdbid.length === 4) {
                if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(this);}
                this.unfilteredChains = null;
                this.PDBparsing = false;
                loadAlignmentViewer(vm.fasta_data);
                this.chains = null;
                this.chainid = [];
                this.hide_chains = true;
                generateChainsFromLiteMol(`https://coords.litemol.org/${pdbid.toLowerCase()}/assembly?id=1&lowPrecisionCoords=1&encoding=BCIF`, "unfilteredChains");
                ajax(`https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/${pdbid.toLowerCase()}`).then(struc_data => {
                    if(vm.unfilteredChains){return;}
                    vm.unfilteredChains = struc_data[pdbid.toLowerCase()];
                }).catch(error => {
                    console.log(error);
                    var elt = document.querySelector("#onFailedChains");
                    if (error.status == 404){
                        elt.innerHTML  = "Couldn't find this PDB on EBI!<br/>Try a different PDB ID."
                    } else if (error.status == 0){
                        elt.innerHTML  = "It looks like PDBe is down! Running alternative chain parser..."
                    } else {
                        elt.innerHTML  = "Problem with parsing the chains! Try a different PDB ID."
                    }
                })
            }
        },
        showAlignment(aln_id, taxid, type_tree) {
            cleanupOnNewAlignment(this, "Loading alignment...");
            this.chainid = [];
            if (type_tree == "orth"){
                var url = `/ortholog-aln-api/${aln_id}/${taxid}`}
            if (type_tree == "para"){
                var url = '/paralog-aln-api/'+aln_id.split(',')[1]}
            if (type_tree == "upload" && !this.uploadSession&& this.cdhitSelectedOpt != "untrunc"){
                var url = '/custom-aln-data'
            }
            if (type_tree == "upload" && !this.uploadSession && this.cdhitSelectedOpt == "untrunc"){
                var url = '/custom-aln-data-nocdhit'
            }
            if (this.uploadSession){
                this.$nextTick(function(){
                    loadAlignmentViewer (vm.fasta_data);
                });
                return;
            }
            ajax(url).then(fasta => {
                if (fasta['TwinCons']){
                    this.custom_aln_twc_flag = fasta['TwinCons'];
                    fetchTWCdata(fasta['Alignment']);
                }
                this.cdHITReport = fasta["cdHitReport"]
                if (this.cdHITReport){
                    let cdNums = this.cdHITReport.split(/comparing sequences from.*\n/)[1].split(/\n/)[1].split(/ +/);
                    this.cdHITnums = [cdNums[1], cdNums[3]];
                }
                this.fastaSeqNames = fasta['Sequence names'];
                window.aaFreqs = fasta['AA frequencies'];
                var barColors = Array(aaFreqs.length).fill('#808080');
                window.barColors = barColors;
                this.fasta_data = fasta['Alignment'];
                this.aa_properties = calculateFrequencyData(fasta['AA frequencies']);
                loadAlignmentViewer (fasta['Alignment']);
            })
        }, showTopologyViewer (pdbid, chainid, fasta){
            this.topology_loaded = false;
            if (chainid.length > 1){this.chainid = chainid[0];}
            const topview_item = document.getElementById("topview");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "")}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            let temp = this.chains.filter(obj => {
                return obj["value"] == chainid;
            })[0];
            let ebi_sequence = temp["sequence"];
            let startIndex = temp["startIndex"];
            let stopIndex = temp["endIndex"];
            let struc_id = `${pdbid}-${temp["entityID"]}-${chainid}`
            if (!this.uploadSession){
                getStructMappingAndTWC (fasta, struc_id, startIndex, stopIndex, ebi_sequence, this);
            }
            loadAlignmentViewer (vm.fasta_data);
            var rna_url = `https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/${pdbid}/chain/${chainid}`
            ajax(rna_url).then(data => {
                if(vm.topology_loaded){return;}
                var entityid = data[pdblower]["molecules"][0].entity_id;
                var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}"></pdb-rna-viewer>`
                document.getElementById('topview').innerHTML = topology_viewer;
                window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
            });
        }, showPDBViewer(pdbid, chainid, entityid){
            showPDBHelper(pdbid, chainid, entityid)
        }, populateECODranges(pdbid, chainid) {
            populateECODranges(pdbid, chainid);
        }, calculateProteinContacts(pdbid, chainid) {
            vm.mapped_aa_contacts_mods = new Map();
            var url = `protein-contacts/${pdbid}/${chainid}`
            ajax(url).then(data => {
                calculateModifiedResidues(pdbid, chainid, this.entityID)
                if(data) {
                    vm.protein_contacts = data;
                    var newContactMap;
                    var filtered_chains = vm.protein_chains.filter(e => e.value in data);

                    vm.protein_chains = filtered_chains;
                    var i = 1.0;
                    var colorMap = new Map();
                    vm.selectSections_proteins = new Map();
                    vm.mapped_aa_contacts_mods.set("Protein Contacts", [])
                    for (var val in vm.protein_contacts) {
                        vm.selectSections_proteins.set(val, [])
                        var color = interpolateLinearly(i/filtered_chains.length, aaColorData.get("Protein contacts")[0])
                        var rgbColor = "rgb(" + color[0][0] + "," + color[0][1] + "," + color[0][2] + ")";
                        colorMap.set(val, rgbColor);
                        //newContactMap.set(vm.protein_contacts, aaColorData.get("Shannon entropy")[0][1]
                        i = i+1;
                        for (var j in vm.protein_contacts[val]) {
                            vm.selectSections_proteins.get(val).push({
                                entity_id: "" + this.entityID,
                                residue_number: vm.protein_contacts[val][j],
                                color: color[1],
                                sideChain: false,
                            });
                            var protein_name = vm.protein_chains.filter(e => e.value == val)[0].text
                            vm.mapped_aa_contacts_mods.get("Protein Contacts").push([vm.protein_contacts[val][j], protein_name])
                        }
                    }                    
                    vm.proteinColorMap = colorMap;
                }
                }).catch(error => {
                    console.log(error)
                })  
        },
        postStructureData(pdbid, chainid) {
            const topview_item = document.getElementById("topview");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading Structure Data ", true)}
            let tempEntities = this.chains.filter(obj => {
                return obj["value"] == chainid;
            });
            let entities = [];
            tempEntities.forEach(function(ent){
                entities.push({ entityID: ent["entityID"], chainID: ent["value"] })
            });
            this.entityID = tempEntities[0]["entityID"];
            postCIFdata(pdbid, entities);
        },downloadAlignmentImage() {
            downloadAlignmentImage(document.querySelector('#alnDiv'));
        },downloadAlignmentData() {
            downloadAlignmentData(vm.fasta_data);
        },downloadCSVData() {
            downloadPyMOLscript();
            downloadCSVData();
        },getExampleFile(url, name) {
            getExampleFile(url, name);
        },cleanFilter(checked_filter, masking_range){
            cleanFilter(checked_filter, masking_range);
        },handleMaskingRanges(mask_range){
            handleMaskingRanges(mask_range);
        },cleanSelection(checked_selection, filter_range){
            cleanSelection(checked_selection, filter_range);
        },cleanCustomMap(checked_customMap){
            cleanCustomMap(checked_customMap);
        },handleCustomMappingData(){
            handleCustomMappingData();
        },populatePDBs(alndata){
            populatePDBs(alndata);
        },populatePDBsFromCustomAln(firstSeq){
            if (this.guideOff){
                populatePDBsFromCustomAln(firstSeq);
            }
        },listSecondaryStructures() {
            listSecondaryStructures();
        },downloadFreqsData(){
            let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
            let anchor = document.createElement('a');
            anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(vm.freqCSV);
            anchor.target = '_blank';
            anchor.download = `PVfreqData-${month}-${date}-${year}.csv`;
            anchor.click();
        },flushDjangoSession(){
            ajaxProper({
                    url: `/flush-session`,
                    type: `GET`,
                    dataType: `text`,
                }).then(response => {
                if (response == 'Success!'){
                    console.log("Session flushed successfully!") 
                }
            })
        }, uploadCustomPDB(){
            uploadCustomPDB();
        }, updateMolStarWithRibosome(checkRibo){
            if(checkRibo&&viewerInstance&&this.pdbid&&this.entityID){
                this.completeRiboContext = false;
                viewerInstance.visual.update({
                    moleculeId: this.pdbid, 
                    assemblyId: '1',
                    bgColor: {r:255,g:255,b:255},
                });
                viewerInstance.events.loadComplete.subscribe(function (e) {
                    let prom = viewerInstance.visual.select({ 
                        data: [{entity_id: `${vm.entityID}` }], 
                        nonSelectedColor: {r:180, g:180, b:180} 
                    });
                    prom.then(function(v){
                        viewerInstance.visual.focus([{ entity_id: `${vm.entityID}` }]);
                        if(viewerInstanceTop&&vm.selected_property){
                            viewerInstanceTop.pluginInstance.displayDomain();
                        }
                    })
                });
            }
            if (!checkRibo&&viewerInstance&&this.pdbid&&this.entityID){
                this.showPDBViewer(this.pdbid, this.chainid[0], this.entityID);
                viewerInstance.events.loadComplete.subscribe(function (e) {
                    if(viewerInstanceTop&&vm.selected_property){
                        viewerInstanceTop.pluginInstance.displayDomain();
                    }
                });
            }
        }
    }, 
    mounted() {
        addFooterImages("footerDiv");
        this.schemesMgr = new schemes();
    },
    created() {
        $(window).bind('beforeunload', function(){
            vm.flushDjangoSession();
        });
    },
}
</script>