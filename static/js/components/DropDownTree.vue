<template>
    <div id="phylo_tree_dropdown">
        <div class="left-sidebar">
            <div id="tree_type" class="btn-group btn-group-toggle" data-toggle="buttons">
                
                <label class="btn btn-outline-dark" style="margin: 0 1% 0 0;width:50%;" for="orthologs" >
                    <input type="radio" id="orthologs" value="orth" v-model="type_tree" v-on:input="cleanTreeOpts()" checked>
                    DESIRE
                </label>
                <!--<label class="btn btn-outline-dark" for="paralogs">
                    <input type="radio" id="paralogs" value="para" v-model="type_tree" v-on:input="cleanTreeOpts()">
                    Paralogs
                </label>-->
                <label class="btn btn-outline-dark" style="margin: 0 0 0 1%;width:50%;" for="upload">
                    <input type="radio" id="upload" value="upload" v-model="type_tree" v-on:input="cleanTreeOpts()">
                    User upload
                </label>
            </div>
            <div id="treeselect" v-if="type_tree=='para'|type_tree=='orth'">
            <treeselect ref="treeselect"
              :load-options="loadOptions"
              v-model="tax_id" 
              v-on:input="loadData(tax_id, type_tree)"
              placeholder="Select a phylogenetic group"
              no-children-text="Loading... or no children"
              :multiple="true" 
              :options="options" 
              :flat="true"
              :limit="3"
              :default-expand-level="1"
              >Loading phylogenetic tree <img src="static/img/loading.gif" alt="BLASTing available PDBs" style="height:25px;">
              </treeselect>
            </div>
            <div v-else>
                Select an alignment file:
                <p><input id="inputUploadFasta" class="btn btn-outline-dark" type = "file" accept=".fasta,.fas,.fa" ref="custom_aln_file" v-on:change="handleFileUpload()"/></p>
                <p><button id="uploadShowFasta" class="btn btn-outline-dark" v-on:click="submitCustomAlignment()" v-if="file&&type_tree=='upload'">Upload the alignment</button></p>
                <p><button id="downloadExampleFasta" class="btn btn-outline-dark" v-on:click="getExampleFile(`static/alignments/EFTU_example.fas`, `rv3ExampleAlignment.fas`)" v-if="!file&&type_tree=='upload'">Download example alignment</button></p>
            </div>
            <p>
                <select class="btn btn-outline-dark dropdown-toggle" id="selectaln" v-if="tax_id" v-model="alnobj">
                    <option v-if="tax_id" :value="null" selected disabled hidden>Select an alignment</option>
                    <option v-if="tax_id" v-for="aln in alignments" v-bind:value="{ id: aln.value, text: aln.text }">{{ aln.text }}</option>
                </select>
            </p>

                <span v-if="alnobj&&alnobj!='custom'">Select PDB for mapping:</span>
                <span v-if="alnobj&&alnobj=='custom'">Type PDB for mapping:</span>

            <p>
                <select class="btn btn-outline-dark dropdown-toggle" id="pdb_input" v-if="alnobj&&alnobj!='custom'" v-model="pdbid">
                    <option :value="null" selected disabled hidden>Select structure</option>
                    <option v-for="pdb in pdbs" v-bind:value="pdb.id">{{pdb.name}}</option>
                </select>
                <autocomplete isAsync:true :items="blastPDBresult" v-if="alnobj&&alnobj=='custom'" v-model="pdbid"></autocomplete>
                <div id="blastingPDBsMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln&&fetchingPDBwithCustomAln!='complete'">
                    <b>BLASTing available PDBs</b>
                    <img src="static/img/loading.gif" alt="BLASTing available PDBs" style="height:25px;">
                </div>
                <span id="completeBLASTsMSG" v-if="alnobj&&alnobj=='custom'&&fetchingPDBwithCustomAln=='complete'"><b>Completed BLAST for similar PDBs.</b></span>
                <div v-if="hide_chains" id="onFailedChains">Looking for available polymers <img src='static/img/loading.gif' alt='Searching available polymers' style='height:25px;'></div>
            </p>
            <p><select multiple class="form-control btn-outline-dark" id="polymerSelect" v-bind:style="{ resize: 'both'}"  v-if="chains&&fasta_data&&pdbid||uploadSession" v-model="chainid" >
                <option :value ="null" selected disabled>Select a polymer</option>
                <option v-for="chain in chains" v-bind:value="chain.value" @click="postStructureData(pdbid, chainid); populateECODranges(); showPDBViewer(pdbid, chainid, chain.entityID); ">{{ chain.text }}</option>
            </select></p>
            <div v-if="structure_mapping">
                <select id="downloadDataBtn" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="downloadMapDataOpt" v-if="topology_loaded">
                    <option :value="null" selected disabled>Download mapped data</option>
                    <option value='csv'>As CSV file</option>
                    <option value='pymol'>As PyMOL script</option>
                </select>
            </div>
            <div v-if="topology_loaded">
                <div id="maskingSection"><p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_filter" v-on:change="cleanFilter(checked_filter, masking_range)">
                        Mask/Unmask 2D and 3D residues</label>
                    </div>
                    <span v-if="checked_filter"><b>Input multiple</b> residue ranges to <b>show</b>, separated by semicolon. <br> For example: 1-80;91-111;</span>
                    <input class="input-group-text" v-if="checked_filter" v-model="masking_range" v-on:input="handleMaskingRanges(masking_range)">
                </p></div>
                <p v-if="correct_mask!=true&&masking_range!=null">Incorrect range syntax!</p>
                <div id="filterSection"><p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_selection" v-on:change="cleanSelection(checked_selection, filter_range)">
                        Cut/Uncut 2D and 3D structures</label>
                    </div>
                    <span v-if="checked_selection"><b>Input single</b> residue range to <b>show</b>, ending with semicolon. <br> For example: 1-80;</span>
                    <input class="input-group-text" v-if="checked_selection" v-model="filter_range" v-on:input="handleFilterRange(filter_range)">
                </p></div>
                <div id="domainSelectionSection"><p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_domain" v-on:change="cleanSelection(checked_domain, true)">
                        Select a domain to show</label>
                    </div>
                </p>
                <p><select multiple class="form-control btn-outline-dark" id="domainSelect" v-bind:style="{ resize: 'both'}"  v-if="domain_list&&checked_domain">
                    <option v-for="domain in domain_list" v-bind:value="selected_domain" @click="handleFilterRange(domain.range)">{{ domain.name + ' ' + domain.range }}</option>
                </select></p>
                </div>
                <div id="customDataSection">
                <p><div class="checkbox">
                        <label><input type="checkbox" v-model="checked_customMap" v-on:change="cleanCustomMap(checked_customMap)">
                        Upload custom mapping data</label>
                        <p><input class="btn btn-outline-dark" id="inputUploadCSV" v-if="checked_customMap" type="file" accept=".csv" ref="custom_csv_file" v-on:change="handleCustomMappingData()"/></p>
                        <p v-if="raiseCustomCSVWarn" v-html="raiseCustomCSVWarn"></p>
                        <p><button class="btn btn-outline-dark" id="downloadExampleCSV" v-if="checked_customMap&&csv_data==null" type="button" v-on:click="getExampleFile(`static/alignments/rv3_example_cusom_mapping.csv`, `rv3ExampleCusomMapping.csv`)">
                        Download example mapping data
                        </button></p>
                    </div>
                </p></div>
            </div>
            <p><div v-if="alnobj" class="checkbox" id="showFrequencies">
                <label><input type="checkbox" v-model="checked_propensities" v-on:change="handlePropensities(checked_propensities)">
                Show amino-acid frequencies</label>
                <select class="btn btn-outline-dark dropdown-toggle" id="propensitiesSubstructure" v-if="checked_propensities&&structure_mapping" v-model="property" v-on:change="getPropensities(property); handlePropensities(checked_propensities)">
                    <option :value="null" selected disabled hidden>Select secondary structure</option>
                    <option :value="0">All residues</option>
                    <option v-for="substructure in substructures" v-bind:value="{ id: substructure.value, text: substructure.text, indices: substructure.indices }">{{ substructure.text }}</option>
                </select>
            </div></p>
        </div>
        <div class="alignment_section">
            <div id="alnif" v-if="alnobj">
                <div id="alnMenu" style="display: flex;">
                    <button id="downloadFastaBtn" class="btn btn-outline-dark" style="margin: 0 1%;" v-if="msavWillMount" type="button" v-on:click="downloadAlignmentData()">
                        Download alignment
                    </button>
                    <select id="downloadAlnImageBtn" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="downloadAlignmentOpt" v-if="msavWillMount">
                        <option :value="null" selected disabled>Download alignment image</option>
                        <option value='full'>Full alignment</option>
                        <option value='visible'>Visible alignment</option>
                    </select>
                    <select id="selectColorMappingProps" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="selected_property" v-if="msavWillMount">
                        <option :value="null" selected disabled>Annotation</option>
                        <option v-for="prop in available_properties" >{{ prop.Name }}</option>
                    </select>
                    <select id="selectAlnColorScheme" class="btn btn-outline-dark dropdown-toggle" style="margin: 0 1%;" v-model="colorScheme" v-if="msavWillMount">
                        <option :value="null" selected disabled>Select a colorscheme</option>
                        <option v-for="colorscheme in availColorschemes" >{{ colorscheme }}</option>
                    </select>
                </div>
                <div id="alnDiv">Loading alignment <img src="static/img/loading.gif" alt="Loading alignment" style="height:25px;"></div>
            </div>
        </div>
        <div class="warningSection">
            <div id="warningPoorStructureAln" v-if="poor_structure_map" >
                <b>Warning, poor alignment between the structure and sequences!!!<br/>
                Found {{poor_structure_map}} poorly aligned residues.
                Proceed with caution or try a different structure.</b>
            </div>
        </div>
        <div class="topology_section">
            <span id="topif" v-if="chainid.length>0">
                <div v-if="!topology_loaded">
                    Loading topology viewer and conservation data <img src="static/img/loading.gif" alt="Loading topology viewer" style="height:25px;">
                </div>
                <div id="topview"></div>
            </span>
        </div>
        <div class = "gradient_section" v-if = "topology_loaded">
            <img id = 'gradientSVG' 
                v-for="prop in available_properties" 
                v-if = "selected_property == prop.Name"
                :src="prop.url"
            ><!--
            <object id="gradientSVG"
                v-for="prop in available_properties" 
                v-if = "selected_property == prop.Name"
                :data="prop.url" type="image/svg+xml">
            </object>-->
        </div>
        <div class="molstar_section">
            <span id="molif" v-if="chainid.length>0">
                <div id ="pdbeMolstarView">
                    Loading Molstar Component <img src="static/img/loading.gif" alt="Loading MolStar" style="height:25px;">
                </div>
            </span>
        </div>
        <div class = "propensity_section">
            <span id="propif" v-if="checked_propensities">
                <div id = "total"></div>
            </span>
        </div>
        <footer>
            <div id="footerDiv" style="display: flex;"></div>
        </footer>
    </div>
</template>


<script>
  import {addFooterImages} from './Footer.js'
  import {initialState} from './DropDownTreeVars.js'
  import {AlnViewer} from './AlignmentViewer.js'
  import {customCSVhandler} from './handleCSVdata.js'
  import {updateProperty} from './handleCSVdata.js'
  import {populatePDBsFromCustomAln} from './populatePDBsFromCustomAln.js'
  import ReactDOM, { render } from 'react-dom';
  import React, { Component } from "react";
  import Treeselect from '@riophae/vue-treeselect'
  import Autocomplete from './Autocomplete.vue'
  import { intersection } from 'lodash';
  import {downloadPyMOLscript} from './handlePyMOLrequest.js'
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
        },colorScheme: function (scheme){
            if (window.PVAlnViewer){
                window.PVAlnViewer.setState({colorScheme:scheme});
            }
        },topology_loaded: function(topology_loaded){
            if (window.tempCSVdata!= null && this.topology_loaded){
                vm.csv_data = window.tempCSVdata;
                window.tempCSVdata = null;
                customCSVhandler(vm.csv_data);
            }
            if (this.selected_property){
                recolorTopStar(this.selected_property);
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
            if (!name){return;}
            if(!aaPropertyConstants.has(name)){return;}
            let min = Math.min(...aaPropertyConstants.get(name));
            let max = Math.max(...aaPropertyConstants.get(name));
            let colormapArray = aaColorData.get(name);
            let propData = this.aa_properties.get(name);
            var separatedData = [];
            var updatedBarColors = [];
            if (this.aa_properties.has(name)){
                propData.forEach(function(data, index){
                    separatedData.push([index+1, Number(math.sum(data).toFixed(2))]);
                })
            } else if (name == 'TwinCons'){
                separatedData = this.unmappedTWCdata;
            } else {
                //assume custom data
                if (this.structure_mapping && window.custom_prop){
                    let invertedMap = _.invert(this.structure_mapping);
                    window.custom_prop.get(name).forEach(function(data){
                        separatedData.push([Number(invertedMap[data[0]]), data[1]]);
                    })
                }
            }
            const [rgbMap, MappingData] = parsePVData(separatedData, min, max, colormapArray);
            rgbMap.forEach(function (data){
                updatedBarColors.push(rgbToHex(...data[0]))
            })
            window.barColors = updatedBarColors;
            var alnDiv = document.querySelector('#alnDiv');
            window.msaOptions.colorScheme = this.colorScheme;
            this.aaPos = window.PVAlnViewer.state.aaPos;
            this.seqPos = window.PVAlnViewer.state.seqPos;
            ReactDOM.unmountComponentAtNode(alnDiv);
            this.msavWillMount = null;
            this.$nextTick(function(){
                ReactDOM.render(
                    <AlnViewer ref={(PVAlnViewer) => {window.PVAlnViewer = PVAlnViewer}}/>, alnDiv
                );
            });
            if (this.topology_loaded){
                recolorTopStar(name);
            }
        }
    },methods: {
        handleFileUpload(){
            this.file = this.$refs.custom_aln_file.files[0];
            if (this.tax_id != null){this.tax_id = null;}
        },
        submitCustomAlignment(){
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
                    $.ajax({
                        url: '/custom-aln-data',
                        data: formData,
                        cache: false,
                        contentType: false,
                        processData: false,
                        method: 'POST',
                        type: 'POST', // For jQuery < 1.9
                        success: function(data){
                            cleanupOnNewAlignment(vm, "Loading alignment...");
                            vm.alnobj = "custom";
                            vm.showAlignment(null, null, "upload");
                        },
                        error: function(error) {
                            alert(`${error.responseText}`);
                        }
                    });
                }else{
                    alert("Check the fasta format of the uploaded file!")
                }
            };
            fr.readAsText(this.file)
        },
        cleanTreeOpts() {
            if (this.uploadSession){return;}
            Object.assign(vm.$data, initialState());
            this.type_tree="upload";
            this.topology_loaded=false;
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
        }, loadData (value, type_tree) {
            if (this.uploadSession){return;}
            if (type_tree == "upload"){this.tax_id = null; return;}
            if (value.length == 0){this.tax_id = null; return;}
            cleanupOnNewAlignment(this, "Select new alignment!");
            vm.chainid = [];
            if (this.alnobj != null) {this.alnobj = null;}
            if (type_tree == "orth"){
                this.alignments = null;
                var url = '/desire-api/taxonomic-groups/?format=json&taxgroup_id__in=' + value
                ajax(url).then(data => {
                    loadOrthAlns(data, this);
                });
            }
            if (type_tree == "para"){
                loadParaAlns (value, this)
            }
        }, getPDBchains(pdbid, aln_id) {
            if (pdbid.length === 4) {
                if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(this);}
                this.chains = null;
                this.chainid = [];
                this.hide_chains = true;
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase()).then(struc_data => {
                    var chain_list = struc_data[pdbid.toLowerCase()];
                    if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
                    if (this.type_tree != "upload") {
                        filterAvailablePolymers(chain_list, aln_id, vm);
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
                        var chainsFromBlast = vm.blastMAPresult.get(pdbid);
                        for (let i = 0; i < chain_list.length; i++) {
                            let chain_listI = chain_list[i]
                            if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
                            if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
                            if (typeof(chain_listI.source[0]) === "undefined") {continue;}
                            let intersectedChains = _.intersection(chainsFromBlast, chain_listI["in_chains"]);
                            intersectedChains.forEach(function(chainVal){
                                chain_options.push({
                                    text: `${chainVal} ${chain_listI["molecule_name"][0]}`,
                                    value: chainVal,
                                    sequence: chain_listI["sequence"],
                                    entityID: chain_listI["entity_id"],
                                    startIndex: chain_listI.source[0].mappings[0].start.residue_number
                                });
                            });
                        }
                        if (chain_options.length === 0) {
                            chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
                        }
                        vm.chains = chain_options;
                        this.hide_chains = null;
                    }
                }).catch(error => {
                    var elt = document.querySelector("#onFailedChains");
                    this.pdbid = null;
                    if (error.status == 404){
                        elt.innerHTML  = "Couldn't find this PDB ID!<br/>Try a different PDB ID."
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
            if (type_tree == "upload"){
                var url = '/custom-aln-data'}
            ajax(url).then(fasta => {
                if (fasta['TwinCons']){
                    this.custom_aln_twc_flag = fasta['TwinCons'];
                    fetchTWCdata(fasta['Alignment']);
                }
                this.fastaSeqNames = fasta['Sequence names'];
                window.aaFreqs = fasta['AA frequencies'];
                var barColors = Array(aaFreqs.length).fill('#BA20B7');
                window.barColors = barColors;
                this.fasta_data = fasta['Alignment'];
                this.aa_properties = calculateFrequencyData(fasta['AA frequencies']);
            })
        },  
         showTopologyViewer (pdbid, chainid, fasta){
            this.topology_loaded = false;
            window.filterRange = "-10000,10000";
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(this);}
            if (chainid.length > 1){this.chainid = chainid[0];}
            const topview_item = document.getElementById("topview");
            const molstar_item = document.getElementById("pdbeMolstarView");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "")}
            if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Loading Molstar Component ", true)}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            let temp = this.chains.filter(obj => {
                return obj["value"] == chainid;
            })[0];
            let ebi_sequence = temp["sequence"];
            let startIndex = temp["startIndex"];
            ajax('/mapSeqAln/', {fasta, ebi_sequence, startIndex}).then(struct_mapping=>{
                this.structure_mapping = struct_mapping;
                if (struct_mapping['BadMappingPositions']){this.poor_structure_map = struct_mapping['BadMappingPositions'];}
                var mapped_aa_properties = mapAAProps(this.aa_properties, struct_mapping);
                if (((this.tax_id != null && this.tax_id.length == 2) || (this.custom_aln_twc_flag != null && this.custom_aln_twc_flag == true) || (this.type_tree == 'para'))) {
                    if (vm.unmappedTWCdata) {
                        mapTWCdata(vm.structure_mapping, vm.unmappedTWCdata, mapped_aa_properties);
                    }
                }
                window.mapped_aa_properties = mapped_aa_properties;
                topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
                topviewer.pluginInstance.createDomainDropdown();
            }).catch(error => {
                //var topview = document.querySelector('#topview');
                console.log(error);
                //this.topology_loaded = 'error';
                //topview.innerHTML = "Failed to load the viewer!<br>Try another structure."
            });
            var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
            ajax(topology_url).then(data => {
                var entityid = Object.keys(data[pdblower])[0];
                vm.coil_residues = filterCoilResidues(data[pdblower][entityid][chainid]["coils"])
                vm.helix_residues = filterCoilResidues(data[pdblower][entityid][chainid]["helices"])
                vm.strand_residues = filterCoilResidues(data[pdblower][entityid][chainid]["strands"])
                listSecondaryStructures();
                var mapping = [];
                var range_string = minIndex.concat("-").concat(maxIndex);
                let ebiMappingURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/'+pdbid;
                ajax(ebiMappingURL).then(data=>{
                    var result = [];
                    customFilter(data, result, "chain_id", chainid);
                    result = result[0];
                    
                    if(result != null) {
                        var pdb_start = parseInt(result["start"]["residue_number"]);
                        var pdb_end = parseInt(result["end"]["residue_number"]);
                        var uniprot_start = parseInt(result["unp_start"]);
                        for (let residue_number_str in range_string.split("-")){
                            var residue_number = parseInt(residue_number_str);
                            if(residue_number >= pdb_start && residue_number <= pdb_end){
                                let offset = uniprot_start - pdb_start;
                                mapping.push(residue_number - offset);
                            }else{
                                mapping.push(residue_number);
                            }
                        }
                    }else{
                        console.log("No mapping for pdb "+pdbid+" and chain"+ chainid)
                        mapping = [range_string.split("-")[0],range_string.split("-")[1]];
                    }
                    var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid} filter-range=${mapping}></pdb-topology-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
                }).catch(error => {
                    var topview = document.querySelector('#topview');
                    console.log(error);
                    this.topology_loaded = 'error';
                    topview.innerHTML = "Failed to fetch the secondary structure!<br>Try another structure."
                });
            }).catch(error => {
                var topview = document.querySelector('#topview');
                console.log(error);
                this.topology_loaded = 'error';
                topview.innerHTML = "Failed to load the viewer!<br>Try another structure."
            });
        }, showPDBViewer(pdbid, chainid, entityid){
            if (document.querySelector("pdbe-molstar")) {return;}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            var viewerInstance = new PDBeMolstarPlugin();
            var options = {
                customData: { url: `https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=${entityid}&encoding=bcif`, 
                                format: 'cif', 
                                binary:true },
                hideCanvasControls: ["expand", "selection", " animation"],
                assemblyId: '1',
                hideControls: true,
                subscribeEvents: true,
                bgColor: {r:255,g:255,b:255},
            }
            var viewerContainer = document.getElementById('pdbeMolstarView');
            viewerInstance.render(viewerContainer, options);
            window.viewerInstance = viewerInstance;
            
            document.addEventListener('PDB.topologyViewer.click', (e) => {
                var molstar= viewerInstance;
                var chainId=e.eventData.chainId;
                var entityId=e.eventData.entityId;
                var residueNumber=e.eventData.residueNumber;
                var types=e.eventData.type;
                molstar.visual.select({
                    data:[
                        {
                            entity_id:entityId,
                            start_residue_number:residueNumber,
                            end_residue_number:residueNumber,
                            color:{r:20, y:100, b:200},
                            focus:false
                        },
                    ],
                })
            })
            document.addEventListener('PDB.topologyViewer.mouseover', (e) => {
                var molstar= viewerInstance;
                var chainId=e.eventData.chainId;
                var entityId=e.eventData.entityId;
                var residueNumber=e.eventData.residueNumber;
                var types=e.eventData.type;
                
                molstar.visual.highlight({
                    data:[
                        {
                            entity_id:entityId,
                            start_residue_number:residueNumber,
                            end_residue_number:residueNumber,
                        },
                    ],
                })
            })
            document.addEventListener('PDB.molstar.mouseover', (e) => {
                var eventData = e.eventData;
                let resi_id = eventData.auth_seq_id;
                if(masked_array && masked_array[resi_id] == false) {
                    viewerInstance.plugin.behaviors.interaction.hover._value.current.loci.kind = "empty-loci"
                }
            });
        }, populateECODranges() {
            $.ajax ({
                type: "GET",
                url: "/alignments/authEcodQuery",
                data: {url: `/desire-api/ECOD-domains/?pdb=${vm.pdbid}&chain=${vm.chainid}`},
                success: function (data){
                    vm.domain_list = []
                    for (var i = 0; i < data.results.length; i++) {
                        let re = /\d+-\d+$/;
                        let range_str = re.exec(data.results[i].pdb_range)[0] + ';';
                        vm.domain_list.push({name: data.results[i].x_name + ' ' + data.results[i].f_name, range: range_str});
                    }
                }
            });
        },postStructureData(pdbid, chainid) {
            const topview_item = document.getElementById("topview");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading Structure Data", true)}
            let tempEntities = this.chains.filter(obj => {
                return obj["value"] == chainid;
            });
            let entityIDS = [];
            tempEntities.forEach(function(ent){
                entityIDS.push(ent["entityID"])
            })
            testingCIFParsing(pdbid, entityIDS);
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
        },handleFilterRange(filter_range){
            handleFilterRange(filter_range);
        },handlePropensities(checked_propensities){
            handlePropensities(checked_propensities);
        },populatePDBs(alndata){
            populatePDBs(alndata);
        },populatePDBsFromCustomAln(firstSeq){
            if (this.guideOff){
                populatePDBsFromCustomAln(firstSeq);
            }
        },getPropensities(sequence_indices) {
            getPropensities(sequence_indices);
        },listSecondaryStructures() {
            listSecondaryStructures();
        }
    }, 
    mounted() {
        addFooterImages("footerDiv");
    }
}
</script>