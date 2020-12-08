<template>
    <div id="phylo_tree_dropdown">
        <div class="left-sidebar">
            <p id="tree_type">
                <input type="radio" id="orthologs" value="orth" v-model="type_tree" v-on:input="cleanTreeOpts()" checked>
                <label for="orthologs" >Orthologs</label>
                <input type="radio" id="paralogs" value="para" v-model="type_tree" v-on:input="cleanTreeOpts()">
                <label for="paralogs">Paralogs</label>
                <input type="radio" id="upload" value="upload" v-model="type_tree" v-on:input="cleanTreeOpts()">
                <label for="upload">Upload</label>
            </p>
            <div id="treeselect" v-if="type_tree=='para'|type_tree=='orth'">
            <treeselect ref="treeselect"
              :load-options="loadOptions"
              v-model="tax_id" 
              v-on:input="loadData(tax_id, type_tree)"
              placeholder="Select a group"
              no-children-text="Loading... or no children"
              :multiple="true" 
              :options="options" 
              :flat="true"
              :limit="2"
              :default-expand-level="1"
              >Loading phylogenetic tree...</treeselect>
            </div>
            <div v-else>
                <p>Select alignment file: </p>
                <p><input type = "file" accept=".fasta,.fas,.fa" ref="custom_aln_file" v-on:change="handleFileUpload()"/></p>
                <p><button v-on:click="submitCustomAlignment()">Upload alignment</button></p>
            </div>
            <p>
                <select id="selectaln" v-if="tax_id" v-model="alnobj">
                    <option v-if="tax_id" :value="null" selected disabled hidden>Select an alignment</option>
                    <option v-if="tax_id" v-for="aln in alignments" v-bind:value="{ id: aln.value, text: aln.text }">{{ aln.text }}</option>
                </select>
            </p>
            <p>
                <span v-if="alnobj&&alnobj!='custom'">Selected alignment: {{ alnobj.text }}.<br></span>
                <span v-if="alnobj">Input PDB and polymer for mapping:</span>
            </p>
            <p>
                <select id="pdb_input" v-if="alnobj&&alnobj!='custom'" v-model="pdbid">
                    <option :value="null" selected disabled hidden>Select a pdb</option>
                    <option v-for="pdb in pdbs" v-bind:value="pdb.id">{{pdb.name}}</option>
                </select>
                <input type="text" id="pdb_input_custom" v-if="alnobj&&alnobj=='custom'" v-model="pdbid" maxlength="4"></input>
                <div v-if="hide_chains" id="onFailedChains">Looking for available polymers...</div>
            </p>
            <p>
                <div v-if="alnobj" class="checkbox">
                    <label><input type="checkbox" v-model="checked_propensities" v-on:change="handlePropensities(checked_propensities)">Propensities</label>
                    <select v-if="checked_propensities&&structure_mapping" v-model="property">
                        <option :value="null" selected disabled hidden>Select a substructure</option>
                        <option v-for="substructure in substructures" v-bind:value="{ id: substructure.value, text: substructure.text }">{{ substructure.text }}</option>
                    </select>
                </div>
            </p>
            <p><select id="polymerSelect" v-bind:style="{ resize: 'both'}" multiple v-if="chains&&fasta_data&&pdbid" v-model="chainid" >
                <option :value ="null" selected disabled>Select polymer</option>
                <option v-for="chain in chains" v-bind:value="chain.value" @click="showTopologyViewer(pdbid, chainid, fasta_data); showPDBViewer(pdbid, chainid, chain.entityID)">{{ chain.text }}</option>
            </select></p>
            <div v-if="poor_structure_map">
                <p style="color:#DE3163"><b>Warning!!!<br>
                Poor structure to alignment mapping!<br>
                There were {{poor_structure_map}} poorly mapped residues!<br>
                Proceed with caution or try different structure.</b></p>
            </div>
            <div v-if="structure_mapping">
                <button id="downloadDataBtn" type="button" v-on:click="downloadCSVData()">
                    Download mapped data
                </button>
            </div>
            <div v-if="topology_loaded != 'False'">
                <div id="maskingSection"><p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_filter" v-on:change="cleanFilter(checked_filter, masking_range)">Masking ranges</label>
                    </div>
                    <span v-if="checked_filter">Residue ranges to show, separated by semicolon. <br> For example: 1-80;91-111;</span>
                    <input v-if="checked_filter" v-model="masking_range" v-on:input="handleMaskingRanges(masking_range)">
                </p></div>
                <p v-if="correct_mask!='True'&&masking_range!=null">Incorrect range syntax!</p>
                <div id="filterSection"><p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_selection" v-on:change="cleanSelection(checked_selection, filter_range)">Filter Range</label>
                    </div>
                    <span v-if="checked_selection">Residue range to show </span>
                    <input v-if="checked_selection" v-model="filter_range" v-on:input="handleFilterRange(filter_range)">
                </p></div>
                <div id="customDataSection">
                <p><div class="checkbox">
                        <label><input type="checkbox" v-model="checked_customMap" v-on:change="cleanCustomMap(checked_customMap)">Custom Data</label>
                        <p><input v-if="checked_customMap" type="file" accept=".csv" ref="custom_csv_file" v-on:change="handleCustomMappingData()"/></p>
                    </div>
                </p></div>
            </div>
        </div>
        <div class="alignment_section">
            <div id="alnif" v-if="alnobj">
                <div id="alnMenu" style="display: flex;">
                    <button id="downloadFastaBtn" style="margin: 0 1%;" v-if="colorScheme"  type="button" v-on:click="downloadAlignmentData()">
                        Download alignment
                    </button>
                    <button id="downloadAlnImageBtn" style="margin: 0 1%;" v-if="colorScheme"  type="button" v-on:click="downloadAlignmentImage()">
                        Download alignment image
                    </button>
                    <select id="selectAlnColorScheme" style="margin: 0 1%;" v-model="colorScheme" v-if="colorScheme">
                        <option :value="null" selected disabled>Select a colorscheme</option>
                        <option v-for="colorscheme in availColorschemes" >{{ colorscheme }}</option>
                    </select>
                </div>
                <div id="alnDiv">Loading alignment...</div>
            </div>
        </div>
        <div class="topology_section">
            <span id="topif" v-if="chainid">
                <div id="topview">Loading topology viewer and conservation data...</div>
            </span>
        </div>
        <div class="molstar_section">
            <span id="molif" v-if="chainid">
                <div id ="pdbeMolstarView">Loading Molstar Component...</div>
            </span>
        </div>
        <div class = "propensity_section">
            <span id="propif" v-if="checked_propensities">
                <div id = "total"></div>
            </span>
        </div>
        <footer >Footer</footer>
    </div>
</template>


<script>
  import {AlnViewer} from './AlignmentViewer.js'
  import ReactDOM, { render } from 'react-dom';
  import React, { Component } from "react";
  import Treeselect from '@riophae/vue-treeselect'
  export default {
      // register the component
      components: { Treeselect },
      data() {
        return {
            tax_id: null,
            alnobj: null,
            options: null,
            alignments: null,
            pdbs: [
                {id: "4v9d", name: "4V9D E. coli"},
                {id: "4v6u", name: "4V6U P. furiosus"},
                {id: "4ug0", name: "4UG0 H. sapiens"},
                ],
            availColorschemes: [
                "buried","cinema","clustal","clustal2","helix","lesk","mae","nucleotide","purine","strand","taylor","turn","zappo",
                ],
            pdbid: null,
            chains: null,
            chainid: null,
            fasta_data: null,
            fastaSeqNames: null,
            colorScheme: null,
            hide_chains: null,
            type_tree: "orth",
            aa_properties: null,
            structure_mapping: null,
            poor_structure_map: null,
            file: null,
            custom_aln_twc_flag: null,
            topology_loaded: 'False',
            twc_loaded: false,
            masking_range: null,
            filter_range: null,
            correct_mask: false,
            checked_filter: false,
            checked_selection: false,
            checked_customMap: false,
            csv_data: null,
            checked_propensities: false,
            coil_residues: null,
            helix_residues: null,
            strand_residues: null,
            substructures: null,
            property: null,
        }
      },
      watch: {
        csv_data: function (csv_data) {
            var topviewer = document.getElementById("PdbeTopViewer");
            var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
            if (csv_data == null){
                selectBoxEle.removeChild(selectBoxEle.childNodes[selectBoxEle.options.length-1]);
                topviewer.pluginInstance.resetDisplay();
                return;
            }
            let custom_data = csv_data.split('\n').map(function(e){
                return e.split(',').map(Number);
            })
            if (custom_data[custom_data.length-1] == 0){custom_data.splice(-1,1)}
            if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
                let vals = custom_data.map(function(v){ return v[1] });
                let indexes = custom_data.map(function(v){ return v[0] });
                window.aaColorData.set("CustomData", [viridis]);
                window.aaPropertyConstants.set("CustomData", [Math.min(...vals), Math.max(...vals)]);
                let coilsOutOfCustom = this.coil_residues.filter(value => !indexes.includes(value));
                window.coilsOutOfCustom = coilsOutOfCustom;
                var custom_prop = new Map();
                custom_prop.set("CustomData", custom_data);
                topviewer.pluginInstance.getAnnotationFromRibovision(custom_prop);
                window.custom_prop = custom_prop;
                var custom_option = document.createElement("option");
                custom_option.setAttribute("value", selectBoxEle.options.length);
                custom_option.appendChild(document.createTextNode("Custom Data"));
                selectBoxEle.appendChild(custom_option);
                if(this.correct_mask == 'True') {
                    var j = topviewer.pluginInstance.domainTypes.length-1;
                    colorResidue(j, window.masked_array);
                }
            }
        },alnobj: function (data){
            this.populatePDBs(data);
            this.showAlignment(data.id, vm.tax_id, vm.type_tree);
        },pdbid: function (pdbid){
            this.getPDBchains(pdbid, vm.alnobj.id);
        },colorScheme: function (scheme){
            if (window.AlnViewer){
                window.AlnViewer.setState({colorScheme:scheme});
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
            cleanupOnNewAlignment(this, "Select new alignment!");
            [this.options, this.tax_id, this.alnobj] = [null, null, null];
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
            if (type_tree == "upload"){this.tax_id = null; return;}
            if (value.length == 0){this.tax_id = null; return;}
            cleanupOnNewAlignment(this, "Select new alignment!");
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
                this.chains = null
                this.hide_chains = true
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase()).then(struc_data => {
                    var chain_list = struc_data[pdbid.toLowerCase()];
                    if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
                    if (this.type_tree != "upload") {
                        filterAvailablePolymers(chain_list, aln_id, vm);
                    } else {
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
            if (type_tree == "orth"){
                var url = `/ortholog-aln-api/${aln_id}/${taxid}`}
            if (type_tree == "para"){
                var url = '/paralog-aln-api/'+aln_id.split(',')[1]}
            if (type_tree == "upload"){
                var url = '/custom-aln-data'}
            ajax(url).then(fasta => {
                if (fasta['TwinCons'] != null){
                    this.custom_aln_twc_flag = fasta['TwinCons']
                }
                vm.fastaSeqNames = fasta['Sequence names'];
                window.aaFreqs = fasta['AA frequencies'];
                var main_elmnt = document.querySelector(".alignment_section");
                var msaHeight = main_elmnt.offsetHeight * 0.8;
                if (msaHeight > 17*(vm.fastaSeqNames.length+2)){
                    var alnifEle = document.querySelector('#alnif');
                    alnifEle.style.paddingTop="10%";
                    msaHeight = 17*(vm.fastaSeqNames.length+2);
                }
                let seqsForMSAViewer = parseFastaSeqForMSAViewer(fasta['Alignment']);
                var msaOptions = {
                    sequences: seqsForMSAViewer,
                    colorScheme: vm.colorScheme,
                    height: msaHeight,
                    width: main_elmnt.offsetWidth * 0.7,
                    tileHeight: 17,
                    tileWidth: 17,
                };
                window.msaOptions = msaOptions;
                ReactDOM.render(
                    <AlnViewer ref={(AlnViewer) => {window.AlnViewer = AlnViewer}}/>,
                    document.getElementById('alnDiv')
                  );
                this.fasta_data = fasta['Alignment'];
                this.aa_properties = calculateFrequencyData(fasta['AA frequencies']);
            })
        }, showTopologyViewer (pdbid, chainid, fasta){
            window.filterRange = "-10000,10000";
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(this);}
            if (chainid.length > 1){this.chainid = chainid[0];}
            const topview_item = document.getElementById("topview");
            const molstar_item = document.getElementById("pdbeMolstarView");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading topology viewer and conservation data...")}
            if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Loading Molstar Component...")}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            let temp = this.chains.filter(obj => {
                return obj["value"] == chainid;
            })[0];
            let ebi_sequence = temp["sequence"];
            let startIndex = temp["startIndex"];
            // let ebi_sequence = this.chains[0]["sequence"];
            ajax('/mapSeqAln/', {fasta, ebi_sequence, startIndex}).then(struct_mapping=>{
                this.structure_mapping = struct_mapping;
                if (struct_mapping['BadMappingPositions']){this.poor_structure_map = struct_mapping['BadMappingPositions'];}
                var mapped_aa_properties = mapAAProps(this.aa_properties, struct_mapping);
                if ((this.tax_id != null && this.tax_id.length == 2) || (this.custom_aln_twc_flag != null && this.custom_aln_twc_flag == true) || (this.type_tree == 'para')) {
                    ajax('/twc-api/', {fasta}).then(twcDataUnmapped => {
                        const build_mapped_props = function(mapped_props, twcDataUnmapped, structure_mapping){
                            mapped_props.set("TwinCons", [])
                            for (let i = 0; i < twcDataUnmapped.length; i++) {
                                let mappedI0 = structure_mapping[twcDataUnmapped[i][0]];
                                if (mappedI0) {
                                    mapped_props.get("TwinCons").push([mappedI0, twcDataUnmapped[i][1]]);
                                }
                            }
                            return mapped_props;
                        }
                        var topviewer = document.getElementById("PdbeTopViewer");
                        mapped_aa_properties = build_mapped_props(mapped_aa_properties, twcDataUnmapped, this.structure_mapping);
                        window.mapped_aa_properties = mapped_aa_properties;
                        if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
                            var empty_props = new Map();
                            let twc_props = build_mapped_props(empty_props, twcDataUnmapped, this.structure_mapping);
                            topviewer.pluginInstance.getAnnotationFromRibovision(twc_props);
                            var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
                            var twc_option = document.createElement("option");
                            twc_option.setAttribute("value", selectBoxEle.options.length);
                            twc_option.appendChild(document.createTextNode("TwinCons"));
                            selectBoxEle.appendChild(twc_option);
                        }
                    })
                }
                window.mapped_aa_properties = mapped_aa_properties;
                var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
                ajax(topology_url).then(data => {
                    var entityid = Object.keys(data[pdblower])[0];
                    vm.coil_residues = filterCoilResidues(data[pdblower][entityid][chainid]["coils"])
                    vm.helix_residues = filterCoilResidues(data[pdblower][entityid][chainid]["helices"])
                    vm.strand_residues = filterCoilResidues(data[pdblower][entityid][chainid]["strands"])
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

                        let data_string = JSON.stringify(Array.from(mapped_aa_properties.entries())).replaceAll(",[[", ":").replaceAll("]],",";").replaceAll("],[",",");
                        let formatted_data_string = data_string.replaceAll("[","").replaceAll("]","").replaceAll("\"","");
                        var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${formatted_data_string} filter-range=${mapping}></pdb-topology-viewer>`
                        document.getElementById('topview').innerHTML = topology_viewer;
                        window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
                        this.topology_loaded = 'True';
                    })
                })
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
        },downloadAlignmentImage() {
            downloadAlignmentImage(document.querySelector('#alnDiv'));
        },downloadAlignmentData() {
            downloadAlignmentData(vm.fasta_data);
        },downloadCSVData() {
            downloadCSVData();
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
        }
    }
}
</script>