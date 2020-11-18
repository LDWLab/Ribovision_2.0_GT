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
            <div v-if="type_tree=='para'|type_tree=='orth'">
            <treeselect 
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
                <p><button v-on:click="submitCustomAlignment()">Upload alignment</button>{% csrf_token %}</p>
            </div>
            <p>
                <select id="selectaln" v-if="tax_id" v-model="alnobj" @change="showAlignment(alnobj.id, tax_id, type_tree)" >
                    <option v-if="tax_id" :value="null" selected disabled hidden>Select an alignment</option>
                    <option v-if="tax_id" v-for="aln in alignments" v-bind:value="{ id: aln.value, text: aln.text }">{{ aln.text }}</option>
                </select>
            </p>
            <p>
                <span v-if="alnobj&&alnobj!='custom'">Selected alignment: {{ alnobj.text }}.<br></span>
                <span v-if="alnobj">Input PDB and polymer for mapping:</span>
            </p>
            <p><input id="pdb_input" v-if="alnobj" v-model="pdbid" v-on:input="getPDBchains(pdbid, alnobj.id)" placeholder="4v9d" maxlength="4">
                <div v-if="alnobj" class="checkbox">
                    <label><input type="checkbox" v-model="checked_propensities" v-on:change="handlePropensities(checked_propensities)">Propensities</label>
                </div>
            </p>
            <p v-if="hide_chains">Looking for available polymers...</p>
            <p><select v-bind:style="{ resize: 'both'}" multiple v-if="chains&&fasta_data" v-model="chainid" >
                <option :value ="null" selected disabled>Select polymer</option>
                <option v-for="chain in chains" v-bind:value="chain.value" @click="showTopologyViewer(pdbid, chainid, fasta_data); showPDBViewer(pdbid, chainid, chain.entityID)">{{ chain.text }}</option>
            </select></p>
            <div v-if="chainid">
                <button type="button" v-on:click="downloadCSVData()">
                    Download mapped data
                </button>
            </div>
            <div v-if="topology_loaded != 'False'">
                <p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_filter" v-on:change="cleanFilter(checked_filter, masking_range)">Masking ranges</label>
                    </div>
                    <span v-if="checked_filter">Residue ranges to show, separated by semicolon. <br> For example: 1-80;91-111;</span>
                    <input v-if="checked_filter" v-model="masking_range" v-on:input="handleMaskingRanges(masking_range)">
                </p>
                <p v-if="correct_mask!='True'&&masking_range!=null">Incorrect range syntax!</p>
                <p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_selection" v-on:change="cleanSelection(checked_selection, filter_range)">Filter Range</label>
                    </div>
                    <span v-if="checked_selection">Residue range to show </span>
                    <input v-if="checked_selection" v-model="filter_range" v-on:input="handleFilterRange(filter_range)">
                </p>
                <p>
                    <div class="checkbox">
                        <label><input type="checkbox" v-model="checked_customMap" v-on:change="cleanCustomMap(checked_customMap)">Custom Data</label>
                        <p><input v-if="checked_customMap" type="file" accept=".csv" ref="custom_csv_file" v-on:change="handleCustomMappingData()"/></p>
                    </div>
                </p>
            </div>
        </div>
        <div class="alignment_section">
            <p></p>
            <span id="alnif" v-if="alnobj">
                <div id="alnDiv">Loading alignment...</div>
                <div id="testaln"></div>
            </span><br>
            <span v-if="aln_meta_data">Selected residue: {{ aln_meta_data }}</span>
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
        <footer >Footer</footer>
    </div>
</template>


<script>
  import {Tooltip} from './Tooltip.js'
  import { actions, 
    MSAViewer, 
    SequenceViewer,
    Labels,
    msaConnect,
    withPositionStore, } from '@plotly/react-msa-viewer';
  import React, { Component } from "react";
  import ReactDOM, { render } from 'react-dom';
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
            pdbid: null,
            chains: null,
            chainid: null,
            aln_meta_data: null,
            fasta_data: null,
            fastaSeqNames: null,
            hide_chains: null,
            type_tree: "orth",
            aa_properties: null,
            structure_mapping: null,
            file: null,
            custom_aln_twc_flag: null,
            topology_loaded: 'False',
            twc_loaded: false,
            masking_range: null,
            filter_range: null,
            correct_mask: false,
            coil_residues: null,
            checked_filter: false,
            checked_selection: false,
            checked_customMap: false,
            csv_data: null,
            checked_propensities: false,
            helix_residues: null,
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
        },
    },methods: {
        handleFileUpload(){
            this.file = this.$refs.custom_aln_file.files[0];
            if (this.tax_id != null){this.tax_id = null;}
        },
        submitCustomAlignment(){
            let formData = new FormData();
            formData.append('custom_aln_file', this.file)
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
                    if (data["results"].length === 2) {
                        function getObjIntersection(o1, o2) {
                            return Object.keys(o1).filter({}.hasOwnProperty.bind(o2));
                        }
                        var alns_first_tax = Object.fromEntries(data["results"][0]["alignment_ids"]);
                        var alns_second_tax = Object.fromEntries(data["results"][1]["alignment_ids"]);
                        var aln_indexes = getObjIntersection(alns_first_tax, alns_second_tax);
                        var fpa = []
                        aln_indexes.forEach(function(alnk) {
                            fpa.push(Array(Number(alnk), alns_first_tax[alnk]))
                        });
                    } else {
                        var fpa = data["results"][0]["alignment_ids"]
                    }
                    var fpa_viz = [];
                    fpa.forEach(function(fkey) {
                        fpa_viz.push({
                            text: fkey[1],
                            value: fkey[0]
                        });
                    });
                    this.alignments = fpa_viz
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
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase())
                    .then(struc_data => {
                        var chain_list = struc_data[pdbid.toLowerCase()];
                        if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
                        if (this.type_tree != "upload") {
                            filterAvailablePolymers(chain_list, aln_id, vm);
                            this.hide_chains = null;
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
                        alert("Problem with parsing the chains:\n" + error)
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
                window.main_elmnt = main_elmnt;
                let seqsForMSAViewer = parseFastaSeqForMSAViewer(fasta['Alignment']);
                var msaOptions = {
                    sequences: seqsForMSAViewer,
                    colorScheme: "clustal2",
                    height: main_elmnt.offsetHeight * 0.9,
                    width: main_elmnt.offsetWidth * 0.75,
                    tileHeight: 18,
                    tileWidth: 18,
                    overflow: "auto",
                };
                window.msaOptions = msaOptions;
                ReactDOM.render(
                    <RV3AlnViewer ref={(RV3AlnViewer) => {window.RV3AlnViewer = RV3AlnViewer}}/>,
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
                    var mapping = [];
                    var range_string = minIndex.concat("-").concat(maxIndex);
                    GetRangeMapping(pdbid, chainid, range_string, mapping);
                    let data_string = JSON.stringify(Array.from(mapped_aa_properties.entries())).replaceAll(",[[", ":").replaceAll("]],",";").replaceAll("],[",",");
                    let formatted_data_string = data_string.replaceAll("[","").replaceAll("]","").replaceAll("\"","");
                    var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${formatted_data_string} filter-range=${mapping}></pdb-topology-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
                    this.topology_loaded = 'True';
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
        }
    }
}

class RV3AlnViewer extends Component {
  state = { 
      tileWidth: 18,
      tileHeight: 18,
      aaPos: 0,
      seqPos: 0,
      width: main_elmnt.offsetWidth * 0.7,
      height: main_elmnt.offsetHeight * 0.9,
      highlight: null 
  };
  handleResize = () => {
      this.setState({
          width: main_elmnt.offsetWidth * 0.7,
          height: main_elmnt.offsetHeight * 0.9
      });
      //var style = document.querySelector('[data="rv3_style"]');
      //style.innerHTML = ".slider::-webkit-slider-thumb { width: "+main_elmnt.offsetWidth*0.05+"px}"
  };
  componentDidMount() {
      window.addEventListener("resize", this.handleResize);
      //var style = document.querySelector('[data="rv3_style"]');
      //style.innerHTML = ".slider::-webkit-slider-thumb { width: "+main_elmnt.offsetWidth*0.05+"px}"
  };
  componentWillUnmount() {
      window.removeEventListener("resize", this.handleResize);
  };
  onResidueMouseEnter = e => {
      if (vm.topology_loaded == 'True'){
          let resiPos = vm.structure_mapping[e.position];
          if (resiPos !== undefined){
              viewerInstanceTop.pluginInstance.highlight(resiPos, resiPos);
              viewerInstance.visual.highlight({
                  data:[{
                          entity_id:vm.entityId,
                          start_residue_number:resiPos,
                          end_residue_number:resiPos,
                      },],
              });
          }
      }
      if (!window.ajaxRun){
          window.ajaxRun = true;
          if (e.position !== undefined){
              registerHoverResiData(e, this, vm);
          }
      } else {
          this.setState({ fold: undefined, phase: undefined });
      }
   };
  onResidueMouseLeave = e => {
      if (vm.topology_loaded == 'True'){
          viewerInstanceTop.pluginInstance.clearHighlight();
          viewerInstance.visual.clearHighlight();
      }
      this.setState({ fold: undefined, phase: undefined });
  };
  highlightRegion = () => {
      const highlight = {
          sequences: {
            from: 0,
            to: 2
          },
          residues: {
            from: 2,
            to: 13
          }
        };
      this.setState({ highlight });
  };
  removeHighlightRegion = () => {
      this.setState({ highlight: null });
  };                
  render() {
      const xPos = this.state.tileWidth * (this.state.aaPos - 1);
      const yPos = this.state.tileHeight * (this.state.seqPos - 1);
      const maxXpos = window.aaFreqs.length - Math.round(((main_elmnt.offsetWidth * 0.7)/this.state.tileWidth))+2;
      const maxYpos = vm.fastaSeqNames.length - Math.round(((main_elmnt.offsetHeight * 0.9)/this.state.tileHeight))+2;
      return (
      <div style={{ display: "flex" }}>
          <div>
            <input
              style = {{ 
                  width: main_elmnt.offsetWidth * 0.7+"px",
                  position: "relative",
                  left: main_elmnt.offsetWidth * 0.2+"px"
              }}
              type="range"
              min="0"
              max={maxXpos}
              value={this.state.aaPos}
              onChange={(evt) => this.setState({ aaPos: evt.target.value })}
              class="slider"
              id="xPosSlider"
              />
          <MSAViewer 
          {...msaOptions}
          ref={(ref) => (this.el = ref)}
          highlight={this.state.highlight}
          width={this.state.width}
          height={this.state.height}
          tileWidth={this.state.tileWidth}
          tileHeight={this.state.tileHeight}
          position={{ xPos, yPos }}
          >
          <div style={{ position: "relative", display: "flex"}}>
          <Labels style={{
              width: main_elmnt.offsetWidth * 0.2
              }}/>
          <div>
              <SequenceViewer
                onResidueMouseEnter={this.onResidueMouseEnter}
                onResidueMouseLeave={this.onResidueMouseLeave}
              />
              {this.state.fold && (
                <div
                  style={{
                    position: "absolute",
                    opacity: 0.8,
                    ...this.state.tooltipPosition,
                  }}
                >
                  <Tooltip>
                    Fold: {this.state.fold} <br></br>
                    Phase: {this.state.phase}
                  </Tooltip>
                </div>
              )}
              </div>
          </div>
          </MSAViewer>
      </div>
      <input
          style={{ 
              width: main_elmnt.offsetHeight*0.9+"px",
          }}
          type="range"
          min="0"
          max={maxYpos}
          value={this.state.seqPos}
          onChange={(evt) => this.setState({ seqPos: evt.target.value })}
          class="slider"
          id="yPosSlider"
          />
      </div>
      );
  }
}
</script>