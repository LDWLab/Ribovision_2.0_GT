function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}
var csrftoken = getCookie('csrftoken');    function csrfSafeMethod(method) {
    // these HTTP methods do not require CSRF protection
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}
$.ajaxSetup({
    beforeSend: function(xhr, settings) {
        if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});

function ajax(url, optional_data='') {
    if (optional_data != ''){
        //var el = document.getElementsByName("csrfmiddlewaretoken");
        //csrf_value = Cookies.get('csrftoken');
        //csrf_value = el[0].getAttribute("value");
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: 'POST',
                dataType: "json",
                data: optional_data,
                headers: {'X-CSRFToken': csrftoken},
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    }else{
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: 'GET',
                dataType: "json",
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    }
}

var pushChainData = function(temp_arr, chain_listI){
    temp_arr.push({
        text: chain_listI["molecule_name"][0],
        value: chain_listI["in_chains"][0],
        sequence: chain_listI["sequence"],
        entityID: chain_listI["entity_id"],
        startIndex: chain_listI.source[0].mappings[0].start.residue_number
    })
    return temp_arr;
}

var filterAvailablePolymers = function(chain_list, aln_id, vueObj) {
    let temp_arr = [];
    let url = `/desire-api/alignments/${aln_id}/?format=json`;
    ajax(url).then( aln_data => {
        for (let i = 0; i < chain_list.length; i++) {
            let chain_listI = chain_list[i]
            if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
            if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
            for (let ix =0; ix < aln_data["polymers"].length; ix++){
                if (aln_data["polymers"][ix]["genedescription"].trim() == chain_list[i]["molecule_name"][0]){
                    temp_arr = pushChainData(temp_arr, chain_listI);
                }
            }
        }
    // console.log("___" + temp_arr[temp_arr.length - 1]["sequence"] + "___");
    chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
    if (chain_options.length === 0) {
        chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
    }
    vueObj.chains = chain_options;
    });
}

var create_deleted_element = function (parent_id, child_id, child_text) {
    const parent = document.getElementById(parent_id);
    const child_elt = document.createElement("div");
    const childText = document.createTextNode(child_text);
    child_elt.setAttribute("id", child_id);
    child_elt.setAttribute("id", child_id);
    child_elt.appendChild(childText);
    parent.appendChild(child_elt);
}

var cleanupOnNewAlignment = function (vueObj, aln_text='') {
    const menu_item = document.querySelector(".smenubar");
    const aln_item = document.getElementById("alnDiv");
    const topview_item = document.getElementById("topview");
    const molstar_item = document.getElementById("pdbeMolstarView");
    const pdb_input = document.getElementById("pdb_input");
    if (menu_item) {menu_item.remove();}
    if (aln_text != ''){
        vueObj.custom_aln_twc_flag == null;
        window.mapped_aa_properties == null;
        if (pdb_input) {
            if (pdb_input.getAttribute("value") != ""){vueObj.pdbid = null;}
        }
        if (vueObj.chains) {vueObj.chains = null;}
        if (vueObj.aln_meta_data) {vueObj.aln_meta_data = null;}
        if (vueObj.fasta_data) {vueObj.fasta_data = null;}
        if (vueObj.frequency_data) {vueObj.frequency_data = null;}
        if (vueObj.topology_loaded) {vueObj.topology_loaded = 'False';}
        if (aln_item) {aln_item.remove(); create_deleted_element("alnif", "alnDiv", aln_text)}
    }
    if (window.masked_array.length > 0) {window.masked_array = [];}
    if (vueObj.masking_range) {vueObj.masking_range = null;}
    //if (vueObj.chainid) {vueObj.chainid = null;}
    if (vueObj.checked_filter) {vueObj.checked_filter = false;}
    if (vueObj.checked_customMap) {vueObj.checked_customMap = false;}
    if (vueObj.csv_data) {vueObj.csv_data = null;}
    if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Select new chain!")}
    if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Select new structure!")}
}

var loadParaOptions = function (action, callback, vm) {
    if (action === "LOAD_ROOT_OPTIONS"){
        ajax('/alignments/showStrucTaxonomy').then(data =>{
            data.isDisabled = true,
            vm.options = [data];
            callback();
        }).catch(error => {
            console.log(error)
        })
    }
}

var loadParaAlns = function (value, vm) {
    vm.alignments = null;
    ajax('/alignments/fold-api/'+value).then(data=>{
        var fpa = data["Folds to polymers to alignments"]
        var fpa_viz = [];
        Object.keys(fpa).forEach(fkey => {
            Object.keys(fpa[fkey]).forEach(pkey => {
                fpa[fkey][pkey].forEach(function (akey){
                    fpa_viz.push({
                        text:  'Alignment '.concat(akey[1],'; fold ',fkey),
                        value: fkey.concat(',',akey)
                    });
                });
            });
        });
        var temp_arr = fpa_viz
        fpa_viz = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
        vm.alignments = fpa_viz
    });
}

var calculateFrequencyData = function (frequencies){
    const multiplyvector = function (a,b){
        return a.map((e,i) => e * b[i]);
    }
    let aaPropertiesData = new Map([
                            ["Charge",[0,0,-1,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0]],
                            ["Hydropathy",[1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]],
                            ["Hydrophobicity",[0.02,0.77,-1.04,-1.14,1.35,-0.80,0.26,1.81,-0.41,1.14,1,-0.77,-0.09,-1.10,-0.42,-0.97,-0.77,1.13,1.71,1.11]],
                            ["Polarity",[0,1.48,49.7,49.9,0.35,0,51.6,0.13,49.5,0.13,1.43,3.38,1.58,3.53,52,1.67,1.66,0.13,2.1,1.61]],
                            ["Mutability",[100,44,86,77,51,50,91,103,72,54,93,104,58,84,83,117,107,98,25,50]],
                            ["Shannon entropy",[0.000000000000001,4.321928094887363]],
                            ["TwinCons",[-2.935,12.065]]
                        ]);
    let aaColorData = new Map([
                            ["Charge",[Blues, Reds]],
                            ["Hydropathy",[Blues, Reds]],
                            ["Hydrophobicity",[Reds, Blues]],
                            ["Polarity",[viridis]],
                            ["Mutability",[viridis]],
                            ["Shannon entropy",[plasma]],
                            ["TwinCons",[RdPu, YlGn]],
                        ]);
    window.aaColorData = aaColorData;
    window.aaPropertyConstants = aaPropertiesData;
    window.selectSections_RV1 = new Map();
    let outPropertyPosition = new Map();
    aaPropertiesData.forEach(function (data, property_name){
        if (property_name == "TwinCons"){return;}
        let const_data = data
        outPropertyPosition.set(property_name, [])
        frequencies.forEach(function (col_frequency) {
            if (property_name == "Shannon entropy"){
                const_data = new Array;
                col_frequency.forEach( function (single_freq){
                    if (single_freq == 0){
                        const_data.push(0)
                    }else{
                        const_data.push(Math.log2(single_freq)*-1)
                    }
                });
            }
            outPropertyPosition.get(property_name).push(multiplyvector(const_data, col_frequency));
        });
    });
    return outPropertyPosition;
}

var mapAAProps = function (aa_properties, mapping){
    let outPropertyMappedPosition = new Map();
    aa_properties.forEach(function (data, property_name){
        outPropertyMappedPosition.set(property_name, [])
        data.forEach(function (data, aln_ix) {
            let mappedI0 = mapping[aln_ix+1];
            if (mappedI0) {
                outPropertyMappedPosition.get(property_name).push([mappedI0, Number(math.sum(data).toFixed(2))]);
            }
        });
    });
    return outPropertyMappedPosition;
}

var filterCoilResidues = function (coil_data){
    const range = (start, stop, step) => Array.from({ length: (stop - start) / step + 1}, (_, i) => start + (i * step));
    let coilResidues = [];
    coil_data.forEach(function (coilRange){
        if (coilRange.start < coilRange.stop){
            coilResidues.push(range(coilRange.start, coilRange.stop, 1))
        }
    })
    return coilResidues.flat()
}

var generateCSVstring = function (mapped_data){
    let properties = Array.from(mapped_data.keys());
    let csv = 'Index,'
    csv += properties.join(',');
    csv += '\n';
    let csv_ix = [];
    
    mapped_data.get(properties[0]).forEach((datapoint) =>{
        csv_ix.push([datapoint[0]]);
    })

    properties.forEach((prop) => {
        let ix = 0;
        mapped_data.get(prop).forEach((datapoint) =>{
            csv_ix[ix].push(datapoint[1]);
            ix += 1;
        })
    })

    csv_ix.forEach((row) => {
        csv += row.join(',');
        csv += '\n';
    })

    return csv;
}

var masked_array = [];
Vue.component('treeselect', VueTreeselect.Treeselect, )

var vm = new Vue({
    el: '#phylo_tree_dropdown',
    delimiters: ['[[', ']]'],
    data: {
        tax_id: null,
        alnobj: null,
        options: null,
        alignments: null,
        pdbid: null,
        chains: null,
        chainid: null,
        aln_meta_data: null,
        fasta_data: null,
        hide_chains: null,
        type_tree: "orth",
        aa_properties: null,
        structure_mapping: null,
        file: null,
        custom_aln_twc_flag: null,
        topology_loaded: 'False',
        masking_range: null,
        correct_mask: false,
        coil_residues: null,
        checked_filter: false,
        checked_customMap: false,
        custom_prop: null,
        csv_data: null,
        checked_propensities: false,
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
                var masked_array = window.masked_array;
                var j = topviewer.pluginInstance.domainTypes.length-1;
                var f = 0;
                while(f < topviewer.pluginInstance.domainTypes[4].data.length) {
                    if(!masked_array[f] && topviewer.pluginInstance.domainTypes[j].data[f]) {
                        topviewer.pluginInstance.domainTypes[j].data[f].color = "rgb(255,255,255)";
                        topviewer.pluginInstance.domainTypes[j].data[f].tooltipMsg = "NaN";                   
                        selectSections_RV1.get(topviewer.pluginInstance.domainTypes[j].label)[f].color = {r: 255, g: 255, b: 255};

                    } if(!masked_array[f] && vm.coil_residues.includes(f) && topviewer.pluginInstance.domainTypes[j].data[f]) {
                        topviewer.pluginInstance.domainTypes[j].data[f].color = "rgb(0,0,0)";
                        topviewer.pluginInstance.domainTypes[j].data[f].tooltipMsg = "NaN";
                    }
                        
                    f++;
                }
            }
        },
    },
    methods: {
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
            cleanupOnNewAlignment(vm, "Select new alignment!");
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
                loadParaOptions(action, callback, vm);
            }
        }, loadData (value, type_tree) {
            if (type_tree == "upload"){this.tax_id = null; return;}
            if (value.length == 0){this.tax_id = null; return;}
            cleanupOnNewAlignment(vm, "Select new alignment!");
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
                loadParaAlns (value, vm)
            }
        }, getPDBchains(pdbid, aln_id) {
            if (pdbid.length === 4) {
                if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(vm);}
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
            cleanupOnNewAlignment(vm, "Loading alignment...");
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
                ReactDOM.render(
                    React.createElement(ReactMSAViewer.MSAViewer, options),
                    document.getElementById('testaln')
                  );
                this.fasta_data = fasta['Alignment'];
                this.aa_properties = calculateFrequencyData(fasta['AA frequencies'])
                var main_elmnt = document.querySelector(".alignment_section")
                var opts = {
                    el: document.getElementById("alnDiv"),
                    seqs: msa.io.fasta.parse(fasta['Alignment']),
                    colorscheme: {
                        scheme: "clustal2",
                    },
                    //columns: {
                    //    hidden: fasta['Gap-only columns'] // hidden columns
                    //},
                    zoomer: {
                        // general
                        alignmentWidth: main_elmnt.offsetWidth * 0.75,
                        alignmentHeight: main_elmnt.offsetHeight * 0.9,
                        columnWidth: 15,
                        rowHeight: 15,
                        labelNameLength: main_elmnt.offsetWidth * 0.18,
                        autoResize: false, // only for the width
                    },
                    conf: {
                        registerMouseHover: true,
                        registerMouseClicks: true,
                    },
                    // smaller menu for JSBin
                    menu: "small",
                    bootstrapMenu: true
                };
                var m = new msa.msa(opts);
                m.render();
                m.g.on("residue:click", function(data) {
                    vm.aln_meta_data = null;
                    const strainQuery = '&res__poldata__strain__strain=';
                    var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(data["rowPos"]) + 1)}&aln=${aln_id}${strainQuery}${fasta['Sequence names'][Number(data["seqId"])]}`
                    ajax(url).then(alnpos_data => {
                        ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
                            vm.aln_meta_data = resiData;
                        });
                    }).catch(error => {
                        vm.aln_meta_data = null;
                        console.log("No residue with alignment position: " + data["rowPos"] + ". In alignment " + aln_id + ". Of species " + fasta['Sequence names'][Number(data["seqId"])]);
                        console.log(error);
                    })
                });
                m.g.on("residue:mousein", function(data) {
                    console.log(data)
                });
            })
        }, showTopologyViewer (pdbid, chainid, fasta){
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(vm);}
            if (chainid.length > 1){this.chainid = chainid[0];}
            const topview_item = document.getElementById("topview");
            const molstar_item = document.getElementById("pdbeMolstarView");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading topology viewer and conservation data...")}
            if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Loading Molstar Component...")}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            let temp = vm.chains.filter(obj => {
                return obj["value"] == chainid;
            })[0];
            let ebi_sequence = temp["sequence"];
            let startIndex = temp["startIndex"];
            // let ebi_sequence = vm.chains[0]["sequence"];
            ajax('/mapSeqAln/', optional_data={fasta, ebi_sequence, startIndex}).then(struct_mapping=>{
                this.structure_mapping = struct_mapping;
                var mapped_aa_properties = mapAAProps(this.aa_properties, struct_mapping);
                if ((this.tax_id != null && this.tax_id.length == 2) || (this.custom_aln_twc_flag != null && this.custom_aln_twc_flag == true) || (this.type_tree == 'para')) {
                    ajax('/twc-api/', optional_data={fasta}).then(twcDataUnmapped => {
                        const build_mapped_props = function(mapped_props, twcDataUnmapped, structure_mapping){
                            mapped_props.set("TwinCons", [])
                            for (i = 0; i < twcDataUnmapped.length; i++) {
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
                    var mapping = [];
                    var range_string = minIndex.concat("-").concat(maxIndex);
                    GetRangeMapping(pdbid, chainid, range_string, mapping);
                    let data_string = JSON.stringify(Array.from(mapped_aa_properties.entries())).replaceAll(",[[", ":").replaceAll("]],",";").replaceAll("],[",",");
                    let formatted_data_string = data_string.replaceAll("[","").replaceAll("]","").replaceAll("\"","");
                    var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${formatted_data_string} filter-range=${mapping}></pdb-topology-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                    this.topology_loaded = 'True';
                })
            });
        }, showPDBViewer(pdbid, chainid, entityid){
            if (document.querySelector("pdbe-molstar")) {return;}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
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
        },cleanFilter(checked_filter){
            if (checked_filter){return;}
            if (this.masking_range == null){return;}
            window.masked_array = [];
            this.masking_range = null;
            var topviewer = document.getElementById("PdbeTopViewer");
            topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
            var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
            topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
            window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
        },isCorrectMask(mask_range){
            window.masking_range_array = null;
            if (mask_range.match(/^(\d+-\d+;)+$/)) {
                var temp_array = mask_range.split(';').join('-').split('-');
                temp_array = temp_array.slice(0, -1)
                var i = 0;
                var isCorrect = true;
                while(i < temp_array.length) {
                    if(i % 2 == 0) {
                        if(Number(temp_array[i]) > Number(temp_array[i + 1])) {
                            isCorrect = false;
                        }
                    }
                    i = i + 1;
                }
                window.masking_range_array = temp_array;
            }
            if(isCorrect) {
                this.correct_mask = 'True';
            }
            return isCorrect;
        },
        initializeMaskedArray() {
            var topviewer = document.getElementById("PdbeTopViewer");
            var masked_array = [];
            var j = 0;
            while(j < mapped_aa_properties.get(topviewer.pluginInstance.domainTypes[4].label).length) {
                masked_array[j] = false;
                i = 0;
                while(i < window.masking_range_array.length && !masked_array[j]) {
                    if(j >= window.masking_range_array[i] && j <= window.masking_range_array[i + 1]) {
                        masked_array[j] = true;
                    }
                    i = i+2;
                }
                j = j+1;
            }
            window.masked_array = masked_array;
            return masked_array;
        },
        handleMaskingRanges(mask_range){
            window.masking_range = mask_range;
            window.masking_range_array = null;
            if (this.isCorrectMask(mask_range)) {   
                var topviewer = document.getElementById("PdbeTopViewer");
                topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
                if(window.custom_prop) {
                    topviewer.pluginInstance.getAnnotationFromRibovision(custom_prop); 
                }
                var masked_array = this.initializeMaskedArray();          
                var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;

                var j = 4;
                while(j < topviewer.pluginInstance.domainTypes.length) {
                    var f = 0;
                    while(f < topviewer.pluginInstance.domainTypes[4].data.length) {
                        if(!masked_array[f] && topviewer.pluginInstance.domainTypes[j].data[f]) {
                            topviewer.pluginInstance.domainTypes[j].data[f].color = "rgb(255,255,255)";
                            topviewer.pluginInstance.domainTypes[j].data[f].tooltipMsg = "NaN";                   
                            selectSections_RV1.get(topviewer.pluginInstance.domainTypes[j].label)[f].color = {r: 255, g: 255, b: 255};

                        } if(!masked_array[f] && vm.coil_residues.includes(f) && topviewer.pluginInstance.domainTypes[j].data[f]) {
                            topviewer.pluginInstance.domainTypes[j].data[f].color = "rgb(0,0,0)";
                            topviewer.pluginInstance.domainTypes[j].data[f].tooltipMsg = "NaN";
                        }
                            
                        f++;
                    }
                    j++;
                    window.masked_array = masked_array;
                    this.correct_mask = 'True';
                }
                let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex]
                
                if (selectedData.data){
                    topviewer.pluginInstance.updateTheme(selectedData.data); 
                    window.viewerInstance.visual.select({data: selectSections_RV1.get(selectedData.label), nonSelectedColor: {r:255,g:255,b:255}});
                    }
                }
        }, cleanCustomMap(checked_customMap){
            if (checked_customMap){return;}
            var topviewer = document.getElementById("PdbeTopViewer");
            topviewer.pluginInstance.domainTypes = topviewer.pluginInstance.domainTypes.filter(obj => {return obj.label !== "CustomData"})
            window.coilsOutOfCustom = null;
            this.csv_data = null;
        }, handleCustomMappingData(){
            const readFile = function (fileInput) {
                var reader = new FileReader();
                reader.onload = function () {
                    vm.csv_data = reader.result;
                };
                reader.readAsBinaryString(fileInput);
            };
            readFile(this.$refs.custom_csv_file.files[0]);

        }, downloadCSVData() {

            csv = generateCSVstring(mapped_aa_properties);

            let anchor = document.createElement('a');
            anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
            anchor.target = '_blank';
            anchor.download = 'rv3data.csv';
            anchor.click();

        }, handlePropensities(checked_propensities){
            if (checked_propensities){
                console.log("Checked")
            }else{
                console.log("UnChecked")
            }
            
        }
    }
})

var options = {
    sequences: [
      {
        name: "seq.1",
        sequence: "MEEPQSDPSIEP-PLSQETFSDLWKLLPENNVLSPLPS-QA-VDDLMLSPDDLAQWLTED"
      },
      {
        name: "seq.2",
        sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
      },
      {
        name: "seq.3",
        sequence: "MEEPQSDLSIEL-PLSQETFSDLWKLLPPNNVLSTLPS-SDSIEE-LFLSENVAGWLEDP"
      },
    ],
    colorScheme: "zappo",
   };
  