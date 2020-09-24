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
                    temp_arr.push({
                        text: chain_listI["molecule_name"][0],
                        value: chain_listI["in_chains"][0],
                        sequence: chain_listI["sequence"],
                        startIndex: chain_listI.source[0].mappings[0].start.residue_number
                    })
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
        if (pdb_input) {
            if (pdb_input.getAttribute("value") != ""){vueObj.pdbid = null;}
        }
        if (vueObj.chains) {vueObj.chains = null;}
        if (vueObj.aln_meta_data) {vueObj.aln_meta_data = null;}
        if (vueObj.fasta_data) {vueObj.fasta_data = null;}
        if (vueObj.frequency_data) {vueObj.frequency_data = null;}
        if (aln_item) {aln_item.remove(); create_deleted_element("alnif", "alnDiv", aln_text)}
    }
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
    var aaProperties = ["Charge","Hydropathy","Hydrophobicity","Polarity","Mutability"]
    let aaPropertiesData = new Map([
                            ["Charge",[0,0,-1,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0]],
                            ["Hydropathy",[1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]],
                            ["Hydrophobicity",[0.02,0.77,-1.04,-1.14,1.35,-0.80,0.26,1.81,-0.41,1.14,1,-0.77,-0.09,-1.10,-0.42,-0.97,-0.77,1.13,1.71,1.11]],
                            ["Polarity",[0,1.48,49.7,49.9,0.35,0,51.6,0.13,49.5,0.13,1.43,3.38,1.58,3.53,52,1.67,1.66,0.13,2.1,1.61]],
                            ["Mutability",[100,44,86,77,51,50,91,103,72,54,93,104,58,84,83,117,107,98,25,50]]
                        ]);
    let outPropertyPosition = new Map();
    aaProperties.forEach(function (prop){
        outPropertyPosition.set(prop, [])
        frequencies.forEach(function (item) {
            outPropertyPosition.get(prop).push(multiplyvector(aaPropertiesData.get(prop), item));
        });
    });
    return outPropertyPosition;
}

var mapAAProps = function (aa_properties, mapping){
    const aaProperties = ["Charge","Hydropathy","Hydrophobicity","Polarity","Mutability"]
    let outPropertyMappedPosition = new Map();
    aaProperties.forEach(function (prop){
        outPropertyMappedPosition.set(prop, [])
        let currProp = aa_properties.get(prop)
        currProp.forEach(function (data, aln_ix) {
            let mappedI0 = mapping[aln_ix+1];
            if (mappedI0) {
                outPropertyMappedPosition.get(prop).push([mappedI0, Number(math.sum(data).toFixed(2))]);
            }
        });
    });
    return outPropertyMappedPosition;
}


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
        custom_aln_file: null
    },
    methods: {
        handleFileUpload(event){
            this.custom_aln_file = this.$refs.file.custom_aln_file[0];
            ajax('/custom-aln-data', optional_data=event.target.files)
        },
        cleanTreeOpts() {
            this.options = null;
            this.tax_id = null;
            cleanupOnNewAlignment(vm, "Select new alignment!");
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
        }, loadData: function(value, type_tree) {
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
                this.chains = null
                this.hide_chains = true
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase())
                    .then(struc_data => {
                        var chain_list = struc_data[pdbid.toLowerCase()];
                        if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
                        filterAvailablePolymers(chain_list, aln_id, vm);
                        this.hide_chains = null;
                    }).catch(error => {
                        alert("No such pdb id: " + pdbid + ".", error)
                    })
            }
        },
        showAlignment(aln_id, taxid, type_tree) {
            cleanupOnNewAlignment(vm, "Loading alignment...");
            if (type_tree == "orth"){
                var url = `/ortholog-aln-api/${aln_id}/${taxid}`}
            if (type_tree == "para"){
                var url = '/paralog-aln-api/'+aln_id.split(',')[1]}
            ajax(url).then(fasta => {
                this.fasta_data = fasta[0];
                this.aa_properties = calculateFrequencyData(fasta[3])
                var main_elmnt = document.querySelector(".alignment_section")
                var opts = {
                    el: document.getElementById("alnDiv"),
                    seqs: msa.io.fasta.parse(fasta[0]),
                    colorscheme: {
                        scheme: "clustal2",
                    },
                    //columns: {
                    //    hidden: fasta[2] // hidden columns
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
                    var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(data["rowPos"]) + 1)}&aln=${aln_id}${strainQuery}${fasta[1][Number(data["seqId"])]}`
                    ajax(url).then(alnpos_data => {
                        ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
                            vm.aln_meta_data = resiData;
                        });
                    }).catch(error => {
                        vm.aln_meta_data = null;
                        console.log("No residue with alignment position: " + data["rowPos"] + ". In alignment " + aln_id + ". Of species " + fasta[1][Number(data["seqId"])]);
                        console.log(error);
                    })
                });
                m.g.on("residue:mousein", function(data) {
                    console.log(data)
                });
            })
        }, showTopologyViewer (pdbid, chainid, fasta){
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(vm);}
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
                if (this.tax_id.length == 2) {
                    ajax('/twc-api/', optional_data={fasta}).then(twcDataUnmapped => {
                        mapped_aa_properties.set("TwinCons", [])
                        for (i = 0; i < twcDataUnmapped.length; i++) {
                            let mappedI0 = this.structure_mapping[twcDataUnmapped[i][0]];
                            if (mappedI0) {
                                mapped_aa_properties.get("TwinCons").push([mappedI0, twcDataUnmapped[i][1]]);
                            }
                        }
                    })
                }
                window.mapped_aa_properties = mapped_aa_properties;
                var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
                ajax(topology_url).then(data => {
                    var entityid = Object.keys(data[pdblower])[0];
                    var mapping = []
                    var range_string = minIndex.concat("-").concat(maxIndex)
                    GetRangeMapping(pdbid, chainid, range_string, mapping)
                    let data_string = JSON.stringify(Array.from(mapped_aa_properties.entries())).replaceAll(",[[", ":").replaceAll("]],",";").replaceAll("],[",",")
                    let formatted_data_string = data_string.replaceAll("[","").replaceAll("]","").replaceAll("\"","")
                    var topology_viewer = `<pdb-topology-viewer entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${formatted_data_string} filter-range=${mapping}></pdb-topology-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                })
            });
        }, showPDBViewer(pdbid, chainid){
            if (document.querySelector("pdbe-molstar")) {return;}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            console.log('PDBV');
            var coordinates_url=`https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=27&encoding=bcif`;
            var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
            console.log(coordinates_url);

            ajax(topology_url).then (data => {
                var entityid = Object.keys(data[pdblower])[0];
                var coordinates_url=`https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?${entityid}&encoding=bcif`;
                console.log(entityid);
                var mapping = []
                var range_string = minIndex.concat("-").concat(maxIndex)
                GetRangeMapping(pdbid, chainid, range_string, mapping)
                console.log(mapping)

                var PDBMolstar_viewer = `<pdbe-molstar id="PdbeMolstarComponent" molecule-id="1cbs" hide-controls="true" subscribe-events="true" ></pdbe-molstar>`
                document.getElementById('pdbeMolstarView').innerHTML = PDBMolstar_viewer;
                var PdbeMolstarComponent = document.getElementById('PdbeMolstarComponent');
                var viewerInstance2 = PdbeMolstarComponent.viewerInstance;

                viewerInstance2.visual.update({
                    customData: { url: `https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=${entityid}&encoding=bcif`, format: 'cif', binary:true },
                        hideCanvasControls: ["expand", "selection", " animation"],
                        assemblyId: '1',                    
                        hideControls: true,                   
                        subscribeEvents: true
                });

                document.addEventListener('PDB.topologyViewer.click', (e) => {
                    var pdbeMolstar=document.getElementById("PdbeMolstarComponent")
                    var molstar= pdbeMolstar.viewerInstance;                            
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
                    var pdbeMolstar=document.getElementById("PdbeMolstarComponent")
                    var molstar= pdbeMolstar.viewerInstance;                            
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
            });
        }
    }
})
