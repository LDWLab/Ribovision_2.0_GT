import {loadAlignmentViewer} from './loadAlignmentViewer.js'
import {ajaxProper} from './ajaxProper.js'

export function getStructMappingAndTWC (fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj){
    vm.sequence = ebi_sequence;
    if (vm.fasta_data){
        let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
        vm.fasta_data = cleanFasta;
    };
    ajax('/mapSeqAln/', {fasta, struc_id}).then(structMappingAndData=>{
        var struct_mapping = structMappingAndData["structureMapping"];
        var largestKey = Math.max(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        var smallestKey = Math.min(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        if ((largestKey != stopIndex || smallestKey != startIndex)&&ebi_sequence){
            ajax('/mapSeqAlnOrig/', {fasta, ebi_sequence, startIndex:1}).then(origStructMappingAndData=>{
                var orig_struct_mapping = origStructMappingAndData["structureMapping"];
                if (structMappingAndData["gapsInStruc"]&&structMappingAndData["gapsInStruc"].length > 0){
                    structMappingAndData["gapsInStruc"].forEach(function(gapTup){
                        let lowMiss = Number(_.invert(struct_mapping)[gapTup[0]]);
                        let topMiss = Number(_.invert(struct_mapping)[gapTup[1]]);
                        if (topMiss-lowMiss>1){
                            for (let i = lowMiss+1; i < topMiss; i++){
                                delete orig_struct_mapping[i];
                            }
                        }
                    })
                }
                assignColorsAndStrucMappings(vueObj, origStructMappingAndData);
            })
        } else {
            assignColorsAndStrucMappings(vueObj, structMappingAndData);
            
        }
    }).catch(error => {
        vueObj.topology_loaded = 'error';
        console.log(error);
        var topview = document.querySelector('#topview');
        topview.innerHTML = "Failed to load the alignment-structure mapping!<br>Try another structure."
    });
}

var assignColorsAndStrucMappings = function (vueObj, struct_mapping){
    vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];
    vueObj.fasta_data = struct_mapping["amendedAln"];
    vueObj.structure_mapping = struct_mapping["structureMapping"];
    loadAlignmentViewer (vueObj.fasta_data);

    var mapped_n_properties = mapNProps(vueObj.n_properties, vueObj.structure_mapping);
    if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
        if (vueObj.unmappedTWCdata) {
            mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_n_properties);
        }
    }
    window.mapped_n_properties = mapped_n_properties;

    vm.sequence=vueObj.fasta_data.split(' ')[1];
    let sequence2=vm.sequence.replaceAll(/-|\n/g, "");
    vm.sequence3=sequence2.substring("sequence".length, sequence2.indexOf(">"));
    vm.sequence4=vm.sequence3;
    //tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
    delayedMapping();
    //retry(delayedMapping, 10, 1000);
}

var delayedMapping = function (){
    
    if ( typeof viewerInstanceTop === 'undefined' || viewerInstanceTop === null ){
        tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
    } else {
        viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_n_properties);
        }
    
}

function retry (fn, maxAttempts = 1, delay = 0, attempts = 0) {
    return Promise.resolve()
      .then(sleeper(delay)).then(fn)
      .catch(err => {
        if (attempts < maxAttempts) {
          return retry (fn, maxAttempts, delay, attempts + 1)
        }
        var topview = document.querySelector('#topview');
        vm.topology_loaded = 'error';
        topview.innerHTML = "EBI topology diagram is taking too long!<br>Trying to generate topology from custom mode..."
        tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
        throw err
      })
  }

var tryCustomTopology = function (pdbid, entityid, chainid){
    vm.topology_loaded = false;

    vm.getR2DT(vm.sequence4);
    vm.URL = `r2dt/${vm.sequence3}`
    var postTopologyURL = `r2dt/${vm.sequence3}/`
    pdbid='cust'; 
    var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api=true subscribe-events=true></pdb-rna-viewer>`
    document.getElementById('topview').innerHTML = topology_viewer;
    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
    ajaxProper({
        url: postTopologyURL,
        type: 'POST',
        dataType: 'json'
    }).then (parsedResponse => {
        if (parsedResponse == "Topology Success!"){
            //var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid} pvapi="true" filter-range=1,100000></pdb-topology-viewer>`
            //document.getElementById('topview').innerHTML = topology_viewer;
            //window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
            console.log("Topology Success");

            
        }
    }).catch(error => {
        var topview = document.querySelector('#topview');
        vm.topology_loaded = 'error';
        topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
        console.log(error.responseText);
    });
}

function sleeper(ms) {
    return function(x) {
        return new Promise(resolve => setTimeout(() => resolve(x), ms));
    };
}

