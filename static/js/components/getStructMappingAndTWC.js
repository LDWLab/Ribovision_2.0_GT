import {loadAlignmentViewer} from './loadAlignmentViewer.js'
import {ajaxProper} from './ajaxProper.js'

export function getStructMappingAndTWC (fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj){
    vm.sequence = ebi_sequence;
    if (vm.fasta_data){
        let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
        vm.fasta_data = cleanFasta;
    };
    const postData = {
        fasta,
        struc_id,
        "cif_mode_flag" : vm.user_uploaded_cif_flag
    };
    ajax('/mapSeqAln/', postData).then(structMappingAndData=>{
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

    var mapped_aa_properties = mapAAProps(vueObj.aa_properties, vueObj.structure_mapping);
    if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
        if (vueObj.unmappedTWCdata) {
            mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_aa_properties);
        }
    }
    window.mapped_aa_properties = mapped_aa_properties;

    vm.sequence=vueObj.fasta_data.split(' ')[vm.user_uploaded_cif_flag === null || vm.user_uploaded_cif_flag ? 1 : 2];
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
        viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties);
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
    //console.log("TCT_2", pdbid, vm.sequence4); 
    //vm.getR2DT(vm.sequence4);
    //vm.URL = `r2dt/${vm.sequence3}`
    //var postTopologyURL = `r2dt/${vm.sequence3}/`
    console.log('eid1',entityid);

    pdbid='cust';

    async function getRNAChain(pdbid) {
      try {
        const returnedObject = await ajax(`full-RNA-seq/${pdbid}/`);
        console.log('RNA_full_sequence_pdb', pdbid);
        const result = returnedObject["RNAseq"];
        console.log('RNA_full_sequence', result);
        return result;
      } catch (error) {
        console.error(error);
        return null;
      }
    }
    
    async function RNAseqCall(pdbid) {
      const RNA_full_sequence = await getRNAChain(pdbid);
      console.log('RNA_full_sequence2', RNA_full_sequence);
      return RNA_full_sequence;
    }
    
    RNAseqCall(pdbid).then(seq1 => {
      console.log('RNA_full_sequence3', seq1);
      vm.getR2DT(seq1);
      vm.URL = `r2dt/${seq1}/${entityid}`
      var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api="true" ></pdb-rna-viewer>` 
      document.getElementById('topview').innerHTML = topology_viewer; 
      window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
      ajaxProper({
        //url: postTopologyURL,
        url: vm.URL,
        type: 'POST',
        dataType: 'json',
        postData : {
            cif_mode_flag : vm.user_uploaded_cif_flag
        }
      }).then (parsedResponse => {
        if (parsedResponse == "Topology Success!"){
            console.log("Topology Success!");          
        }
      }).catch(error => {
        var topview = document.querySelector('#topview');
        vm.topology_loaded = 'error';
        topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
        console.log(error.responseText);
        });
    }).catch(error => {
      console.error(error);
    });
}
    //vm.getR2DT(vm.pdbSeq);
    //vm.URL = `r2dt/${vm.pdbSeq}/`
    //vm.getR2DT(seq1);
    //vm.URL = `r2dt/${seq1}/`
    
    //var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api="true" ></pdb-rna-viewer>` 
    //document.getElementById('topview').innerHTML = topology_viewer; 
    //window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
   // ajaxProper({
        //url: postTopologyURL,
        //url: vm.URL,
       //type: 'POST',
        //dataType: 'json'
    //}).then (parsedResponse => {
        //if (parsedResponse == "Topology Success!"){
            //console.log("Topology Success!");          
        //}
    //}).catch(error => {
        //var topview = document.querySelector('#topview');
        //vm.topology_loaded = 'error';
        //topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
        //console.log(error.responseText);
    //});
//}

function sleeper(ms) {
    return function(x) {
        return new Promise(resolve => setTimeout(() => resolve(x), ms));
    };
}
