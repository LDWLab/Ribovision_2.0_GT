import {loadAlignmentViewer} from './loadAlignmentViewer.js'
export function getStructMappingAndTWC (fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj){
    ajax('/mapSeqAln/', {fasta, struc_id}).then(struct_mapping=>{
        var largestKey = Math.max(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        var smallestKey = Math.min(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        if ((largestKey != stopIndex || smallestKey != startIndex)&&ebi_sequence){
            ajax('/mapSeqAlnOrig/', {fasta, ebi_sequence, startIndex:1}).then(orig_struct_mapping=>{
                if (struct_mapping["gapsInStruc"].length > 0){
                    struct_mapping["gapsInStruc"].forEach(function(gapTup){
                        let lowMiss = Number(_.invert(struct_mapping)[gapTup[0]]);
                        let topMiss = Number(_.invert(struct_mapping)[gapTup[1]]);
                        if (topMiss-lowMiss>1){
                            for (let i = lowMiss+1; i < topMiss; i++){
                                delete orig_struct_mapping[i];
                            }
                        }
                    })
                }
                assignColorsAndStrucMappings(vueObj, orig_struct_mapping);
            })
        } else {
            assignColorsAndStrucMappings(vueObj, struct_mapping)
        }
    }).catch(error => {
        var topview = document.querySelector('#topview');
        console.log(error);
        vueObj.topology_loaded = 'error';
        topview.innerHTML = "Failed to load the alignment-structure mapping!<br>Try another structure."
    });
}

var assignColorsAndStrucMappings = function (vueObj, struct_mapping){
    vueObj.structure_mapping = struct_mapping;
    vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];
    vueObj.fasta_data = struct_mapping["amendedAln"]
    loadAlignmentViewer (struct_mapping["amendedAln"]);
    var mapped_aa_properties = mapAAProps(vueObj.aa_properties, struct_mapping);
    if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
        if (vueObj.unmappedTWCdata) {
            mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_aa_properties);
        }
    }
    window.mapped_aa_properties = mapped_aa_properties;
    retry(delayedMapping, 10, 1000);
}

var delayedMapping = function (){
    viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
    viewerInstanceTop.pluginInstance.createDomainDropdown();
}

function retry (fn, maxAttempts = 1, delay = 0, attempts = 0) {
    return Promise.resolve()
      .then(fn)
      .catch(err => {
        if (attempts < maxAttempts) {
          return retry (fn, maxAttempts, delay, attempts + 1)
        }
        throw err
      })
  }