import { loadAlignmentViewer } from './loadAlignmentViewer.js'
import { ajaxProper } from './ajaxProper.js'


const typeMappings = {
  'aes': { '1a': '1', '2': '1', '10': '1', '13': '1', '7': '1', '17': '1', '22': '1', '0': '2', '7b': '2', '9': '2', '11': '2', '16': '2', '20': '2', '18': '2', '3a': '2', '3b': '3', '5a': '3', '4': '3', '15': '3', '24': '3', '6': '3', '23': '3', '27': '3', '1': '4', '4a': '4', '7a': '4', '8': '4', '14a': '4', '25': '4', '26': '4', '21': '4', '1b': '5', '3': '5', '5': '5', '12': '5', '14': '5', '19': '5', '20a': '5' },
  'AES': { '1.0': '1', '12.0': '1', '19.0': '1', '22.0': '1', '7.0': '1', '34.0': '1', '35.0': '1', '49.0': '1', '28.0': '1', '20.0': '1', '52.0': '1', '58.0': '1', '37.0': '1', '53': '2', '47': '3', '8': '2', '54': '2', '10': '2', '16': '2', '32': '2', '5': '2', '12a': '2', '18': '2', '50': '2', '4a': '3', '21': '3', '4': '3', '26': '3', '14': '3', '25': '3', '33': '3', '38': '3', '11': '3', '51': '3', '39': '3', '40': '3', '6': '4', '6a': '4', '55': '4', '45': '4', '10a': '4', '23': '4', '31': '4', '56': '4', '15': '4', '2': '4', '43': '4', '46': '4', '30': '4', '59': '4', '41': '5', '24': '5', '17': '5', '48': '5', '9': '5', '32a': '5', '57': '5', '15a': '5', '29': '5', '44': '5', '36': '5', '3': '5' },
  'helix': { '2': '1', '6': '1', '9_3a': '1', '9_3b': '1', '13': '1', '17': '1', '21_6d': '1', '24': '1', '26a': '1', '33': '1', '36': '1', '40': '1', '42': '1', '44': '1', '41_10': '1', '3': '2', '6a': '2', '10': '2', '11': '2', '14': '2', '15': '2', '21_6b': '2', '22': '2', '27': '2', '28': '2', '31': '2', '32': '2', '37': '2', '38': '2', '41': '2', '45': '2', '4': '3', '7': '3', '9': '3', '19': '3', '21': '3', '21_3c': '3', '23': '3', '26': '3', '29': '3', '34': '3', '39_9': '3', '44_12': '3', '1': '4', '5': '4', '8': '4', '12': '4', '16': '4', '18': '4', '20': '4', '21_6a': '4', '23a': '4', '25': '4', '26_7': '4', '30': '4', '35': '4', '39': '4', '43': '4' },
  'Helix': { '5S1': '1', '5S2': '3', '5S3': '2', '5S5': '2', '5S4': '4', '5': '1', '9': '1', '13': '1', '18': '1', '21': '1', '23': '1', '25': '1', '25_7b': '1', '27': '1', '30': '1', '35a': '1', '36': '1', '38a': '1', '41': '1', '43a': '1', '47': '1', '49': '1', '51': '1', '55': '1', '58': '1', '64': '1', '78': '1', '79': '1', '74': '1', '82': '1', '86': '1', '90': '1', '94': '1', '98_39b': '1', '63_27': '1', '79_31a': '1', '2': '2', '6': '2', '10': '2', '14': '2', '19': '2', '25_7a': '2', '28': '2', '33': '2', '38': '2', '40': '2', '25a': '2', '42': '2', '52': '2', '49a': '2', '54_20a': '2', '59': '2', '63': '2', '61': '2', '67': '2', '75': '2', '79_31': '2', '83': '2', '88': '2', '92': '2', '97': '2', '99': '2', '31_9': '2', '3': '3', '7': '3', '11': '3', '15': '3', '19a': '3', '22': '3', '25_7d': '3', '26': '3', '32': '3', '34': '3', '39': '3', '44': '3', '45': '3', '48': '3', '49b': '3', '53': '3', '56': '3', '63a': '3', '63_27b': '4', '66': '3', '69': '3', '70??': '3', '76': '3', '79_31c': '3', '80': '3', '84': '3', '91': '3', '93': '3', '95': '3', '98': '3', '101': '3', '9_3': '3', '4': '4', '8': '4', '10_4': '4', '12': '4', '16': '4', '20': '4', '24': '4', '25_7c': '4', '31': '4', '35': '4', '37': '4', '38_12': '4', '43': '4', '46': '4', '50': '4', '52_19': '4', '54': '4', '57': '4', '60': '4', '62': '4', '68': '4', '73': '4', '77': '4', '79_31b': '4', '81': '4', '85': '4', '87': '4', '89': '4', '96': '4', '100': '4', '98_39a': '4' }
}
async function getBanName(pdbId, PchainId) {
  try {
    const apiUrl = `https://api.ribosome.xyz/neo4j/get_banclass_for_chain/?pdbid=${pdbId}&auth_asym_id=${PchainId}&format=json`
    return await (await fetch(apiUrl)).json();
  } catch (e) {
    //console.log(`Ban naming is not available!`, e);
    return void 0;
  };
}

function sleeper(ms) {
  return function (x) {
    return new Promise(resolve => setTimeout(() => resolve(x), ms));
  };
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function colorStructure(fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj, structMappingAndData, full_sequence_from_pdb) {
  const fix_colors = require('./graphColorPrediction.js');
  
  let struct_mapping = structMappingAndData["structureMapping"];

  vm.struct_to_alignment_mapping = Object.fromEntries(Object.entries(struct_mapping).map(([key, value]) => [value, key]));
  let associatedDataMappedPerType3D = {};

  while((Object.keys(vm.associatedDataCache).length == 0) && (vm.cifPdbMode == null)){
    await sleep(5000);
    console.log("waiting on associatedDataCache");
  }
  let associatedDataCache = vm.associatedDataCache;  

  // for (let [alignmentIndexAsString, structureIndex] of Object.entries(struct_mapping)) {
  //   let alignmentIndex = Number.parseInt(alignmentIndexAsString);
  //   // todo fix this else block, it goes there and we can't see the helix
  //   if (alignmentIndex in associatedDataCache) {
  //     for (let { type, value } of associatedDataCache[alignmentIndex]) {
  //       if (!(type in associatedDataMappedPerType3D)) {
  //         associatedDataMappedPerType3D[type] = [];
  //       }

  //       if (type in typeMappings) {
  //         let selectedDataDict = typeMappings[type] || {};
  //         value = selectedDataDict[value] || 0;
  //       } else {
  //         if (value.length === 0) {
  //           value = "0";
  //         }
  //         value = Number.parseInt(value);
  //       }
  //       let associatedDataI = [
  //         structureIndex,
  //         value
  //       ];
  //       associatedDataMappedPerType3D[type].push(associatedDataI);
  //     }
  //   }
  // }
  
  // vm.AD_headers = [];
  // vm.associatedDataMappedPerType_3D = associatedDataMappedPerType3D;
  // vm.associatedDataMappedPerType_3D = fix_colors(
  //   viewerInstanceTop.viewInstance.uiTemplateService.apiData.sequence,
  //   viewerInstanceTop.viewInstance.uiTemplateService.baseStrs.get('cWW')[1],
  //   associatedDataMappedPerType3D
  // );
  
  var largestKey = Math.max(...Object.values(struct_mapping).filter(a => typeof (a) == "number"))
  var smallestKey = Math.min(...Object.values(struct_mapping).filter(a => typeof (a) == "number"))
  if ((largestKey != stopIndex || smallestKey != startIndex) && ebi_sequence) {
    ajax('/mapSeqAlnOrig/', { fasta, ebi_sequence, startIndex: 1 }).then(origStructMappingAndData => {
      var orig_struct_mapping = origStructMappingAndData["structureMapping"];
      
      vm.st_mapping2D = orig_struct_mapping;
      vm.st_mapping3D = struct_mapping;

      vm.struct_to_alignment_mapping = Object.fromEntries(Object.entries(orig_struct_mapping).map(([key, value]) => [value, key]));
      let associatedDataMappedPerType2D = {};
      
      for (let [alignmentIndexAsString, structureIndex] of Object.entries(orig_struct_mapping)) {
        let alignmentIndex = Number.parseInt(alignmentIndexAsString);
        // associatedDataCache has correct boundaries
        // let associatedDataCache = vm.associatedDataCache;
        if (alignmentIndex in associatedDataCache) {
          for (let { type, value } of associatedDataCache[alignmentIndex]) {
            if (!(type in associatedDataMappedPerType2D)) {
              associatedDataMappedPerType2D[type] = [];
            }

            if (type in typeMappings) {
              let selectedDataDict = typeMappings[type] || {};
              value = selectedDataDict[value] || 0;
            } else {
              if (value.length === 0) {
                value = "0";
              }
              value = Number.parseInt(value);
            }
            let associatedDataI = [
              structureIndex,
              value
            ];
            associatedDataMappedPerType2D[type].push(associatedDataI);
          }
        }
      }
      // todo: insert color patching function here.  
      // associatedDataMappedPerType.helix[0]=[6, '1']
      
      vm.AD_headers = [];
      // console.log(associatedDataMappedPerType);
      // vm.associatedDataMappedPerType = associatedDataMappedPerType;
      // console.log('associatedDataMappedPerType2D', associatedDataMappedPerType2D);
      vm.associatedDataMappedPerType_2D = fix_colors(
        viewerInstanceTop.viewInstance.uiTemplateService.apiData.sequence,
        viewerInstanceTop.viewInstance.uiTemplateService.baseStrs.get('cWW')[1],
        associatedDataMappedPerType2D
      );

      // map 3d and 2d to be the same 
      vm.mapping3D_2D = {};

      for (let [k,v] of Object.entries(vm.st_mapping3D)){
          if (Object.keys(vm.st_mapping2D).includes(k)){
            vm.mapping3D_2D[v] = vm.st_mapping2D[k];
          }
      }

      vm.colorsMap2D = {};
      for(let [type, list] of Object.entries(vm.associatedDataMappedPerType_2D)){
        let colorMap = {};
        
        for(let [k, v] of Object.values(list)){ 
          colorMap[k] = v;
        }
        vm.colorsMap2D[type] = colorMap;
      
      }

      vm.associatedDataMappedPerType_3D = {};
      for(let [type, list] of Object.entries(vm.associatedDataMappedPerType_2D)){
        let colorMap = vm.colorsMap2D[type];
        let newColors = [];
        
        console.log("Doing: ", type);

        for(let k of Object.values(vm.st_mapping3D)){
          newColors.push([k, colorMap[vm.mapping3D_2D[k]]])
        }

        vm.associatedDataMappedPerType_3D[type] = newColors;

      }

      // vm.associatedDataMappedPerType_3D = fix_colors(
      //   viewerInstanceTop.viewInstance.uiTemplateService.apiData.sequence,
      //   viewerInstanceTop.viewInstance.uiTemplateService.baseStrs.get('cWW')[1],
      //   vm.associatedDataMappedPerType_3D
      // );


      

      if (structMappingAndData["gapsInStruc"] && structMappingAndData["gapsInStruc"].length > 0) {
        structMappingAndData["gapsInStruc"].forEach(function (gapTup) {
          let lowMiss = Number(_.invert(struct_mapping)[gapTup[0]]);
          let topMiss = Number(_.invert(struct_mapping)[gapTup[1]]);
          if (topMiss - lowMiss > 1) {
            for (let i = lowMiss + 1; i < topMiss; i++) {
              delete orig_struct_mapping[i];
            }
          }
        })
      }
      assignColorsAndStrucMappings(vueObj, origStructMappingAndData, structMappingAndData);
    })
  } else {
    assignColorsAndStrucMappings(vueObj, structMappingAndData,null);

  }
}

export function getStructMappingAndTWC(fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj, full_sequence_from_pdb = "") {
  
  vm.structFailed = false
  vm.sequence = ebi_sequence;

  if (vm.fasta_data) {
    let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
    vm.fasta_data = cleanFasta;
  };

  const postData = {
    fasta,
    struc_id,
    "cif_mode_flag": vm.user_uploaded_cif_flag,
    //hardcoded_structure: full_sequence_from_pdb
    hardcoded_structure: vm.customFullSequence
  };
  
  ajax('/mapSeqAln/', postData).then(x => {colorStructure(fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj, x, full_sequence_from_pdb)}).catch(error => {
    vueObj.topology_loaded = 'error';
    console.log(error);
    var topview = document.querySelector('#topview');
    topview.innerHTML = "Failed to load the alignment-structure mapping!<br>Try another structure."
  });
  
}

var assignColorsAndStrucMappings = function (vueObj, struct_mapping, struct_mapping3D) {
  vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];
  vueObj.fasta_data = struct_mapping["amendedAln"];
  vueObj.structure_mapping = struct_mapping["structureMapping"];

  //vueObj.poor_structure_map3D = struct_mapping3D['BadMappingPositions'];
  //vueObj.fasta_data3D = struct_mapping3D["amendedAln"];
  if(struct_mapping3D) {
    vueObj.structure_mapping3D = struct_mapping3D["structureMapping"];
  }

  loadAlignmentViewer(vueObj.fasta_data);
  //loadAlignmentViewer(vueObj.fasta_data3D);

  let mapped_aa_properties = mapAAProps(vueObj.aa_properties, vueObj.structure_mapping);
  let mapped_aa_properties3D = null
  if(struct_mapping3D) {
    mapped_aa_properties3D = mapAAProps(vueObj.aa_properties, vueObj.structure_mapping3D);
}
  
  if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
    if (vueObj.unmappedTWCdata) {
      mapTWCdata(vueObj.structure_mapping, vueObj.structure_mapping3D, vueObj.unmappedTWCdata, mapped_aa_properties, mapped_aa_properties3D);
      // mapTWCdata(vueObj.structure_mapping3D, vueObj.unmappedTWCdata, mapped_aa_properties3D);
    } else {
      window.mapped_aa_properties = mapped_aa_properties;
      window.mapped_aa_properties3D = mapped_aa_properties3D;
    }
  } else {
      window.mapped_aa_properties = mapped_aa_properties;
      window.mapped_aa_properties3D = mapped_aa_properties3D;
  }
  //console.log(vueObj.fasta_data);
  vm.sequence = vueObj.fasta_data.split(' ')[vm.user_uploaded_cif_flag === null || vm.user_uploaded_cif_flag ? 1 : 2];
  let sequence2 = vm.sequence.replaceAll(/-|\n/g, "");
  vm.sequence3 = sequence2.substring("sequence".length, sequence2.indexOf(">"));
  vm.sequence4 = vm.sequence3;
  //tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
  delayedMapping();
  //retry(delayedMapping, 10, 1000);
}

var delayedMapping = function () {
  //console.log("delayed mapping")
  if (typeof viewerInstanceTop === 'undefined' || viewerInstanceTop === null || vm.structFailed) {
    tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
  } else {
    if (vm.type_tree != "upload") {
      if (vm.topology_loaded) {
        try {
          var topviewer = document.getElementById("PdbeTopViewer");
          
          for (const [type, _] of Object.entries(vm.associatedDataMappedPerType_2D)) {
            let associatedDataMappedPerTypeI_2D = vm.associatedDataMappedPerType_2D[type];
            let associatedDataMappedPerTypeI_3D = vm.associatedDataMappedPerType_3D[type];

            associatedDataMappedPerTypeI_2D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            associatedDataMappedPerTypeI_3D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            const AD_header = type;
            // const ADDataArray = associatedDataMappedPerTypeI;
            vm.AD_headers.push(AD_header);
            mapAssociatedData(associatedDataMappedPerTypeI_2D, associatedDataMappedPerTypeI_3D, AD_header, topviewer);
          }
        } catch (error) {
          console.log(error)
          console.log("Mapping associated data failed")
        }
      }
      //viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties);  

      else {
        // console.log('mapped_aa_properties', mapped_aa_properties);
        viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties, mapped_aa_properties3D);
        // viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties3D);
        setTimeout(delayedMapping, 500);
      }
    }
  }
}

function retry(fn, maxAttempts = 1, delay = 0, attempts = 0) {
  return Promise.resolve()
    .then(sleeper(delay)).then(fn)
    .catch(err => {
      if (attempts < maxAttempts) {
        return retry(fn, maxAttempts, delay, attempts + 1)
      }
      var topview = document.querySelector('#topview');
      vm.topology_loaded = 'error';
      topview.innerHTML = "EBI topology diagram is taking too long!<br>Trying to generate topology from custom mode..."
      tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
      throw err
    })
}

var tryCustomTopology = function (pdbid, entityid, chainid) {
  vm.topology_loaded = false;
  //console.log("TCT_2", pdbid, vm.sequence4); 
  //vm.getR2DT(vm.sequence4);
  //vm.URL = `r2dt/${vm.sequence3}`
  //var postTopologyURL = `r2dt/${vm.sequence3}/`
  //console.log('eid1',entityid);

  //pdbid='cust';
  if (pdbid == '' || pdbid == 'cust') {
    pdbid = 'cust'
  }

  async function getRNAChain(pdbid) {
    try {
      const returnedObject = await ajax(`full-RNA-seq/${pdbid}/${chainid}`);
      //console.log('RNA_full_sequence_cif', pdbid);
      const result = returnedObject["RNAseq"];

      //console.log('RNA_full_sequence', result);
      return result;
    } catch (error) {
      console.error(error);
      return null;
    }
  }


  async function RNAseqCall(pdbid) {
    if (vm.user_uploaded_cif_flag == true) {
      const RNA_full_sequence = await getRNAChain(pdbid);
      //console.log('RNA_full_sequence2c', RNA_full_sequence);
      return RNA_full_sequence;
    }
    //console.log('cf',vm.user_uploaded_cif_flag)
    if (vm.user_uploaded_cif_flag == false) {
      const RNA_full_sequence = vm.customFullSequence;
      //console.log('RNA_full_sequence2p', vm.customFullSequence, RNA_full_sequence);
      return RNA_full_sequence;
    }

    //return RNA_full_sequence;
    //return vm.customFullSequence;
  }

  RNAseqCall(pdbid).then(seq1 => {
    //
    if (pdbid == "cust") {
      vm.sequence_for_r2dt = seq1
    }
    else
      vm.sequence_for_r2dt = vm.sequence3;
    //console.log('RNA_full_sequence3', seq1);
    vm.getR2DT(seq1);
    vm.URL = `r2dt/${entityid}`
    var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api="true"></pdb-rna-viewer>`
    document.getElementById('topview').innerHTML = topology_viewer;
    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");

    function success(parsedResponse) {
      if (parsedResponse == "Topology Success!") {
        //console.log("Topology Success!");          
      }
      // const banName = getBanName(pdbid, 'H')
      vm.json_structures_from_r2dt = parsedResponse;
      window.viewerInstanceTop.viewInstance.uiTemplateService.render(parsedResponse.RNA_2D_json, parsedResponse.RNA_BP_json, parsedResponse.RNA_BP_json, undefined, window.viewerInstanceTop.viewInstance);
      if (vm.structFailed) {
        vm.AD_headers = [];
        var topviewer = document.getElementById("PdbeTopViewer");
        try {
          
          for (const [type, _] of Object.entries(vm.associatedDataMappedPerType_2D)) {
            let associatedDataMappedPerTypeI_2D = vm.associatedDataMappedPerType_2D[type];
            let associatedDataMappedPerTypeI_3D = vm.associatedDataMappedPerType_3D[type];

            associatedDataMappedPerTypeI_2D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            associatedDataMappedPerTypeI_3D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            const AD_header = type;
            // const ADDataArray = associatedDataMappedPerTypeI;
            vm.AD_headers.push(AD_header);
            mapAssociatedData(associatedDataMappedPerTypeI_2D, associatedDataMappedPerTypeI_3D, AD_header, topviewer);
          }
        } catch (error) {
          console.log("Mapping associated data failed")
        }
      }
    }
    function handle_error(error) {
      var topview = document.querySelector('#topview');
      vm.topology_loaded = 'error';
      topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
      console.log(error.responseText);
    }
    call_r2dt("POST", success, handle_error);
  }).catch(error => {
    console.error(error);
  });
}

export function call_r2dt(request_method, success = function () {/* Do nothing. */ }, handle_error = function () {/* Do nothing. */ }) {
  let keys_as_string = JSON.stringify({
    cif_mode_flag: vm.user_uploaded_cif_flag,
    cif_file_path: vm.cif_file_path
  });
  let all_lines = [
    keys_as_string,
    "\n",
    ...vm.sequence_for_r2dt.split("\n")
  ];
  let sequence_file = new File(all_lines, "my_sequence.txt", {
    type: "text/plain"
  });
  const formData = new FormData();
  formData.append("custom_seq_file", sequence_file);
  $.ajax({
    url: vm.URL,
    data: formData,
    cache: false,
    contentType: false,
    processData: false,
    method: request_method,
    dataType: 'json',
    type: request_method, // For jQuery < 1.9
    success,
    handle_error
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

