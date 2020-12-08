var registerHoverResiData = function (e, tooltipObj){
  if (vm.type_tree == 'upload'){
      //Figure out how to do the hover in this case;
      return;
  }
  const strainQuery = '&res__poldata__strain__strain=';
  var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(e.position) + 1)}&aln=${vm.alnobj.id}${strainQuery}${vm.fastaSeqNames[Number(e.i)]}`
  ajax(url).then(alnpos_data => {
    var alnViewCanvasEle = document.querySelector("#alnDiv canvas:nth-of-type(1)");
    var alnViewLabelsEle = document.querySelector("#alnViewerLabels");
    let boundLabelBox = alnViewLabelsEle.getBoundingClientRect();
    let boundingBox = absolutePosition(alnViewCanvasEle);
    let relativeBox = alnViewCanvasEle.getBoundingClientRect();
    if (alnpos_data.count != 0){
        ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
            if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
              let tooltipPosition = {
                top: mousePos.y-boundingBox.top+15 +"px",
                left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+8 +"px",
              };
              if (resiData["Structural fold"][0] !== undefined && resiData["Associated data"][0] !== undefined){
                  tooltipObj.setState({
                  fold: resiData["Structural fold"][0][1],
                  phase: resiData["Associated data"][0][1],
                  tooltipPosition,
                });
              }else{
                  tooltipObj.setState({
                  fold: 'NA',
                  phase: 'NA',
                  tooltipPosition,
                });
              }
            }
            window.ajaxRun = false;
        });
    }else{
        if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
            let tooltipPosition = {
                top: mousePos.y-boundingBox.top+15 +"px",
                left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+5 +"px",
            };
            window.ajaxRun = false;
            tooltipObj.setState({
                fold: 'NA',
                phase: 'NA',
                tooltipPosition,
            });
        }
    }
  }).catch(error => {
    window.ajaxRun = false;
    console.log(error);
  })
  return true;
};

function handleMaskingRanges(mask_range){
  vm.masking_range = mask_range;
  window.masking_range_array = null;
  if (isCorrectMask(mask_range)) {   
      var topviewer = document.getElementById("PdbeTopViewer");
      topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
      if(window.custom_prop) {
          topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop); 
      }
      window.masked_array = initializeMaskedArray();          
      var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;

      var index = 4;
      while(index < topviewer.pluginInstance.domainTypes.length) {
          colorResidue(index, window.masked_array);
          index++;
      }
      let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex]
      
      if (selectedData.data){
          topviewer.pluginInstance.updateTheme(selectedData.data); 
          window.viewerInstance.visual.select({data: selectSections_RV1.get(selectedData.label), nonSelectedColor: {r:255,g:255,b:255}});
          }
      vm.correct_mask = 'True';
  } else {
      vm.correct_mask = 'False';
  }
};
function handleFilterRange(filter_range) {
  const temp_array = filter_range.split('-');
  if (filter_range.match(/^\d+-\d+/) && Number(temp_array[0]) < Number(temp_array[1])) {
      vm.filter_range = filter_range;
      window.filterRange = temp_array.join(",");
      var topviewer = document.getElementById("PdbeTopViewer");
      var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
      topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
      viewerInstance.visual.update({
          customData: {
              url: `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/residueRange?entityId=${topviewer.entityId}&range=${filter_range}&encoding=bcif`,
              format: 'cif',
              binary:true },
          assemblyId: '1',
          subscribeEvents: true,
          bgColor: {r:255,g:255,b:255},
      });
      viewerInstance.events.loadComplete.subscribe(() => { 
          let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
          if(selectSections_RV1.get(selectedData.label)) {
              let select_sections = selectSections_RV1.get(selectedData.label).slice(Number(temp_array[0]), Number(temp_array[1])+1);
              window.viewerInstance.visual.select({
              data: select_sections,
              nonSelectedColor: {r:255,g:255,b:255}});
          }
          //var selectedDomain = topviewer.pluginInstance.domainTypes[selectedIndex];
          //topviewer.updateTheme(selectedDomain.data);
       });
       topviewer.pluginInstance.initPainting(window.select_sections)
       let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
       topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
       topviewer.pluginInstance.updateTheme(selectedData.data); 
  }else{
      //
  }
};

function colorResidue(index, masked_array) {
  var topviewer = document.getElementById("PdbeTopViewer");
  var f = 0;
  while(f < topviewer.pluginInstance.domainTypes[4].data.length) {
      if(!masked_array[f] && topviewer.pluginInstance.domainTypes[index].data[f]) {
          topviewer.pluginInstance.domainTypes[index].data[f].color = "rgb(255,255,255)";
          topviewer.pluginInstance.domainTypes[index].data[f].tooltipMsg = "NaN";                   
          selectSections_RV1.get(topviewer.pluginInstance.domainTypes[index].label)[f].color = {r: 255, g: 255, b: 255};

      } if(!masked_array[f] && vm.coil_residues.includes(f) && topviewer.pluginInstance.domainTypes[index].data[f]) {
          topviewer.pluginInstance.domainTypes[index].data[f].color = "rgb(0,0,0)";
          topviewer.pluginInstance.domainTypes[index].data[f].tooltipMsg = "NaN";
      }                        
      f++;
  }
};
function cleanCustomMap(checked_customMap){
  if (checked_customMap){return;}
  var topviewer = document.getElementById("PdbeTopViewer");
  topviewer.pluginInstance.domainTypes = topviewer.pluginInstance.domainTypes.filter(obj => {return obj.label !== "CustomData"})
  window.coilsOutOfCustom = null;
  //window.custom_prop = null;
  vm.csv_data = null;
};
function handleCustomMappingData(){
  const readFile = function (fileInput) {
      var reader = new FileReader();
      reader.onload = function () {
          vm.csv_data = reader.result;
      };
      reader.readAsBinaryString(fileInput);
  };
  readFile(vm.$refs.custom_csv_file.files[0]);

};

function cleanFilter(checked_filter, masking_range){
  if (checked_filter){return;}
  if (masking_range == null){return;}
  window.masked_array = [];
  vm.masking_range = null;
  var topviewer = document.getElementById("PdbeTopViewer");
  topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
  if(window.custom_prop) {
      topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
  }
  var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
  topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
  window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
};
function cleanSelection(checked_selection, filter_range){
  if (checked_selection){return;}
  if (filter_range == null){return;}
  var topviewer = document.getElementById("PdbeTopViewer");
  var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
  vm.filter_range = null;
  window.filterRange = "-10000,10000";
  topviewer.pluginInstance.initPainting();
  viewerInstance.visual.update({
      customData: {
          url: `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/chains?entityId=${topviewer.entityId}&encoding=bcif`,
          format: 'cif',
          binary:true },
      assemblyId: '1',
      subscribeEvents: true,
      bgColor: {r:255,g:255,b:255},});
  topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
  if(window.custom_prop) {
      topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
  }
  topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
  window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
};

function parseConsecutiveIndices(structureTypeString, structureList, indicesList) {
    let structureIndex = 0;
    if (indicesList.length == 0) {
        return;
    }
    let previousIndex = indicesList[0];
    let currentStructure = [previousIndex];
    for (let i = 1; i < indicesList.length; i++) {
        let currentIndex = indicesList[i];
        if (currentIndex == previousIndex + 1) {
            currentStructure.push(currentIndex);
        } else {
            let structureObject = {};
            structureObject.text = structureTypeString + " #" + structureIndex;
            structureObject.value = structureIndex;
            structureObject.indices = currentStructure;
            structureList.push(structureObject);
            structureIndex++;
            currentStructure = [currentIndex];
        }
        previousIndex = currentIndex;
    }
}

var populatePDBs = function (alndata){
    let alnPolurl = `/desire-api/polymers/?alns_of_polymer=${alndata.id}`
    ajax(alnPolurl).then(polymersForAln => {
        let trueNom = polymersForAln.results[0].nomgd.split('/')[5]
        let url = `/desire-api/old-nomenclatures/?n_b_y_h_a=BAN&nn_fk=${trueNom}`;
        ajax(url).then(oldnomData => {
            if (oldnomData.count == 0){return;}
            let oldName = oldnomData.results[0].old_name.replace(/^(.{2})(0)/,"$1")
            let riboXYZurl = `https://ribosome.xyz:8000/neo4j/gmo_nom_class/?banName=${oldName}&format=json`
            ajax(riboXYZurl).then(data => {
                var pdb_entries = []
                data.forEach(function(entry){
                    let pdb_text = `${entry.parent} ${entry.orgname[0].slice(0,39)}`
                    pdb_entries.push({id: entry.parent.toLowerCase(), name:pdb_text})
                });
                if (pdb_entries.length == 0){return;}
                vm.pdbs.push(...pdb_entries);
            }).catch(error => {
                console.log(error);
            })
        }).catch(error => {
            console.log(error);
        })
    }).catch(error => {
            console.log(error);
    })
}

var customFilter = function (object, result, key, value){
    if(object.hasOwnProperty(key) && object[key] == value)
        result.push(object);
    for(var i=0; i<Object.keys(object).length; i++){
        let nextObj = object[Object.keys(object)[i]];
        if(typeof nextObj == "object" && nextObj != null){
            customFilter(nextObj, result, key, value);
        }
    }
}

function handlePropensities(checked_propensities){
    if (checked_propensities){
        // console.log(document.getElementById("selectaln"));
        var coilsListOfLists = []
        var strandsListOfLists = []
        var helicesListOfLists = []
        // ajax("https://www.ebi.ac.uk/pdbe/api/topology/entry/" + vm.pdbid).then(topology => {
        //     let topology_entries_map = new Map(Object.entries(topology[vm.pdbid]))
        //     let topology_entries_map_iterator = topology_entries_map.entries();
        //     let entry = true;
        //     while (entry) {
        //         entry = topology_entries_map_iterator.next().value;
        //         console.log(entry);
        //         for (let i = 1; i < entry.length(); i++) {
        //             console.log("\t" + entry[i]);
        //         }
        //     }
        // });
        parseConsecutiveIndices("Coil", coilsListOfLists, vm.coil_residues);
        parseConsecutiveIndices("Strand", strandsListOfLists, vm.strand_residues);
        parseConsecutiveIndices("Helix", helicesListOfLists, vm.helix_residues);
        vm.substructures = []
        Array.prototype.push.apply(vm.substructures, coilsListOfLists);
        Array.prototype.push.apply(vm.substructures, strandsListOfLists);
        Array.prototype.push.apply(vm.substructures, helicesListOfLists);
        let sequence_indices = ["1", "2", "3"];
        // let sequence_indices = prompt("Enter sequence indices (comma-separated): ").replace(/^\s+|\s+$/gm,'').split(',')
        for (let i = 0; i < sequence_indices.length; i++) {
            sequence_indices[i] = parseInt(sequence_indices[i])
        }
        var full = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
        let url = null;
        if (vm.structure_mapping) {
            let alignment_indices = []
            let inverse_structure_mapping = {}
            for (var key in vm.structure_mapping) {
                let value = vm.structure_mapping[key]
                inverse_structure_mapping[value] = key
            }
            for (var sequence_index of sequence_indices) {
                alignment_indices.push(inverse_structure_mapping[sequence_index])
            }
            var indices = alignment_indices.join(',')
            url = `/propensity-data/${vm.alnobj.id}/${vm.tax_id.join(',')}`
        } else {
            url = `/propensity-data/${vm.alnobj.id}/${vm.tax_id}`
            indices = '';
        }
        ajax(url, {indices}).then(data => {
            build_propensity_graph(data['amino acid'], full, vm.alnobj.text + ' ' + 'Amino Acid Propensities for ' + vm.tax_id.join(' '), 'total');
        });

    }
}

var build_propensity_graph = function (data, amino_acids, title, div) {
    // reformat data into trace format
    var data2 = {};
    var hover = {};
    // initialize arrays for each amino acid in the JavaScript object
    for (aa of amino_acids) {
        data2[aa] = [];
        hover[aa] = [];
    }
    // for a species, add the propensities for each amino acid to thei
    for (species in data) {
        for (aa of amino_acids) {
            data2[aa].push(data[species][aa])
            hover[aa].push(data[species]['name'].replace(/_/g, ' '))
        }
    };
    // build traces
    traces = []
    for (aa of amino_acids) {
        // styling - can we make this a stripplot?
        var newtrace = {
            y: data2[aa],
            type: 'box',
            boxpoints: 'all',
            jitter: 0.2,
            pointpos: 0,
            name: aa,
            text: hover[aa]};
        traces.push(newtrace)
    };
    
    var layout = {
        title: title,
        xaxis: {title: 'amino acid group'},
        yaxis: {title: 'propensity'},
        hovermode: 'closest',
        hoveron: 'points'};
    Plotly.newPlot(div, traces, layout);
}