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

function isCorrectMask(mask_range){
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
    return isCorrect;
  };

function initializeMaskedArray() {
    var topviewer = document.getElementById("PdbeTopViewer");
    var masked_array = [];
    var j = 0;
    while(j < mapped_aa_properties.get(topviewer.pluginInstance.domainTypes[4].label).length) {
        masked_array[j] = false;
        var i = 0;
        while(i < window.masking_range_array.length && !masked_array[j]) {
            if(j >= window.masking_range_array[i] && j <= window.masking_range_array[i + 1]) {
                masked_array[j] = true;
            }
            i = i+2;
        }
        j = j+1;
    }
    return masked_array;
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

      var index = 1;
      while(index < topviewer.pluginInstance.domainTypes.length) {
          colorResidue(index, window.masked_array);
          index++;
      }
      let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex]
      
      if (selectedData.data){
          topviewer.pluginInstance.updateTheme(selectedData.data); 
          window.viewerInstance.visual.select({data: selectSections_RV1.get(selectedData.label), nonSelectedColor: {r:255,g:255,b:255}});
          }
      vm.correct_mask = true;
  } else {
      vm.correct_mask = false;
  }
};
function handleDomainRange(domain_range) {
    //handleFilterRange(domain_range);
    domain_array = domain_range.split(';');
    if(domain_array.length == 2) {
        vm.masking_range = null;
        vm.checked_filter = false;
        handleFilterRange(domain_range);
    } else {
        var first = domain_array[0].split('-')[0];
        var last = domain_array[domain_array.length - 2].split('-')[1];
        var full_range = first + "-" + last + ";";
        vm.checked_filter = true;
        vm.handleMaskingRanges(domain_range);
        handleFilterRange(full_range);
    }
}
function handleFilterRange(filter_range) {
    if (filter_range.match(/^\d+-\d+;/)) {
        handlePropensities(vm.checked_propensities);
        var filter_range = filter_range.slice(0, -1);
        const temp_array = filter_range.split('-');
        if (Number(temp_array[0]) < Number(temp_array[1])){
            window.filterRange = temp_array.join(",");
            var topviewer = document.getElementById("PdbeTopViewer");
            var selectBoxOut = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
            var selectedIndexOut = indexMatchingText(selectBoxOut.options, vm.selected_property);
            var coordURL = `https://coords.litemol.org/${vm.pdbid.toLowerCase()}/residueRange?entityId=${topviewer.entityId}&authAsymId=${topviewer.chainId}&range=${filter_range}&encoding=bcif`
            //var coordURL = `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/residueRange?entityId=${topviewer.entityId}&range=${filter_range}&encoding=bcif`
            topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
            viewerInstance.visual.update({
                customData: {
                    url: coordURL,
                    format: 'cif',
                    binary:true },
                assemblyId: '1',
                subscribeEvents: true,
                bgColor: {r:255,g:255,b:255},
            });
            viewerInstance.events.loadComplete.subscribe(() => { 
                if(!vm.selected_property){return;}
                let rangeArr = window.filterRange.split(',');
                let selectBox = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
                let selectedIndex = indexMatchingText(selectBox.options, vm.selected_property);
                let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
                if(selectSections_RV1.get(selectedData.label)) {
                    var select_sections = selectSections_RV1.get(selectedData.label).filter(resi3D  => {
                        if (resi3D.start_residue_number >= Number(rangeArr[0]) && resi3D.start_residue_number <= Number(rangeArr[1])){
                            return resi3D;
                        }
                    })
                    window.viewerInstance.visual.select({
                    data: select_sections,
                    nonSelectedColor: {r:255,g:255,b:255}});
                }
                if (selectedIndex > 0){
                    var selectedDomain = topviewer.pluginInstance.domainTypes[selectedIndex];
                    topviewer.pluginInstance.updateTheme(selectedDomain.data);
                }
                selectBox.selectedIndex = selectedIndex;
            });
            topviewer.pluginInstance.alreadyRan = false;
            topviewer.pluginInstance.initPainting(window.select_sections)
            let selectedData = topviewer.pluginInstance.domainTypes[selectedIndexOut];
            topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
            if(selectedIndexOut > 0) {
                topviewer.pluginInstance.updateTheme(selectedData.data);
            }
            if(vm.correct_mask){
                handleMaskingRanges(vm.masking_range)
            }
        }else{
            //Swapped start end
        }
    }else{
        //Incorrect syntax
    }
};

function colorResidue(index, masked_array) {
    viewerInstanceTop.pluginInstance.domainTypes[index].data.forEach(function(resiEntry){
        if (!masked_array[resiEntry.start]){
            resiEntry.color = "rgb(255,255,255)";
            resiEntry.tooltipMsg = "NaN";
        } 
        if (!masked_array[resiEntry.start] && vm.coil_residues.includes(resiEntry.start)){
            resiEntry.color = "rgb(0,0,0)";
            resiEntry.tooltipMsg = "NaN";
        }
    })
    selectSections_RV1.get(viewerInstanceTop.pluginInstance.domainTypes[index].label).forEach(function(resiEntry){
        if (!masked_array[resiEntry.start_residue_number]){
            resiEntry.color = {r: 255, g: 255, b: 255};
        }
    })
};
function clearInputFile(f){
    if(f.value){
        try{
            f.value = ''; //for IE11, latest Chrome/Firefox/Opera...
        }catch(err){ }
        if(f.value){ //for IE5 ~ IE10
            var form = document.createElement('form'),
                parentNode = f.parentNode, ref = f.nextSibling;
            form.appendChild(f);
            form.reset();
            parentNode.insertBefore(f,ref);
        }
    }
}

function cleanCustomMap(checked_customMap){
    if (vm.uploadSession){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    if (!topviewer || !topviewer.pluginInstance.domainTypes){
        if (checked_customMap){return;}
        var sliceAvailProp = Array.prototype.slice.call(vm.available_properties).filter(availProp => {
            return vm.custom_headers.includes(availProp.Name)
        })
        const setSlice = new Set(sliceAvailProp.map(a=>{return a.Name}));
        const newArray = vm.available_properties.filter(obj => !setSlice.has(obj.Name));
        vm.available_properties = newArray;
        return;
    }
    var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    topviewer.pluginInstance.domainTypes = topviewer.pluginInstance.domainTypes.filter(obj => {
        return !vm.custom_headers.includes(obj.label)
    })
    
    var sliceChildren = Array.prototype.slice.call(selectBoxEle.childNodes).filter(optionsNode => {
        return vm.custom_headers.includes(optionsNode.label)
    })
    
    sliceChildren.forEach(function(){
        selectBoxEle.removeChild(selectBoxEle.childNodes[selectBoxEle.options.length-1]);
        vm.available_properties.splice(-1,1)
    })

    if (checked_customMap){return;}
    window.coilsOutOfCustom = null;
    window.custom_prop = null;
    vm.csv_data = null;
    vm.custom_headers = [];
};
function handleCustomMappingData(){
  const readFile = function (fileInput) {
      var reader = new FileReader();
      reader.onload = function () {
          vm.csv_data = reader.result.replace("\u00EF\u00BB\u00BF", '');
      };
      reader.readAsBinaryString(fileInput);
  };
  readFile(vm.$refs.custom_csv_file.files[0]);
};

var displayMappingDataByIndex = function(topviewer, selectedIndex){
    var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    topviewer.pluginInstance.resetTheme();
    topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data);
    window.viewerInstance.visual.select({
        data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), 
        nonSelectedColor: {r:255,g:255,b:255}
    });
    selectBoxEle.selectedIndex = selectedIndex;
    vm.selected_property = topviewer.pluginInstance.domainTypes[selectedIndex].label;
}

var mapCustomMappingData = function(custom_data, custom_data_name, topviewer){
    var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    let vals = custom_data.map(function(v){ return v[1] });
    let indexes = custom_data.map(function(v){ return v[0] });
    window.aaColorData.set(custom_data_name, [viridis]);
    window.aaPropertyConstants.set(custom_data_name, [Math.min(...vals), Math.max(...vals)]);
    let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
    window.coilsOutOfCustom = coilsOutOfCustom;
    var custom_prop = new Map();
    custom_prop.set(custom_data_name, custom_data);
    if (window.custom_prop){
        window.custom_prop.set(custom_data_name, custom_data)
    } else {
        window.custom_prop = custom_prop;
    }
    topviewer.pluginInstance.getAnnotationFromRibovision(custom_prop);
    var custom_option = document.createElement("option");
    custom_option.setAttribute("value", selectBoxEle.options.length);
    custom_option.appendChild(document.createTextNode(custom_data_name));
    selectBoxEle.appendChild(custom_option);
    if (!vm.available_properties.some(prop => prop.Name === custom_data_name)){
        vm.available_properties.push({Name:custom_data_name, url:"static/alignments/svg/Custom.svg"})
    }
    if(vm.correct_mask) {
        var j = topviewer.pluginInstance.domainTypes.length-1;
        colorResidue(j, window.masked_array);
    }
}

var getExampleFile = function(url, name){
    $.ajax({
        url: url,
        type: 'GET',
        dataType: "text",
        success: function(data) {
            let anchor = document.createElement('a');
            anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(data);
            anchor.target = '_blank';
            anchor.download = name;
            anchor.click();
        },
    })
};

function cleanFilter(checked_filter, masking_range){
  if (checked_filter){return;}
  if (masking_range == null){return;}
  window.masked_array = [];
  vm.masking_range = null;
  vm.correct_mask = null;
  var topviewer = document.getElementById("PdbeTopViewer");
  topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
  if(window.custom_prop) {
      topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
  }
  var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
  if (selectedIndex > 0){
    topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
  }
  window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
};
function cleanSelection(checked_selection, filter_range){
  if (checked_selection || filter_range == null || !vm.pdbid){return;}
  var selectBox = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
  var newIndex = indexMatchingText(selectBox.options, vm.selected_property);
  vm.filter_range = null;
  window.filterRange = "-10000,10000";
  viewerInstanceTop.pluginInstance.alreadyRan = false;
  viewerInstanceTop.pluginInstance.initPainting();
  var coordURL = `https://coords.litemol.org/${vm.pdbid.toLowerCase()}/chains?entityId=${viewerInstanceTop.entityId}&authAsymId=${viewerInstanceTop.chainId}&encoding=bcif`;
  //var coordURL = `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/chains?entityId=${topviewer.entityId}&encoding=bcif`;
  viewerInstance.visual.update({
      customData: {
          url: coordURL,
          format: 'cif',
          binary:true },
      assemblyId: '1',
      subscribeEvents: true,
      bgColor: {r:255,g:255,b:255},
    }).finally(response => {
        viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
        if(window.custom_prop) {
            viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
        }
        if(newIndex > 0) {
            viewerInstanceTop.pluginInstance.updateTheme(viewerInstanceTop.pluginInstance.domainTypes[newIndex].data); 
        }
        if (response){
            window.viewerInstance.visual.select({data: selectSections_RV1.get(vm.selected_property), nonSelectedColor: {r:255,g:255,b:255}});
        }
        if(vm.correct_mask) {
            handleMaskingRanges(vm.masking_range)
        }
          handlePropensities(vm.checked_propensities);
    });
};

var populatePDBs = function (alndata){
    if (alndata != null){
        let alnPolurl = `/desire-api/polymers/?alns_of_polymer=${alndata.id}`
        ajax(alnPolurl).then(polymersForAln => {
            let trueNom = polymersForAln.results[0].nomgd.split('/')[5];
            var polNames = polymersForAln.results.map(entry => entry.genedescription.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,''));
            let url = `/desire-api/old-nomenclatures/?n_b_y_h_a=BAN&nn_fk=${trueNom}`;
            ajax(url).then(oldnomData => {
                if (oldnomData.count == 0){return;}
                let oldName = oldnomData.results[0].old_name.replace(/^(.{2})(0)/,"$1")
                let riboXYZurl = `https://api.ribosome.xyz/neo4j/gmo_nom_class/?banName=${oldName}&format=json`
                ajax(riboXYZurl).then(data => {
                    var pdb_entries = []
                    data.forEach(function(entry){
                        let pdb_text = `${entry.parent_rcsb_id} ${entry.rcsb_source_organism_description[0]}`
                        let pdbxDescription = entry.rcsb_pdbx_description.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'')
                        if (polNames.includes(pdbxDescription)){
                            pdb_entries.push({id: entry.parent_rcsb_id.toLowerCase(), name:pdb_text})
                        }
                    });
                    if (pdb_entries.length == 0){return;}
                    vm.pdbs.push(...pdb_entries.sort((a, b) => (a.id > b.id) ? 1 : -1));
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

function getPropensities(property) {
    let indices = null;
    if (vm.structure_mapping && property && property!=0) {
        var sequence_indices = property.indices;
        let alignment_indices = []
        let inverse_structure_mapping = _.invert(vm.structure_mapping);
        for (var sequence_index of sequence_indices) {
            if (inverse_structure_mapping[sequence_index]){
                alignment_indices.push(inverse_structure_mapping[sequence_index])
            }
        }
        indices = alignment_indices.join(',')
    } else {
        indices = '';
    }
    vm.propensity_indices = indices
    vm.fasta_data
}

function handlePropensityIndicesOnTruncatedStructure(indices, startTrunc, endTrunc){
    var newIndices = '';
    var invertedMap = _.invert(vm.structure_mapping);
    if (!indices){
        tempIndices = [];
        vm.all_residues.forEach(function(resi){
            tempIndices.push(invertedMap[resi]);
        })
        indices = tempIndices.join(',');
    }
    indices.split(',').forEach(function(entry){
        if (invertedMap[startTrunc] <= Number(entry) &&  Number(entry) <= invertedMap[endTrunc]){
            newIndices+=`${entry},`;
        }
    })
    return newIndices.slice(0, -1);
}

function handlePropensities(checked_propensities) {
    if (checked_propensities) {
        var title = 'Amino Acid Frequencies'
        if (vm.type_tree == "orth"){
            var title = `${vm.alnobj.text} ${title}`;
        }
        if (vm.property){
            title += ` for ${vm.property.text}`
        }
        let indices = vm.propensity_indices;
        if (vm.selected_domain.length > 0){
            title += `<br>of ECOD domain ${vm.selected_domain[0].name}`;
            var domainIndices = '';
            var invertedMap = _.invert(vm.structure_mapping);
            vm.selected_domain[0].range.split(';').forEach(function(singleRange){
                if (singleRange == ''){return;}
                var startDomain = Number(invertedMap[Number(singleRange.split('-')[0])]);
                var endDomain = Number(invertedMap[Number(singleRange.split('-')[1])]);
                for (let i = startDomain; i <= endDomain; i++){
                    domainIndices += `${i},`;
                }
            });
            if (!indices || indices == ''){
                indices = domainIndices.slice(0, -1);
            } else {
                let tempIndices = _.intersection(indices.split(','),domainIndices.slice(0, -1).split(','));
                indices = tempIndices.join(',')
            }
        }
        if (vm.filter_range){
            var startRange = Number(vm.filter_range.split('-')[0]);
            var endRange = Number(vm.filter_range.split('-')[1].slice(0, -1));
            title += `<br>between positions ${startRange} and ${endRange}`;
            indices = handlePropensityIndicesOnTruncatedStructure(indices, startRange, endRange);
        }
        var full = ['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W'];
        let customFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
        if (indices) {
            ajax("/propensity-data-custom/", {indices, customFasta}).then(data => {
                storeFrequencyData(data['amino acid'],title);
                build_propensity_graph(data['amino acid'], full, title, 'total');
            });
        } else {
            ajax("/propensity-data-custom/", {customFasta}).then(data => {
                storeFrequencyData(data['amino acid'],title);
                build_propensity_graph(data['amino acid'], full, title, 'total');
            });
        }
    }
}

var listSecondaryStructures = function() {
    vm.substructures = []
    let coilObject = {
        text: "Coil residues",
        value: 0,
        indices: vm.coil_residues,
    },
    strandObject = {
        text: "Strand residues",
        value: 1,
        indices: vm.strand_residues,
    },
    helixObject = {
        text: "Helix residues",
        value: 2,
        indices: vm.helix_residues
    };
    Array.prototype.push.apply(vm.substructures, [coilObject, strandObject, helixObject])
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
        xaxis: {title: 'Amino Acid'},
        yaxis: {title: 'Frequency'},
        hovermode: 'closest',
        hoveron: 'points',
    };
    Plotly.newPlot(div, traces, layout);
    var myPlot = document.getElementById(div);
    myPlot.scrollIntoView();
    myPlot.on('plotly_hover', function(data){
        data.points.map(function(d){
            let seqname = d.text.split(' ').slice(1,).join(' ')
            PVAlnViewer.highlightRegion({
                sequences: {from: vm.fastaSeqNames.indexOf(seqname), to: vm.fastaSeqNames.indexOf(seqname)},
                residues: {from: 0, to: vm.fasta_data.split('>')[1].split('\n')[1].length}
            })
        });
        //Plotly.Fx.hover(div,[
        //    { curveNumber:0, pointNumber:1 },
        //    { curveNumber:1, pointNumber:2 },
        //    { curveNumber:2, pointNumber:3 },
        //]);
    }).on('plotly_unhover', function(data){
        //
    });
}

var storeFrequencyData = function (data, title){
    csv_string = '';
    aaNames = [];
    aaList = [];
    once = true;
    for (var key of Object.keys(data)) {
        var entry = data[key];
        if (once){
            once = false;
            aaNames = Object.keys(entry).slice(1);
            csv_string += 'Species\\AA,'+aaNames.join(',')+'\n';
        }
        csv_string += `${entry["name"]},`
        aaNames.forEach(function(aa){
            csv_string += `${entry[aa]},`
        })
        csv_string.slice(0, -1);
        csv_string += '\n';
    }
    vm.freqCSV = `${title.replace('<br>',' ')}\n${csv_string}`;
}