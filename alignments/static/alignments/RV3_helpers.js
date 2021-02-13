var absolutePosition = function (el) {
  var
      found,
      left = 0,
      top = 0,
      width = 0,
      height = 0,
      offsetBase = absolutePosition.offsetBase;
  if (!offsetBase && document.body) {
      offsetBase = absolutePosition.offsetBase = document.createElement('div');
      offsetBase.style.cssText = 'position:absolute;left:0;top:0';
      document.body.appendChild(offsetBase);
  }
  if (el && el.ownerDocument === document && 'getBoundingClientRect' in el && offsetBase) {
      var boundingRect = el.getBoundingClientRect();
      var baseRect = offsetBase.getBoundingClientRect();
      found = true;
      left = boundingRect.left - baseRect.left;
      top = boundingRect.top - baseRect.top;
      width = boundingRect.right - boundingRect.left;
      height = boundingRect.bottom - boundingRect.top;
  }
  return {
      found: found,
      left: left,
      top: top,
      width: width,
      height: height,
      right: left + width,
      bottom: top + height
  };
};

var parseFastaString = function(fastaString){
    let arrayFasta = [];
    let tempFasta = String(fastaString).split('>');
    tempFasta = tempFasta.filter(n => n);
    tempFasta.forEach(seq =>{
        let splitSeq = seq.split(/\n/);
        arrayFasta.push(splitSeq[0]);
        arrayFasta.push(splitSeq.slice(1).join(''))
    });
    return arrayFasta;
}

var validateFasta = function (fasta) {
    //From here https://www.blopig.com/blog/2013/03/a-javascript-function-to-validate-fasta-sequences/
    
    if (!fasta) { // check there is something first of all
        return false;
    }
    
    fastaArr = parseFastaString(fasta);
    var fastaSeqs = '';
    fastaArr.map(function(element, index) {
        if (index % 2 == 1){
            fastaSeqs += fastaArr[index];
        }
    });

    if (!fastaSeqs) { // is it empty whatever we collected ? re-check not efficient 
        return false;
    }

    return /^[-ACDEFGHIKLMNPQRSTUVWY\s]+$/i.test(fastaSeqs);
}

var parseFastaSeqForMSAViewer = function (fasta){
    let outSeqs = [];
    arrayFasta = parseFastaString(fasta);
    arrayFasta.map(function(element, index) {
        if (index % 2 == 0){
            let seqName = element.replaceAll('_', ' ').replaceAll('>', '');
            let seqObj = {'name': seqName, 'sequence': arrayFasta[index+1]}
            outSeqs.push(seqObj);
        }
    });
    return outSeqs;
};

(function() {
  var mousePos;
  document.onmousemove = handleMouseMove;
  function handleMouseMove(event) {
      var eventDoc, doc, body;
      event = event || window.event; // IE-ism
      // If pageX/Y aren't available and clientX/Y are,
      // calculate pageX/Y - logic taken from jQuery.
      // (This is to support old IE)
      if (event.pageX == null && event.clientX != null) {
          eventDoc = (event.target && event.target.ownerDocument) || document;
          doc = eventDoc.documentElement;
          body = eventDoc.body;
          event.pageX = event.clientX +
            (doc && doc.scrollLeft || body && body.scrollLeft || 0) -
            (doc && doc.clientLeft || body && body.clientLeft || 0);
          event.pageY = event.clientY +
            (doc && doc.scrollTop  || body && body.scrollTop  || 0) -
            (doc && doc.clientTop  || body && body.clientTop  || 0 );
      }
      mousePos = {
        x: event.pageX,
        y: event.pageY
    };
    window.mousePos = mousePos;
  }
})();

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

function downloadCSVData() {
  let csv = generateCSVstring(mapped_aa_properties);
  let anchor = document.createElement('a');
  anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
  anchor.target = '_blank';
  anchor.download = 'PVdata.csv';
  anchor.click();
};

var downloadAlignmentData = function(fastaString){
    let anchor = document.createElement('a');
    anchor.href = 'data:text;charset=utf-8,' + encodeURIComponent(fastaString);
    anchor.target = '_blank';
    anchor.download = 'PValignment.fas';
    anchor.click();
}

var downloadAlignmentImage = function(){
    var alnDiv = document.querySelector("#MSAViewer");
    var anchor = document.createElement('a');
    html2canvas(alnDiv).then(canvas => {
        var imageData = canvas.toDataURL("image/png");
        anchor.href = imageData.replace(/^data:image\/png/, "data:application/octet-stream");
        anchor.target = '_blank';
        anchor.download = 'PValignment.png';
        anchor.click();
    })
}

var downloadFullAlignmentImage = function (){
    var handleCanvasErr = function (err, labelsDiv, initialLabelsWidth){
        labelsDiv.style.width = initialLabelsWidth;
        PVAlnViewer.handleResize();
        alert("Couldn't parse the alignment. Check the console for error.")
        console.log(err)
    }
    var alnLength = vm.fasta_data.split('>')[1].split('\n')[1].length;
    PVAlnViewer.setState({
        aaPos: 0,
        seqPos: 0,
        height: (vm.fastaSeqNames.length+2)*17, 
        width: (alnLength+5)*17,
    });
    var labelsDiv = document.querySelector("#alnViewerLabels");
    var longestName = '';
    labelsDiv.firstElementChild.firstElementChild.children.forEach(function (labelNode) {
        if (labelNode.textContent.length > longestName.length) {
            longestName = labelNode.textContent;
          }
    });
    var initialLabelsWidth = labelsDiv.style.width;
    var maxLabelsWidth = getWidthOfText(longestName, 'Arial', '15px');
    labelsDiv.style.width = maxLabelsWidth;
    var anchor = document.createElement('a');
    var alnDiv = document.querySelector("#MSAViewer");
    //Have to do 2 nested html2canvas so that canvas gets rerendered at the 0,0 position.
    html2canvas(alnDiv).then(() => {
        var alnDivDownload = document.querySelector("#MSAViewer");
        html2canvas(alnDivDownload).then(canvas => {
            var imageData = canvas.toDataURL("image/png");
            labelsDiv.style.width = initialLabelsWidth;
            PVAlnViewer.handleResize();
            anchor.href = imageData.replace(/^data:image\/png/, "data:application/octet-stream");
            anchor.target = '_blank';
            anchor.download = 'PValignment.png';
            anchor.click();
        }).catch(err => {
            handleCanvasErr(err, labelsDiv, initialLabelsWidth);
        })
    }).catch(err => {
        handleCanvasErr(err, labelsDiv, initialLabelsWidth);
    });
}

var getWidthOfText = function (txt, fontname, fontsize){
    if(getWidthOfText.c === undefined){
        getWidthOfText.c=document.createElement('canvas');
        getWidthOfText.ctx=getWidthOfText.c.getContext('2d');
    }
    var fontspec = fontsize + ' ' + fontname;
    if(getWidthOfText.ctx.font !== fontspec)
        getWidthOfText.ctx.font = fontspec;
    return getWidthOfText.ctx.measureText(txt).width;
}

var pushChainData = function(temp_arr, chain_listI){
  try{
    temp_arr.push({
        text: chain_listI["molecule_name"][0],
        value: chain_listI["in_chains"][0],
        sequence: chain_listI["sequence"],
        entityID: chain_listI["entity_id"],
        startIndex: chain_listI.source[0].mappings[0].start.residue_number
    })
    }catch(err){console.log(err);}
  return temp_arr;
};

var filterAvailablePolymers = function(chain_list, aln_id, vueObj) {
  let temp_arr = [];
  let url = `/desire-api/alignments/${aln_id}/?format=json`;
  ajax(url).then( aln_data => {
      for (let i = 0; i < chain_list.length; i++) {
          let chain_listI = chain_list[i]
          if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
          if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
          for (let ix =0; ix < aln_data["polymers"].length; ix++){
              let desirePolymerName = aln_data["polymers"][ix]["genedescription"].trim();
              let pdbePolymerNames = chain_list[i]["molecule_name"];
              for (let nameIx =0; nameIx < pdbePolymerNames.length; nameIx++){
                  let pdbeName = pdbePolymerNames[nameIx].replace(/-[A-Z]{1}$/,'');
                  if (pdbeName == desirePolymerName){
                    temp_arr = pushChainData(temp_arr, chain_listI);
                    break;
                  }
              }
          }
      }
  let chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
  if (chain_options.length === 0) {
      var elt = document.querySelector("#onFailedChains");
      elt.innerHTML  = "Couldn't find a matching chain!<br>Try a different PDB ID."
      vueObj.pdbid = null;
      chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
  }else{
    vueObj.hide_chains = null;
  }
  vueObj.chains = chain_options;
  });
};

var create_deleted_element = function (parent_id, child_id, child_text, optionalLoadIMG=null) {
    const parent = document.getElementById(parent_id);
    const child_elt = document.createElement("div");
    const childText = document.createTextNode(child_text);
    child_elt.setAttribute("id", child_id);
    child_elt.setAttribute("id", child_id);
    child_elt.appendChild(childText);
    if (optionalLoadIMG){
        let imgElt = document.createElement("img");
        imgElt.setAttribute("src","static/img/loading.gif");
        imgElt.setAttribute("alt","Loading");
        imgElt.setAttribute("style","height:25px;");
        child_elt.appendChild(imgElt);
    }
    parent.appendChild(child_elt);
};

var cleanupOnNewAlignment = function (vueObj, aln_text='') {
    if (vm.uploadSession){return;}
    const menu_item = document.querySelector(".smenubar");
    const aln_item = document.getElementById("alnDiv");
    const topview_item = document.getElementById("topview");
    const molstar_item = document.getElementById("pdbeMolstarView");
    const pdb_input = document.getElementById("pdb_input");
    if (menu_item) {menu_item.remove();}
    if (aln_text != ''){
        vueObj.custom_aln_twc_flag = null;
        vueObj.pdbs = [
            {id: "4v9d", name: "4V9D E. coli"},
            {id: "4v6u", name: "4V6U P. furiosus"},
            {id: "4ug0", name: "4UG0 H. sapiens"},
        ];
        vueObj.colorScheme = 'clustal2';
        vueObj.aaPos = 0;
        vueObj.seqPos = 0;
        vueObj.msavWillMount = null;
        vueObj.unmappedTWCdata = null;
        if (pdb_input) {
            if (pdb_input.getAttribute("value") != ""){vueObj.pdbid = null;}
        }
        if (vueObj.chains) {vueObj.chains = null;}
        if (vueObj.aln_meta_data) {vueObj.aln_meta_data = null;}
        if (vueObj.fasta_data) {vueObj.fasta_data = null;}
        if (vueObj.fastaSeqNames) {vueObj.fastaSeqNames = null;}
        if (vueObj.frequency_data) {vueObj.frequency_data = null;}
        if (aln_item) {aln_item.remove(); create_deleted_element("alnif", "alnDiv", aln_text, true)}
    }
    window.mapped_aa_properties = null;
    vueObj.checked_propensities = null;
    vueObj.structure_mapping = null;
    vueObj.poor_structure_map = null;
    vueObj.selected_property = null;
    window.ajaxRun = false;
    window.custom_prop = null;
    if (vueObj.topology_loaded) {vueObj.topology_loaded = false;}
    if (vueObj.raiseCustomCSVWarn) {vueObj.raiseCustomCSVWarn = null;}
    if (window.masked_array.length > 0) {window.masked_array = [];}
    if (vueObj.masking_range) {vueObj.masking_range = null;}
    if (vueObj.checked_filter) {vueObj.checked_filter = false;}
    if (vueObj.checked_selection) {vueObj.checked_selection = false;}
    if (vueObj.checked_customMap) {vueObj.checked_customMap = false;}
    if (vueObj.csv_data) {vueObj.csv_data = null;}
    if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Select new chain!")}
    if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Select new structure!")}
};

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
};

var intersection = function () {
    var result = [];
    var lists;
    if(arguments.length === 1) {
        lists = arguments[0];
    } else {
        lists = arguments;
    }
    for(var i = 0; i < lists.length; i++) {
        var currentList = lists[i];
        for(var y = 0; y < currentList.length; y++) {
            var currentValue = currentList[y];
            if(result.indexOf(currentValue) === -1) {
                if(lists.filter(function(obj) { return obj.indexOf(currentValue) == -1 }).length == 0) {
                    result.push(currentValue);
            }
          }
        }
    }
    return result;
}

var objectify = function (array){
    return array.reduce(function(p, c) {
         p[c[0]] = [c[1], c[2]];
         return p;
    }, {});
}

var loadOrthAlns = function(data, vm){
    if (data["results"].length === 1) {
        var fpa = data["results"][0]["alignment_ids"]
    } else if (data["results"].length > 1) {
        var fpa = [];
        var alnid_maps = [];
        var alnid_keys = [];
        data["results"].forEach(function(tax_result){
            var temp_map = objectify(tax_result["alignment_ids"])
            alnid_maps.push(temp_map);
            alnid_keys.push(Object.keys(temp_map));
        });
        var filtered_keys = intersection(alnid_keys);
        filtered_keys.forEach(function(alnk){
            fpa.push(Array(Number(alnk), alnid_maps[0][alnk][0], alnid_maps[0][alnk][1]))
        })
    } else {
        var fpa = [null, 'No alignments found', "PROMALS3D"]
    }
    var fpa_viz = [];
    fpa.forEach(function(fkey) {
        if (fkey[2] == "PROMALS3D"){
            fpa_viz.push({
                text: fkey[1],
                value: fkey[0]
            });
        }
    });
    vm.alignments = fpa_viz
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
};

var setGlobalProperties = function(){
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
        ["Mutability",[plasma]],
        ["Shannon entropy",[plasma]],
        ["TwinCons",[Reds, Greens]],
        //["TwinCons",[Reds, Blues]],
        //["TwinCons",[RdPu, YlGn]],
    ]);
    window.aaColorData = aaColorData;
    window.aaPropertyConstants = aaPropertiesData;
    window.selectSections_RV1 = new Map();
    return aaPropertiesData;
}

var calculateFrequencyData = function (frequencies){
  const multiplyvector = function (a,b){
      return a.map((e,i) => e * b[i]);
  }
  aaPropertiesData = setGlobalProperties();
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
};

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
};

var filterCoilResidues = function (coil_data){
  const range = (start, stop, step) => Array.from({ length: (stop - start) / step + 1}, (_, i) => start + (i * step));
  let coilResidues = [];
  coil_data.forEach(function (coilRange){
      if (coilRange.start < coilRange.stop){
          coilResidues.push(range(coilRange.start, coilRange.stop, 1))
      }
  })
  return coilResidues.flat()
};

var generateCSVstring = function (mapped_data){
  let properties = Array.from(mapped_data.keys());
  let csv = 'Index,Alignment index,'
  csv += properties.join(',');
  csv += '\n';
  let csv_ix = [];
  
  mapped_data.get(properties[0]).forEach((datapoint) =>{
      let alnIx = _.invert(vm.structure_mapping)[datapoint[0]];
      csv_ix.push([datapoint[0], alnIx]);
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
};

function replacer(key, value) {
    const originalObject = this[key];
    if(originalObject instanceof Map) {
      return {
        dataType: 'Map',
        value: Array.from(originalObject.entries()), // or with spread: value: [...originalObject]
      };
    } else {
      return value;
    }
  }

function reviver(key, value) {
    if(typeof value === 'object' && value !== null) {
      if (value.dataType === 'Map') {
        return new Map(value.value);
      }
    }
    return value;
  }

var parsePVData = function (separatedData, lowVal, highVal, colormapArray, masking=null) {
        let TWCData = new Map();
        let TWCrgbMap = new Map();    
        separatedData.forEach(function (item, index) {
            let parsedItem = item[0];
            //if(!masking || masking[index]) {
                let itemValue = item[1];
                TWCData.set(parsedItem, itemValue);
                if (colormapArray.length === 1) {
                    let newValue = itemValue - lowVal;
                    TWCrgbMap.set(parsedItem, interpolateLinearly(newValue/(highVal - lowVal), colormapArray[0]));
                }
                else {
                    if (itemValue === 'NA'){
                        TWCrgbMap.set(parsedItem, [[192, 192, 192], {r:192, g:192, b:192}]);
                    } else if (itemValue < 0){
                        TWCrgbMap.set(parsedItem, interpolateLinearly(itemValue/lowVal, colormapArray[0]));
                    } else {
                        TWCrgbMap.set(parsedItem, interpolateLinearly(itemValue/highVal, colormapArray[1]));
                    }
                }
            //}
            /*else {
                TWCrgbMap.set(parsedItem, [[255, 255, 255], {r:0, g:0, b:0, a:.4}]);
                TWCData.set(parsedItem, null);
            }*/
        });
        return [TWCrgbMap, TWCData];
    }

var indexMatchingText = function(ele, text) {
    for (var i=0; i<ele.length;i++) {
        if (ele[i].childNodes[0].nodeValue === text){
            return i;
        }
    }
    return undefined;
}

function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
  }
var rgbToHex = function(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

var build_mapped_props = function(mapped_props, twcDataUnmapped, structure_mapping){
    mapped_props.set("TwinCons", [])
    for (let i = 0; i < twcDataUnmapped.length; i++) {
        let mappedI0 = structure_mapping[twcDataUnmapped[i][0]];
        if (mappedI0) {
            mapped_props.get("TwinCons").push([mappedI0, twcDataUnmapped[i][1]]);
        }
    }
    return mapped_props;
}

var mapTWCdata = function (structMap, twcDataUnmapped, mapped_aa_properties){
    var topviewer = document.getElementById("PdbeTopViewer");
    mapped_aa_properties = build_mapped_props(mapped_aa_properties, twcDataUnmapped, structMap);
    window.mapped_aa_properties = mapped_aa_properties;
    if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
        var empty_props = new Map();
        let twc_props = build_mapped_props(empty_props, twcDataUnmapped, structMap);
        topviewer.pluginInstance.getAnnotationFromRibovision(twc_props);
        var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
        var twc_option = document.createElement("option");
        twc_option.setAttribute("value", selectBoxEle.options.length);
        twc_option.appendChild(document.createTextNode("TwinCons"));
        selectBoxEle.appendChild(twc_option);
    }
}

var fetchTWCdata = function (fasta){
    ajax('/twc-api/', {fasta}).then(twcDataUnmapped => {
        vm.unmappedTWCdata = twcDataUnmapped;
        var settedProps = new Set(vm.available_properties.map(a=>{return a.Name}))
        if (!settedProps.has("TwinCons")){
            vm.available_properties.unshift({Name:"TwinCons", url:"static/alignments/svg/TwinCons.svg"});
        }
    })
}

var recolorTopStar = function (name){
    var selectBox = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    var newIndex = indexMatchingText(selectBox.options, name)
    selectBox.selectedIndex = newIndex; 
    viewerInstanceTop.pluginInstance.displayDomain();
}

var masked_array = [];
