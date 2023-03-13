var annotationArraySE = [];
var annotationArrayTWC = [];
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
    let tempFasta = String(fastaString).split('\n>');
    tempFasta[0] = tempFasta[0].slice(1);
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
        alert("Empty file was uploaded!");
        return false;
    }
    
    fastaArr = parseFastaString(fasta);
    var fastaSeqs = '';
    var nameSeqs = '';
    var badName = false
    if (fastaArr.length > 4000){
        alert("Fasta file has over 2000 sequences! We currently do not support that many sequences.");
        return false;
    }

    if (fastaArr[1].length*fastaArr.length > 2000000){
        alert("Fasta file has over 1000000 letters! We currently do not support such big files.");
        return false;
    }

    fastaArr.map(function(element, index) {
        if (index % 2 == 1){
            fastaSeqs += fastaArr[index];
        } else {
            if (fastaArr[index].includes('>')){
                badName = '>';
            }
            if (fastaArr[index].includes('Structure sequence')){
                badName = 'struct';
            }
        }
    });

    if (!fastaSeqs) { // is it empty whatever we collected ? re-check not efficient 
        alert("No sequences were found in the file!");
        return false;
    }

    if (badName == '>'){
        alert("The character > should appear only once in sequence headers!");
        return false;
    } else if (badName == 'struct'){
        alert("Structure sequence is a protected sequence id! ProteoVision uses it to append the structure-derived sequence!");
        return false;
    }

    if (!/^[-ACDEFGHIKLMNPQRSTUVWYX\s]+$/i.test(fastaSeqs)){
        alert("Found non-standard characters in the sequences!");
        return false;
    }

    return true;
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

function downloadCSVData() {
  let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
  combined_map = new Map([...mapped_n_properties, ...vm.mapped_n_contacts_mods]);
  let csv = generateCSVstring(combined_map);
  let anchor = document.createElement('a');
  anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
  anchor.target = '_blank';
  anchor.download = `PVData-${month}-${date}-${year}.csv`;
  anchor.click();
};

var downloadAlignmentData = function(fastaString){
    let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
    let anchor = document.createElement('a');
    anchor.href = 'data:text;charset=utf-8,' + encodeURIComponent(fastaString);
    anchor.target = '_blank';
    anchor.download = `PValignment-${month}-${date}-${year}.fas`;
    anchor.click();
}

var downloadAlignmentImage = function(){
    var [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
    var alnDiv = document.querySelector("#MSAViewer");
    var styleHeight = Number(alnDiv.firstElementChild.style.height.replace("px",""))
    alnDiv.firstElementChild.style.height = styleHeight + 50;
    var anchor = document.createElement('a');
    html2canvas(alnDiv).then(canvas => {
        var imageData = canvas.toDataURL("image/png");
        anchor.href = imageData.replace(/^data:image\/png/, "data:application/octet-stream");
        anchor.target = '_blank';
        anchor.download = `PValignment-${month}-${date}-${year}.png`;
        anchor.click();
    })
}

var downloadFullAlignmentImage = function (){
    var [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
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
            anchor.download = `PValignment-full-${month}-${date}-${year}.png`;
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
}

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
            {id: "7k00", name: "7K00 E. coli"},
            {id: "4v6u", name: "4V6U P. furiosus"},
            {id: "4v6x", name: "4V6X H. sapiens"},
        ];
        vueObj.colorScheme = 'nucleotide';
        vueObj.fetchUNtruncatedAln = false;
        vueObj.cdHITReport = false;
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
    window.mapped_n_properties = null;
    vueObj.checkedRNA = false,
    vueObj.customPDBid = null,
    vueObj.pdbStart = null,
    vueObj.pdbEnd = null,
    vueObj.pdbSeq = null,
    vueObj.customPDBsuccess = null,
    vueObj.PDBparsing = false;
    vueObj.entityID = null,
    vueObj.unfilteredChains = null,
    vueObj.hide_chains = null,
    vueObj.all_residues = null;
    vueObj.coil_residues = null;
    vueObj.helix_residues = null;
    vueObj.strand_residues = null;
    vueObj.selected_property = null;
    vueObj.structure_mapping = null;
    vueObj.poor_structure_map = null;
    vueObj.freqCSV = null;
    window.ajaxRun = false;
    window.custom_prop = null;
    if (vueObj.fasta_data) {vueObj.fasta_data = vueObj.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");}
    if (vueObj.topology_loaded) {vueObj.topology_loaded = false;}
    if (vueObj.raiseCustomCSVWarn) {vueObj.raiseCustomCSVWarn = null;}
    if (window.masked_array.length > 0) {window.masked_array = [];}
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

var pushChainData = function(temp_arr, chain_listI){
    try{
      temp_arr.push({
          text: chain_listI["molecule_name"][0],
          value: chain_listI["in_chains"][0],
          sequence: chain_listI["sequence"],
          entityID: chain_listI["entity_id"],
          startIndex: chain_listI.source[0].mappings[0].start.residue_number,
          endIndex: chain_listI.source[0].mappings[0].end.residue_number
      })
      }catch(err){console.log(err);}
    return temp_arr;
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
        if (fkey[2] == "GSD_LSD_rRNA"){
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
    let nPropertiesData = new Map([
        ["Shannon entropy",[0.000000000000001,2.0]],
        ["TwinCons",[-2.25,6.75]]
    ]);
    let nColorData = new Map([
        ["Shannon entropy",[plasma]],
        ["Protein contacts",[rainbow]],
        //["TwinCons",[Reds, Greens]],
        ["TwinCons",[Reds, Blues]],
        //["TwinCons",[RdPu, YlGn]],
    ]);
    window.nColorData = nColorData;
    window.nPropertyConstants = nPropertiesData;
    window.selectSections_RV1 = new Map();
    return nPropertiesData;
}

var calculateFrequencyData = function (frequencies){
  const multiplyvector = function (a,b){
      return a.map((e,i) => e * b[i]);
  }
  nPropertiesData = setGlobalProperties();
  let outPropertyPosition = new Map();
  nPropertiesData.forEach(function (data, property_name){
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

var mapNProps = function (n_properties, mapping){
  let outPropertyMappedPosition = new Map();
  n_properties.forEach(function (data, property_name){
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
  let csv_map = new Map()
  mapped_data.get(properties[0]).forEach((datapoint) =>{
      let alnIx = _.invert(vm.structure_mapping)[datapoint[0]];
      //csv_ix.push([datapoint[0], alnIx]);
      csv_map.set(datapoint[0], [alnIx])
  })
  i = 1
  properties.forEach((prop) => {
      //let ix = 0;
      mapped_data.get(prop).forEach((datapoint) =>{
          if (csv_map.has(datapoint[0])) {
            if (prop != "Protein Contacts" && prop != "Modified Residues") {
                csv_map.get(datapoint[0]).push(datapoint[1])
            } else if (prop == "Protein Contacts" || prop == "Modified Residues") {
                if (csv_map.get(datapoint[0]).length == i + 1){
                    csv_map.get(datapoint[0])[i] = csv_map.get(datapoint[0])[i] + " " + datapoint[1]
                }                
                else {
                    while (csv_map.get(datapoint[0]).length < i) {
                        csv_map.get(datapoint[0]).push(" ")
                    }
                    csv_map.get(datapoint[0]).push(datapoint[1])
                }
            }
        }
          //ix += 1;
    })
      i += 1
  })

  /*csv_ix.forEach((row) => {
      csv += row.join(',');
      csv += '\n';
  })*/
  csv_map.forEach((value, row) => {
    csv += row + ',' + value.join(',');
    csv += '\n';
})
  return csv;
};
var unSelectNucleotide = function(event, pdbId, label_seq_id, isUnobserved) {
    event.stopImmediatePropagation();
    this.clearHighlight(pdbId);
    const ttEle = document.getElementById(`${pdbId}-rnaTopologyTooltip`);
    ttEle.style.display = 'none';

    if(!isUnobserved) {
        const evData = { pdbId, label_seq_id }
        const textElement = document.querySelector(`.rnaview_${pdbId}_${label_seq_id}`);
        CustomEvents.dispatchCustomEvent(this.pdbevents['PDB.RNA.viewer.mouseout'], evData, textElement);
    }
};
var clearHighlight = function(pdbId) {
    var selected = 5;
    document.querySelector(`svg.rnaTopoSvg`).getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${selected}`)[0].setAttribute("fill","323232");
//document.querySelector(`.rnaTopoSvgHighlight_${pdbId}`)!.innerHTML = "";
};
var getEntropyAnnotations = function (separatedData, lowVal, highVal, chainid) {
    annotationArraySE.length=0;
    for (var i = 1; i < 101; i++) {
        annotationArraySE.push({"annotation":i,"ids":[]})
    }
    separatedData.forEach(function (item, index) {
        let parsedItem = item[0];
        let itemValue = item[1];
        let newValue = itemValue - lowVal;
        let normalizedVal = Math.round(newValue/(highVal - lowVal) * 99)
        annotationArraySE[normalizedVal].ids.push(chainid + " " + parsedItem)
    })
    return annotationArraySE;
};

var getTWCAnnotations = function (separatedData, lowVal, highVal, chainid) {
    annotationArrayTWC.length=0;
    for (var i = 1; i < 101; i++) {
        annotationArrayTWC.push({"annotation":i,"ids":[]})
    }
    separatedData.forEach(function (item, index) {
        let parsedItem = item[0];
        let itemValue = item[1];
        
        let newValue = itemValue - lowVal;
        let normalizedVal = Math.round(newValue/(highVal - lowVal) * 99)

        /*if (itemValue < 0){
            console.log('IV_neg', itemValue);
            let normalizedVal = Math.round(itemValue/(lowVal) * 99);
        } 
        else if (itemValue > 0){
            let normalizedVal = Math.round(itemValue/(highVal) * 99);
        }*/
        annotationArrayTWC[normalizedVal].ids.push(chainid + " " + parsedItem)
    })
    return annotationArrayTWC;
}
var getAnnotationArray = function() {
    return {'SE':annotationArraySE,'TWC':annotationArrayTWC };
}   
var parsePVData = function (separatedData, lowVal, highVal, colormapArray, masking=null) {
    let TWCData = new Map();
    let TWCrgbMap = new Map(); 
    let TWCrgbMapPalette = new Map(); 
    let TWCrgbPalette=[];
    let IL=[];
    let ILN=[];
    for (var i = 0; i < 75; i++) {
        TWCrgbMapPalette.set(i, interpolateLinearly((75-i)/75, colormapArray[1]));
        IL=interpolateLinearly((75-i)/75, colormapArray[1]);
        TWCrgbPalette.push(IL[0]);
    };
    for (var i = 0; i < 25; i++) {
        TWCrgbMapPalette.set(i, interpolateLinearly(i/25, colormapArray[0]));
        ILN=interpolateLinearly(i/25, colormapArray[0]);
        TWCrgbPalette.push(ILN[0]);
    };
    separatedData.forEach(function (item, index) {
        let parsedItem = item[0];
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
var hexToRgb = function (hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
      r: parseInt(result[1], 16),
      g: parseInt(result[2], 16),
      b: parseInt(result[3], 16)
    } : null;
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

var mapTWCdata = function (structMap, twcDataUnmapped, mapped_n_properties){
    var topviewer = document.getElementById("PdbeTopViewer");
    mapped_n_properties = build_mapped_props(mapped_n_properties, twcDataUnmapped, structMap);
    window.mapped_n_properties = mapped_n_properties;
    if (topviewer != null && topviewer.viewInstance.uiTemplateService.domainTypes != undefined){
        var empty_props = new Map();
        let twc_props = build_mapped_props(empty_props, twcDataUnmapped, structMap);
        topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(twc_props);
        //var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
        //var twc_option = document.createElement("option");
        //twc_option.setAttribute("value", selectBoxEle.options.length);
        //twc_option.appendChild(document.createTextNode("TwinCons"));
        //selectBoxEle.appendChild(twc_option);
    }
}
var showPDBHelper = function(pdbid, chainid, entityid) {
    const molstar_item = document.getElementById("pdbeMolstarView");
    if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Loading Molstar Component ", true)}
    var pdblower = pdbid.toLocaleLowerCase();
    if (pdbid == "cust"){
        var coordURL = `/custom-struc-data/${pdblower}-${entityid}-${chainid}`;
        var binaryCif = false;
        var structFormat = "cif";
    }else{
        //var coordURL = `https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=${entityid}&encoding=bcif`
        //var coordURL = `https://coords.litemol.org/${pdblower}/chains?entityId=${entityid}&authAsymId=${chainid}&encoding=bcif`;
        var coordURL = `https://models.rcsb.org/v1/${pdblower}/atoms?label_entity_id=${entityid}&encoding=bcif`
        var binaryCif = true;
        var structFormat = "bcif";
    }
    window.pdblower = pdblower;
    var viewerInstance = new PDBeMolstarPlugin();
    vm.viewer_options = {
        customData: { url: coordURL,
                        format: structFormat, 
                        binary: binaryCif },
        hideCanvasControls: ["expand", "selection", " animation"],
        assemblyId: '1',
        hideControls: true,
        subscribeEvents: true,
        bgColor: {r:255,g:255,b:255},
    }
    var viewerContainer = document.getElementById('pdbeMolstarView');
    viewerInstance.render(viewerContainer, vm.viewer_options);
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
                    residue_number:residueNumber,
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
                    residue_number:residueNumber,
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
var drawCircle = function (pdbId, i, color){
    const circle = document.querySelector(`svg.rnaTopoSvg`).getElementsByClassName(`circle_${pdbId}_${i}`)[0]
    const nucleotide = document.querySelector(`svg.rnaTopoSvg`).getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${i}`)[0]
    const BBox = nucleotide.getBBox()
    const nx = (BBox.x + BBox.width/2)
    const ny = (BBox.y + BBox.height/2)
    circle.setAttribute("cx", nx)
    circle.setAttribute("cy", ny)
    circle.setAttribute("stroke", `${color}`);
    circle.setAttribute("fill", `${color}`);
    circle.style.display = "block";
}
var calculateModifiedResidues = function(pdbid, chainid, entityid) {
    var url = `modified-residues/${pdbid}/${chainid}`
    ajax(url).then(data => {
        let offset = 0
        let modifiedData = new Map()
        let modifications = []
        for (let val in data.Modified) {
            if(modifications.indexOf(data.Modified[val][0]) < 0) {
                modifications.push(data.Modified[val][0])
                modifiedData.set(data.Modified[val][0], [])
            }
            index = data.Modified[val][1] - offset
            modifiedData.get(data.Modified[val][0]).push(index)
            offset += 4
        }
        vm.modified_residues = modifiedData
        var i = 1.0;
        var colorMap = new Map();
        vm.selectSections_modified = new Map();
        vm.mapped_n_contacts_mods.set("Modified Residues", [])
        for (var val of modifications) {
            vm.selectSections_modified.set(val, [])
            //Need to add modifications color scheme, using PC for now
            var color = interpolateLinearly(i/modifications.length, nColorData.get("Protein contacts")[0])
            var rgbColor = "rgb(" + color[0][0] + "," + color[0][1] + "," + color[0][2] + ")";
            colorMap.set(val, rgbColor);
            //newContactMap.set(vm.protein_contacts, nColorData.get("Shannon entropy")[0][1]
            i = i+1;
            for (var j of vm.modified_residues.get(val)) {
                vm.selectSections_modified.get(val).push({
                    entity_id: "" + entityid,
                    residue_number: j, 
                    color: color[1],
                    sideChain: false,
                });
                vm.mapped_n_contacts_mods.get("Modified Residues").push([j, val])
            }
        }                 
        vm.modifiedColorMap = colorMap;
        if(data.Modified.length > 0) {
            vm.modified = true
        }
        //viewerInstanceTop.viewInstance.uiTemplateService.colorMap(); 
    });
}
var showContactsHelper = function(entityid) {
    var protein_data = new Map();
    protein_data.set("contacts", [])
    protein_data.get("contacts").push({entity_id: entityid, focus: true})
    for (let val in vm.pchainid) {
        var chain = vm.pchainid[val];
        for (let entry in vm.selectSections_proteins.get(chain)) {
            protein_data.get("contacts").push(vm.selectSections_proteins.get(chain)[entry])
        }
    }
    const mapSort1 = protein_data.get("contacts").sort((a, b) => a.residue_number - b.residue_number);
    viewerInstance.visual.select({
        data: mapSort1, 
        nonSelectedColor: {r:255,g:255,b:255}
        }).catch(err => {
            console.log(err);
            vm.$nextTick(function(){
                viewerInstance.visual.select({
                    data: mapSort1,
                    nonSelectedColor: {r:255,g:255,b:255}
                })
            })
        })
}
const sleep = (delay) => new Promise (( resolve) => setTimeout (resolve, delay))
var showProteins3D = function() {
    showPDBHelper(vm.pdbid, vm.chainid, vm.entityID)
    
    const showProteins = async () => {
        await sleep (5000)
        colorData = []
        for (let val in vm.pchainid) {
            auth_id = vm.pchainid[val]
            chain = vm.protein_chains.filter(e => e.value == auth_id)[0]
            eID = chain.entityID
            data = {url: `https://www.ebi.ac.uk/pdbe/model-server/v1/${vm.pdbid}/atoms?label_entity_id=${eID}&auth_asym_id=${auth_id}&encoding=bcif`, format: 'cif', binary:true, bgColor: {r: 255, g: 255, b: 255}}
            viewerInstance.visual.update({customData: data, bgColor: {r: 255, g: 255, b: 255}}, false)
            color = vm.proteinColorMap.get(auth_id)
            split_color = color.split('(')[1].split(')')[0].split(',')
            formatted_color = {r: split_color[0], g: split_color[1], b: split_color[2]}
            colorData.push(formatted_color)
        }
        await sleep (10000)
        viewerInstance.visual.colorByChain(colorData)
    }
    showProteins()
}
var showModificationsAndContactsHelper = function(entityid) {
    if (vm.pchainid.length > 0){
        showProteins3D()
    }
    var modified_data = new Map();
    modified_data.set("mods", [])
    modified_data.get("mods").push({entity_id: entityid, focus: true})
    for (let val in vm.pchainid) {
        var chain = vm.pchainid[val];
        for (let entry in vm.selectSections_proteins.get(chain)) {
            modified_data.get("mods").push(vm.selectSections_proteins.get(chain)[entry])
        }
    }
    for (let val of vm.modifications) {
        for (let entry of vm.selectSections_modified.get(val)) {
            modified_data.get("mods").push(entry)
        }
    }
    /*window.viewerInstance.visual.select({
        data: [],
        nonSelectedColor: {r:255,g:255,b:255}
    })*/
    const mapSort1 = modified_data.get("mods").sort((a, b) => a.residue_number - b.residue_number);
    const selectColors = async() => {
        if(vm.pchainid.length > 0) {
            await sleep(5000)
        }
        viewerInstance.visual.select({
            data: mapSort1, 
            nonSelectedColor: {r:255,g:255,b:255}
            }).catch(err => {
                console.log(err);
                vm.$nextTick(function(){
                    viewerInstance.visual.select({
                        data: mapSort1,
                        nonSelectedColor: {r:255,g:255,b:255}
                    })
                })
            })
    }
    selectColors()
        
}
var recolorTopStar = function (name){
    var selectBox = viewerInstanceTop.viewInstance.targetEle.querySelector('.mappingSelectbox');
    var newIndex = indexMatchingText(selectBox.options, name);
    //var selectedDomain = viewerInstanceTop.viewInstance.uiTemplateService.domainTypes[newIndex];
    selectBox.selectedIndex = newIndex; 
    //if(selectSections_RV1.get(name).length < 1600) {
    /*
        viewerInstance.visual.select({
            data: selectSections_RV1.get(name), 
            nonSelectedColor: {r:255,g:255,b:255}
        }).catch(err => {
            console.log(err);
            vm.$nextTick(function(){
                viewerInstance.visual.select({
                    data: selectSections_RV1.get(name), 
                    nonSelectedColor: {r:255,g:255,b:255}
                })
            })
        }) 
        */
    //}
    if(name == "Shannon entropy") {
        viewerInstance.visual.clearSelection();
        viewerInstance.coloring.shannonEntropy({ sequence: true, het: false, keepStyle: true });
    }   else if(name == "TwinCons") {
        viewerInstance.visual.clearSelection();
        viewerInstance.coloring.twinCons({ sequence: true, het: false, keepStyle: true });
    } else if(name == "Select data") {
        viewerInstance.visual.reset({ theme: true })
    }
    
    viewerInstanceTop.viewInstance.uiTemplateService.colorMap(); 
}

var masked_array = [];
