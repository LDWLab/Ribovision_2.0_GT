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
  let csv = generateCSVstring(mapped_aa_properties);
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
    window.mapped_aa_properties = null;
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
    vueObj.checked_propensities = null;
    vueObj.domain_or_selection = null;
    vueObj.checked_domain = null;
    vueObj.selected_domain = [];
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
    let aaPropertiesData = new Map([
        ["Shannon entropy",[0.000000000000001,2.0]],
        ["TwinCons",[-2.935,12.065]]
    ]);
    let aaColorData = new Map([
        ["Shannon entropy",[plasma]],
        ["Protein contacts",[rainbow]],
        //["TwinCons",[Reds, Greens]],
        ["TwinCons",[Reds, Blues]],
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
}
var clearHighlight = function(pdbId) {
        var selected = 5;
        document.querySelector(`svg.rnaTopoSvg`).getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${selected}`)[0].setAttribute("fill","323232");
    //document.querySelector(`.rnaTopoSvgHighlight_${pdbId}`)!.innerHTML = "";
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

var mapTWCdata = function (structMap, twcDataUnmapped, mapped_aa_properties){
    var topviewer = document.getElementById("PdbeTopViewer");
    mapped_aa_properties = build_mapped_props(mapped_aa_properties, twcDataUnmapped, structMap);
    window.mapped_aa_properties = mapped_aa_properties;
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
    /*window.viewerInstance.visual.select({
        data: [],
        nonSelectedColor: {r:255,g:255,b:255}
    })*/
    const mapSort1 = protein_data.get("contacts").sort((a, b) => a.start_residue_number - b.start_residue_number);
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
var recolorTopStar = function (name){
    var selectBox = viewerInstanceTop.viewInstance.targetEle.querySelector('.mappingSelectbox');
    var newIndex = indexMatchingText(selectBox.options, name);
    //var selectedDomain = viewerInstanceTop.viewInstance.uiTemplateService.domainTypes[newIndex];
    selectBox.selectedIndex = newIndex; 
    viewerInstanceTop.viewInstance.uiTemplateService.colorMap(); 
    if(selectSections_RV1.get(name).length < 1600) {
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
    }
}

var masked_array = [];
