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

var parseFastaSeqForMSAViewer = function (fasta){
  let outSeqs = [];
  let arrayFasta = String(fasta).split('\n').slice(0, -1);
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

function downloadCSVData() {

  let csv = generateCSVstring(mapped_aa_properties);

  let anchor = document.createElement('a');
  anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
  anchor.target = '_blank';
  anchor.download = 'rv3data.csv';
  anchor.click();

};

function handlePropensities(checked_propensities){
  if (checked_propensities){
      console.log("Checked")
  }else{
      console.log("UnChecked")
  }
  
};


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
};
var csrftoken = getCookie('csrftoken');    function csrfSafeMethod(method) {
  // these HTTP methods do not require CSRF protection
  return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
};
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
};

var pushChainData = function(temp_arr, chain_listI){
  temp_arr.push({
      text: chain_listI["molecule_name"][0],
      value: chain_listI["in_chains"][0],
      sequence: chain_listI["sequence"],
      entityID: chain_listI["entity_id"],
      startIndex: chain_listI.source[0].mappings[0].start.residue_number
  })
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
              if (aln_data["polymers"][ix]["genedescription"].trim() == chain_list[i]["molecule_name"][0]){
                  temp_arr = pushChainData(temp_arr, chain_listI);
              }
          }
      }
  // console.log("___" + temp_arr[temp_arr.length - 1]["sequence"] + "___");
  let chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
  if (chain_options.length === 0) {
      chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
  }
  vueObj.chains = chain_options;
  });
};

var create_deleted_element = function (parent_id, child_id, child_text) {
  const parent = document.getElementById(parent_id);
  const child_elt = document.createElement("div");
  const childText = document.createTextNode(child_text);
  child_elt.setAttribute("id", child_id);
  child_elt.setAttribute("id", child_id);
  child_elt.appendChild(childText);
  parent.appendChild(child_elt);
};

var cleanupOnNewAlignment = function (vueObj, aln_text='') {
  const menu_item = document.querySelector(".smenubar");
  const aln_item = document.getElementById("alnDiv");
  const topview_item = document.getElementById("topview");
  const molstar_item = document.getElementById("pdbeMolstarView");
  const pdb_input = document.getElementById("pdb_input");
  if (menu_item) {menu_item.remove();}
  if (aln_text != ''){
      vueObj.custom_aln_twc_flag == null;
      window.mapped_aa_properties == null;
      if (pdb_input) {
          if (pdb_input.getAttribute("value") != ""){vueObj.pdbid = null;}
      }
      if (vueObj.chains) {vueObj.chains = null;}
      if (vueObj.aln_meta_data) {vueObj.aln_meta_data = null;}
      if (vueObj.fasta_data) {vueObj.fasta_data = null;}
      if (vueObj.fastaSeqNames) {vueObj.fastaSeqNames = null;}
      if (vueObj.frequency_data) {vueObj.frequency_data = null;}
      if (vueObj.topology_loaded) {vueObj.topology_loaded = 'False';}
      if (aln_item) {aln_item.remove(); create_deleted_element("alnif", "alnDiv", aln_text)}
  }
  window.ajaxRun = false;
  if (window.masked_array.length > 0) {window.masked_array = [];}
  if (vueObj.masking_range) {vueObj.masking_range = null;}
  //if (vueObj.chainid) {vueObj.chainid = null;}
  if (vueObj.checked_filter) {vueObj.checked_filter = false;}
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

var calculateFrequencyData = function (frequencies){
  const multiplyvector = function (a,b){
      return a.map((e,i) => e * b[i]);
  }
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
                          ["Mutability",[viridis]],
                          ["Shannon entropy",[plasma]],
                          ["TwinCons",[RdPu, YlGn]],
                      ]);
  window.aaColorData = aaColorData;
  window.aaPropertyConstants = aaPropertiesData;
  window.selectSections_RV1 = new Map();
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
  let csv = 'Index,'
  csv += properties.join(',');
  csv += '\n';
  let csv_ix = [];
  
  mapped_data.get(properties[0]).forEach((datapoint) =>{
      csv_ix.push([datapoint[0]]);
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

var masked_array = [];
