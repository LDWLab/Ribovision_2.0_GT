var registerHoverResiData = function (e, tooltipObj){
  const strainQuery = '&res__poldata__strain__strain=';
  var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(e.position) + 1)}&aln=${vm.alnobj.id}${strainQuery}${vm.fastaSeqNames[Number(e.i)]}`
  ajax(url).then(alnpos_data => {
      var alnViewCanvasEle = document.querySelector("#alnDiv canvas:nth-of-type(1)");
      var alnViewLabelsEle = document.querySelector("#alnDiv div:nth-of-type(2)");
      let boundLabelBox = alnViewLabelsEle.getBoundingClientRect();
      let boundingBox = absolutePosition(alnViewCanvasEle);
      let relativeBox = alnViewCanvasEle.getBoundingClientRect();
      if (alnpos_data.count != 0){
      ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
          if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
              let tooltipPosition = {
                top: mousePos.y-boundingBox.top+5 +"px",
                left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+5 +"px",
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
                  top: mousePos.y-boundingBox.top+5 +"px",
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
          subscribeEvents: true
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
      subscribeEvents: true});
  topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
  if(window.custom_prop) {
      topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
  }
  topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
  window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
};
