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
    if (!topviewer || !topviewer.viewInstance.uiTemplateService.domainTypes){
        if (checked_customMap){return;}
        var sliceAvailProp = Array.prototype.slice.call(vm.available_properties).filter(availProp => {
            return vm.custom_headers.includes(availProp.Name)
        })
        const setSlice = new Set(sliceAvailProp.map(a=>{return a.Name}));
        const newArray = vm.available_properties.filter(obj => !setSlice.has(obj.Name));
        vm.available_properties = newArray;
        return;
    }
    //var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
    var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    topviewer.viewInstance.uiTemplateService.domainTypes = topviewer.viewInstance.uiTemplateService.domainTypes.filter(obj => {
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
var mapCustomMappingData = function(custom_data, custom_data_name, topviewer){
    
    //var selectBoxEle = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
    var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    
    let vals = custom_data.map(function(v){ return v[1] });
    let indexes = custom_data.map(function(v){ return v[0] });
    window.nColorData.set(custom_data_name, [viridis]);
    window.nPropertyConstants.set(custom_data_name, [Math.min(...vals), Math.max(...vals)]);
    //let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
    //window.coilsOutOfCustom = coilsOutOfCustom;
    var custom_prop = new Map();
    custom_prop.set(custom_data_name, custom_data);
    if (window.custom_prop){
        window.custom_prop.set(custom_data_name, custom_data)
    } else {
        window.custom_prop = custom_prop;
    }
    topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(custom_prop);
    var custom_option = document.createElement("option");
    custom_option.setAttribute("value", selectBoxEle.options.length);
    custom_option.appendChild(document.createTextNode(custom_data_name));
    selectBoxEle.appendChild(custom_option);
    if (!vm.available_properties.some(prop => prop.Name === custom_data_name)){
        vm.available_properties.push({Name:custom_data_name, url:"static/alignments/svg/Custom.svg"})
    }
    if(vm.correct_mask) {
        var j = topviewer.viewInstance.uiTemplateService.domainTypes.length-1;
        colorResidue(j, window.masked_array);
    }
}
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
    //var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //topviewer.pluginInstance.resetTheme();
    //topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data);
    window.viewerInstance.visual.select({
        data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), 
        nonSelectedColor: {r:255,g:255,b:255}
    });
    selectBoxEle.selectedIndex = selectedIndex;
    vm.selected_property = topviewer.pluginInstance.domainTypes[selectedIndex].label;
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
                let riboXYZurl = `https://ribosome.xyz:8000/neo4j/gmo_nom_class/?banName=${oldName}&format=json`
                ajax(riboXYZurl).then(data => {
                    var pdb_entries = []
                    data.forEach(function(entry){
                        let pdb_text = `${entry.parent} ${entry.orgname[0].slice(0,39)}`
                        let pdbxDescription = entry.protein.rcsb_pdbx_description.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'')
                        if (polNames.includes(pdbxDescription)){
                            pdb_entries.push({id: entry.parent.toLowerCase(), name:pdb_text})
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