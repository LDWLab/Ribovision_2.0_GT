export function customCSVhandler(csv_data) {
    if (vm.uploadSession){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    cleanCustomMap(vm.checked_customMap);
    vm.raiseCustomCSVWarn = null;
    vm.custom_headers = [];
    if (csv_data == null){
        topviewer.pluginInstance.resetTheme();
        window.viewerInstance.visual.select({data: null, nonSelectedColor: {r:255,g:255,b:255}});
        return;
    }
    var csvArray = csv_data.split('\n');
    var custom_header = csvArray.shift().replace(/^\s+|\s+$/g, '').split(',');
    var headerLength = custom_header.length;
    if (csvArray[csvArray.length-1] == 0 || csvArray[csvArray.length-1] == ''){csvArray.splice(-1,1)}
    if (custom_header[0] != 'Index'){
        vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>No defined index header!'
        return;
    }
    if (custom_header.length < 2 || custom_header[custom_header.length-1] == ''){
          vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>Bad data header definition!'
          return;
    }
    if (custom_header.some(r=> Array.from(mapped_aa_properties.keys()).includes(r))){
        vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>Header definitions exist in RV3!<br/>Use different headers.'
        return;
    }
    if (new Set(custom_header).size !== custom_header.length){
        vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>Duplicate header definitions.'
        return;
    }
    let customDataObj = csvArray.map(function(e){
        let stringDat = e.replace(/^\s+|\s+$/g, '')
        let currentDat = stringDat.split(',');
        if (stringDat[stringDat.length-1] != ',' && currentDat.length == headerLength){
            return { ix: currentDat[0], data: currentDat.slice(1) }
        } else {
            return 'MISMATCH'
        }
    })
    if (customDataObj.includes('MISMATCH')){
        vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>Mismatch between header and data!'
        return;
    }
    let customDataArrays = [];
    let customDataNames = [];
    var colIndex = 0;
    while (colIndex < headerLength-1) {
        let tempArr = [];
        customDataObj.forEach((row) => 
            tempArr.push([Number(row.ix), Number(row.data[colIndex])])
        )
        customDataNames.push()
        customDataArrays.push(tempArr);
        colIndex += 1;
    }

    if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
        for (let ix = 0; ix < customDataArrays.length; ix++) {
            vm.custom_headers.push(custom_header[ix+1]);
            mapCustomMappingData(customDataArrays[ix], custom_header[ix+1], topviewer);
        }
        displayMappingDataByIndex(topviewer, topviewer.pluginInstance.domainTypes.length-1);
    }
}
