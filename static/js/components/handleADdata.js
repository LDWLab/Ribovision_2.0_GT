export function customADhandler(ad_data) {
    if (vm.uploadSession){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    cleanCustomMap(vm.checked_customMap);
    vm.raiseCustomADWarn = null;
    vm.custom_headers = [];
    if (ad_data == null){
        topviewer.pluginInstance.resetTheme();
        window.viewerInstance.visual.select({data: null, nonSelectedColor: {r:255,g:255,b:255}});
        return;
    }

    var adArray = ad_data.split('\n');
    var ad_header = adArray.shift().replace(/^\s+|\s+$/g, '').split(',');
    var headerLength = ad_header.length;
    if (adArray[adArray.length-1] == 0 || adArray[adArray.length-1] == ''){adArray.splice(-1,1)}
    if (custom_header[0] != 'Index'){
        vm.raiseCustomADWarn = 'Bad AD format:<br/>No defined index header!'
        return;
    }
    if (custom_header.length < 2 || custom_header[custom_header.length-1] == ''){
          vm.raiseCustomADWarn = 'Bad AD format:<br/>Bad data header definition!'
          return;
    }
    if (custom_header.some(r=> Array.from(mapped_aa_properties.keys()).includes(r))){
        vm.raiseCustomADWarn = 'Bad CSV format:<br/>Header definitions exist in RV3!<br/>Use different headers.'
        return;
    }
    if (new Set(custom_header).size !== custom_header.length){
        vm.raiseCustomCSVWarn = 'Bad CSV format:<br/>Duplicate header definitions.'
        return;
    }
    let associatedDataObj = adArray.map(function(e){
        let stringDat = e.replace(/^\s+|\s+$/g, '')
        let currentDat = stringDat.split(',');
        if (stringDat[stringDat.length-1] != ',' && currentDat.length == headerLength){
            return { ix: currentDat[0], data: currentDat.slice(1) }
        } else {
            return 'MISMATCH'
        }
    })
    if (associatedDataObj.includes('MISMATCH')){
        vm.raiseCustomADWarn = 'Bad CSV format:<br/>Mismatch between header and data!'
        return;
    }
    let ass0ciatedDataArrays = [];
    let associatedDataNames = [];
    var colIndex = 0;
    while (colIndex < headerLength-1) {
        let tempArr = [];
        associatedDataObj.forEach((row) => 
            tempArr.push([Number(row.ix), Number(row.data[colIndex])])
        )
        associatedDataNames.push()
        associatedDataArrays.push(tempArr);
        colIndex += 1;
    }

    if (topviewer != null && topviewer.viewInstance.uiTemplateService.domainTypes != undefined){
        for (let ix = 0; ix < associatedDataArrays.length; ix++) {
            vm.associated_headers.push(associated_header[ix+1]);
            mapCustomMappingData(associatedDataArrays[ix], associated_header[ix+1], topviewer);
        }
        displayMappingDataByIndex(topviewer, topviewer.viewInstance.uiTemplateService.domainTypes.length-1);
    }
}
