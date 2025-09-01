export function customADhandler(ad_data) {
    if (vm.uploadSession){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    cleanCustomMap(vm.checked_customMap);
    vm.raiseCustomADWarn = null;
    vm.associated_headers = [];
    if (ad_data == null){
        topviewer.pluginInstance.resetTheme();
        window.viewerInstance.visual.select({data: null, nonSelectedColor: {r:255,g:255,b:255}});
        return;
    }

    var adArray = ad_data.split('\n');
    var ad_header = adArray.shift().replace(/^\s+|\s+$/g, '').split(',');
    var headerLength = ad_header.length;
    if (adArray[adArray.length-1] == 0 || adArray[adArray.length-1] == ''){adArray.splice(-1,1)}
    if (ad_header[0] != 'Index'){
        vm.raiseCustomADWarn = 'Bad AD format:<br/>No defined index header!'
        return;
    }
    if (ad_header.length < 2 || ad_header[ad_header.length-1] == ''){
          vm.raiseCustomADWarn = 'Bad AD format:<br/>Bad data header definition!'
        return;
    }
    if (ad_header.some(r=> Array.from(mapped_aa_properties.keys()).includes(r))){
        vm.raiseCustomADWarn = 'Bad AD format:<br/>Header definitions exist in RV3!<br/>Use different headers.'
        return;
    }
    if (new Set(ad_header).size !== ad_header.length){
        vm.raiseCustomADWarn = 'Bad AD format:<br/>Duplicate header definitions.'
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
        vm.raiseCustomADWarn = 'Bad AD format:<br/>Mismatch between header and data!'
        return;
    }
    let associatedDataArrays = [];
    let associatedDataNames = [];
    var colIndex = 0;
    while (colIndex < headerLength-1) {
        let tempArr = [];
        associatedDataObj.forEach((row) => 
            tempArr.push([Number(row.ix), Number(row.data[colIndex])])
        )
        associatedDataNames.push(ad_header[colIndex+1]);
        associatedDataArrays.push(tempArr);
        colIndex += 1;
    }

    if (topviewer != null && topviewer.viewInstance.uiTemplateService.domainTypes != undefined){
        for (let ix = 0; ix < associatedDataArrays.length; ix++) {
            vm.associated_headers.push(ad_header[ix+1]);
            
            // Analyze data characteristics to determine appropriate colormap
            let dataValues = associatedDataArrays[ix].map(row => row[1]);
            let minVal = Math.min(...dataValues);
            let maxVal = Math.max(...dataValues);
            let dataRange = maxVal - minVal;
            let hasNegativeValues = dataValues.some(val => val < 0);
            let hasPositiveValues = dataValues.some(val => val > 0);
            let isDiverging = hasNegativeValues && hasPositiveValues;
            
            // Select appropriate colormap based on data characteristics
            let selectedColormap;
            if (isDiverging) {
                // Use cool-warm diverging colormap for data that spans negative and positive values
                selectedColormap = window.coolwarm || window.RdBu; // fallback to RdBu if coolwarm not available
            } else {
                // Use jet colormap for continuous data
                selectedColormap = window.jet || window.viridis; // fallback to viridis if jet not available
            }
            
            // Pass the selected colormap to the mapping function
            mapCustomMappingData(associatedDataArrays[ix], ad_header[ix+1], topviewer, selectedColormap);
        }
        displayMappingDataByIndex(topviewer, topviewer.viewInstance.uiTemplateService.domainTypes.length-1);
    }
}
