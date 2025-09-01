// Import colormaps
import { jet, coolwarm, viridis } from '../../alignments/colormaps.js';

export function customCSVhandler(csv_data) {
    if (vm.uploadSession){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    cleanCustomMap(vm.checked_customMap);
    vm.raiseCustomCSVWarn = null;
    vm.custom_headers = [];
    vm.AD_headers = [];
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
        customDataNames.push(custom_header[colIndex+1]);
        customDataArrays.push(tempArr);
        colIndex += 1;
    }
    
    if (topviewer != null && topviewer.viewInstance.uiTemplateService.domainTypes != undefined){
        for (let ix = 0; ix < customDataArrays.length; ix++) {
            vm.custom_headers.push(custom_header[ix+1]);
            
            // Analyze data characteristics to determine appropriate colormap
            let dataValues = customDataArrays[ix].map(row => row[1]);
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
            mapCustomMappingData(customDataArrays[ix], custom_header[ix+1], topviewer, selectedColormap);
        }
        //displayMappingDataByIndex(topviewer, topviewer.viewInstance.uiTemplateService.domainTypes.length-1);
    }
}
