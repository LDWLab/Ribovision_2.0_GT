import { forEach } from "lodash";

var buildPyMOLscriptASstring = function () {
    if (vm.uploadSession){return;}
    var pmlString = new String("\
    set hash_max, 500\n\
    set cartoon_loop_radius,0.4\n\
    set cartoon_tube_radius,1\n\
    set cartoon_ladder_radius,0.6\n\
    set cartoon_oval_length,1.4\n\
    set cartoon_oval_width,0.6\n\
    set ray_opaque_background, off\n\
    bg_color black\n\
    set ray_trace_mode,1\n\
    set ray_shadows,0\n");
    var coloringString = new String();

    var filteredChain = vm.chains.filter(function(itm){
        return vm.chainid[0].indexOf(itm.value) > -1;
    });
    
    var chainText = filteredChain[0].text.replace(/ /g, '_').replace(/'/g, '').replace(/\W/g, '')

    vm.chains.text;
    pmlString += `fetch ${vm.pdbid}, async=0\n`;
    pmlString += `disable all\n`;
    window.selectSections_RV1.forEach(function(val,key){
        var objName = `${key}_${chainText}`.replace(/ /g,'_')
        pmlString += `create ${objName}, ${vm.pdbid} and chain ${vm.chainid[0]}\n`;
        val.shift();
        let colorMap = {};
        val.forEach(function(resiProps){
            let hexColor = rgbToHex(resiProps.color.r,resiProps.color.g,resiProps.color.b);
            //coloringString += `color ${hexColor}, ${objName} and resi ${resiProps.start_residue_number}\n`
            if (colorMap.hasOwnProperty(hexColor)) {
                
                colorMap[hexColor].push(resiProps.residue_number);
            }
            else{
                colorMap[hexColor] = [resiProps.residue_number]

            }
            // coloringString += `color ${hexColor}, ${objName} and resi ${resiProps.residue_number}\n`
        });

        for(let [hexColor, values] of Object.entries(colorMap)){
            // let stringDef = `color ${hexColor}, ${objName}`;
            let auxString = `color ${hexColor}, ${objName} and (`;
            // coloringString += 

            // Find the locations where there are jumps in indices
            let continous = [];
            let uniqueValues = [...new Set(values)];

            if (uniqueValues.length == 0){
                continue;
            }

            // 43, 265, 701
            for (let i = 0; i < uniqueValues.length; i++) {
                if(continous.length == 0){
                    continous.push(uniqueValues[i]);
                }
                else if(uniqueValues[i] - 1 == continous[continous.length-1]){
                    continous.push(uniqueValues[i]);
                }
                else{
                    if(auxString.slice(-2) != " ("){
                        auxString += ' or '; 
                    }

                    if(continous.length == 1){
                        // coloringString += `${stringDef} and resi ${continous[0]}\n`;
                        auxString += `resi ${continous[0]}`;
                    }
                    else{
                        auxString += `resi ${continous[0]}-${continous[continous.length-1]}`;
                    }

                    continous = [uniqueValues[i]];
                }
            }

            if(continous.length != 0){
                if(auxString.slice(-2) != " ("){
                    auxString += ' or '; 
                }
            
                if(continous.length == 1){
                    // coloringString += `${stringDef} and resi ${continous[0]}\n`;
                    auxString += `resi ${continous[0]}`;
                }
                else{
                    auxString += `resi ${continous[0]}-${continous[continous.length-1]}`;
                }
                
            }
            coloringString += (auxString + ")\n");

        }

        console.log(key, JSON.stringify(colorMap));
    })


    pmlString += coloringString;
    
    return pmlString;
}

var downloadPyMOLscript = function () {
    let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
    var PyMOLstring = buildPyMOLscriptASstring();
    if (vm.uploadSession){return;}
    let anchor = document.createElement('a');
    anchor.href = 'data:text;charset=utf-8,' + encodeURIComponent(PyMOLstring);
    anchor.target = '_blank';
    anchor.download = `PVPyMOL_${vm.pdbid}-${month}-${date}-${year}.pml`;
    anchor.click();
}

function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
  }

var rgbToHex = function(r, g, b) {
    return "0x" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

export { buildPyMOLscriptASstring, downloadPyMOLscript }