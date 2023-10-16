import { forEach } from "lodash";
import JSZip from "jszip";
import { saveAs } from "file-saver";
var buildPyMOLscriptASstring = function (pdb_filename) {
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
    console.log('pdb_data',vm.pdbdata);
    var filteredChain = vm.chains.filter(function(itm){
        return vm.chainid[0].indexOf(itm.value) > -1;
    });
    
    var chainText = filteredChain[0].text.replace(/ /g, '_').replace(/'/g, '').replace(/\W/g, '')
    
    vm.chains.text;
    pmlString += `load ./${pdb_filename}\n`;
    //pmlString += vm.pdbdata;
    pmlString += `disable all\n`;
    
    window.selectSections_RV1.forEach(function(val,key){
        var objName = `${key}_${chainText}`.replace(/ /g,'_')
        pmlString += `create ${objName}, cust and chain ${vm.chainid[0]}\n`;
        val.shift();
        val.forEach(function(resiProps){
            let hexColor = rgbToHex(resiProps.color.r,resiProps.color.g,resiProps.color.b);
            //coloringString += `color ${hexColor}, ${objName} and resi ${resiProps.start_residue_number}\n`
            coloringString += `color ${hexColor}, ${objName} and resi ${resiProps.residue_number}\n`
        })
    })
    pmlString += coloringString;
    return pmlString;
}

var downloadPyMOLcustomscript = function () {
    let [month, date, year] = new Date().toLocaleDateString("en-US").split("/");
    let zip = new JSZip();  
    if (vm.user_uploaded_cif_flag == false){
        let pml_filename = `RV_PyMOL_cust-${month}-${date}-${year}.pml`;
        let pdb_filename = `cust.pdb`;
        let pdbString = vm.pdbdata;
        var PyMOLstring = buildPyMOLscriptASstring(pdb_filename);
        zip.file(pml_filename, PyMOLstring);
        zip.file(pdb_filename, pdbString);  
    }
    

    if (vm.user_uploaded_cif_flag == true){
        let pml_filename = `RV_PyMOL_cust-${month}-${date}-${year}.pml`;
        let cif_filename = `cust.cif`;
        let cifString = vm.cifdata;
        var PyMOLstring = buildPyMOLscriptASstring(cif_filename);
        zip.file(pml_filename, PyMOLstring);
        zip.file(cif_filename, cifString);  
    }

    

    if (vm.uploadSession){return;}
    //let anchor = document.createElement('a');
    //anchor.href = 'data:text;charset=utf-8,' + encodeURIComponent(PyMOLstring);
    //anchor.target = '_blank';
    //anchor.download = `PVPyMOL_${vm.pdbid}-${month}-${date}-${year}.pml`;
    //anchor.click();
    zip
        .generateAsync({ type: "blob" })
        .then(function (content) {
        saveAs(content, `PVPyMOL_${vm.pdbid}-${month}-${date}-${year}.zip`);
    });   
}


function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
  }

var rgbToHex = function(r, g, b) {
    return "0x" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

export { buildPyMOLscriptASstring, downloadPyMOLcustomscript }