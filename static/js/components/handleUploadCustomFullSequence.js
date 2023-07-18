import {ajaxProper} from './ajaxProper.js'
import parsePdb from 'parse-pdb'
import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function uploadCustomFullSequence(hardcoded_structure = ""){
    if (vm.$refs.customFullSequence.files.length == 0){return;}
    var fr = new FileReader();
    fr.onload = function() {
        vm.customFullSequence = fr.result;
    }
    fr.readAsText(vm.$refs.customFullSequence.files[0]);
}