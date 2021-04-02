import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function loadViewersWithCustomUploadStructure(){

    var pdbid = vm.customPDBid.split('-')[0];
    var chainid = vm.customPDBid.split('-')[2];
    var entityid = vm.customPDBid.split('-')[1];
    //topviewer with custom url
    getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, 1, 100000, null, vm)
    vm.showPDBViewer(pdbid, chainid, entityid);

}

