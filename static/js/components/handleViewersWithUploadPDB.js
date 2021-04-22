

export function loadViewersWithCustomUploadStructure(){

    var pdbid = vm.customPDBid.split('-')[0];
    var chainid = vm.customPDBid.split('-')[2];
    var entityid = vm.customPDBid.split('-')[1];
    vm.entityID = entityid;
    vm.chainid = [chainid];
    vm.showPDBViewer(pdbid, chainid, entityid);

}

