import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function loadViewersWithCustomUploadStructure(){

    var pdbid = vm.customPDBid.split('-')[0];
    var chainid = vm.customPDBid.split('-')[2];
    var entityid = vm.customPDBid.split('-')[1];
    vm.entityID = entityid;
    vm.chainid = [chainid];
    var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid} pvapi="true" filter-range=1,100000></pdb-topology-viewer>`
    document.getElementById('topview').innerHTML = topology_viewer;
    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
    getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, 1, 100000, null, vm);
    vm.showPDBViewer(pdbid, chainid, entityid);

}

