import {ajaxProper} from './ajaxProper.js'
import parsePdb from 'parse-pdb'
import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function uploadCustomPDB(){
    console.log('uCPDB',  vm.$refs.customPDBfile);
    if (vm.$refs.customPDBfile.files.length == 0){return;}
    vm.PDBparsing = true;
    vm.customPDBsuccess = null;
    vm.customPDBid = null;
    submitCustomPDB(vm.$refs.customPDBfile.files[0]);
    clearInputFile(document.getElementById('uploadCustomPDB'));
}

function submitCustomPDB(file){
    var fr = new FileReader();
    fr.onload = function(){
        if (validatePDB(fr.result)){
            checkAndPopulateChains(fr.result).then (chainID => {
                if (chainID.length > 1){
                    alert("Detected multiple chains! Currently supports only single chain PDBs!");
                    return;
                } else{
                    vm.customPDBid = `cust-1-${chainID[0]}`;
                    postPDBdata("cust", { entityID: "1", chainID: chainID[0], stringData: fr.result });
                }
            }).catch(error => {
                console.log(error)
                alert("Check the PDB format of the uploaded file!")
            })
        }else{
            alert("Check the PDB format of the uploaded file!")
        }
    };
    fr.readAsText(file)
}

function checkAndPopulateChains(strucString){
    //check for chains. For now allow only single chain PDBs
    return new Promise((resolve, reject) => {
        try {
            let parsed = parsePdb(strucString);
            let chains = [];
            parsed.atoms.forEach(function(atom){
                chains.push(atom.chainID)
            });
            var uniqueChains = chains.filter((v, i, a) => a.indexOf(v) === i);
            resolve(uniqueChains)
        }catch(error){
            reject(error)
        }
    })
}

function validatePDB(strucString){
    try {
        let parsed = parsePdb(strucString);
        vm.pdbStart = parsed.atoms[0].resSeq;
        vm.pdbEnd = parsed.atoms[parsed.atoms.length-1].resSeq
        vm.pdbSeq = getSeqFromAtoms(parsed.atoms)
        return true;
    }catch(error){
        console.log(error);
        return false;
    }
}

var getSeqFromAtoms = function (atoms){
    let set  = new Set(atoms.map(item => [item.resSeq, item.resName]).map(JSON.stringify));
    let seqData = Array.from(set).map(JSON.parse);
    var sequence = '';
    seqData.forEach((a)=>{
        try{
            var oneLetter = threeLetterToOne[a[1]];
        }catch{
            var oneLetter = 'X';
        }
        sequence += oneLetter;
    })
    return sequence;
}

var threeLetterToOne = {
    ADE: 'A',
    GUA: 'G',
    CYT: 'C',
    URA: 'U'
   
}

function postPDBdata (pdbID, entities){
    vm.postedPDBEntities = false;
    let parseURL = `custom-struc-data/${pdbID}`;
    var stringEntities = JSON.stringify(entities); 
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        vm.PDBparsing = false;
        if (parsedResponse == "Success!"){
            console.log("Posted PDB data successfully!");
            vm.customPDBsuccess = true;
            getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd, null, vm);
            ajaxProper({
                url: postTopologyURL,
                type: 'POST',
                dataType: 'json'
            }).then (parsedResponse => {
                if (parsedResponse == "Success!"){ 
                    var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rvapi="true"></pdb-rna-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
                }
            }).catch(error => {
                var topview = document.querySelector('#topview');
                vm.topology_loaded = 'error';
                topview.innerHTML = "Failed to generate topology from the structure file!!! <br>Try different PDB."
                console.log(error.responseText);
            });
        }
    }).catch(error => {
        vm.PDBparsing = 'error';
        console.log(error.responseText);
    });
}