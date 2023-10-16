import {ajaxProper} from './ajaxProper.js'
import parsePdb from 'parse-pdb'
import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function uploadCustomPDB(full_sequence_from_pdb = ""){
    if (!vm.customFullSequence) {
        return;
    }
    vm.user_uploaded_cif_flag = false;
    //console.log('uCPDB',  vm.$refs.customPDBfile);
    //console.log('FS', vm.customFullSequence);
    if (vm.$refs.customPDBfile.files.length == 0){return;}
    vm.PDBparsing = true;
    vm.customPDBsuccess = null;
    vm.customPDBid = null;
    
    submitCustomPDB(vm.$refs.customPDBfile.files[0], vm.customFullSequence);
    clearInputFile(document.getElementById('uploadCustomPDB'));
}

function submitCustomPDB(file, full_sequence_from_pdb = ""){
    var fr = new FileReader();
    postFullSeq(full_sequence_from_pdb);
    console.log('FSFPDB', full_sequence_from_pdb )
    fr.onload = function(){
        if (validatePDB(fr.result)){
            checkAndPopulateChains(fr.result).then (chainID => {
                if (chainID.length > 1){
                    alert("Detected multiple chains! Currently supports only single chain PDBs!");
                    return;
                } else{
                    vm.customPDBid = `cust-1-${chainID[0]}`;
                    vm.chains = [{value: chainID[0], text: chainID[0]}];
                    //console.log(fr.result);
                    vm.pdbdata = fr.result;
                    vm.pdbcust=true;
                    postPDBdata("cust", { entityID: "1", chainID: chainID[0], stringData: fr.result, full_sequence_from_pdb });
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

function postPDBdata (pdbID, entities, full_sequence_from_pdb = ""){
    vm.postedPDBEntities = false;
    let parseURL = `custom-struc-data-pdb/${pdbID}`;
    var stringEntities = JSON.stringify(entities); 
    console.log("PostPDB_FS", full_sequence_from_pdb);
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        vm.PDBparsing = false;
        if (parsedResponse == "Success!") {
            vm.customPDBsuccess = true;
            getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd, null, vm, full_sequence_from_pdb);
        } else if ("successFlag" in parsedResponse && parsedResponse.successFlag) {
            vm.cif_file_path = parsedResponse.cif_file_path;
            vm.customPDBsuccess = true;
            getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd, null, vm, full_sequence_from_pdb);
        }
    }).catch(error => {
        vm.PDBparsing = 'error';
        console.log(error.responseText);
    });
}

function postFullSeq (FullSeq){
    vm.postedFullSeq = false;
    let URL = `custom-struc-full-seq/`;
    var stringSeq = JSON.stringify(FullSeq); 
    ajaxProper({
        url: URL,
        type: 'POST',
        dataType: 'json',
        postData: {"sequence": stringSeq}
    }).then (parsedResponse => {
        vm.PDBparsing = false;
        if (parsedResponse == "Success!"){
            vm.postedFullSeq = true;
            
        }
    }).catch(error => {
        vm.PDBparsing = 'error';
        console.log(error.responseText);
    });
}
