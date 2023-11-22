import {ajaxProper} from './ajaxProper.js'
import parsePdb from 'parse-pdb'
import parseMmcif from 'parse-mmcif'
import {getStructMappingAndTWC} from './getStructMappingAndTWC.js'

export function uploadCustomPDB(){
    console.log("uploadCustomPDB()");
    vm.user_uploaded_cif_flag = false;
    console.log('uCPDB',  vm.$refs.customPDBfile);
    if (vm.$refs.customPDBfile.files.length == 0){return;}
    vm.PDBparsing = true;
    vm.customPDBsuccess = null;
    vm.customPDBid = null;
    //vm.customChain='a';
    //submitCustomPDB(vm.$refs.customPDBfile.files[0]);
    submitCustomCIF(vm.$refs.customPDBfile.files[0]);
    clearInputFile(document.getElementById('uploadCustomPDB'));
}

export function uploadCustomCIF() {
    console.log("uploadCustomCIF()");
    vm.user_uploaded_cif_flag = true;
    if (vm.$refs.customCIFfile.files.length === 0) {
        return;
    }
    const entity_id = Number.parseInt(vm.$refs.entity_id.value);
    if (Number.isNaN(entity_id)) {
        return;
    }
    vm.PDBparsing = true;
    vm.customPDBsuccess = null;
    vm.customPDBid = null;
    vm.customChain='a';
    vm.customEntity=vm.$refs.entity_id.value;
    submitCustomCIF(vm.$refs.customCIFfile.files[0], entity_id);
    clearInputFile(document.getElementById('uploadCustomCIF'));
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
                    vm.customPDBid = `cust-1-${vm.customChain}`;
                    postPDBdata("cust", { entityID: "1", chainID: vm.customChain, stringData: fr.result });
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

function submitCustomCIF(file, entity_id = -1){
    console.log("entity_id", entity_id);
    var fr = new FileReader();
    fr.onload = function(){
        if (validateCIF(fr.result)){
            checkAndPopulateChainsCIF(fr.result).then (chainID => {
                if (chainID.length > 1000){
                    alert("Detected multiple chains! Currently supports only single chain CIFs!");
                    return;
                } else{
                    vm.customPDBid = `cust-${entity_id}-${vm.customChain}`;
                    postPDBdata("cust", { entityID: `${entity_id}`, chainID: vm.customChain, stringData: fr.result });
                    vm.chains = [{value: vm.customChain, text: vm.customChain}];
                    vm.cifdata = fr.result;
                    console.log(fr.result);
                    
                    vm.cifcust=true;
                    //postPDBdata("cust", { entityID: "67", chainID: vm.customChain, stringData: str(fr) });
                }
            }).catch(error => {
                console.log(error)
                alert("Check the CIF format of the uploaded file!!")
            })
        }else{
            alert("Check the CIF format of the uploaded file!!!")
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

function checkAndPopulateChainsCIF(strucString){
    //check for chains. For now allow only single chain PDBs
    return new Promise((resolve, reject) => {
        try {
            let parsed = parseMmcif(strucString);
            console.log('parsed1', parsed);
            let chains = [];
            parsed.atoms.forEach(function(atom){
                chains.push(atom.auth_asym_id)

            });
            //console.log('chains', chains)
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
        //console.log('parsed', parsed)
        vm.pdbStart = parsed.atoms[0].resSeq;
        vm.pdbEnd = parsed.atoms[parsed.atoms.length-1].resSeq
        vm.pdbSeq = getSeqFromAtoms(parsed.atoms)
        return true;
    }catch(error){
        console.log(error);
        return false;
    }
}

function validateCIF(strucString){
    var RNA_chain=vm.customChain;
    try {
        let parsed = parseMmcif(strucString);
        let RNA_atoms=[];
        parsed.atoms.forEach(function(atom){
            if  (atom.label_entity_id == vm.customEntity){
            RNA_atoms.push(atom);
            vm.customChain=atom.auth_asym_id;
            }
        });

        //console.log('RNA_atoms', RNA_atoms)
        vm.pdbStart = RNA_atoms[0].auth_seq_id;
        //vm.pdbEnd = parsed.atoms[parsed.atoms.length-1].auth_seq_id;
        vm.pdbEnd = RNA_atoms[RNA_atoms.length-1].auth_seq_id;
        //vm.pdbSeq = getSeqFromAtomsCIF(parsed.atoms);
        //vm.pdbSeq = getSeqFromAtomsCIF(RNA_atoms);
        //console.log(vm.pdbSeq);
        //vm.pdbSeq = 'GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACXAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCUU';
        console.log('PDB_SEQ',vm.pdbSeq);
        vm.pdbSeq = getSeqFromAtomsCIF(RNA_atoms);
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
var getSeqFromAtomsCIF = function (atoms){
    let set  = new Set(atoms.map(item => [item.auth_seq_id, item.auth_comp_id]).map(JSON.stringify));
    let seqData = Array.from(set).map(JSON.parse);
    //console.log('seqData',set, seqData);
    var sequence = '';
    seqData.forEach((a)=>{
        try{
            //var oneLetter = threeLetterToOne[a[1]];
            //console.log('letter', a[1]), threeLetterToOne[a[1]];
            var oneLetter = a[1];
        }catch{
            var oneLetter = 'X';
        }
        sequence += oneLetter;
    })
    console.log('sequence', sequence);
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
    let parseURL = `custom-struc-data-cif/${pdbID}`;
    var stringEntities = JSON.stringify(entities); 
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        vm.PDBparsing = false;
        if (parsedResponse == "Success!") {
            vm.customPDBsuccess = true;
            getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd, null, vm);
        } else if ("successFlag" in parsedResponse && parsedResponse.successFlag) {
            vm.cif_file_path = parsedResponse.cif_file_path;
            vm.customPDBsuccess = true;
            getStructMappingAndTWC (vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd, null, vm);
        }
    }).catch(error => {
        vm.PDBparsing = 'error';
        console.log(error.responseText);
    });
}
/*
function postPDBdata (pdbID, entities){
    vm.postedPDBEntities = false;
    let parseURL = `custom-struc-data/${pdbID}`;
    //let postTopologyURL = `proOrigamiPOSTTopology/${pdbID}-${entities.entityID}-${entities.chainID}`;
    var stringEntities = JSON.stringify(entities); 
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        vm.PDBparsing = false;
        if (parsedResponse == "Success!"){
            console.log("CIF data posted successfully!");
            vm.customPDBsuccess = true;
            console.log('TWC', vm.fasta_data, vm.customPDBid, vm.pdbStart, vm.pdbEnd)
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
*/
