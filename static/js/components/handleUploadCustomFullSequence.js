import {uploadCustomPDB} from './handleUploadPDB_RNA.js'

export function uploadCustomFullSequence(){
    if (vm.$refs.customFullSequence.files.length == 0){return;}
    var fr = new FileReader();
    fr.onload = function() {
        let fasta_content = fr.result;
        vm.customFullSequence = fasta_content.split("\n").filter((line) => !/^>/.test(line)).join("");
        uploadCustomPDB();
    }
    fr.readAsText(vm.$refs.customFullSequence.files[0]);
}