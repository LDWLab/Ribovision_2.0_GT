export function uploadCustomFullSequence(){
    if (vm.$refs.customFullSequence.files.length == 0){return;}
    var fr = new FileReader();
    fr.onload = function() {
        vm.customFullSequence = fr.result;
    }
    fr.readAsText(vm.$refs.customFullSequence.files[0]);
}