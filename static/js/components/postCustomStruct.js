import {ajaxProper} from './ajaxProper.js'
export function postCIFdata (pdbID, entities){
    vm.postedPDBEntities = false;
    let parseURL = `custom-struc-data/${pdbID}`;
    var stringEntities = JSON.stringify(entities); 
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        if (parsedResponse == "Success!"){
            vm.postedPDBEntities = true;
        }
    }).catch(error => {
        var topview = document.querySelector('#topview');
        console.log(error);
        vm.topology_loaded = 'error';
        topview.innerHTML = "Failed to POST structure to our server!<br>Try refreshing the webpage."
    });
}