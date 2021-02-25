import {ajaxProper} from './ajaxProper.js'
export function postCIFdata (pdbID, entityIDS){
    vm.postedPDBEntities = false;
    let parseURL = `custom-struc-data/${pdbID}`;
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entityIDS": entityIDS}
    }).then (parsedResponse => {
        if (parsedResponse == "Success!"){
            vm.postedPDBEntities = true;
        }
    }).catch(error => {
        var topview = document.querySelector('#topview');
        console.log(error);
        this.topology_loaded = 'error';
        topview.innerHTML = "Failed to POST structure to our server!<br>Try refreshing the webpage."
    });
}