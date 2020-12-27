import {ajaxProper} from './ajaxProper.js'
export function testingCIFParsing (pdbID, entityIDS){
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
        console.log(error);
    });
}