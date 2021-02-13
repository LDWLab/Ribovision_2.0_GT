export function testingCIFParsing (){
    let ebiURL = 'https://www.ebi.ac.uk/pdbe/coordinates/4v9d/chains?entityId=28';
    ajaxProper({
        url: ebiURL,
        type: 'GET',
        dataType: 'text',
    }).then(strucData => {
        let parseURL = `custom-struc-data/4V9D_28`
        ajaxProper({
            url: parseURL,
            type: 'POST',
            dataType: 'text',
            postData: {"custom_structure": strucData}
        }).then (parsedResponse => {
            console.log(parsedResponse);
        }).catch(error => {
            console.log(error);
        });
    }).catch(error => {
        console.log(error);
    });
}

var ajaxProper = function ({url, type, dataType, postData}={}){
    if (type === 'POST'){
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: type,
                dataType: dataType,
                data: postData,
                headers: {'X-CSRFToken': csrftoken},
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    } else if (type === 'GET'){
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: type,
                dataType: dataType,
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    }
}