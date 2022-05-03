export function ajaxProper ({url, type, dataType, postData}={}){
    if (type === 'POST'){
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: type,
                dataType: dataType,
                data: postData,
                traditional: true,
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