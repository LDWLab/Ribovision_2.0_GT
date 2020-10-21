import $ from jquery

function ajax(url) {
    return new Promise((resolve, reject) => {
        $.ajax({
            url: url,
            type: 'GET',
            dataType: "json",
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

ajax("bla")