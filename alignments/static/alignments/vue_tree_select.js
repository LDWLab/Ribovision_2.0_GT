
function ajax(url) {
	return new Promise((resolve, reject) => {
		$.ajax({
			url: url,
			type: 'GET',
			dataType: "json",
			success: function(data) {
				resolve([data])
			},
			error: function (error) {
				console.log(`Error ${error}`);
				reject(error)
			}
		})
	})
}

const iterate = (obj, output_obj) => {
    Object.keys(obj).forEach(key => {
    console.log('key: '+ key + ', value: '+obj[key]);
    if (typeof obj[key] === 'object Object') {
            iterate(obj[key])
        }
    })
}


Vue.component('treeselect', VueTreeselect.Treeselect,)

new Vue({
	el: '#phylo_tree_dropdown',
	delimiters: ['[[',']]'],
	data: {
		value: null,
		valuef: null,
		options: null,
		items: null
	},
	methods: {
		limiter(e) {
			if(e.length > 1) {
				alert('You can only select one group!')
				e.pop()
			}
		},
		loadOptions({action, callback}) {
			if (action === "LOAD_ROOT_OPTIONS"){
				ajax('/alignments/showStrucTaxonomy').then(data =>{
					data[0].isDisabled = true,
					this.options = data;
					callback();
				}).catch(error => {
					console.log(error)
				})
			}
		},
		loadData:function(value) {
			this.items = null;
			ajax('/alignments/fold-api/'+value)
			.then(data=>{
				//Fix up our data
				var fpa = data[0]["Folds to polymers to alignments"]
				var fpa_viz = [];
				Object.keys(fpa).forEach(fkey => {
					Object.keys(fpa[fkey]).forEach(pkey => {
						fpa[fkey][pkey].forEach(function (akey){
							fpa_viz.push({
								text:  akey[1],
								value: pkey.concat(',',akey)
							});
						});
					});
				});
				var temp_arr = fpa_viz
				fpa_viz = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
				this.items = fpa_viz
			});
		}
	}
})
