
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

Vue.component('treeselect', VueTreeselect.Treeselect,{
	methods:{
		
	}
})

new Vue({
	el: '#phylo_tree_dropdown',
	data: {
		value: null,
		options: null
	},
	methods: {
		limiter(e) {
			if(e.length > 2) {
				alert('You can only select two groups!')
				e.pop()
			}
		},
		loadOptions({action, callback}) {
			if (action === "LOAD_ROOT_OPTIONS"){
				ajax('/alignments/showStrucTaxonomy').then(data =>{
					this.options = data;
					callback();
				}).catch(error => {
					console.log(error)
				})
			}
		}
	}
})
