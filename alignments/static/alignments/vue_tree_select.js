
$.ajax({
	url: '/alignments/showStrucTaxonomy',
	type: "GET",
	dataType: "json",
	success: function (data) {
		console.log(data);
		callback(data)
	},
	error: function (error) {
		console.log(`Error ${error}`);
	}
});

Vue.component('treeselect', VueTreeselect.Treeselect)

function callback(tree){
	new Vue({
		el: '#phylo_tree_dropdown',
		data: {
			value: null,
			options:[
				tree
			]
		},
		methods: {
			limiter(e) {
				if(e.length > 2) {
					alert('You can only select two groups!')
					e.pop()
				}
			},
		}
	})
}
