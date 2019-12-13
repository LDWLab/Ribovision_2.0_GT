 var tree;
 $.ajax({
     url: '/alignments/showTaxonomy',
     type: "GET",
 	dataType: "json",
     success: function (data) {
 		console.log(data);
 		tree = data;
     },
     error: function (error) {
         console.log(`Error ${error}`);
	 },
	 async: false		//Essentially means everything after this will wait for the ajax call to happen,
 });					//which is precisely what we wanted to do. 
						//We need to be careful how long the python function takes to construct the tree.
						//In principle not a good way to do it, but it works for now.

 console.log('Tree is in the JS global scope: ', tree);

Vue.component('tree-menu', { 
	delimiters: ['[[',']]'],
	template: '#tree-menu',
	props: [ 'nodes', 'label', 'depth', 'taxID' ],
	data() {
	   return {
		 showChildren: false
	   }
	},
	computed: {
	  iconClasses() {
		return {
		  'fa-plus-square-o': !this.showChildren,
		  'fa-minus-square-o': this.showChildren
		}
	  },
	  labelClasses() {
		return { 'has-children': this.nodes }
	  },
	  indent() {
		return { transform: `translate(${this.depth * 50}px)` }
	  }
	},
	methods: {
		toggleChildren() {
		this.showChildren = !this.showChildren;
		// if (this.showChildren) {
			var button = document.getElementById("getAlignment");
			//var newValue = "{% url 'alignments:detail' " + this.taxID + " %}"
			var newValue = this.taxID;
			button.setAttribute("value2", newValue);
		// }
		}
	}
});
  
  new Vue({
	el: '#app',
	data: {
	  tree
	}
  })