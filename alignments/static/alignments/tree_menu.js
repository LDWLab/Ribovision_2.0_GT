 var tree;
 $.ajax({
     url: '/alignments/showTaxonomy',
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
			var
				value1 = button.getAttribute("value2"),
				value2 = this.taxID;
			button.setAttribute("value1", value1);
			button.setAttribute("value2", value2);
		// }
		}
	}
});
function callback(tree){
  new Vue({
	el: '#app',
	data: {
	  tree
	}
  })
}