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

Vue.component('tree-menu', { 
	delimiters: ['[[',']]'],
	template: '#tree-menu',
	props: [ 'children', 'label', 'depth', 'id' ],
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
		return { 'has-children': this.children }
	  },
	  indent() {
		return { transform: `translate(${this.depth * 50}px)` }
	  }
	},
	methods: {
		toggleChildren() {
		this.showChildren = !this.showChildren;
		// if (this.showChildren) {
			//var newValue = "{% url 'alignments:detail' " + this.taxID + " %}"
			setTaxID2(getTaxID1());
			setTaxID1(this.id);
			prepareMethod1();

			// var
			// 	form = document.getElementById("mainForm"),
			// 	taxID1 = form.getAttribute("taxID2"),
			// 	taxID2 = this.taxID;
			// form.setAttribute("taxID1", taxID1);
			// form.setAttribute("taxID2", taxID2);
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