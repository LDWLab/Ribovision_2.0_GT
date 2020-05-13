function get_taxonomy_api(taxID) {
	$.ajax({
		url: '/alignments/showTaxonomy-api/'+taxID,
		type: "GET",
		dataType: "json",
		success: function (data) {
			callback(data);
		},
		error: function (error) {
			console.log(`Error ${error}`);
		}
	});
}

function callback(json_data){
	if (json_data.nodes.length > 0){
		var index=0;
		var catOptions = "<option value="+json_data.taxID+">Select taxonomic group</option>";
		var grouplabel = json_data.label.concat('-').concat(json_data.level)
		catOptions += "<optgroup class=parent_level label="+grouplabel+" value="+json_data.taxID+" selected=selected>";
		while (index < json_data.nodes.length) { 
			current_taxid = json_data.nodes[index][1]
			current_label = json_data.nodes[index][0]
			catOptions+= "<option onclick value="+current_taxid+">"+current_label+"</option>"
			index++; 
		}
		catOptions += "</optgroup>"
		catOptions += "<option onclick value="+json_data.parent+">Go back</option>";
		document.getElementById("dropdown").innerHTML = catOptions;
	}
}

window.onload = function(){
	get_taxonomy_api(0)
};