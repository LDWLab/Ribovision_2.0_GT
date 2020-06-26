function customFilter(object, result, key, value){
	if(object.hasOwnProperty(key) && object[key] == value)
		result.push(object);
	for(var i=0; i<Object.keys(object).length; i++){
		if(typeof object[Object.keys(object)[i]] == "object"){
			customFilter(object[Object.keys(object)[i]], result, key, value);
		}
	}
}
function GetRangeMapping(pdbid, chainid, range, mapping){
	console.log(range)
	$.ajax({
	url: 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/'+pdbid,
	type: "GET",
	dataType: "json",
	success: function (data) {
		console.log(data);
		var result = []
		customFilter(data, result, "chain_id", chainid);
		result = result[0]
		if(result != null) {
			var pdb_start = parseInt(result["start"]["residue_number"]);
			var pdb_end = parseInt(result["end"]["residue_number"]);
			var uniprot_start = parseInt(result["unp_start"]);
			for (residue_number_str of range.split("-")){
				var residue_number = parseInt(residue_number_str);
				if(residue_number >= pdb_start && residue_number <= pdb_end){
					offset = uniprot_start - pdb_start;
					mapping.push(residue_number - offset);
				}else{
					mapping.push(residue_number);
				}
			}
			return mapping
		}else{
			console.log("No mapping for pdb "+pdbid+" and chain"+ chainid)
			return [range.split("-")[0],range.split("-")[1]]
		}
	},
	error: function (error) {
		console.log(`Error ${error}`);},
		 async: false
	});
}