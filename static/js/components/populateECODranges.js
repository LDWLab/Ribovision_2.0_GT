export function populateECODranges (pdbid, chainid) {
	vm.domain_list = [];
	$.ajax ({
		type: "GET",
		url: `/desire-api/ECOD-domains/?pdb=${pdbid}&chain=${chainid[0]}`,
		//url: "/alignments/authEcodQuery",
		//data: {url: `/desire-api/ECOD-domains/?pdb=${pdbid}&chain=${chainid[0]}`},
		success: function (data){
			for (var i = 0; i < data.results.length; i++) {
				if (data.results[i].chain == chainid[0]){
					var range_str = null;
					var fetchedRange = data.results[i].pdb_range;
					if (fetchedRange.split(';').length > 1){
						//add handling of multiple ranges with masking the positions between them
						let rangeArray = fetchedRange.split(';');
						range_str = '';
						let re = /\d+-\d+/
						for(var j in rangeArray) {
                            range_str += re.exec(rangeArray[j]) + ';'
                        }
					}
					if(!range_str) {
						let re = /\d+-\d+$/;
						range_str = re.exec(fetchedRange)[0] + ';';
					}
					let xName = data.results[i].x_name.replace("NO_X_NAME", "");
					let fName = data.results[i].f_name.replace("NO_F_NAME", "");
					let domName = `${xName} ${fName}`.trim()
					console.log(range_str)
					vm.domain_list.push({name: domName, range: range_str});
				}
			}
			if (vm.domain_list.length == 0){
				vm.domain_list.push({name: "No ECOD match!", range: null});
			}
		}, error: function(xhr, status, error){
			vm.domain_list.push({name: "Failed getting ECOD domains!", range: null});
			var errorMessage = xhr.status + ': ' + xhr.statusText
			console.log('Error - ' + errorMessage);
		}
	});
};