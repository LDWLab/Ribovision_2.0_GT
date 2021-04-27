export function populatePDBsFromCustomAln (firstSeq) {
	vm.fetchingPDBwithCustomAln = true;
	$.ajax({
		url: "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run",
		beforeSend: function(xhr) { 
		  xhr.setRequestHeader('Content-Type','application/x-www-form-urlencoded'); 
		  xhr.setRequestHeader('Accept','text/plain'); 
		},
		type: 'POST',
		contentType: 'application/json',
		data: `email=anton.petrov%40biology.gatech.edu&program=blastp&stype=protein&sequence=${firstSeq}&database=pdb`,
		success: function (data){
			repeatingFunc(data, firstSeq.length);
		},
		error: function(error){
			vm.fetchingPDBwithCustomAln = 'error';
			console.log(error);
		}
	});
};

var repeatingFunc = function(jobid, qLength) {
	ebiAjax(`https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/${jobid}`).then(jobStatus=>{
		if (jobStatus == 'RUNNING'){
			setTimeout(repeatingFunc, 2500, jobid, qLength);
		} else if (jobStatus == 'FINISHED'){
			fetchBLASTresult(jobid, qLength);
		} else {
			vm.fetchingPDBwithCustomAln = 'error';
			console.log(`Error with blast job! Job finished with status ${jobStatus}`);
		}
	}).catch(error => {
		vm.fetchingPDBwithCustomAln = 'error';
		console.log(`Error with status blast job! Error: ${error}`);
	})
}

var fetchBLASTresult = function (jobID, qLength){
	ebiAjax(`https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/${jobID}/out?format=10`).then(csvResult=>{
		let csvArr = csvResult.split(/\n/g);
		var filteredPDBs = new Map();
		var tempPDB = new Array();
		csvArr.forEach(function(blastRow){
			let blastEval = Number(blastRow.split(',')[10]);
			let blastLength = Number(blastRow.split(',')[3]);
			let blastCoverage = blastLength/qLength;
			if (blastEval < 0.00001 && blastCoverage > 0.75){
				var pdb = blastRow.split(',')[1].split(/:|_/)[1];
				var chain = blastRow.split(',')[1].split(/:|_/)[2];
				if (filteredPDBs.has(pdb)){
					filteredPDBs.get(pdb).push(chain);
				} else {
					tempPDB.push({pdb:pdb, coverage:blastCoverage, eval:blastEval})
					filteredPDBs.set(pdb, [chain]);
				}
			}
		})
		if (tempPDB.length > 0){
			tempPDB.sort((a, b) => parseFloat(b.blastCoverage) - parseFloat(a.blastCoverage) || parseFloat(a.blastEval) - parseFloat(b.blastEval));
			vm.blastMAPresult = filteredPDBs;
			fetchAndParsePDBnames(constructRCSBGraphQuery(tempPDB), tempPDB);
		} else {
			vm.fetchingPDBwithCustomAln = 'none';
		}
	}).catch(error => {
		vm.fetchingPDBwithCustomAln = 'error';
		console.log(`Error with result blast job! Error: ${error}`);
	})
}

var ebiAjax = function (url){
	return new Promise((resolve, reject) => {
		$.ajax({
			url: url,
			type: 'GET',
			dataType: "text",
			success: function(data) {
				resolve(data);
			},
			error: function(error) {
				console.log(`Error ${error}`);
				vm.fetchingPDBwithCustomAln = 'error';
				reject(error);
			}
		})
	})
}

var constructRCSBGraphQuery = function (pdblist){
	var queryString = '{polymer_entities(entity_ids: ["';
	pdblist.forEach(function(pdb){
		queryString += `${pdb.pdb}_1","`
	});
	var query = queryString.slice(0,-1);
	query += ']) {rcsb_entity_source_organism {ncbi_scientific_name}}}'
	return encodeURIComponent(query);
}

var fetchAndParsePDBnames = function(query, tempPDB){
	$.ajax({
		url: `https://data.rcsb.org/graphql?query=${query}`,
		type: 'GET',
		contentType: 'json',
		success: function (data){
			var tempNames = [];
			data.data.polymer_entities.forEach(function(pol){
				if (pol.rcsb_entity_source_organism && pol.rcsb_entity_source_organism[0]){
					tempNames.push(pol.rcsb_entity_source_organism[0].ncbi_scientific_name);
				} else {
					tempNames.push('');
				}
			})
			combineAndAssignNames(tempNames, tempPDB);
		},
		error: function(error){
			var tempNames = [];
			tempPDB.forEach(function(pol){
				tempNames.push('');
			});
			combineAndAssignNames(tempNames, tempPDB);
			console.log(error);
		}
	});
}

var combineAndAssignNames = function(tempNames, tempPDB){
	var namedPDBs = tempPDB.map(function (entry, index){
		return { id:entry.pdb, name: `${entry.pdb} ${tempNames[index]}`, stats: `E-value: ${entry.eval}\nCoverage: ${entry.coverage.toFixed(2)}`}
	 });
	vm.blastPDBresult.push(...namedPDBs);
	vm.fetchingPDBwithCustomAln = 'complete';
}

//qaccver		EMBOSS_001,
//saccver		PDB:4V8S_BC,
//pident		48.507,
//length		402,
//mismatch		187,
//gapopen		4,
//qstart		1,
//qend			401,
//sstart		8,
//send			390,
//evalue		3.73e-115,
//bitscore		344
//https://www.ebi.ac.uk/Tools/common/tools/help/index.html?tool=ncbiblast
//https://data.rcsb.org/#gql-api
