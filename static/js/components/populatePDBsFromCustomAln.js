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
			setTimeout(repeatingFunc(jobid, qLength), 1000);
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
					tempPDB.push(pdb)
					filteredPDBs.set(pdb, [chain]);
				}
			}
		})
		if (tempPDB.length > 0){
			vm.blastPDBresult.push(...tempPDB.sort());
			vm.blastMAPresult = filteredPDBs;
			vm.fetchingPDBwithCustomAln = 'complete';
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


