export function readLoadRV3State (fileInput) {
	var reader = new FileReader();
	reader.onload = function () {
		var uploadedState = JSON.parse(reader.result, reviver);
		vm.uploadSession = true;
		Object.assign(vm.$data, uploadedState);
		let aaPropertiesData = setGlobalProperties();
		window.selectSections_RV1 = uploadedState["window.selectSections_RV1"];
		vm.$nextTick(function(){
			vm.chains = uploadedState.chains;
			vm.$nextTick(function(){
				vm.fasta_data = uploadedState.fasta_data;
				vm.alnobj = uploadedState.alnobj;
				vm.$nextTick(function(){
					var polSele = document.querySelector("#polymerSelect");
					polSele.forEach(function(childElt){
						if (childElt.value == uploadedState.chainid){
							childElt.click();
						}
					})
					window.selectSections_RV1 = uploadedState["window.selectSections_RV1"];
					if (uploadedState.checked_propensities){
						vm.checked_propensities = uploadedState.checked_propensities;
						handlePropensities(vm.checked_propensities);
					}
					if (uploadedState.checked_customMap){
						vm.checked_customMap = uploadedState.checked_customMap;
						window.tempCSVdata = uploadedState.csv_data;
					}
					vm.uploadSession = false;
				});
			});
		});
	};
	reader.readAsBinaryString(fileInput);
};