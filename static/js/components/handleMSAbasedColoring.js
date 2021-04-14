export function colorByMSAColorScheme(scheme, vm) {
	if (vm.selected_property){
		vm.selected_property = null;
		viewerInstanceTop.pluginInstance.resetDisplay();
		//Reset the 3D colors?
	}
	if (vm.chains.length > 0){
		var tempEntity = vm.chains.filter(obj => {
			return obj["entityID"] == vm.entityID;
		});
	} else {
		//assume customPDB
		var tempEntity = [{ entityID: vm.entityID, startIndex: vm.pdbStart, endIndex: vm.pdbEnd, sequence: vm.pdbSeq }]
		return;
	}
	var currScheme = vm.schemesMgr.getScheme(scheme);
	var colorData2D = [];
	var colorData3D = [{entity_id: `${vm.entityID}`, focus: true}];
	for (var i=1; i < tempEntity[0].sequence.length; i++) {
		if (i < tempEntity[0].startIndex || i > tempEntity[0].endIndex){
			colorData2D.push({color: 'rgb(180,180,180)', start: i, end: i,
			});
			colorData3D.push({
				entity_id: `${vm.entityID}`, start_residue_number: i, end_residue_number: i,
				color: {r:180, g:180, b:180}
			})
		} else {
			let rgbColors = hexToRgb(currScheme.getColor(tempEntity[0].sequence[i-1]));
			colorData2D.push({
				start: i, end: i, color: `rgb(${rgbColors.r},${rgbColors.g},${rgbColors.b})`,
			});
			colorData3D.push({
				entity_id: `${vm.entityID}`, start_residue_number: i, end_residue_number: i,
				color: {r: rgbColors.r, g: rgbColors.g, b: rgbColors.b}
			})
		}
	}
	vm.colorSchemeData = [colorData2D,colorData3D];
	viewerInstanceTop.pluginInstance.updateTheme(colorData2D);
	viewerInstance.visual.select({ data: colorData3D, nonSelectedColor: {r:180,g:180,b:180}})
}