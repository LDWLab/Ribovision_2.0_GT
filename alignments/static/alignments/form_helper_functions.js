function setTaxID1(tax1) {
	taxID1 = tax1;
}

function setTaxID2(tax2) {
	taxID2 = tax2;
}

function setminIndex(minin){
	window.minIndex = minin;
}

function setmaxIndex(maxin){
	window.maxIndex = maxin;
}

function getTaxID1() {
	return taxID1;
}

function getTaxID2() {
	return taxID2;
}

function prepareMethod1() {
	visualizerURLSuffix = pdb + "/" + chain + "/" + alignmentName + "/" + taxID1 + "/" + taxID2;
	if (typeof minIndex !== 'undefined' && typeof maxIndex !== 'undefined') {
		visualizerURLSuffix += "/" +minIndex+ "/" +maxIndex
	}
	visualizerURL = alignmentName + "/" + taxID1 + "/" + taxID2 + "/" + pdb
	mainForm.action = "twincons/" + visualizerURLSuffix;
	submitButton.disabled = !(alignmentName && taxID1 && taxID2 && chain);
	visualizer.disabled = submitButton.disabled;
}

function prepareMethod2() {
	visualizerURLSuffix = pdb + "/" + chain;
	visualizerURL = alignmentName + "/" + taxID1 + "/" + taxID2 + "/" + pdb
	mainForm.action = "upload/twincons/" + visualizerURLSuffix;
	submitButton.disabled = !(pdb && chain && pdb && (textbox.disabled || fileSelector.disabled));
	visualizer.disabled = submitButton.disabled;
}

//testing here (still slight problems when preparing with ranges)
function prepareMethod(upload = new Boolean(false)) {
	
	if (upload){
		visualizerURLSuffix = pdb + "/" + chain;
		if (typeof minIndex !== 'undefined' && typeof maxIndex !== 'undefined') {
			visualizerURLSuffix += "/" +minIndex+ "/" +maxIndex
		}
		mainForm.action = "upload/twincons/" + visualizerURLSuffix;
		submitButton.disabled = !(pdb && chain && pdb && (textbox.disabled || fileSelector.disabled));
	}else{
		visualizerURLSuffix = alignmentName + "/" + taxID1 + "/" + taxID2 + "/" + pdb;
		if (typeof minIndex !== 'undefined' && typeof maxIndex !== 'undefined') {
			visualizerURLSuffix += "/" +minIndex+ "/" +maxIndex
		}
		mainForm.action = "twincons/" + visualizerURLSuffix;
		submitButton.disabled = !(alignmentName && taxID1 && taxID2);
		
	}
	visualizer.disabled = submitButton.disabled;
}