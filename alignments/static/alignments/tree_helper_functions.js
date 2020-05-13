function setTaxID1(tax1) {
	taxID1 = tax1;
}

function setTaxID2(tax2) {
	taxID2 = tax2;
}

function getTaxID1() {
	return taxID1;
}

function getTaxID2() {
	return taxID2;
}

function prepareMethod1() {
	visualizerURLSuffix = alignmentName + "/" + taxID1 + "/" + taxID2 + "/" + pdb;
	mainForm.action = "twincons/" + visualizerURLSuffix;
	submitButton.disabled = !(alignmentName && taxID1 && taxID2);
	visualizer.disabled = submitButton.disabled;
}

function prepareMethod2() {
	visualizerURLSuffix = pdb + "/" + chain;
	mainForm.action = "twincons/" + visualizerURLSuffix;
	submitButton.disabled = !(pdb && chain && pdb && (textbox.disabled || fileSelector.disabled));
	visualizer.disabled = submitButton.disabled;
}