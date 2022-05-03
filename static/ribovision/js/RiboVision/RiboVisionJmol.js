/* RiboVision 1.2 script library RiboVision.js 6:42 PM 04/25/2015 Chad R. Bernier


based on:
 *
 * Copyright (C) 2012,2013  RiboEvo, Georgia Institute of Technology, apollo.chemistry.gatech.edu
 *
 * Contact: Bernier.C.R@gatech.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *  02111-1307  USA.
 */

// for documentation see apollo.chemistry.gatech.edu/RiboVision/Documentation
//This doesn't exist and this probably won't be the final license.

// Initialize Jmol
//jmolInitialize("./jmol");
//Jmol.debugCode = true
/*
var myJmol1;
var myInfo1 = {
        height: '100%',
        width: '100%',
        jarFile: "JmolApplet.jar",
        jarPath: '..',
        j2sPath: "j2s",
        use: 'HTML5',
	console: "myJmol1_infodiv",
        debug: false
};*/

var JmolInfo = {
	addSelectionOptions: false,
	color: "#FFFFFF",
	debug: false,
	defaultModel: "",
	height: "100%",
	j2sPath: "./jmol/j2s", 
	isSigned: true,
	jarFile: "JmolAppletSigned0.jar",
	jarPath: "./jmol/java",
	//j2sPath: "jmol/jsmol/j2s",
	memoryLimit: 1024,
	readyFunction: null,
	script: null,
	serverURL: "./jmol/jsmol.php",
	src: null,
	use: "Java HTML5",
	width: "100%"
};	

function init3D(){
	Jmol.setDocument(0);
	Jmol.getApplet("myJmol", JmolInfo);
	$("#the3DpanelDiv").html(Jmol.getAppletHtml("myJmol", JmolInfo));

}
function waitFor3Dinit(dataStructure){
    if(typeof myJmol !== "undefined"){
        //variable exists, do what you want
		if (dataStructure & dataStructure.StructureName != undefined){
			load3Dstructure(dataStructure.JmolState);
		}
    }
    else{
        setTimeout(function(){
            waitFor3Dinit(dataStructure);
        },250);
    }
}
function waitFor3Dload(){
    if(typeof myJmol !== "undefined"){
        //variable exists, do what you want
		init3dStructures();
    }
    else{
        setTimeout(function(){
            waitFor3Dload();
        },250);
    }
}
function init3dStructures() {
	//var jscript = "display " + rvDataSets[speciesIndex].SpeciesEntry.Jmol_Model_Num_rRNA + ".1";
	//Jmol.script(myJmol, jscript);
	updateModel();
	
	// if (rvDataSets[1]) {
		// var rna_chains = ":" + rvDataSets[0].SpeciesEntry.New_RNA_Chains.replace(/:/,' or :') + " or :" + rvDataSets[1].SpeciesEntry.New_RNA_Chains.replace(/:/,' or :');
	// } else {
		// var rna_chains = ":" + rvDataSets[0].SpeciesEntry.New_RNA_Chains.replace(/:/,' or :');
	// }
	
	// var sele="rna and (" + rna_chains + ")";
	// Struct.setSelection(sele);
	// Struct.centerView(true);
	
}
function load3Dstructure(JmolState){
	if(myJmol!=null){
		Jmol.script(myJmol, "script states/" + JmolState + ".spt");
	}
}

function update3DProteinsLow(newcolor){
	var JscriptP = "set hideNotSelected false;";
	for (var i = 0; i < newcolor.length; i++) {
		JscriptP += "select (" + (rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rProtein) + ".1 and :" + 
		rvDataSets[0].SpeciesEntry.RNA_Chains_rProtein[1][rvDataSets[0].SpeciesEntry.RNA_Chains_rProtein[2].indexOf(seleProt[i])] +
		"); color Cartoon opaque [" + newcolor[i].replace("#", "x") + "];";
	}
	//JscriptP+="display " + (rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rRNA ) + ".1, " + (rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rProtein ) + ".1;";
	//jmolScript(JscriptP);
	Jmol.script(myJmol, JscriptP);
}

function colorMappingLoop3DLow(changeProteins){
	var Jscript = "display (selected), (";
	var JscriptP = "set hideNotSelected false;";
	//Jscript += rvds.SpeciesEntry.Jmol_Model_Num_rProtein + ".1 and (";
	//Jscript += " or ";
	
	for (var i = 0; i < changeProteins.length; i++) {
		Jscript += ":" + changeProteins[i].foundProt + " or ";
		JscriptP += "select (" +  ":" + changeProteins[i].foundProt + 
			"); color Cartoon opaque [" + changeProteins[i].newcolor.replace("#", "x") + "];spacefill off;";
	}
	Jscript=Jscript.slice(0,-4);
	Jscript += ")";
		
	Jmol.script(myJmol, Jscript);
	Jmol.script(myJmol, JscriptP);
				
}

function update3Dcolors() {
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	var script = "set hideNotSelected false;";
	var r0,
	r1,
	curr_chain,
	curr_color,
	compare_color,
	n,
	m;
	$.each(rvDataSets, function (index, rvds) {

	
		if (rvds.Residues[0] == undefined){return};
		//r0=rvds.Residues[0].resNum.replace(/[^:]*:/g,"");
		r0 = rvds.Residues[0].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
		curr_chain = rvds.Residues[0].ChainID;
		var targetLayer=rvds.getLinkedLayer();
		//rvds.Residues[0].CurrentData=targetLayer.Data[0];

		curr_color = colorNameToHex(targetLayer.dataLayerColors[0]);
		
		if (!curr_color || curr_color === '#000000') {
			curr_color = '#858585';
		}
		for (var i = 1; i < rvds.Residues.length; i++) {
			var residue = rvds.Residues[i];
			var residueLast = rvds.Residues[i - 1];
			var residueLastColor = targetLayer.dataLayerColors[i - 1];
			//rvds.Residues[i].CurrentData=targetLayer.Data[i];
			
			if (!residueLastColor){
				residueLastColor = '#858585';
			}
			if (residue.ChainID != "") {
				if (curr_chain == "") {
					curr_chain = residue.ChainID;
					curr_color = colorNameToHex(targetLayer.dataLayerColors[i]);
					if (!curr_color || curr_color === '#000000') {
						curr_color = '#858585';
					}
					r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, ""); ;
				} else if (residue.ChainID == null) {
					curr_chain = residue.ChainID;
					curr_color = colorNameToHex(targetLayer.dataLayerColors[i]);
					if (!curr_color || curr_color === '#000000') {
						curr_color = '#858585';
					}
					r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, ""); ;
				} else {
					if (!targetLayer.dataLayerColors[i]){
						compare_color = '#858585';
					} else {
						compare_color = colorNameToHex(targetLayer.dataLayerColors[i]);
					}
					if (((compare_color != colorNameToHex(residueLastColor)) || (curr_chain != residue.ChainID)) || (i == (rvds.Residues.length - 1))) {
						r1 = residueLast.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, ""); ;
						n = r1.match(/[A-z]/g);
						if (n != undefined) {
							r1 = r1.replace(n, "^" + n);
						}
						//if (colorNameToHex(residueLastColor).indexOf("#") == -1) {
							//script += "select " + (SubunitNames.indexOf(rvds.SpeciesEntry.Subunit) + 1) + ".1 and :" + curr_chain + " and (" + r0 + " - " + r1 + "); color Cartoon opaque [x" + curr_color + "]; ";
							//script += "select " + rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and " + r0 + " - " + r1 + ":" + curr_chain + "; color Cartoon opaque [x" + curr_color + "]; ";
							//script += "select " + rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and " + r0 + " - " + r1 + ":" + curr_chain + "; color opaque [x" + curr_color + "]; ";

						//} else {
							//script += "select " + rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and " + r0 + " - " + r1 + ":" + curr_chain + "; color Cartoon opaque [" + curr_color.replace("#", "x") + "]; ";
							//script += "select " + rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and " + r0 + " - " + r1 + ":" + curr_chain + "; color opaque [" + curr_color.replace("#", "x") + "]; ";
							script += "select " + r0 + " - " + r1 + ":" + curr_chain + "; color Cartoon opaque [" + curr_color.replace("#", "x") + "]; ";
							script += "select " + r0 + " - " + r1 + ":" + curr_chain + "; color opaque [" + curr_color.replace("#", "x") + "]; ";
						//}
						r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, ""); ;
						m = r0.match(/[A-z]/g);
						if (m != undefined) {
							r0 = r0.replace(m, "^" + m);
						}
						if (residue.ChainID != "") {
							curr_chain = residue.ChainID;
						}
						curr_color = colorNameToHex(targetLayer.dataLayerColors[i]);
						if (!curr_color || curr_color === '#000000') {
							curr_color = '#858585';
						}
					}
				}
			}
		}
		//if (colorNameToHex(residueLastColor).indexOf("#") == -1) {
			//script += "select " + (rvds.SpeciesEntry.Jmol_Model_Num_rRNA) + ".1 and "  + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color Cartoon opaque [x" + curr_color + "]; ";
			//script += "select " + (rvds.SpeciesEntry.Jmol_Model_Num_rRNA) + ".1 and "  + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color opaque [x" + curr_color + "]; ";
		//} else {
			//script += "select " + (rvds.SpeciesEntry.Jmol_Model_Num_rRNA) + ".1 and "  + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color Cartoon opaque [" + curr_color.replace("#", "x") + "]; ";
			//script += "select " + (rvds.SpeciesEntry.Jmol_Model_Num_rRNA) + ".1 and "  + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color opaque [" + curr_color.replace("#", "x") + "]; ";
			script += "select " + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color Cartoon opaque [" + curr_color.replace("#", "x") + "]; ";
			script += "select " + r0 + " - " + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + ":" + curr_chain + "; color opaque [" + curr_color.replace("#", "x") + "]; ";
		//}
		//updateSelectionDiv();
		//jmolScript(script);
	})
	Jmol.script(myJmol, script);
}

//////////////////////////////// Jmol Functions ////////////////////////////////
function updateModel() {
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	var n;
	var script;
	$.each(rvDataSets, function (index, rvds) {
		//Come back and support multiple selections?
		var targetSelection = rvds.getSelection($('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name'));
		if (targetSelection.Residues.length > 0){
			if (typeof script == 'undefined'){
				//script='set hideNotSelected true;select (';
			}
			//script += rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and (";
			for (var i = 0; i < targetSelection.Residues.length; i++) {
				if (targetSelection.Residues[i].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "") != null) {
					if (i > 0 || (index > 0 && i > 0)) {
						script += " or ";
					};
					n = targetSelection.Residues[i].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "").match(/[A-z]/g);
					if (n != null) {
						r1 = targetSelection.Residues[i].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "").replace(n, "^" + n);
					} else {
						r1 = targetSelection.Residues[i].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
					}
					script += r1 + ":" + targetSelection.Residues[i].ChainID;
				}
			}
			if (index !== rvDataSets.length - 1){
				script += ") or ";
			}
		} else {
			refreshModel();
		}
	});
	if (typeof script != 'undefined'){
		script += "));center selected;";
		Jmol.script(myJmol, script);
	}
}

function refreshModel() {
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	var script= "set hideNotSelected true;select (";
	$.each(rvDataSets, function (index, rvds) {
		if (index > 0 ){
			script +=" or ";
		}
		//script += rvds.SpeciesEntry.Jmol_Model_Num_rRNA + ".1 and (";
		script += "(";
		if (rvds.SpeciesEntry.RNA_Chains){
			var psplit=rvds.SpeciesEntry.RNA_Chains.split(";");
			for (var ii = 0; ii < psplit.length; ii++) {
				script += ":" + psplit[ii];
				if (ii < (psplit.length - 1)) {
					script += " or ";
				}
			}
			
		}
		script += ")";
	});
	script += "); center selected;";
	Jmol.script(myJmol, script);
}

function resetColorState() {
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	clearColor(false);
	Jmol.script(myJmol, "script states/" + rvDataSets[0].SpeciesEntry.Jmol_Script);
	var jscript = "display " + rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rRNA + ".1";
	
	Jmol.script(myJmol, jscript);
	//commandSelect();
	updateModel();
}

function save3dImgJmol() {
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	var jmlImgB64 = Jmol.getPropertyAsString(myJmol,'image');
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "saveJmolImg.php");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("name", "content");
	hiddenField.setAttribute("value", jmlImgB64);
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
}

function save3dImg() {
	AgreeFunction = save3dImgJmol();
	checkSavePrivacyStatus();
}

function ColorProteins3D(ColorProteins){
	if($('input[name="jp"][value=off]').is(':checked')){
		return;
	}
	if (rvDataSets[0].Residues[0] == undefined){return};
	
	var script = "set hideNotSelected false;";
	$.each(ColorProteins, function (index,value){
		var ressplit = value.ResNum.split("_");
		if (ressplit[0] !== "undefined"){
			if (colorNameToHex(value.Color).indexOf("#") == -1) {
				script += "select " + (rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rProtein) + ".1 and :" + ressplit[0] + " and " + ressplit[1].replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + "; color Cartoon opaque [x" + value.Color + "]; ";
			} else {
				script += "select " + (rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rProtein) + ".1 and :" + ressplit[0] + " and " + ressplit[1].replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + "; color Cartoon opaque [" + value.Color.replace("#", "x") + "]; ";
			}
		}
	});
	
	Jmol.script(myJmol, script);
}
///////////////////////////////////////////////////////////////////////////////