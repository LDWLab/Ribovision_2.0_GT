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



///////////////////////// New Experimental Section  ???? //////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////// Label Functions ///////////////////////////////////
function initLabels(speciesSplit,customResidues) {
	$.each(speciesSplit, function (speciesIndex, species) {
		if (rvDataSets[speciesIndex] != undefined){
			rvDataSets[speciesIndex].addLabels([], []);
			if (species != "None" && species != "custom") {
				var TextLabels=[];
				var LineLabels=[];
				// $.postJSON('api/RiboVision/v1.0/textLabels', {
					// TextLabels : species
				// }, function (textlabels) {
					// TextLabels = textlabels;
					// rvDataSets[speciesIndex].clearCanvas("labels");
					// rvDataSets[speciesIndex].addLabels(TextLabels, LineLabels);
					// rvDataSets[speciesIndex].drawLabels("labels");
				// })
				// $.getJSON('api/RiboVision/v1.0/lineLabels', {
					// LineLabels : rvDataSets[speciesIndex].SpeciesEntry.LineLabels
				// }, function (linelabels) {
					// LineLabels = linelabels;
					// rvDataSets[speciesIndex].clearCanvas("labels");
					// rvDataSets[speciesIndex].addLabels(TextLabels, LineLabels);
					// rvDataSets[speciesIndex].drawLabels("labels");
				// })
				
				$.ajax({
					url: 'api/RiboVision/v1.0/textLabels',
					type: 'POST',
					contentType: 'application/json', 
					accept: 'application/json',
					data: JSON.stringify(species),
					cache: false,
					success: function(textlabels) {
						TextLabels = textlabels;
						rvDataSets[speciesIndex].clearCanvas("labels");
						rvDataSets[speciesIndex].addLabels(TextLabels, LineLabels);
						rvDataSets[speciesIndex].drawLabels("labels");
					}
				})
				$.ajax({
					url: 'api/RiboVision/v1.0/lineLabels',
					type: 'POST',
					dataType: 'json',
					data: JSON.stringify(species),
					cache: false,
					success: function(linelabels) {
						LineLabels = linelabels;
						rvDataSets[speciesIndex].clearCanvas("labels");
						rvDataSets[speciesIndex].addLabels(TextLabels, LineLabels);
						rvDataSets[speciesIndex].drawLabels("labels");
					}
				});;
				
			} else {
				rvDataSets[speciesIndex].clearCanvas("labels");
				if (customResidues){
					var customLabels=processCustomLabels(customResidues);
					rvDataSets[speciesIndex].addLabels(customLabels.TextLabels, customLabels.LineLabels);
					rvDataSets[speciesIndex].drawLabels("labels");
					console.log("Text:", rvDataSets[speciesIndex].rvTextLabels, "Line:", rvDataSets[speciesIndex].rvLineLabels)
				} else {
					console.log("Cleaning unnecessarily")
					rvDataSets[speciesIndex].addLabels([], []);
				}
			}
		}
	});
}

function processCustomLabels(customResidues){
	//alert("custom labels");
	var customLabels=[];
	customLabels["TextLabels"]=[];
	customLabels["LineLabels"]=[];

	
	var labeledResidues = $.grep(customResidues, function(e){ return e.LabelX != ''; });
			
	$.each(labeledResidues, function (index, data) {
		customLabels.TextLabels[index]={};
		customLabels.TextLabels[index].Fill=data.LabelColor
		customLabels.TextLabels[index].Font="MyriadPro-Regular";
		customLabels.TextLabels[index].FontSize=parseFloat(data.LabelFontSize);
		customLabels.TextLabels[index].LabelText=data.LabelSymbol
		customLabels.TextLabels[index].X=parseFloat(data.LabelX);
		customLabels.TextLabels[index].Y=parseFloat(data.LabelY);
		customLabels.TextLabels[index].id=index;
		customLabels.LineLabels[index]={};
		customLabels.LineLabels[index].id=index;
		customLabels.LineLabels[index].Fill="";
		customLabels.LineLabels[index].Stroke=data.LineColor;
		customLabels.LineLabels[index].StrokeLineJoin="round";
		customLabels.LineLabels[index].StrokeMiterLimit=parseFloat("10.000");
		customLabels.LineLabels[index].StrokeWidth=parseFloat(data.LineThickness);
		customLabels.LineLabels[index].X1=parseFloat(data.LineX1);
		customLabels.LineLabels[index].X2=parseFloat(data.LineX2);
		customLabels.LineLabels[index].Y1=parseFloat(data.LineY1);
		customLabels.LineLabels[index].Y2=parseFloat(data.LineY1);	
	});	
	return customLabels;
}
//function drawLabels(){}
///////////////////////////////////////////////////////////////////////////////


////////////////////////// Window Functions ///////////////////////////////////
function pan(dx,dy){
	rvViews[0].panXY(dx,dy);
}

function resizeElements(noDraw) {
	var MajorBorderSize = 1;
	var width = $(window).width();
	var height = $(window).height();
	var MainMenuFrac = 0.65;
	//ToolBar
	$(".toolBarBtn").css("height",0.04*height);
	$(".toolBarBtn").css("width",0.04*height);
	$("#toolBar").css("width",1.1*$(".toolBarBtn").first().outerHeight() + 4);
	var toolBarWidth = $("#toolBar").outerWidth();
	$("#toolBar").css("left",width-toolBarWidth);
	//Menu
	$("#menu").css('height', 0.93 * height);
	$("#menu").css('width', 0.16 * width);
	var xcorr = $("#menu").outerWidth();
	
	//Top Menu Section
	$("#topMenu").css('left', xcorr );
	$("#topMenu").css('top', 0);
	if($('input[name="nl"][value=on]').is(':checked')){
		$("#topMenu").show();
		$("#topMenu").outerHeight(TopDivide * height);
		$("#topMenu").css('width', width - xcorr - toolBarWidth);
		//var ycorr = $("#topMenu").outerHeight();
	} else {
		$("#topMenu").hide();
		$("#topMenu").outerHeight(0);
		$("#topMenu").css('width', 0);
	}
	var ycorr = $("#topMenu").outerHeight();
	if($('input[name="3dp"][value=on]').is(':checked')){
		var lp = (width - xcorr - toolBarWidth) * PanelDivide;
	} else {
		var lp = (width - xcorr - toolBarWidth);
	}
	var rp = (width - xcorr - toolBarWidth) - lp;
	var s = (height - ycorr);
	
	//SiteInfo
	$("#SiteInfo").css('width', xcorr);
	
	//ExportData
	$("#ExportData").css('width', xcorr);

	//MainMenu
	$("#MainMenu").css('width', MainMenuFrac * xcorr); //65%, make room for new MiniLayer
	$("#MainMenu").css('height', (0.93 * height) - parseFloat($("#SiteInfo").css('height')) - parseFloat($("#ExportData").css('height')));
	$("#MainMenu").css('top', parseFloat($("#SiteInfo").css('height')));
	
	//LinkSection
	$("#LinkSection").css('width', (1 - MainMenuFrac) * xcorr);
	$("#LinkSection").css('height',parseFloat($("#LinkSection .ui-widget-header3").css('height')) + parseFloat($(".miniLayerName").first().css('height')) + 13);
	$("#LinkSection").css('top', (0.93 * height) - parseFloat($("#ExportData").css('height')) - parseFloat($("#LinkSection").css('height')));
	
	//MiniLayer
	$("#MiniLayer").css('width', (1 - MainMenuFrac) * xcorr); //35%, for new MiniLayer
	$("#MiniLayer").css('height', (0.93 * height) - parseFloat($("#SiteInfo").css('height')) - parseFloat($("#ExportData").css('height')) - parseFloat($("#LinkSection").css('height')));
	$("#MiniLayer").css('top', parseFloat($("#SiteInfo").css('height')));
	
	//SideBarAccordian
	$("#SideBarAccordian").accordion("refresh");
	
	//Canvas Section
	$("#canvasDiv").css('height', s);
	$("#canvasDiv").css('width', lp);
	$("#canvasDiv").css('left', xcorr - 1);
	$("#canvasDiv").css('top', ycorr - 1);
	
	$("#canvasDiv canvas").attr({
		width : lp - 2 * MajorBorderSize,
		height : s - 2 * MajorBorderSize
	});
	$("#canvasDiv canvas").css('top', 0);
	$("#canvasDiv canvas").css('left', 0);
	
	//Navigator Section
	$("#navigator").css('left', xcorr);
	$("#navigator").css('top', ycorr + parseFloat($("#canvaslabel").css("height")));
	
	//3D Section
	$("#the3DpanelDiv").css('height', s);
	$("#the3DpanelDiv").css('width', rp + 1);
	$("#the3DpanelDiv").css('left', xcorr + parseFloat($("#canvasDiv").css('width')) - 1);
	$("#the3DpanelDiv").css('top', ycorr - 1);
	$("#the3DpanelDiv canvas").attr({
		width : rp + 1,
		height : s
	});
	resize3D()
	//Jmol.resizeApplet(myJmol,[(rp - 2 * MajorBorderSize),(s - 2 * MajorBorderSize)]);
	//$("#myJmol_object").css('height', );
	//$("#myJmol_object").css('width', );
	//$("#myJmol_object").css('top', 0);
	//$("#myJmol_object").css('left', 0);
	
	// Layer Panel
	$( "#LayerDialog" ).dialog( "option", "height", s - 2 * MajorBorderSize );	
	$( "#LayerDialog" ).dialog("widget").position({
		my: "right top",
		at: "right top",
		of: $( "#canvasDiv" )
	});
	// Color Panel
	$( "#ColorDialog" ).dialog( "option", "height", 0.75 * s - 2 * MajorBorderSize );	
	$( "#ColorDialog" ).dialog("widget").position({
		my: "right top",
		at: "right top",
		of: $( "#canvasDiv" )
	});
	
	//LogoDiv
	$("#LogoDiv").css('width', xcorr);
	$("#LogoDiv").css('height', 0.07 * height);
	$("#LogoDiv").css('top', $("#menu").css('height'));
	
	if (rvViews[0]){
		rvViews[0].width = rvDataSets[0].HighlightLayer.Canvas.width;
		rvViews[0].height = rvDataSets[0].HighlightLayer.Canvas.height;
		rvViews[0].clientWidth = rvDataSets[0].HighlightLayer.Canvas.clientWidth;
		rvViews[0].clientHeight = rvDataSets[0].HighlightLayer.Canvas.clientHeight;
		rvViews[0].windowHeight=height;
		rvViews[0].xcorr=xcorr;
		rvViews[0].ycorr=ycorr;
		if (noDraw!==true){
			$.each(rvDataSets, function (index, value) {
				value.refreshCanvases();
			});
		}
		drawNavLine();
	}
}

function resetView() {
	rvViews[0].centerZoom(rvViews[0].defaultScale);
	rvViews[0].panXY(rvViews[0].defaultX,rvViews[0].defaultY);
	$.each(rvDataSets, function (index, value) {
		value.refreshCanvases();
	});
	$("#slider").slider("value", 20);
}
///////////////////////////////////////////////////////////////////////////////


/////////////////////////// Selection Functions ///////////////////////////////
function dragHandle(event) {
	rvViews[0].pan(event);
}

function getSelectedLine(event){
	var nx = (event.clientX - rvViews[0].x - rvViews[0].xcorr + 2) / rvViews[0].scale; //subtract 250 for the menu width
	var ny = (event.clientY - rvViews[0].y - rvViews[0].ycorr + 2) / rvViews[0].scale; //subtract 80 for the info height
	//var zoomEnabled = $('input[name="za"][value=on]').is(':checked');
	if(ActiveBasePairSet != undefined){
		for (var i = 0; i < ActiveBasePairSet.length; i++) {
			var residue_i = MainResidueMap[ActiveBasePairSet[i].residue_i];
			var residue_j = MainResidueMap[ActiveBasePairSet[i].residue_j];
			
			var jdist = Math.sqrt(((nx - residue_i.X)*(nx - residue_i.X) + (ny - residue_i.Y)*(ny - residue_i.Y)));
			var kdist = Math.sqrt(((nx - residue_j.X)*(nx - residue_j.X) + (ny - residue_j.Y)*(ny - residue_j.Y)));
			var jkdist = Math.sqrt(((residue_i.X - residue_j.X)*(residue_i.X - residue_j.X) + (residue_i.Y - residue_j.Y)*(residue_i.Y - residue_j.Y)));

			if(zoomEnabled){
				var jkdist = Math.sqrt(((residue_i.X - residue_j.X)*(residue_i.X - residue_j.X) + (residue_i.Y - residue_j.Y)*(residue_i.Y - residue_j.Y)));
				if((150 - rvViews[0].scale*23) > jkdist){
					continue;
				}
				if(( (residue_i.X*rvViews[0].scale+rvViews[0].x < 0) || (residue_i.X*rvViews[0].scale+rvViews[0].x > rvViews[0].clientWidth) || (residue_i.Y*rvViews[0].scale+rvViews[0].y < 0) ||  (residue_i.Y*rvViews[0].scale+rvViews[0].y > rvViews[0].clientHeight))
					&& ( (residue_j.X*rvViews[0].scale+rvViews[0].x < 0) || (residue_j.X*rvViews[0].scale+rvViews[0].x > rvViews[0].clientWidth) || (residue_j.Y*rvViews[0].scale+rvViews[0].y < 0) ||  (residue_j.Y*rvViews[0].scale+rvViews[0].y > rvViews[0].clientHeight)) )  {
					continue;
				}
			}
			if( (jdist+kdist - jkdist) < .03){
				return i;
			}
		}
	}
	return -1;
}

function getSelected(event) {
	var nx = (event.clientX - rvViews[0].x - rvViews[0].xcorr + 2) / rvViews[0].scale; //subtract 250 for the menu width
	var ny = (event.clientY - rvViews[0].y - rvViews[0].ycorr + 2) / rvViews[0].scale; //subtract 80 for the info height
	
	if (ResiduePositions[0] != undefined) {
		for (var j = 0; j < rvDataSets.length; j++){
			if (ResiduePositions[j]){ //Need to check this, because this function could run before the second set is put in here.
				for (var i = 0; i < ResiduePositions[j].length; i++) {
					var res = ResiduePositions[j][i];
					if (((nx > res.X) ? nx - res.X : res.X - nx) + ((ny > res.Y) ? ny - res.Y : res.Y - ny) < 2) {
						return [i,j];
					}
				}
			}
		}
		return [-1,-1];
	} else {
		return [-1,-1];
	}
}

function expandSelection(command, SelectionName,SpeciesIndex) {
	if (!SelectionName){
		SelectionName = "Main";
	}
	var rvds=rvDataSets[SpeciesIndex];
	var targetSelection=rvds.getSelection(SelectionName);
	for (var i = 0; i < command.length; i++) {
		var com = command[i];
		if (!com){return;};
		var comsplit = com.split(":");
		if (comsplit.length > 1) {
			//: detected, so either a single residue or a range with molecule names
			var index = comsplit[1].indexOf("-");
			if (index != -1) {
				//Range detected
				
				var start = comsplit[0] + ":" + comsplit[1].substring(1, index);
				var end = comsplit[0] + ":" + comsplit[1].substring(index + 1, comsplit[1].length - 1);
				
				if (start && end) {
					var start_ind = MainResidueMap[start].index;
					var end_ind = MainResidueMap[end].index;
					
					for (var j = start_ind; j <= end_ind; j++) {
						var res = rvds.Residues[j];
						if (res != undefined){
							targetSelection.Residues.push(res);
						}
					}
				}
			} else {
				// Single residue detected
				var check_residue = MainResidueMap[ comsplit[0] + ":" + comsplit[1]];
				if (!check_residue){
					// Assume protein.
					var chainID = rvds.SpeciesEntry.RNA_Chains_rProtein[rvds.SpeciesEntry.Molecule_Names_rProtein.indexOf(comsplit[0])];
					var aloneRes = chainID + "_" + comsplit[1];
					// Skip the residue list check for now
					//var alone_ind = rvDataSets[0].ResidueList_rProtein.indexOf(aloneRes);
					//if (alone_ind >=0){
						//var targetSelection=rvDataSets[0].getSelection(SelectionName);
					targetSelection.Residues_rProtein.push(aloneRes);
					return aloneRes;
				} else if (check_residue && check_residue.rvds_index == SpeciesIndex){
					// RNA in this dataset
					targetSelection.Residues.push(rvds.Residues[check_residue.index]);
				} else {
					// do nothing
				}
			}
		} else if (comsplit[0] != "") {
			//It is either a basic range selection or a regular expression search. Go by first number/letter for now.
			if (comsplit[0].search(/\d/) > -1){
				var index = comsplit[0].indexOf("-");
				if (index != -1) {
					//assume first moleculename
					var molName = rvds.SpeciesEntry.Moleclue_Names[0];
					var start = molName + ":" + comsplit[0].substring(0, index);
					var end = molName + ":" + comsplit[0].substring(index + 1, comsplit[0].length);
					
					if (start && end) {
						var start_ind = MainResidueMap[start].index;
						var end_ind = MainResidueMap[end].index;
						for (var j = start_ind; j <= end_ind; j++) {
							var targetSelection=rvds.getSelection(SelectionName);
							targetSelection.Residues.push(rvds.Residues[j]);
						}
					}
				} else {
					var molName = rvds.SpeciesEntry.Moleclue_Names[0];
					var aloneRes = molName + ":" + comsplit[0];
					var alone_ind = MainResidueMap[aloneRes].index;
					if (alone_ind >=0){
						var targetSelection=rvds.getSelection(SelectionName);
						targetSelection.Residues.push(rvds.Residues[alone_ind]);
					}
				}
			} else {
				var re = new RegExp(comsplit[0],"g");
				while ((match = re.exec(rvds.SequenceList)) != null) {
					var start_ind = match.index;
					var end_ind = match.index + match[0].length - 1;
					for (var j = start_ind; j <= end_ind; j++) {
						var targetSelection=rvds.getSelection(SelectionName);
						targetSelection.Residues.push(rvds.Residues[j]);
					}
				}
			}
		}
	}
}

function commandSelect(command,SeleName) {
	if (!command) {
		var command = document.getElementById("commandline").value;
	}
	if (!SeleName) {
		var SeleName = $('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name');
	}
	command = command.split(";");
	
	$.each(rvDataSets, function(index,rvds){
		expandSelection(command,SeleName,index);
		updateSelectionDiv(SeleName,index);
		rvds.drawResidues("residues");
		rvds.drawSelection("selected");
		rvds.drawBasePairs("lines");
	});
	drawNavLine();
}

function selectResidue(event) {
	if (drag) {
		var curX = (event.clientX - rvViews[0].xcorr + 2 - rvViews[0].x) / rvViews[0].scale;
		var curY = (event.clientY - rvViews[0].ycorr + 2 - rvViews[0].y) / rvViews[0].scale;
		if (rvDataSets[0].Residues != undefined) {
			for (var SpeciesIndex = 0; SpeciesIndex < rvDataSets.length; SpeciesIndex++){
				var targetSelection = rvDataSets[SpeciesIndex].getSelection($('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name'));
				for (var i = 0; i < ResiduePositions[SpeciesIndex].length; i++) {
					if (rvViews[0].startX <= ResiduePositions[SpeciesIndex][i].X && ResiduePositions[SpeciesIndex][i].X <= curX && rvViews[0].startY <= ResiduePositions[SpeciesIndex][i].Y && ResiduePositions[SpeciesIndex][i].Y <= curY) {
						targetSelection.Residues.push(rvDataSets[SpeciesIndex].Residues[i]);
					}
					if (rvViews[0].startX >= ResiduePositions[SpeciesIndex][i].X && ResiduePositions[SpeciesIndex][i].X >= curX && rvViews[0].startY <= ResiduePositions[SpeciesIndex][i].Y && ResiduePositions[SpeciesIndex][i].Y <= curY) {
						targetSelection.Residues.push(rvDataSets[SpeciesIndex].Residues[i]);
					}
					if (rvViews[0].startX <= ResiduePositions[SpeciesIndex][i].X && ResiduePositions[SpeciesIndex][i].X <= curX && rvViews[0].startY >= ResiduePositions[SpeciesIndex][i].Y && ResiduePositions[SpeciesIndex][i].Y >= curY) {
						targetSelection.Residues.push(rvDataSets[SpeciesIndex].Residues[i]);
					}
					if (rvViews[0].startX >= ResiduePositions[SpeciesIndex][i].X && ResiduePositions[SpeciesIndex][i].X >= curX && rvViews[0].startY >= ResiduePositions[SpeciesIndex][i].Y && ResiduePositions[SpeciesIndex][i].Y >= curY) {
						targetSelection.Residues.push(rvDataSets[SpeciesIndex].Residues[i]);
					}
				}
				//Unselect code
				var sel = getSelected(event);
				if (sel[0] >= 0) {
					var res = rvDataSets[sel[1]].Residues[sel[0]];
					var result = $.grep(targetSelection.Residues, function(e){ return e.map_Index == res.map_Index; });
					//alert(result[0].resNum);
					
					if (result[0]) {
						targetSelection.Residues = $.grep(targetSelection.Residues, function(e){ return e.map_Index !== result[0].map_Index; });
					} else {
						targetSelection.Residues.push(res);
					}
				}
				rvDataSets[SpeciesIndex].drawResidues("residues");
				rvDataSets[SpeciesIndex].drawSelection("selected");
				rvDataSets[SpeciesIndex].drawBasePairs("lines");
				
				//drawLabels();
				updateSelectionDiv(targetSelection.Name,SpeciesIndex);
				drawNavLine();
				drag = false;
			}
		} else {
			drag = false;
		}
	}
	$("#canvasDiv").off("mouseup", selectResidue);			
	$("#canvasDiv").on("mousemove", mouseMoveFunction);

}

function updateSelectionDiv(SeleName,SpeciesIndex) {
	if (updateSelectionDiv.text == undefined){
		updateSelectionDiv.text=[''];
	};
	
	if (!SeleName){
		SeleName = "Main";
	}
	updateSelectionDiv.text[SpeciesIndex] = "";
	var targetSelection = rvDataSets[SpeciesIndex].getSelection(SeleName);
	for (var i = 0; i < targetSelection.Residues.length; i++) {
		res = targetSelection.Residues[i];
		//text = text + rvDataSets[0].SpeciesEntry.Molecule_Names[rvDataSets[0].SpeciesEntry.RNA_Chains.indexOf(res.ChainID)] + ":" + res.resNum.replace(/[^:]*:/g, "") + "( " + res.CurrentData + " ); ";
		//text[index] = text + rvds.SpeciesEntry.Molecule_Names[rvds.SpeciesEntry.RNA_Chains.indexOf(res.ChainID)] + ":" + res.resNum.replace(/[^:]*:/g, "") + "; ";
		if (res != undefined){
			updateSelectionDiv.text[SpeciesIndex] = updateSelectionDiv.text[SpeciesIndex] + res.resNum + "; ";
		}
	}
	//$("#selectDiv").html(text)
	$("[name=" + SeleName + "]").find(".selectionContent").find("[name=selectDiv]").text(updateSelectionDiv.text[SpeciesIndex]);
}

function clearSelection(AllFlag) {
	$.each(rvDataSets, function (SpeciesIndex,rvds){
		if (AllFlag){
			$.each(rvds.Selections, function(index,value){
				value.Residues = [];
				updateSelectionDiv(value.Name,SpeciesIndex);
			});
		} else {
			var targetSelection = rvds.getSelection($('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name'));
			targetSelection.Residues = []
			updateSelectionDiv(targetSelection.Name,SpeciesIndex);
			updateModel();
		}
		rvds.drawResidues("residues");
		rvds.drawSelection("selected");
		rvds.drawBasePairs("lines");
		drawNavLine();
	});
}

function selectionBox(event) {
	rvViews[0].startX = (event.clientX - rvViews[0].x - $("#menu").width()) / rvViews[0].scale;
	rvViews[0].startY = (event.clientY - rvViews[0].y - $("#topMenu").height()) / rvViews[0].scale;
	drag = true;
}

function clearLineSelection(event) {
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////// Color Functions /////////////////////////////////////
function colorResidue(event) {
	var sel = getSelected(event);
	if (sel[0] >=0) {
		var targetLayer=rvDataSets[sel[1]].getSelectedLayer();
		var color = colorNameToHex($("#MainColor").val());
		targetLayer.dataLayerColors[sel[0]]=color;	
		switch (targetLayer.Type){
			case "residues" : 
				var res = rvDataSets[sel[1]].Residues[sel[0]]
				res.color = color;
				rvDataSets[sel[1]].drawResidues("residues");
				break;
			case "circles" :
				rvDataSets[sel[1]].refreshResiduesExpanded(targetLayer.LayerName);
				break;
			case "contour" :
				rvDataSets[sel[1]].refreshResiduesExpanded(targetLayer.LayerName);
				break;
			default:
		}
		$.each(rvDataSets, function (index, rvds) {
			rvds.refreshCanvases();
		});
		//drawLabels();
		update3Dcolors();
	}
}

function colorLine(event) {
	var seleLine = getSelectedLine(event);
	if(seleLine >=0 ){
		var	targetLayer=rvDataSets[0].getLayerByType("lines");
		var color = colorNameToHex($("#LineColor").val());
		var value=targetLayer[0].Data[seleLine];
		
		var grd = targetLayer[0].CanvasContext.createLinearGradient(MainResidueMap[value.residue_i].X, MainResidueMap[value.residue_i].Y, MainResidueMap[value.residue_j].X, MainResidueMap[value.residue_j].Y);
		grd.addColorStop(0, "rgba(" + h2d(color.slice(1, 3)) + "," + h2d(color.slice(3, 5)) + "," + h2d(color.slice(5)) + "," + value.opacity + ")");
		grd.addColorStop(1, "rgba(" + h2d(color.slice(1, 3)) + "," + h2d(color.slice(3, 5)) + "," + h2d(color.slice(5)) + "," + value.opacity + ")");
								
		targetLayer[0].Data[seleLine]["color"] = grd;		
		ActiveBasePairSet[seleLine]["color"]=color;
		ActiveBasePairSet[seleLine]["color_hex"]=color;
		$.each(rvDataSets, function (index, rvds) {
			rvds.refreshCanvases();
		});
	}
}

function clearColor(update3D) {
	var	targetLayer=rvDataSets[0].getLayerByType("residues");
	if (arguments.length < 1) {
		var update3D = true;
	}
	for (var i = 0; i < rvDataSets[0].Residues.length; i++) {
		rvDataSets[0].Residues[i].color = "#000000";
		targetLayer[0].dataLayerColors[i]= undefined;
		//targetLayer[0].Data[i] = rvDataSets[0].Residues[i].map_Index;
		targetLayer[0].Data[i] = undefined;
	}
	rvDataSets[0].drawResidues("residues");
	//drawLabels();
	
	if (update3D) {
		update3Dcolors();
		updateModel();
	}
	
}

function colorSelection() {
	$.each(rvDataSets, function (index, rvds) {
		var targetLayer=rvds.getSelectedLayer();
		var color = colorNameToHex($("#MainColor").val());
		var targetSelection = rvds.getSelection($('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name'));
		for (var i = 0; i < targetSelection.Residues.length; i++) {
			targetLayer.dataLayerColors[targetSelection.Residues[i].map_Index - 1] = color;
		}
		rvds.drawDataCircles(targetLayer.LayerName,undefined,undefined,true);
		rvds.drawResidues(targetLayer.LayerName,undefined,undefined,true);
		rvds.refreshCanvases();
	});
	update3Dcolors();
}


function colorProcess(DataInput, indexMode,targetLayer,colors,SwitchPoint,SkipDraw) {
	var color_data = new Array();
	var color_data_IN = new Array();
	var data = new Array();
	var DataPoints = 0;
	if (DataInput.IncludeData){
		color_data_IN = DataInput.IncludeData;
		if (DataInput.ExtraData){
			color_data_IN = color_data_IN.concat(DataInput.ExtraData);
		}
		data = DataInput.IncludeData;
	} else {
		color_data_IN = DataInput;
		data = DataInput;
	}
	$.each(color_data_IN, function (index, value) {
		var f = parseFloat(value);
		if (!isNaN(f)){
			color_data.push(f);
		}
	});
	$.each(data, function (index, value) {
		var f = parseFloat(value);
		if (isNaN(f)){
			data[index]=undefined;
		} else {
			data[index]=f;
		}
	});
	
	var min = Math.min.apply(Math, color_data);
	var max = Math.max.apply(Math, color_data);
	var range = max - min;
	
	if (arguments.length < 3 || targetLayer==undefined) {
		var targetLayer = rvDataSets[0].getSelectedLayer();
	}
	
	if (!targetLayer) {
		$("#dialog-selection-warning p").text("Please select a valid layer and try again.");
		$("#dialog-selection-warning").dialog("open");
		return;
	}
	
	if (indexMode == "1") {
		var dataIndices = data;
	} else {
		var dataIndices = new Array;
		for (var i = 0; i < data.length; i++) {
			if (SwitchPoint !=undefined){
				dataIndices[i] = data[i] >= SwitchPoint ? 1 : 0;
			} else {
				dataIndices[i] = Math.round((data[i] - min) / range * (colors.length - 1));
			}
		}
	}
	
	if (!SkipDraw){
		targetLayer.Data = data;
		switch (targetLayer.Type) {
			case "circles":
				rvDataSets[targetLayer.SetNumber].drawDataCircles(targetLayer.LayerName, dataIndices, colors);
				break;
			case "residues":
				rvDataSets[targetLayer.SetNumber].drawResidues(targetLayer.LayerName, dataIndices, colors);
				break;
			case "contour":
				rvDataSets[targetLayer.SetNumber].drawContourLines(targetLayer.LayerName, dataIndices, colors);
				break;
			default:
				$( "#dialog-layer-type-error" ).dialog("open")
		}
		//update3Dcolors();
	}
	return dataIndices
}

function colorMappingLoop(targetLayer, seleProt, seleProtNames, OverRideColors) {
	if (arguments.length >= 4) {
		var colors2 = OverRideColors;
	} else {
		var colors2 = ColorLists.Rainbow1;
	}
	$("#ProtList").multiselect("widget").find(".ui-multiselect-checkboxes").find("span").css("color","black"); // Reset Protein colors to black.

	//Interaction Part
	var array_of_checked_values = $("#PrimaryInteractionList").multiselect("getChecked").map(function(){
	   return this.value;	
	}).get();
	var interactionchoice = array_of_checked_values[0];
	//var interactionchoice = $('#PrimaryInteractionList').val();
	var p = interactionchoice.indexOf("_NPN");
	if ( p >=0 ){
		ActiveBasePairSet=[];
		$.each(rvDataSets, function (index, rvds){
			//rvds.BasePairs = [];
			rvds.clearCanvas("lines");
		});
	}
	
	
	//Canvas Part
	colorMappingLoopCanvas(targetLayer,seleProt,seleProtNames,colors2);
	
	//Jmol Part
	colorMappingLoop3D(seleProt,colors2);
	
	//Menu Part
	for (var i = 0; i < seleProt.length; i++) {
		if (seleProt.length > 1) {
			var val = Math.round((i / (seleProt.length - 1)) * (colors2.length - 1));
		} else {
			var val = 0;
		}	

		var newcolor = (val < 0 || val >= colors2.length) ? "#000000" : colors2[val];
		var ProtName = $.grep($("#ProtList").multiselect("getChecked"), function(e) {
			return e.value == seleProt[i];
		});
		$(ProtName).next().css("color",newcolor);
		if (p > 0) {
			appendBasePairs(interactionchoice, seleProt[i]);
		}
	}	
	drawNavLine();
}

function colorMappingLoopCanvas(targetLayer,seleProt,seleProtNames,colors2){
	//Rest Canvases Part
	if (targetLayer){
		var setNumber=targetLayer.SetNumber;
		//var targetLayer = rvDataSets[0].getSelectedLayer();
		if (targetLayer.Type === "circles"){
			targetLayer.clearCanvas();
			targetLayer.clearData();
		}
		if (targetLayer.Type === "contour"){
			targetLayer.clearCanvas();
			targetLayer.clearData();
		}
		if (targetLayer.Type === "residues"){
			//targetLayer.clearCanvas();
			clearColor(false);
			targetLayer.clearData();
		}
		targetLayer.Data = new Array;
		for (var j = 0; j < rvDataSets[setNumber].Residues.length; j++) {
			targetLayer.Data[j] = " ";
		}
	}
	
	for (var i = 0; i < seleProt.length; i++) {
		if (seleProt.length > 1) {
			var val = Math.round((i / (seleProt.length - 1)) * (colors2.length - 1));
		} else {
			var val = 0;
		}	

		var newcolor = (val < 0 || val >= colors2.length) ? "#000000" : colors2[val];
		
		if (targetLayer){
			var dataIndices = new Array;
			for (var jj = 0; jj < rvDataSets[setNumber].Residues.length; jj++) {
				if (rvDataSets[setNumber].Residues[jj][seleProt[i]] && rvDataSets[setNumber].Residues[jj][seleProt[i]] >0){
					dataIndices[jj] = rvDataSets[setNumber].Residues[jj][seleProt[i]];
					targetLayer.Data[jj] = targetLayer.Data[jj] + seleProtNames[i] + " ";
				}
			}
			rvDataSets[setNumber].drawDataCircles(targetLayer.LayerName, dataIndices, ["#000000", newcolor], true);
			rvDataSets[setNumber].drawContourLines(targetLayer.LayerName, dataIndices, ["#000000", newcolor], true);
			rvDataSets[setNumber].drawResidues(targetLayer.LayerName, dataIndices, ["#000000", newcolor], true);
		}
	}
}

function colorMappingLoop3D(seleProt,colors2){
	var numProteinsList=[];
	var changeProteins=[];
	
	$.each(rvDataSets, function (index, rvds){
		numProteinsList[index]=0;

		var foundProt;
		for (var i = 0; i < seleProt.length; i++) {
			if (seleProt.length > 1) {
				var val = Math.round((i / (seleProt.length - 1)) * (colors2.length - 1));
			} else {
				var val = 0;
			}	

			var newcolor = (val < 0 || val >= colors2.length) ? "#000000" : colors2[val];	
			foundProt=rvds.SpeciesEntry.RNA_Chains_rProtein[rvds.SpeciesEntry.internal_protein_names.indexOf(seleProt[i])];
			if (foundProt){
				changeProteins.push({"foundProt" : foundProt,"newcolor" : newcolor});
			}
		}
		

	});
	/*Jscript += ");";
	if($('input[name="3dp"][value=on]').is(':checked')){
		//update3Dcolors();
		//refreshModel();
	}*/
	colorMappingLoop3DLow(changeProteins);
}


function update3DProteins(seleProt, OverRideColors) {
	if($('input[name="3dp"][value=off]').is(':checked')){
		return;
	}
	if (arguments.length >= 2) {
		var colors2 = OverRideColors;
	} else {
		var colors2 = ColorLists.Rainbow1;
	}
	var newcolor=[];
		
	for (var i = 0; i < seleProt.length; i++) {
		if (seleProt.length > 1) {
			var val = Math.round((i / (seleProt.length - 1)) * (colors2.length - 1));
		} else {
			var val = 0;
		}
		newcolor[i] = (val < 0 || val >= colors2.length) ? "#000000" : colors2[val];
	}
	
	update3DProteinsLow(newcolor);
}

function colorMapping(targetLayer, colName, colorlist = "Rainbow1" , indexMode = false, rePlaceData) {
	if (!targetLayer) {
		$("#dialog-selection-warning p").text("Please select a valid layer and try again.");
		$("#dialog-selection-warning").dialog("open");
		return;
	}
	var colors = ColorLists[colorlist];
	if(rvDataSets[targetLayer.SetNumber].Residues[0][colName] !=undefined){
		var data = new Array;
		for (var j = 0; j < rvDataSets[targetLayer.SetNumber].Residues.length; j++) {
			data[j] = rvDataSets[targetLayer.SetNumber].Residues[j][colName];
		}
		colorProcess(data, indexMode,targetLayer,colors);
	} else {
		$.ajax({
			url: 'api/RiboVision/v1.0/fetchStructData',
			type: 'POST',
			contentType: 'application/json', 
			accept: 'application/json',
			data: JSON.stringify([rvDataSets[targetLayer.SetNumber].SpeciesEntry.SS_Table, structureName[0].StructureName,colName]),
			cache: false,
			success: function(newdata) {
				var data = new Array;
				for (var j = 0; j < rvDataSets[targetLayer.SetNumber].Residues.length; j++) {
					rvDataSets[targetLayer.SetNumber].Residues[j][colName] = newdata[j].Value;
					data[j] = newdata[j].Value;
				}
				//alert("Database");
				colorProcess(data, indexMode,targetLayer,colors);
			}
		})
	}	
}

function colorNameToHex(color,prefix='#',nullcolor=false) {
	var colors = SupportedColors //Global Variable to save time.
	if (color) {
		var newcolorH = color.match(/#[\dABCDEFabcdef]{6,6}$/);
		if ((newcolorH  !=null) && newcolorH[0].length === 7){
			return newcolorH[0].replace('#',prefix);
		} else if (typeof colors[color.toLowerCase().replace(/\s+/g, '')] != 'undefined'){
			return prefix + colors[color.toLowerCase().replace(/\s+/g, '')];
		} else {
			console.log('Unrecognized color "' + color + '"');
			//return false;
			return prefix + '868686';
		}
	} else {
		return nullcolor;
	}
}
///////////////////////////////////////////////////////////////////////////////


/////////////////////// Residue Functions /////////////////////////////////////
//function refreshResiduesExpanded(targetLayer){}

//function drawResiduesExpanded(targetLayer,dataIndices,ColorArray){}

//function drawResidues(){}
///////////////////////////////////////////////////////////////////////////////


//////////////////////// Interaction Functions ////////////////////////////////
//function drawBasePairs(){}

function appendBasePairs(BasePairTable, colName) {
	var p = BasePairTable.indexOf("_NPN");
	if (p < 0) {
		$.getJSON('api/RiboVision/v1.0/basePairs', {
			BasePairs : BasePairTable
		}, function (basePairs2) {
			ActiveBasePairSet=ActiveBasePairSet.concat(basePairs2);
			//rvDataSets[0].BasePairs = rvDataSets[0].BasePairs.concat(basePairs2);
			FullBasePairSet = ActiveBasePairSet;
			//rvDataSets[0].drawBasePairs("lines");
		});
	} else {
		//var dd = document.getElementById("ProtList");
		//var colName = dd.options[dd.selectedIndex].value;
		$.getJSON('api/RiboVision/v1.0/basePairs', {
			ProtBasePairs : BasePairTable,
			ProtChain : colName
		}, function (basePairs2) {
			ActiveBasePairSet=ActiveBasePairSet.concat(basePairs2);
			//rvDataSets[0].BasePairs = rvDataSets[0].BasePairs.concat(basePairs2);
			FullBasePairSet = ActiveBasePairSet;
			//rvDataSets[0].drawBasePairs("lines");
		});
	}
}

function refreshBasePairs(BasePairTable) {
	FullBasePairSet=[];
			
	if (BasePairTable == "clear_lines") {
		$.each(rvDataSets, function(index, rvds){
			rvds.clearCanvas("lines");
		})
		
	} else if (BasePairTable == "proteins") {
		var array_of_checked_values = $("#ProtList").multiselect("getChecked").map(function () {
				return this.value;
			}).get();
		var array_of_checked_titles = $("#ProtList").multiselect("getChecked").map(function () {
				return this.title;
		}).get();
		var ims = document.getElementById("SecondaryInteractionList");
		ims.options.length = 0;
		ims.options[0] = new Option("NPN", "NPN");
		ims.options[0].setAttribute("selected", "selected");
		$("#SecondaryInteractionList").multiselect("refresh");
		colorMappingLoop(undefined,array_of_checked_values,array_of_checked_titles);
	} else {
		$.ajax({
			url: 'api/RiboVision/v1.0/fetchInteractions',
			type: 'POST',
			contentType: 'application/json', 
			accept: 'application/json',
			data: JSON.stringify([structureName[0].StructureName,BasePairTable]),
			cache: false,
			success: function(basePairs2) {
				$.each(basePairs2, function (ind, item) {
					item.lineWidth = 0.75;
					item.opacity = 0.5;
					item.color_hex = '#231F20';
					//item.residue_i=rvDataSets[index].SpeciesEntry.Molecule_Names[rvDataSets[index].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[index].Residues[item.resIndex1].ChainID)] + ":" + rvDataSets[index].Residues[item.resIndex1].resNum.replace(/[^:]*:/g, "");
					//item.residue_j=rvDataSets[index].SpeciesEntry.Molecule_Names[rvDataSets[index].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[index].Residues[item.resIndex2].ChainID)] + ":" + rvDataSets[index].Residues[item.resIndex2].resNum.replace(/[^:]*:/g, "");
				});
				FullBasePairSet=FullBasePairSet.concat(basePairs2);	
				ActiveBasePairSet = FullBasePairSet;
				$.each(rvDataSets, function(index, rvds){
					rvds.drawBasePairs("lines");
				})
				// Set interaction submenu to allow for subsets of these basepairs to be displayed. 
				// For now, let's set a global variable to store the whole table, so that it doesn't have to be refetched everytime a subset is chosen. 
				// This will get better when I revamp who BasePair interactions work
	
				//rvDataSets[index].FullBasePairSet = rvDataSets[index].BasePairs;
				var BP_Type = [];	
				$.each(FullBasePairSet, function (ind, item) {
					BP_Type.push(item.bp_type);
				});
				BP_TypeU = $.grep(BP_Type, function (v, k) {
					return $.inArray(v, BP_Type) === k;
				});
				
				var ims = document.getElementById("SecondaryInteractionList");
				ims.options.length = 0;
				$.each(BP_TypeU, function (ind, item) {
					ims.options[ind] = new Option(item, item);
					ims.options[ind].setAttribute("selected", "selected");
				});	
				$("#SecondaryInteractionList").multiselect("refresh");
				
			}
		})
	}
}

function filterBasePairs(IncludeTypes){
	ActiveBasePairSet=[];
	$.each(FullBasePairSet, function (index, value){
		if ($.inArray(value.bp_type,IncludeTypes) >= 0){
			ActiveBasePairSet.push(value);
		}
	});
	$.each(rvDataSets, function (index, rvds){
		rvds.drawBasePairs("lines");
	});
}
///////////////////////////////////////////////////////////////////////////////


//////////////////////////// Mouse Functions //////////////////////////////////
function mouseEventFunction(event) {
	var BaseViewMode = $('input[name="bv"][value=on]').is(':checked');
	$("#ResidueTip").tooltip("close");
	$("#InteractionTip").tooltip("close");
	if (event.handleObj.origType == "mousedown" && !BaseViewMode) {
		if (onebuttonmode == "select" || (event.which == 3 && event.altKey == false) || (event.which == 1 && event.shiftKey == true)) {
			$("#canvasDiv").off("mousemove", dragHandle);
			$("#canvasDiv").off("mousemove", mouseMoveFunction);
			selectionBox(event);
			$("#canvasDiv").on("mousemove", dragSelBox);
			$("#canvasDiv").on("mouseup", selectResidue);
		} else if (onebuttonmode == "selectL" || (event.which == 3 && event.altKey == true )) {
			$("#canvasDiv").off("mousemove", dragHandle);
			//$("#canvasDiv").off("mousemove", mouseMoveFunction);
		} else if (onebuttonmode == "color" || (event.which == 2 && event.altKey == false) || (event.which == 1 && event.ctrlKey == true && event.altKey == false)) {
			$("#canvasDiv").off("mousemove", dragHandle);
			colorResidue(event);
		} else if (onebuttonmode == "colorL" || (event.which == 2 && event.altKey == true) || (event.which == 1 && event.ctrlKey == true && event.altKey == true)) {
			$("#canvasDiv").off("mousemove", dragHandle);
			colorLine(event);
		} else {
			rvViews[0].lastX = event.clientX;
			rvViews[0].lastY = event.clientY;
			$("#canvasDiv").on("mousemove", dragHandle);
		}
	} else if (event.handleObj.origType == "mousedown" && BaseViewMode) {
		$("#canvasDiv").off("mousemove", dragHandle);
		BaseViewCenter(event);
	}
	if (event.handleObj.origType == "mouseup") {
		$("#canvasDiv").off("mousemove", dragHandle);
		$("#canvasDiv").off("mousemove", dragSelBox);
		$.each(rvDataSets, function (index, rvds) {
			rvds.HighlightLayer.clearCanvas();
		});
	}
}

function mouseMoveFunction(event){
	//mouseMoveFunction.count = mouseMoveFunction.count || 1 // mouseMoveFunction.count is undefined at first 
	//console.log(mouseMoveFunction.count++);
	var sel = getSelected(event);
	if (sel[1] >=0 && rvDataSets[sel[1]].Font_Size_Canvas){
		var Font_Size_Canvas=rvDataSets[sel[1]].Font_Size_Canvas;
	} else {
		var Font_Size_Canvas=3.1;
	}
	if (sel[1] >=0 && rvDataSets[sel[1]].Circle_Radius){
		var Circle_Radius=rvDataSets[sel[1]].Circle_Radius;
	} else {
		var Circle_Radius=1.7;
	}
	
	$.each(rvDataSets, function (index, rvds) {
		rvds.HighlightLayer.clearCanvas();
	});
	$("#ResidueTip").tooltip("close");
	$("#InteractionTip").tooltip("close");
	switch (onebuttonmode){
		case "selectL":
			return;
			break;
		case "colorL":
			return;
			break;
		case "select":
			return;
			break;
		case "move":
			if (event.altKey == true && event.ctrlKey == false){
				var seleLine = getSelectedLine(event);
				if(seleLine >=0 ){
					var j = ActiveBasePairSet[seleLine].residue_i;
					var k = ActiveBasePairSet[seleLine].residue_j;
					var residue_i = MainResidueMap[j];
					var residue_j = MainResidueMap[k];
					
					rvDataSets[0].HighlightLayer.CanvasContext.strokeStyle = "#6666ff";
					rvDataSets[0].HighlightLayer.CanvasContext.beginPath();
					rvDataSets[0].HighlightLayer.CanvasContext.moveTo(residue_i.X, residue_i.Y);
					rvDataSets[0].HighlightLayer.CanvasContext.lineTo(residue_j.X, residue_j.Y);
					rvDataSets[0].HighlightLayer.CanvasContext.closePath();
					rvDataSets[0].HighlightLayer.CanvasContext.stroke();
					if(isRTon){
						if(typeof movewaitL != 'undefined'){
							clearTimeout(movewaitL);
						}
						movewaitL = setTimeout(function(){
							createInfoWindow(seleLine,"lines");
							$("#InteractionTip").css("bottom",rvViews[0].windowHeight - event.clientY);
							$("#InteractionTip").css("left",event.clientX);
							$("#InteractionTip").tooltip("open");},300);
					}
				}	
			} else if (event.ctrlKey == true && event.altKey == false) {
				//var sel = getSelected(event);
				if (sel[0] >=0) {
					var targetLayer=rvDataSets[sel[1]].getSelectedLayer();
					switch (targetLayer.Type){
						case "residues" : 
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = "#000000";
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.font = Font_Size_Canvas + 'pt "Myriad Pro", Calibri, Arial';
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.textBaseline = "middle";
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.textAlign = "center";
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillStyle = colorNameToHex($("#MainColor").val());
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillText(rvDataSets[sel[1]].Residues[sel[0]].resName, ResiduePositions[sel[1]][sel[0]].X, ResiduePositions[sel[1]][sel[0]].Y);
							break;
						case "circles" :
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.arc(ResiduePositions[sel[1]][sel[0]].X, ResiduePositions[sel[1]][sel[0]].Y, (targetLayer.ScaleFactor * Circle_Radius), 0, 2 * Math.PI, false);
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = colorNameToHex($("#MainColor").val());
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
							if (targetLayer.Filled) {
								rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillStyle = colorNameToHex($("#MainColor").val());
								rvDataSets[sel[1]].HighlightLayer.CanvasContext.fill();
							}
							break;
						case "contour" :
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineJoin = "round";  
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.moveTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X1 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y1 - .3);
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X2 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y2 - .3);
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X3 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y3 - .3);
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = colorNameToHex($("#MainColor").val());
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.0 * 1.9 * Circle_Radius;					
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
							break;
						default :
					}
				}
			} else if ( event.ctrlKey == true && event.altKey == true) {
				var seleLine = getSelectedLine(event);
				if(seleLine >=0 ){
					var j = ActiveBasePairSet[seleLine].residue_i;
					var k = ActiveBasePairSet[seleLine].residue_i;
					var residue_i = MainResidueMap[j];
					var residue_j = MainResidueMap[k];
					rvDataSets[0].HighlightLayer.CanvasContext.strokeStyle = colorNameToHex($("#LineColor").val());
					rvDataSets[0].HighlightLayer.CanvasContext.beginPath();
					rvDataSets[0].HighlightLayer.CanvasContext.moveTo(residue_i.X, residue_i.Y);
					rvDataSets[0].HighlightLayer.CanvasContext.lineTo(residue_j.X, residue_j.Y);
					rvDataSets[0].HighlightLayer.CanvasContext.closePath();
					rvDataSets[0].HighlightLayer.CanvasContext.stroke();
					if(isRTon){
						if(typeof movewaitL != 'undefined'){
							clearTimeout(movewaitL);
						}
						movewaitL = setTimeout(function(){
							createInfoWindow(seleLine,"lines");
							$("#InteractionTip").css("bottom",rvViews[0].windowHeight - event.clientY);
							$("#InteractionTip").css("left",event.clientX);
							$("#InteractionTip").tooltip("open");},300);
					}
				}	
			} else {
				//var sel = getSelected(event);
				if (sel[0] >=0){
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.arc(ResiduePositions[sel[1]][sel[0]].X, ResiduePositions[sel[1]][sel[0]].Y, 1.176 * Circle_Radius, 0, 2 * Math.PI, false);
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = "#6666ff";
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineWidth=Circle_Radius/1.7;
					rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
					if(isRTon){
						if(typeof movewait != 'undefined'){
							clearTimeout(movewait);
						}
						movewait = setTimeout(function(){
							createInfoWindow(sel,"residue");
							$("#ResidueTip").css("bottom",rvViews[0].windowHeight - event.clientY);
							$("#ResidueTip").css("left",event.clientX);
							$("#ResidueTip").tooltip("open");},100);
					}
				}
			}
			break;
		case "color":
			//var sel = getSelected(event);
			if (sel[0] != -1) {
				var targetLayer=rvDataSets[sel[1]].getSelectedLayer();
				switch (targetLayer.Type){
					case "residues" : 
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = "#000000";
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.font = Font_Size_Canvas + 'pt "Myriad Pro", Calibri, Arial';
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.textBaseline = "middle";
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.textAlign = "center";
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillStyle = colorNameToHex($("#MainColor").val());
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillText(rvDataSets[sel[1]].Residues[sel[0]].resName, ResiduePositions[sel[1]][sel[0]].X, ResiduePositions[sel[1]][sel[0]].Y);
						break;
					case "circles" :
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.arc(ResiduePositions[sel[1]][sel[0]].X, rResiduePositions[sel[1]][sel[0]].Y, (targetLayer.ScaleFactor * Circle_Radius), 0, 2 * Math.PI, false);
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = colorNameToHex($("#MainColor").val());
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
						if (targetLayer.Filled) {
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.fillStyle = colorNameToHex($("#MainColor").val());
							rvDataSets[sel[1]].HighlightLayer.CanvasContext.fill();
						}
						break;
					case "contour" :
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineJoin = "round";  
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.moveTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X1 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y1 - .3);
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X2 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y2 - .3);
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineTo(rvDataSets[sel[1]].ContourLinePoints[sel[0]].X3 - .05, rvDataSets[sel[1]].ContourLinePoints[sel[0]].Y3 - .3);
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = colorNameToHex($("#MainColor").val());
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.0 * 1.9 * Circle_Radius;					
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
						rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
						break;
					default :
				}
			}
			break;
		default: 
	}
}

function dragSelBox(event){
	$.each(rvDataSets, function (index, rvds) {
		rvds.HighlightLayer.clearCanvas();
	});
	rvViews[0].drag(event);
}

function mouseWheelFunction(event,delta){
	rvViews[0].zoom(event, delta);
	
	var sel = getSelected(event);
	$.each(rvDataSets, function (index, rvds) {
		rvds.HighlightLayer.clearCanvas();
	});
	
	if (sel[0] == -1) {
		//document.getElementById("currentDiv").innerHTML = "<br/>";
	} else {
		rvDataSets[sel[1]].HighlightLayer.CanvasContext.beginPath();
		rvDataSets[sel[1]].HighlightLayer.CanvasContext.arc(ResiduePositions[sel[1]][sel[0]].X, ResiduePositions[sel[1]][sel[0]].Y, 2, 0, 2 * Math.PI, false);
		rvDataSets[sel[1]].HighlightLayer.CanvasContext.closePath();
		rvDataSets[sel[1]].HighlightLayer.CanvasContext.strokeStyle = "#6666ff";
		rvDataSets[sel[1]].HighlightLayer.CanvasContext.stroke();
	}
	if (drag) {
		rvViews[0].drag(event);
	}
	return false;
}

function colorLineSelection(event) {
}

///////For popup window////
function createInfoWindow(Sele,InfoMode){
	if (InfoMode === "residue"){
		addPopUpWindowResidue(Sele);
		$("#ResidueTip").tooltip("option","content",$("#residuetip").html());
	} else {
		addPopUpWindowLine(Sele);
		$("#InteractionTip").tooltip("option","content",$("#interactiontip").html());
	}
}
///////////////////////////////////////////////////////////////////////////////


// Misc Functions

function BaseViewCenter(event){
	if($('input[name="3dp"][value=off]').is(':checked')){
		return;
	}
	var sel = getSelected(event);
	if (sel[0] != -1) {
		var res = rvDataSets[sel[1]].Residues[sel[0]];
		var script = "center " + (rvDataSets[sel[1]].SpeciesEntry.Jmol_Model_Num_rRNA) + ".1 and " + res.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') +":" + res.ChainID;
		Jmol.script(myJmol, script);
	}
}
function modeSelect(mode) {
	onebuttonmode = mode;
}
function ProcessBubble(ui,targetLayerName){
	switch ($(ui).parent().attr("id")) {
		case "AlnBubbles" :
			$.each(rvDataSets, function(index, value) {
				var targetLayer = rvDataSets[index].getLayer(targetLayerName);
				colorMapping(targetLayer,ui.data("colName"),ui.data("OverRideColors"),ui.data("indexMode"),ui.data("rePlaceData"));
				drawNavLine();
			});
			break;
		case "StructDataBubbles" :
			$.each(rvDataSets, function(index, value) {
				var targetLayer = rvDataSets[index].getLayer(targetLayerName);
				updateStructData(ui,targetLayer);
			});
			break;
		case "ProteinBubbles" : 
			var array_of_checked_values = $("#ProtList").multiselect("getChecked").map(function () {
					return this.value;
				}).get();
			var array_of_checked_titles = $("#ProtList").multiselect("getChecked").map(function () {
				return this.title;
			}).get();
			
			$.each(rvDataSets, function(index, value) {
				var targetLayer = rvDataSets[index].getLayer(targetLayerName);
				colorMappingLoop(targetLayer,array_of_checked_values,array_of_checked_titles);
			});
			break;
		case "CustomDataBubbles" :
			//CustomDataProcess will use the first selection. I make it a rule that when dragging a custom data bubble,
			// a new selection is created for this purpose. 
			$.each(rvDataSets, function(index, value) {
				rvDataSets[index].addSelection();
			});
			customDataProcess($(ui[0]).attr("filename"),targetLayerName);;
			break;
		default :
			//debugger;
			//alert("other");
	}
	
}
function updateStructData(ui,targetLayer) {
	colorMapping(targetLayer,ui.data("ColName"),ui.data("OverRideColors"),ui.data("indexMode"),ui.data("rePlaceData"));
	drawNavLine();
}

function openRvState() {
	var PrivacyStatus = get_cookie("privacy_status_data");
	var PrivacyString = "This feature does not currently upload any data to our server. We don't have a privacy policy at this time"
		 + " because one isn't needed. We can not see these data you are about to graph. Click \"I agree\" to acknowledge acceptance of our policy.";
	 
	 var FileReaderFile = $("#files2")[0].files; // FileList object
	 
	 AgreeFunction = function () {
		for (var i = 0; i < FileReaderFile.length; i++) {
			reader = new FileReader();
			switch (FileReaderFile[i].name.split('.').pop().toLowerCase()){
				case "zip":
					//alert("zip");
					reader.readAsBinaryString(FileReaderFile[i]);
					reader.onload = function (){
						//var loader = new ZipLoader(reader.result);
						var loader = new ZipLoader(reader.result);
						//var data = loader.load('localfile.zip://Ribovision_State.rvs.txt');
						alert("Unfinished zip support. Please unzip your file.");
					}
					break;	
				default:
					reader.readAsText(FileReaderFile[i]);
					reader.onload = function () {
						var rvSaveState = JSON.parse(reader.result);
						rvDataSets[0]=rvDataSets[0].fromJSON(rvSaveState["RvDS"]);
						// Re stringify a few things for compatibility / symmetry with local storage
						rvSaveState["rvLayers"] = JSON.stringify(rvDataSets[0].Layers);
						rvSaveState["rvSelections"] = JSON.stringify(rvDataSets[0].Selections);
						rvSaveState["rvLastSpecies"] = rvDataSets[0].Name;
						if($('input[name="3dp"][value=on]').is(':checked')){
							Jmol.script(myJmol, "script states/" + rvDataSets[0].SpeciesEntry.Jmol_Script);
							var jscript = "display " + rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rRNA + ".1";
							Jmol.script(myJmol, jscript);
						}				
						processRvState(rvSaveState);
						updateModel();
						update3Dcolors();
						if($('input[name="3dp"][value=on]').is(':checked')){
							var a = rvSaveState.rvJmolOrientation.match(/reset[^\n]+/);
							Jmol.script(myJmol, a[0]);
						}
				}
					break;
			}
			
		}
	 };
	 
	 if (PrivacyStatus != "Agreed") {
		$("#Privacy-confirm").text(PrivacyString);
		CurrPrivacyCookie = "privacy_status_data";
		$("#Privacy-confirm").dialog('open');
	} else {
		AgreeFunction();
	}
}
//This is custom structure mode
function ImportStructureFileSelect(event) {
	
	var FileReaderFile = event.target.files; // FileList object
	for (var i = 0; i < FileReaderFile.length; i++) {
		var reader = new FileReader();
		reader.readAsText(FileReaderFile[i]);
		
		reader.onload = function () {
			var result = reader.result;
			var customStructure=$.csv.toObjects(result);
			var customkeys = Object.keys(customStructure[0]);
			//Input checks
			
			//Load Species
			var d=$('#speciesList').iosMenu().menu()[0];
			if( typeof d.species_array == 'undefined' ) {
				d.species_array = ['',''];
			}
			
			// New two structure mode. Custom mode put structure in first slot.
			//No two structure custom mode yet
			loadSpecies("custom",customStructure);
			
			
			
			
		}
	}
}

//This is custom data mode
function ImportDataFileSelect(event) {
	var PrivacyStatus = get_cookie("privacy_status_data");
	var PrivacyString = "This feature does not currently upload any data to our server. We don't have a privacy policy at this time"
		 + " because one isn't needed. We can not see these data you are about to graph. Click \"I agree\" to acknowledge acceptance of our policy.";
	var NewData;
	var command;
	var ColorList=[];
	var ColorListU;
	var DataList=[];
	var DataListU;
	var ColorGrad=[];
	
	var FileReaderFile = event.target.files; // FileList object
	$("#CustomDataBubbles").find(".dataBubble").remove();
	AgreeFunction = function (event) {
		for (var i = 0; i < FileReaderFile.length; i++) {
			var reader = new FileReader();
			reader.readAsText(FileReaderFile[i]);
			$("#CustomDataBubbles").append($('<h3 class="dataBubble ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom" style="text-align:center;padding:0.2em">')
				.text("User Data").attr('name',"CustomData").attr('title',"Custom"));
			$("#CustomDataBubbles").sortable({
				update : function (event, ui) {
				},
				items : ".dataBubble"
			});
			reader.onload = function () {
				// Normalize new lines
				var result = reader.result.replace(/[\r|[\r\n]]/g, "\n"); 
				$.each(rvDataSets, function (SpeciesIndex,rvds) {
					//Process File
					rvds.addCustomData($.csv.toObjects(result));
					if(rvds.CustomData.length >0){
						var customkeys = Object.keys(rvds.CustomData[0]);
					} else {
						var customkeys =[];
					}
					if ($.inArray("DataDescription", customkeys) >= 0) {
						$("#ImportDataFileDiv").find(".DataDescription").html(rvds.CustomData[0]["DataDescription"]);
						$("#CustomDataBubbles").find(".dataBubble").attr("title",rvds.CustomData[0]["DataDescription"].replace(/(<([^>]+)>)/ig,""));
					} else {
						$("#ImportDataFileDiv").find(".DataDescription").html("Data Description is missing.");
					}
					////clicky.log(window.location.pathname + window.location.hash,'User Data Import');
					$("#CustomDataBubbles").find(".dataBubble").attr("FileName",FileReaderFile[0].name);
					
					//Add Custom Lines support here. 
					if($.inArray("Residue_i", customkeys) >= 0) {
						if ($.inArray("Residue_j", customkeys) >= 0){
							var FullBasePairSet =[];
							var targetLayer = rvds.getLayerByType("lines");
							if ($.inArray("ColorCol", customkeys) >= 0) {
								$(".oneLayerGroup[name=" + targetLayer[0].LayerName + "]").find(".layerContent").find("div[name=llm]").find("select").multiselect("widget").find(":radio:eq(1)").each(function(){
									this.click();
								});
								var processColor=true;
							} else {
								var processColor=false;
							}
							if ($.inArray("Opacity", customkeys) >= 0) {
								var processOpacity=true;
							} else {
								var processOpacity=false;
							}
							if ($.inArray("LineWidth", customkeys) >= 0) {
								var processLineWidth=true;
							} else {
								var processLineWidth=false;
							}
							$.each(rvds.CustomData, function (index,value){
								if (processColor){
									var color = colorNameToHex(value.ColorCol);
								} else {
									var color = colorNameToHex("#231F20");
								}
								if (processOpacity){
									var Opacity = value.Opacity;
								} else {
									var Opacity = 0.5;
								}
								if (processLineWidth){
									var LineWidth = value.LineWidth;
								} else {
									var LineWidth = 1.0;
								}
								var grd = targetLayer[0].CanvasContext.createLinearGradient(MainResidueMap[value.Residue_i].X, MainResidueMap[value.Residue_i].Y, MainResidueMap[value.Residue_j].X, MainResidueMap[value.Residue_j].Y);
								grd.addColorStop(0, "rgba(" + h2d(color.slice(1, 3)) + "," + h2d(color.slice(3, 5)) + "," + h2d(color.slice(5)) + "," + Opacity + ")");
								grd.addColorStop(1, "rgba(" + h2d(color.slice(1, 3)) + "," + h2d(color.slice(3, 5)) + "," + h2d(color.slice(5)) + "," + Opacity + ")");
								
								FullBasePairSet.push({
									bp_type: value.Int_Type,
									color: grd,
									color_hex: color,
									opacity: Opacity,
									lineWidth: LineWidth,
									id: (index + 1).toString(),
									pairIndex: (index + 1).toString(),
									//resIndex1 : j,
									//resIndex2 : k
									residue_i: value.Residue_i,
									residue_j: value.Residue_j
								});
							});
							//rvds.BasePairs=FullBasePairSet;
							//rvds.FullBasePairSet=FullBasePairSet;
							
							targetLayer[0].DataLabel = FileReaderFile[0].name;
							$("[name=" + targetLayer[0].LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(targetLayer[0].DataLabel);
							rvds.drawBasePairs("lines",null);
						} else {
							alert("Expected Residue_j column. Please check input.");
						}
					}
				});
			};
		}
	};
	
	if (PrivacyStatus != "Agreed") {
		$("#Privacy-confirm").text(PrivacyString);
		CurrPrivacyCookie = "privacy_status_data";
		$("#Privacy-confirm").dialog('open');
	} else {
		AgreeFunction(event);
	}
}
	

function customDataProcess(data_label,targetLayerName){
	var colors=[];
	$.each(rvDataSets, function (index, rvds) {
		var targetLayer = rvds.getLayer(targetLayerName);
		targetLayer.DataLabel = data_label
		var NewData;
	
		$("[name=" + targetLayer.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text("User File:").append($("<br>")).append(targetLayer.DataLabel);
		targetLayer.clearData();
		
		if(rvDataSets[targetLayer.SetNumber].CustomData.length >0){
			var customkeys = Object.keys(rvds.CustomData[0]);
		} else {
			var customkeys =[];
		}

		NewData = CustomDataExpand(targetLayer);
		targetLayer.Data = NewData.IncludeData;
		var targetSelection = rvds.Selections[0];
		SelectionMenu(targetSelection);
		RefreshSelectionMenu();
		//Make new selection invisible. 
		$(".oneSelectionGroup[name=" + targetSelection.Name +"]").find(".checkBoxDIV-S").find(".visibilityCheckImg").attr("value","invisible");
		$(".oneSelectionGroup[name=" + targetSelection.Name +"]").find(".checkBoxDIV-S").find(".visibilityCheckImg").attr("src","/static/ribovision/images/invisible.png");
		rvds.drawSelection("selected");
		
		if ($.inArray("ColorPalette", customkeys) >= 0 || $.inArray("TwoColorMode", customkeys) >= 0){
			//console.log(rvDataSets[0].CustomData[0].ColorPalette);
			$.each(rvds.CustomData, function (index,value){
				if (value.ColorPalette && value.ColorPalette !="") {
					colors.push(value.ColorPalette);
				} else if (value.TwoColorMode && value.TwoColorMode !="") {
					colors.push(value.TwoColorMode);
				}else {
					return;
				}
			});
			//Assume if length one, user mean to name a predefined gradient.
			if (colors.length == 1){
				colors = window[colors[0]];
			}
			//console.log(colors);
		} else {
			colors = ColorLists.Rainbow1;
		}
		
		if (targetLayer.Type === "selected"){

		} else {
			if ($.inArray("FontWeight", customkeys) >= 0) {
				$.each(NewData.Weight, function (index,value){
					if (value != undefined) {
						rvds.Residues[index]["font-weight"]=value
					} else {
						rvds.Residues[index]["font-weight"]="normal";
					}
				});
				rvds.drawResidues("residues");
				rvds.refreshResiduesExpanded(targetLayer.LayerName);
				update3Dcolors();
			}
			if ($.inArray("TwoColorMode", customkeys) >= 0 & $.inArray("SwitchPoint", customkeys) >= 0) {
				var SwitchPoint=parseFloat(rvds.CustomData[0]["SwitchPoint"]);
			} else if ($.inArray("TwoColorMode", customkeys) >= 0){
				var SwitchPoint=0;
			} else {
				var SwitchPoint=undefined;
			}
			if ($.inArray("ColorCol", customkeys) >= 0) {
				rvds.drawResidues("residues");
				rvds.refreshResiduesExpanded(targetLayer.LayerName);
				update3Dcolors();
			} else if ($.inArray("DataCol", customkeys) >= 0) {
				colorProcess(NewData,undefined,targetLayer,colors,SwitchPoint);
			} else if ($.inArray("FontWeight", customkeys) >= 0){
				//Do nothing, maybe need more here later;
			} else {
				//alert("No recognized columns found. Please check input.");
			}
		}

		updateSelectionDiv(targetSelection.Name,targetLayer.SetNumber);
		drawNavLine();
	})
		
	//Finishes doing Nucleotides, move on to rProteins
	CustomProcessProteins(colors);
}

function CustomProcessProteins(colors){
	//var rProtein=undefined;

	$.each(rvDataSets, function(index,rvds){
		var NewData = [];
		var ColorProteins=new Array;
		for (var ii = 0; ii < rvds.CustomData.length; ii++) {
			var customkeys = Object.keys(rvds.CustomData[ii]);
			// Assume proteins are single. Come back and add range support sometime.
			//var command = rvds.CustomData[ii]["resNum"];

			//var rProtein = expandSelection([command], targetSelection.Name,0);

			var targetSelection = rvds.Selections[0];
			var rProtein = rvds.CustomData[ii]["resNum"];
			if (rProtein){
				if ($.inArray("DataCol", customkeys) >= 0) {
					ColorProteins.push({ResNum : rProtein, Color : undefined});
					if (isNaN(parseFloat(rvds.CustomData[ii]["DataCol"]))){
						NewData.push(rvds.CustomData[ii]["DataCol"]);
					} else {
						NewData.push(parseFloat(rvds.CustomData[ii]["DataCol"]));
					}
				} else if ($.inArray("ColorCol", customkeys) >= 0) {
					ColorProteins.push({ResNum : rProtein, Color : rvds.CustomData[ii]["ColorCol"]});
				}
			}
		}
		var dataIndices = colorProcess(NewData,undefined,undefined,colors,undefined,true);
		$.each(dataIndices, function (index,value){
			ColorProteins[index]["Color"] = colors[value];
		});
		rvds.ColorProteins = ColorProteins;
	})
	ColorProteins3D();
}

function ColorProteinsPyMOL(){
	var script = "";
	
	$.each(rvDataSets, function(index,rvds){
		var ColorProteins = rvds.ColorProteins;
		if (rvds.Residues[0] == undefined){return};
		
		
		
		// Protein Section
		for (var jj = 0; jj < rvds.SpeciesEntry.Molecule_Names_rProtein.length; jj++) {
			script += "copy " + rvds.SpeciesEntry.Species_Abr + "_rp_" 
				+ rvds.SpeciesEntry.Molecule_Names_rProtein[jj].replace(/\(/g, "_").replace(/\)/g, "") 
				+ "_custom, " + rvds.SpeciesEntry.Species_Abr + "_rp_" 
				+ rvds.SpeciesEntry.Molecule_Names_rProtein[jj].replace(/\(/g, "_").replace(/\)/g, "") 
				+ "\n";
		}
		
		$.each(ColorProteins, function (index,value){
			var ressplit = value.ResNum.split(":");
			if (ressplit[0] !== "undefined"){
				var curr_color = value.Color;
				//var h = rvds.SpeciesEntry.Molecule_Names_rProtein.indexOf(ressplit[0]);
				script += "color " + curr_color.replace("#", "0x") + ", " 
					+ rvds.SpeciesEntry.Species_Abr + "_rp_" 
					+ ressplit[0].replace(/\(/g,"_").replace(/\)/g,"") 
					+ "_custom" 
					+ " and resi " + ressplit[1].replace(/[^:]*:/g, "").replace(/[^:]*:/g, '') + "\n";
			}
		});
	})

	script += "\ndisable *rp*\n";
	return script;
}

function resetFileInput($element) {
	var clone = $element.clone(false, false);
	$element.replaceWith(clone);
	alert(42);
}

function CustomDataExpand(targetLayer){
	var SeleLen = 0;
	var NewData = [];
	var FontWeight = [];
	$.each(rvDataSets[targetLayer.SetNumber].Residues, function (index,value){
		NewData[index]=undefined;
		FontWeight[index]=undefined;
	});

	var ExtraData = [];
	if(rvDataSets[targetLayer.SetNumber].CustomData.length >0){
		var customkeys = Object.keys(rvDataSets[targetLayer.SetNumber].CustomData[0]);
	} else {
		var customkeys =[];
	}
	if($.inArray("resNum", customkeys) >= 0){
		for (var ii = 0; ii < rvDataSets[targetLayer.SetNumber].CustomData.length; ii++) {
			var command = rvDataSets[targetLayer.SetNumber].CustomData[ii]["resNum"].split(";");
			var targetSelection = rvDataSets[targetLayer.SetNumber].Selections[0];
			expandSelection(command, targetSelection.Name,targetLayer.SetNumber);
			var l = targetSelection.Residues.length;
			if (l == SeleLen){
				if ($.inArray("DataCol", customkeys) >= 0) {
					if (isNaN(parseFloat(rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"]))){
						ExtraData.push(rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"]);
					} else {
						ExtraData.push(parseFloat(rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"]));
					}
				}
			} else {
				for (var iii = SeleLen; iii < l; iii++) {
					if (targetSelection.Residues[iii].resNum >= 0) {
						//var ressplit = targetSelection.Residues[iii].resNum.split(":");
						//var ResName = rvDataSets[targetLayer.SetNumber].SpeciesEntry.RNA_Chains[rvDataSets[targetLayer.SetNumber].SpeciesEntry.Molecule_Names.indexOf(ressplit[0])] + "_" + ressplit[1];	
						var ResName = targetSelection.Residues[iii].uResName;
					} else {
						alert("this mode is being deprecated. This shouldn't happen any more");
					//	var chainID =  targetSelection.Residues[iii].ChainID;
					//	var ResName = chainID + "_" + targetSelection.Residues[iii].resNum;
					}
					var k = MainResidueMap[ResName].index;
					
					if ($.inArray("DataCol", customkeys) >= 0) {
						if (isNaN(parseFloat(rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"]))){
							NewData[k] = rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"];
						} else {
							NewData[k] = parseFloat(rvDataSets[targetLayer.SetNumber].CustomData[ii]["DataCol"]);
						}
					}
					if ($.inArray("ColorCol", customkeys) >= 0) {
						targetLayer.dataLayerColors[k] = colorNameToHex(rvDataSets[targetLayer.SetNumber].CustomData[ii]["ColorCol"]);
					}
					if ($.inArray("FontWeight", customkeys) >= 0) {
						FontWeight[k] = rvDataSets[targetLayer.SetNumber].CustomData[ii]["FontWeight"];
					}
					SeleLen = l;
				}
			}
			
		}
	}
	return {IncludeData : NewData,ExtraData : ExtraData, Weight : FontWeight}
}
///////////////////////////////////////////////////////////////////////////////


////////////////////////////// Cookie Functions ///////////////////////////////
function set_cookie(cookie_name, cookie_value, lifespan_in_days, valid_domain) {
	// http://www.thesitewizard.com/javascripts/cookies.shtml
	var domain_string = valid_domain ? ("; domain=" + valid_domain) : '';
	document.cookie = cookie_name +
		"=" + encodeURIComponent(cookie_value) +
		"; max-age=" + 60 * 60 *
		24 * lifespan_in_days +
		"; path=/" + domain_string;
}

function get_cookie(cookie_name) {
	// http://www.thesitewizard.com/javascripts/cookies.shtml
	var cookie_string = document.cookie;
	if (cookie_string.length != 0) {
		var pattern = cookie_name + "=[^;]*";
		var patt = new RegExp(pattern, "g");
		var cookie_value = cookie_string.match(patt);
		if (cookie_value == null) {
			return "";
		} else {
			var p = cookie_value[0].split("=");
			return decodeURIComponent(p[1]);
		}
	}
}

function checkSavePrivacyStatus() {
	var PrivacyStatus = get_cookie("privacy_status_text");
	var PrivacyString = "This feature requires uploading data to our server. Our privacy policy at this time is to not look at these files without permission."
		 + " We will set up autodelete for these files soon. Click \"I agree\" to acknowledge acceptance of our policy.";
	
	if (PrivacyStatus != "Agreed") {
		$("#Privacy-confirm").text(PrivacyString);
		CurrPrivacyCookie = "privacy_status_text";
		$("#Privacy-confirm").dialog('open');
	} else {
		AgreeFunction();
		//clicky.log(window.location.pathname + window.location.hash,'User Download','download');
		Histats_variables.push("SaveSomething","Yes");
		Histats.track_event('b');

	}
}
///////////////////////////////////////////////////////////////////////////////


//////////////////////////////// Save Functions ///////////////////////////////
function saveNavLine() {
	var tmp  = document.getElementById("NavLineDiv");
	var svg = tmp.getElementsByTagName("svg")[0];
	// Extract the data as SVG text string
	var svg_xml = (new XMLSerializer).serializeToString(svg);
	
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/save1D");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'svg' : svg_xml}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
}

function retrieveRvState(filename) {
	SaveStateFileName=filename;
	$.post('retrieveRvState.php', {
		datasetname : SaveStateFileName,
		username : UserName
	}, function (RvSaveState) {
		var rvSaveState = JSON.parse(RvSaveState);
		rvDataSets[0]=rvDataSets[0].fromJSON(rvSaveState["RvDS"]);
		// Re stringify a few things for compatibility / symmetry with local storage
		rvSaveState["rvLayers"] = JSON.stringify(rvDataSets[0].Layers);
		rvSaveState["rvSelections"] = JSON.stringify(rvDataSets[0].Selections);
		rvSaveState["rvLastSpecies"] = rvDataSets[0].Name;
		if($('input[name="3dp"][value=on]').is(':checked')){
			Jmol.script(myJmol, "script states/" + rvDataSets[0].SpeciesEntry.Jmol_Script);
			var jscript = "display " + rvDataSets[0].SpeciesEntry.Jmol_Model_Num_rRNA + ".1";
			Jmol.script(myJmol, jscript);
		}
		processRvState(rvSaveState);
	});
}

function storeRvState(filename){
	SaveStateFileName=filename;
	AgreeFunction = function () {
		var RvSaveState = {};
		RvSaveState["RvDS"] = JSON.stringify(rvDataSets[0]);
		RvSaveState["rvView"] = JSON.stringify(rvViews[0]);
		if($('input[name="3dp"][value=on]').is(':checked')){
			RvSaveState["rvJmolOrientation"] = Jmol.evaluateVar(myJmol,"script('show orientation')");
		}
	
		if($("input[name='PanelSizesCheck']").attr("checked")){
			var po = {
				PanelDivide : PanelDivide,
				TopDivide : TopDivide
			}
			RvSaveState["rvPanelSizes"] = JSON.stringify(po);
		}
		if($("input[name='MouseModeCheck']").attr("checked")){
			RvSaveState["rvMouseMode"] = onebuttonmode;	
		}
		
		$form = $("<form></form>");
		$('body').append($form);
		data = {
			content: JSON.stringify(RvSaveState,null,'\t'),
			datasetname : SaveStateFileName,
			username : UserName
		};
		$.post("storeRvState.php", data, function(d) {});
	}
	checkSavePrivacyStatus(SpeciesIndex);
}
function saveRvState(filename){
	SaveStateFileName=filename;
	AgreeFunction = function () {
		var RvSaveState = {};
		RvSaveState["RvDS"] = JSON.stringify(rvDataSets[0]);
		RvSaveState["rvView"] = JSON.stringify(rvViews[0]);
		if($('input[name="3dp"][value=on]').is(':checked')){
			RvSaveState["rvJmolOrientation"] = Jmol.evaluateVar(myJmol,"script('show orientation')");
		}
		
		if($("input[name='PanelSizesCheck']").attr("checked")){
			var po = {
				PanelDivide : PanelDivide,
				TopDivide : TopDivide
			}
			RvSaveState["rvPanelSizes"] = JSON.stringify(po);
		}
		if($("input[name='MouseModeCheck']").attr("checked")){
			RvSaveState["rvMouseMode"] = onebuttonmode;	
		}
		
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "saveRvState.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", JSON.stringify(RvSaveState,null,'\t'));
		var hiddenField2 = document.createElement("input");
		hiddenField2.setAttribute("type", "hidden");
		hiddenField2.setAttribute("name", "datasetname");
		hiddenField2.setAttribute("value", SaveStateFileName);
		form.appendChild(hiddenField);
		form.appendChild(hiddenField2);
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus(SpeciesIndex);
}
function saveSVG() {
	var CS = canvasToSVG();
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/save2D");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'svg' : CS.SVG, 'ext' : 'svg'}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
}

function saveJPG() {
	var CS = canvasToSVG();
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/save2D");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'svg' : CS.SVG, 'ext' : 'jpg'}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
    //JSON.stringify([600/72*CS.PaperSize[0],600/72*CS.PaperSize[1]]));
}

function savePNG() {
	var CS = canvasToSVG();
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/save2D");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'svg' : CS.SVG, 'ext' : 'png'}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
    //JSON.stringify([600/72*CS.PaperSize[0],600/72*CS.PaperSize[1]]));
}

function savePDF() {
	var CS = canvasToSVG();
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/save2D");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'svg' : CS.SVG, 'ext' : 'pdf'}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
    //JSON.stringify([600/72*CS.PaperSize[0],600/72*CS.PaperSize[1]]));
}

function savePML(){
	var structureName = rvDataSets[0].SpeciesEntry.Species_Abr;
	var dsLayers=[];
	$.each(rvDataSets, function(SpeciesIndex,rvds){
		dsLayers = dsLayers.concat(rvds.getLayerByType(["residues","circles","contour","selected"]));
	});
	
	//Form Submit;
	var form = document.createElement("form");
	form.setAttribute("method", "post");
	form.setAttribute("action", "api/RiboVision/v1.0/savepml");
	form.setAttribute("target", "_blank");
	var hiddenField = document.createElement("input");
	hiddenField.setAttribute("type", "hidden");
	hiddenField.setAttribute("enctype", "text/plain");
	hiddenField.setAttribute("name", "data");
	hiddenField.setAttribute("value", JSON.stringify({'StructureName' : structureName, 'Layers' : dsLayers}));
	form.appendChild(hiddenField);
	document.body.appendChild(form);
	form.submit();
}
/* function savePML() {
	AgreeFunction = function () {
		var script = "";
		var structureName = rvDataSets[0].SpeciesEntry.Species_Abr
		//var PDB_Obj_Names = [];
		//var PDB_files =[];
		
		//Default option
		script += "set bg_rgb, white\n";
		//mmCif File. Assume first and second structure (subunits) come from the same cif file. 
		script += "load " + rvDataSets[0].SpeciesEntry.StructureName + ".cif, " + structureName + "\n";
		//script += "as cartoon, " + PDB_Obj_Names[0] + "\n";
		script += "disable " + structureName + "\n";
		
		$.each(rvDataSets, function(SpeciesIndex,rvds){
			script += "\n";
			// Layers to PyMOL
			var dsLayers = rvds.getLayerByType(["residues","circles","contour","selected"]);
			$.each(dsLayers, function (key, value) {
				script += layerToPML(structureName,value,SpeciesIndex);
			});
			script += "\n";
			
			//Proteins to PyMOL
			script += proteinsToPML(structureName,SpeciesIndex);
			script += "\n";
			
			//Proteins to PyMOL (Custom)
			script += ColorProteinsPyMOL();
			script += "\n";
			
			//Selection to PyMOL
			$.each(rvds.Selections, function (key, value) {
				script += selectionToPML(structureName,value,SpeciesIndex);
			});
			script += "\ndisable RV_Sele_*\n";
			
			//Default zoom
			script += "\nzoom all\n";
		});
		
		//Form Submit;
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "savePML.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", script);
		form.appendChild(hiddenField);
		
		// PDB_filesU = $.grep(PDB_files, function (v, k) {
			// return $.inArray(v, PDB_files) === k;
		// });
		var hiddenField2 = document.createElement("input");
		hiddenField2.setAttribute("type", "hidden");
		hiddenField2.setAttribute("name", "pdbfiles");
		hiddenField2.setAttribute("value",  rvDataSets[0].SpeciesEntry.StructureName + '.cif');
		form.appendChild(hiddenField2);
		
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus();
} */

function layerToPML(PDB_Obj_Names,targetLayer,SpeciesIndex) {
	var PyMOL_obj = [];
	var script = "";
	var r0,r1,curr_chain,curr_color;
	if (rvDataSets[SpeciesIndex].Residues[0] == undefined){return};
	
	for (var j = 0; j < rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names.length; j++) {
		PyMOL_obj[j] = rvDataSets[SpeciesIndex].SpeciesEntry.Species_Abr + "_" + rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[j] + "_" + targetLayer.LayerName;
		script += "create " + PyMOL_obj[j] + ", " + PDB_Obj_Names[0] + " and chain " + rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains[j] + "\n";
		if (targetLayer.Linked){
			script += "enable " + PyMOL_obj[j] + "\n";
		} else {
			script += "disable " + PyMOL_obj[j] + "\n";
		}
	}
	script += "\n";
	
	r0 = rvDataSets[SpeciesIndex].Residues[0].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
	curr_chain = rvDataSets[SpeciesIndex].Residues[0].ChainID;
	curr_color = colorNameToHex(targetLayer.dataLayerColors[0]);
	if (!curr_color || curr_color === '#000000') {
		curr_color = '#858585';
	}
	for (var i = 1; i < rvDataSets[SpeciesIndex].Residues.length; i++) {
		var residue = rvDataSets[SpeciesIndex].Residues[i];
		var residueLast = rvDataSets[SpeciesIndex].Residues[i - 1];
		var residueLastColor = targetLayer.dataLayerColors[i - 1];
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
				r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
			} else if (residue.ChainID == null) {
				curr_chain = residue.ChainID;
				curr_color = colorNameToHex(targetLayer.dataLayerColors[i]);
				if (!curr_color || curr_color === '#000000') {
					curr_color = '#858585';
				}
				r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
			} else {
				if (!targetLayer.dataLayerColors[i]){
					compare_color = '#858585';
				} else {
					compare_color = colorNameToHex(targetLayer.dataLayerColors[i]);
				}
				if (((compare_color != colorNameToHex(residueLastColor)) || (curr_chain != residue.ChainID)) || (i == (rvDataSets[SpeciesIndex].Residues.length - 1))) {
					r1 = residueLast.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, ""); ;
					
					if (r0 === r1){
						if (colorNameToHex(residueLastColor).indexOf("#") == -1) {
							script += "color 0x" + curr_color + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "\n";
						} else {
							script += "color " + curr_color.replace("#", "0x") + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "\n";
						}
					} else {
						if (colorNameToHex(residueLastColor).indexOf("#") == -1) {
							script += "color 0x" + curr_color + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "-" + r1 + "\n";
						} else {
							script += "color " + curr_color.replace("#", "0x") + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "-" + r1 + "\n";
						}
					}
											
					r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
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
	if (colorNameToHex(residueLastColor).indexOf("#") == -1) {
		script += "color 0x" + curr_color + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "-" + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "") + "\n";
	} else {
		script += "color " + curr_color.replace("#", "0x") + ", " + PyMOL_obj[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(curr_chain)] + " and resi " + r0 + "-" + residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "") + "\n";
	}
	script += "\n";
	
	return script;
}

function proteinsToPML(structureName,SpeciesIndex){
	var script = "";
	var curr_color;
	// Protein Section
	for (var jj = 0; jj < rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names_rProtein.length; jj++) {

		script += "create " + rvDataSets[SpeciesIndex].SpeciesEntry.Species_Abr + "_rp_" 
		+ rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names_rProtein[jj].replace(/\(/g, "_").replace(/\)/g, "") 
		+ ", " + structureName + " and chain " + rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains_rProtein[jj].replace(/\(/g, "_").replace(/\)/g, "") 
		+ "\n";
	}
	script += "\ndisable *_rp_*\n";
	script += "color wheat, *_rp_*\n";
	
	var array_of_checked_values = $("#ProtList").multiselect("getChecked").map(function () {
			return this.value;
		}).get();
	
	for (jjj = 0; jjj < array_of_checked_values.length; jjj++) {
		var h = rvDataSets[SpeciesIndex].SpeciesEntry.SubunitProtChains[2].indexOf(array_of_checked_values[jjj]);
		var ProtName = $.grep($("#ProtList").multiselect("getChecked"), function(e) {
			return e.value == array_of_checked_values[jjj];
		});
		if (h >= 0) {
			curr_color = rgb2hex($(ProtName).next().css("color"));
			if(curr_color.indexOf("#") == -1) {
				curr_color = "0x" + curr_color;
			} else {
				curr_color = curr_color.replace("#", "0x");
			}
			script += "color " + curr_color + ", " + rvDataSets[SpeciesIndex].SpeciesEntry.Species_Abr + "_rp_" + rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names_rProtein[h].replace(/\(/g,"_").replace(/\)/g,"") + "\n";
			script += "enable " + rvDataSets[SpeciesIndex].SpeciesEntry.Species_Abr + "_rp_" + rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names_rProtein[h].replace(/\(/g,"_").replace(/\)/g,"") + "\n";
		}
	}
	return script;
}

function selectionToPML(structureName,targetSelection,SpeciesIndex){
	var script = "";
	var PyMOL_obj = "RV_Sele_" + targetSelection.Name;
	var r0,r1,curr_chain;
	var DoneNow=false;
	if (rvDataSets[SpeciesIndex].Residues[0] == undefined){return};
	
	script += "create " + PyMOL_obj + ", " + "resi 0\n";
	
	var SeleResidues=targetSelection.Residues.sort(function (a, b) {;
		return (Number(a.map_Index) - Number(b.map_Index));
	});
	
	if(SeleResidues.length==0){return ''};
	r0 = SeleResidues[0].resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
	curr_chain = SeleResidues[0].ChainID;
	for (var i = 1; i < SeleResidues.length; i++) {
		var residue = SeleResidues[i];
		var residueLast = SeleResidues[i - 1];
		
		if (residue.ChainID != "") {
			if (curr_chain == "") {
				curr_chain = residue.ChainID;
				r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
			} else if (residue.ChainID == null) {
				curr_chain = residue.ChainID;
				r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
			} else {
				if ((residue.map_Index - residueLast.map_Index > 1 ) || (curr_chain != residue.ChainID) || (i == (SeleResidues.length - 1))) {
					if ((i == (SeleResidues.length - 1)) && (curr_chain == residue.ChainID) && (residue.map_Index - residueLast.map_Index == 1 )){
						r1 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
						DoneNow=true;
					} else {
						r1 = residueLast.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
					}
					
					if (r0 === r1){
						script += "create " + PyMOL_obj + ", " + PyMOL_obj + " or (" + structureName + " and chain " + curr_chain + " and resi " + r0 + ")\n";
					} else {
						script += "create " + PyMOL_obj + ", " + PyMOL_obj + " or (" + structureName + " and chain " + curr_chain + " and resi " + r0 + "-" + r1 + ")\n";
					}
											
					r0 = residue.resNum.replace(/[^:]*:/g, "").replace(/[^:]*:/g, "");
					if (residue.ChainID != "") {
						curr_chain = residue.ChainID;
					}
				}
			}
		}
	}
	if(!DoneNow){
		script += "create " + PyMOL_obj + ", " + PyMOL_obj + " or (" + structureName + " and chain " + curr_chain + " and resi " + r0 + ")\n";
	}
	script += "color " + targetSelection.Color.replace("#", "0x") + ", " + PyMOL_obj + "\n\n";
	return script;
}

function canvasToSVG() {
	var ChosenSide;
	var AllMode = $('input[name="savelayers"][value=all]').attr("checked");
	var resMod = $('input[name="ptmod"][value=on]').is(':checked');
	
	//This section is assuming the conditions, 1) only two structures allowed max, 2) only one can be landscape.
	//Condition 1 is active. Condition 2 will be forced soon.
	if (!rvDataSets[1] && rvDataSets[0].SpeciesEntry.Orientation == "landscape"){
		var paperSize=[792,612];
	} else if (!rvDataSets[1] && rvDataSets[0].SpeciesEntry.Orientation != "landscape"){
		var paperSize=[612,792];
	} else {
		if (rvDataSets[0].SpeciesEntry.Orientation == "landscape" || rvDataSets[1].SpeciesEntry.Orientation == "landscape"){
			var paperSize=[1404,792];
		} else {
			var paperSize=[1224,792];
		}
	}
	
	var mapsize = paperSize[0] + " " + paperSize[1];
	var mapsize2 = 'width="' + paperSize[0] +'" height="' + paperSize[1] + '" ';
	
	output = "<?xml version='1.0' encoding='UTF-8'?>\n" +
		'<svg version="1.1" baseProfile="basic" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" ' +
		mapsize2 + 'viewBox="0 0 ' + mapsize + '" xml:space="preserve">\n';
	
	
	$.each(rvDataSets, function (SpeciesIndex, rvds) {
		var pageOffsetX = rvds.PageOffset[0];
		var pageOffsetY = rvds.PageOffset[1];
		
		if (rvds.Font_Size_SVG){
			var Font_Size_SVG=parseFloat(rvds.Font_Size_SVG);
		} else {
			var Font_Size_SVG=3.9;
		}
		
        //Replace with database entries. 
		//This will not work correctly any more, due to the two structure mode. Will need to shift these when loaded second.
		// var elReq = $.ajax({
			// url: "images/" + rvds.Name + "_ExtraLabels.svg",
			// dataType: "text", 
			// cache: false,
			// async: false
		 // });
	 
		// elReq.done(function (data) {
			// output = output + data;
		// });
		output += '<g id="Structure_' + (SpeciesIndex + 1) + '">\n';
		$.each(rvds.Layers, function (index, value) {
			if (rvds.Circle_Radius){
				var radius=parseFloat(rvds.Circle_Radius) * value.ScaleFactor;
			} else {
				var radius= 1.7 * value.ScaleFactor;
			}
			var undefined_color='undefined';
			
			if (AllMode || value.Visible){
				switch (value.Type) {
					case "lines":
						output += '<g id="' + value.LayerName + '">\n';
						$.each(value.Data, function(index,base_pair){
							var residue_i = MainResidueMap[base_pair.residue_i];
							var residue_j = MainResidueMap[base_pair.residue_j];
							output = output + '<line fill="none" stroke="' + base_pair.color_hex + '" stroke-opacity="' + base_pair.opacity + '" stroke-width="' + base_pair.lineWidth + '" stroke-linejoin="round" stroke-miterlimit="10" x1="' + residue_i.X.toFixed(3) + '" y1="' + residue_i.Y.toFixed(3) + '" x2="' + residue_j.X.toFixed(3) + '" y2="' + residue_j.Y.toFixed(3) + '"/>\n';
						});
						output = output + '</g>\n';
						break;
					case "contour" :
						output = output + '<g id="' + value.LayerName + '">\n';
						if (value.OutLineMode){
							output = output + '<g id="' + 'Outline' + '">\n';
							for (var j = 0; j < rvDataSets[SpeciesIndex].ExtraContourLineSegments.length; j++) {
								var ExtraContourLineSegment = rvDataSets[SpeciesIndex].ExtraContourLineSegments[j];
								output = output + '<polyline fill="none" stroke-dasharray="4, 4" stroke-linecap="round" stroke="' + '#000000' + '" stroke-opacity="' + '1' + '" stroke-width="' + '1.25' + 
								'" stroke-linejoin="round" stroke-miterlimit="10" points="' + parseFloat(ExtraContourLineSegment.X1).toFixed(3) + ',' + parseFloat(ExtraContourLineSegment.Y1).toFixed(3)
								+ ' ' + parseFloat(ExtraContourLineSegment.X2).toFixed(3) + ',' + parseFloat(ExtraContourLineSegment.Y2).toFixed(3)
								+ ' "/>\n';
							}
							for (var j = 0; j < rvDataSets[SpeciesIndex].ContourLinePoints.length; j++) {
								var ContourLinePoint = rvDataSets[SpeciesIndex].ContourLinePoints[j];
								if (value.dataLayerColors[j]){
									var PointColor='#000000';
								} else {
									var PointColor=undefined_color;
								}
								
								output = output + '<polyline fill="none" stroke-linecap="round" stroke="' + PointColor + '" stroke-opacity="' + '1' + '" stroke-width="' + value.ScaleFactor * 1.5 * 1.9 * radius + 
								'" stroke-linejoin="round" stroke-miterlimit="10" points="' + parseFloat(ContourLinePoint.X1).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y1).toFixed(3)
								+ ' ' + parseFloat(ContourLinePoint.X2).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y2).toFixed(3)
								+ ' ' + parseFloat(ContourLinePoint.X3).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y3).toFixed(3)
								+ ' "/>\n';
							}
							output = output + '</g>\n';
						}
						output = output + '<g id="' + 'Main_Data' + '">\n';
						for (var j = 0; j < rvDataSets[SpeciesIndex].ContourLinePoints.length; j++) {
							var ContourLinePoint = rvDataSets[SpeciesIndex].ContourLinePoints[j];
							output = output + '<polyline fill="none" stroke-linecap="round" stroke="' + value.dataLayerColors[j] + '" stroke-opacity="' + '1' + '" stroke-width="' + value.ScaleFactor * 1.9 * radius + 
							'" stroke-linejoin="round" stroke-miterlimit="10" points="' + parseFloat(ContourLinePoint.X1).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y1).toFixed(3)
							+ ' ' + parseFloat(ContourLinePoint.X2).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y2).toFixed(3)
							+ ' ' + parseFloat(ContourLinePoint.X3).toFixed(3) + ',' + parseFloat(ContourLinePoint.Y3).toFixed(3)
							+ ' "/>\n';
						}
						output = output + '</g>\n';
						output = output + '</g>\n';
						break;
					case "labels":
						output = output + '<g id="' + value.LayerName + '_Lines">\n';
						for (var iii = 0; iii < rvds.rvLineLabels.length; iii++) {
							var LabelLine = rvds.rvLineLabels[iii];
							output = output + '<line fill="none" stroke="#211E1F" stroke-width="0.25" stroke-linejoin="round" stroke-miterlimit="10" x1="' + (parseFloat(LabelLine.X1) + pageOffsetX).toFixed(3) + '" y1="' + (parseFloat(LabelLine.Y1) + pageOffsetY).toFixed(3) + '" x2="' + (parseFloat(LabelLine.X2) + pageOffsetX).toFixed(3) + '" y2="' + (parseFloat(LabelLine.Y2) + pageOffsetY).toFixed(3) + '"/>\n';
						}
						output = output + '</g>\n';
						
						output = output + '<g id="' + value.LayerName + '_Text">\n';
						for (var ii = 0; ii < rvds.rvTextLabels.length; ii++) {
							var LabelData = rvds.rvTextLabels[ii];
							output = output + '<text transform="matrix(1 0 0 1 ' + (parseFloat(LabelData.X) + pageOffsetX).toFixed(3) + ' ' + (parseFloat(LabelData.Y) + pageOffsetY).toFixed(3) + ')" fill="' + LabelData.Fill + '" font-family="Myriad Pro" font-size="' + LabelData.FontSize + '">' + LabelData.LabelText + '</text>\n';
						}
						output = output + '</g>\n';
						break;
					case "residues":
						output = output + '<g id="' + value.LayerName + '">\n';
						var xcorr = -0.439 * Font_Size_SVG + 0.4346; // magic font corrections.
						var ycorr = 0.2944 * Font_Size_SVG - 0.0033;
						
						// Add PageOffset to corrections for convenience
						xcorr += pageOffsetX;
						ycorr += pageOffsetY;
						
						//console.log(xcorr,ycorr);
						for (var i = 0; i < rvds.Residues.length; i++) {
							var residue = rvds.Residues[i];
							if (resMod){
								var resName = residue.modResName;
								if(residue.modResName == '&'){
									resName = '&amp;';
								}
								if(residue.modResName == '<'){
									resName = '&lt;';
								}
							} else {
								var resName = residue.resName;
							}				
							
							output = output + '<text id="' + residue.resName + "_" + residue.uResName + '" transform="matrix(1 0 0 1 ' + (parseFloat(residue.X) + xcorr).toFixed(3) + ' ' + (parseFloat(residue.Y) + ycorr).toFixed(3) + ')" fill="' + residue.color + '" font-family="Myriad Pro" ' + 'font-weight="' + residue["font-weight"] + '" font-size="' + Font_Size_SVG + '">' + resName + '</text>\n';
						}
						output = output + '</g>\n';
						break;
					case "circles":
						output = output + '<g id="' + value.LayerName + '">\n';
						
						for (var i = 0; i < rvds.Residues.length; i++) {
							var residue = rvds.Residues[i];
							if (residue && value.dataLayerColors[i]) {
								if (value.Filled) {
									output = output + '<circle id="' + residue.uResName + '" fill="' + value.dataLayerColors[i] + '" stroke="' + value.dataLayerColors[i] + '" stroke-width="0.5" stroke-miterlimit="10" cx="' + (parseFloat(residue.X) + pageOffsetX).toFixed(3) + '" cy="' + (parseFloat(residue.Y).toFixed(3)+ pageOffsetY) + '" r="' + radius + '"/>\n';
								} else {
									output = output + '<circle id="' + residue.uResName + '" fill="' + 'none' + '" stroke="' + value.dataLayerColors[i] + '" stroke-width="0.5" stroke-miterlimit="10" cx="' + (parseFloat(residue.X) + pageOffsetX).toFixed(3) + '" cy="' + (parseFloat(residue.Y) + pageOffsetY).toFixed(3) + '" r="' + radius + '"/>\n';
								}
							}
						}
						output = output + '</g>\n';
						break;
					case "selected":
						output = output + '<g id="' + value.LayerName + '">\n';
						var SelectionList =[];
						$('.checkBoxDIV-S').find(".visibilityCheckImg[value=visible]").parent().parent().each(function (index){SelectionList.push($(this).attr("name"))});

						$.each(SelectionList, function (index,SelectionName) {
							var targetSelection = rvds.getSelection(SelectionName);
							output = output + '<g id="' + targetSelection.Name + '">\n';
							$.each(targetSelection.Residues, function (index,residue){
								output = output + '<circle id="' + residue.uResName + '" fill="' + 'none' + '" stroke="' + targetSelection.Color + '" stroke-width="0.5" stroke-miterlimit="10" cx="' + (parseFloat(residue.X) + pageOffsetX).toFixed(3) + '" cy="' + (parseFloat(residue.Y) + pageOffsetY).toFixed(3) + '" r="' + radius + '"/>\n';
							});
							output = output + '</g>\n';
						});
						output = output + '</g>\n';
						break;
						
					default:
						break;
				}
			}
		});
		output = output + '</g>\n';
	});
	output = output + watermark(true);
	output = output + '</svg>';
	return { 'SVG': output, "PaperSize" : paperSize };
}
function computeSeqTable(SpeciesIndex){
	var WholeSet="";
	var OutputString="";
	var WholeSet= WholeSet + "WholeSet,WholeSet\n";
	var WholeSet= WholeSet + "resNum,resName\n";
	$.each(rvDataSets[SpeciesIndex].Residues, function (index,residue) {
		WholeSet+= rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(residue.ChainID)] + ":" + residue.resNum.replace(/[^:]*:/g, "") + "," + residue.resName + "\n";
	});
	var OutputObject = $.csv.toArrays(WholeSet);
	$.each(rvDataSets[SpeciesIndex].Selections, function (index, selection) {
		if (selection.Residues.length >0){
			OutputObject[0].push(selection.Name,selection.Name);
			OutputObject[1].push("resNum","resName");
			$.each(selection.Residues, function (index, residue) {
				OutputObject[2+index].push(rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(residue.ChainID)] + ":" + residue.resNum.replace(/[^:]*:/g, ""),residue.resName);
			});
		}
	});
	
	$.each(OutputObject, function (index, outputline){
		OutputString+=outputline.toString() + "\n";
	});
	
	
	return OutputString;
}
function saveSeqTable(SpeciesIndex){
	AgreeFunction = function (SpeciesIndex) {
		var ST = computeSeqTable(SpeciesIndex);
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "saveSeqTable.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", "\ufeff" + ST);
		form.appendChild(hiddenField);
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus(SpeciesIndex);
}
function computeSeqDataTable(SpeciesIndex){
	var WholeSet="";
	var OutputString="";
	var WholeSet= WholeSet + "WholeSet,WholeSet\n";
	var WholeSet= WholeSet + "resNum,resName\n";
	$.each(rvDataSets[SpeciesIndex].Residues, function (index,residue) {
		WholeSet+= rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(residue.ChainID)] + ":" + residue.resNum.replace(/[^:]*:/g, "") + "," + residue.resName + "\n";
	});
	var OutputObject = $.csv.toArrays(WholeSet);
	$.each(rvDataSets[SpeciesIndex].Layers, function (index, targetLayer) {
		if (targetLayer.Data.length >0){
			if (targetLayer.Type == "lines"){
				/*$.each(targetLayer.Data, function (index, data) {
					OutputObject[2+index].push(data);
				});*/
			} else {
				OutputObject[0].push(targetLayer.LayerName);
				OutputObject[1].push('"' + targetLayer.DataLabel.replace(/,/g,'\\comma\\') + '"'); // Come back and solve this comma,tab nonsense.
				$.each(targetLayer.Data, function (index, data) {
					OutputObject[2+index].push(data);
				});
			}
		}
	});
	
	$.each(OutputObject, function (index, outputline){
		OutputString+=outputline.toString() + "\n";
	});
	
	return OutputString.replace(/,/g,'\t').replace(/\\comma\\/g,',');
}
function computeInteractionDataTable(){
	var WholeSet="";
	var WholeSet= WholeSet + "Residue_i,ResidueName_i,Residue_j,ResidueName_j,Int_Type,ColorCol\n";
	//var WholeSet= WholeSet + "resNum,resName,resNum,resName,bp_type\n";
	
	$.each(rvDataSets[SpeciesIndex].BasePairs, function (index,basepair) {
		var j = basepair.resIndex1;
		var k = basepair.resIndex2;
		
		if (rvDataSets[SpeciesIndex].Residues[j].resNum.indexOf(":") >= 0 ){
			var ResName1 = rvDataSets[SpeciesIndex].Residues[j].resNum;
		} else {
			var ResName1 = rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[SpeciesIndex].Residues[j].ChainID)] +
			":" + rvDataSets[SpeciesIndex].Residues[j].resNum;
		}
		if (rvDataSets[SpeciesIndex].Residues[k].resNum.indexOf(":") >= 0 ){
			var ResName2 = rvDataSets[SpeciesIndex].Residues[k].resNum;
		} else {
			var ResName2 = rvDataSets[SpeciesIndex].SpeciesEntry.Molecule_Names[rvDataSets[SpeciesIndex].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[SpeciesIndex].Residues[k].ChainID)] +
			":" + rvDataSets[SpeciesIndex].Residues[k].resNum;
		}
	
		WholeSet+= ResName1 + "," + rvDataSets[SpeciesIndex].Residues[j].resName + "," + ResName2 + "," + rvDataSets[SpeciesIndex].Residues[k].resName + "," + basepair.bp_type + "," + basepair.color_hex + "\n";
	});
	
	return WholeSet.replace(/,/g,'\t');
}
function saveInteractionDataTable(SpeciesIndex){
	AgreeFunction = function (SpeciesIndex) {
		var SDT = computeInteractionDataTable(SpeciesIndex);
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "saveInteractionDataTable.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", "\ufeff" + SDT);
		form.appendChild(hiddenField);
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus(SpeciesIndex);
}
function saveSeqDataTable(SpeciesIndex){
	AgreeFunction = function (SpeciesIndex) {
		var SDT = computeSeqDataTable(SpeciesIndex);
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "saveSeqDataTable.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", "\ufeff" + SDT);
		form.appendChild(hiddenField);
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus(SpeciesIndex);
}

function saveFasta(){
	AgreeFunction = function () {
		var FA = computeFasta();
		var form = document.createElement("form");
		form.setAttribute("method", "post");
		form.setAttribute("action", "saveFasta.php");
		form.setAttribute("target", "_blank");
		var hiddenField = document.createElement("input");
		hiddenField.setAttribute("type", "hidden");
		hiddenField.setAttribute("name", "content");
		hiddenField.setAttribute("value", FA);
		form.appendChild(hiddenField);
		document.body.appendChild(form);
		form.submit();
	}
	checkSavePrivacyStatus();
}
function computeFasta(){
	var WholeSet="";
	var OutputString="";
	$.each(rvDataSets, function(index, rvds){
		OutputString+=">" + rvds.Name + "\n";
		
		$.each(rvds.Residues, function (index,residue) {
			OutputString+= residue.resName;
		});
		OutputString+="\n\n"
		
		$.each(rvds.Selections, function (index, selection) {
			if (selection.Residues.length >0){
				OutputString+=">" + selection.Name + "\n";
				$.each(selection.Residues, function (index, residue) {
					OutputString+= residue.resName;
				});
				OutputString+="\n\n"
			}
		});
	});
	return OutputString;
}
///////////////////////////////////////////////////////////////////////////////



/////////////////////////////// Load Data Functions ///////////////////////////
function populateInteractionMenu(structureName) {
	$.ajax({
		url: 'api/RiboVision/v1.0/fetchInteractionsMenu',
		type: 'POST',
		contentType: 'application/json', 
		accept: 'application/json',
		data: JSON.stringify([structureName]),
		cache: false,
		success: function(BPList) {
			var il = document.getElementById("PrimaryInteractionList");
			$.each(BPList, function(index, bp_entry){
				il.options[index + 1] = new Option(bp_entry.bp_group, bp_entry.bp_group);
			})
			$("#PrimaryInteractionList").multiselect("refresh");
		}
	})
	
	// $.each(rvDataSets, function(index,rvds){	
		// var il = document.getElementById("PrimaryInteractionList");
		// var BPList = rvds.SpeciesEntry.InterActionMenu.split(";");
		// if (speciesIndex == 0 && BPList[0] != ""){
			// for (var iii = 0; iii < BPList.length; iii++) {
				// var NewBPair = BPList[iii].split(":");
				// if (il.options[iii + 1]) {
					// il.options[iii + 1].value = NewBPair[1] + ';' + il.options[iii + 1].value;
				// } else {
					// il.options[iii + 1] = new Option(NewBPair[0], NewBPair[1]);
				// }
			// }
		// } else if (BPList[0] != ""){
			// for (var iii = 0; iii < BPList.length; iii++) {
				// var NewBPair = BPList[iii].split(":");
				// if (il.options[iii + 1]) {
					// il.options[iii + 1].value = il.options[iii + 1].value + ';' + NewBPair[1];
				// } else {
					// il.options[iii + 1] = new Option(NewBPair[0], NewBPair[1]);
				// }
			// }
		// }
	// })
}

function populateStructDataMenu(structureName) {
	$.ajax({
		url: 'api/RiboVision/v1.0/structdatamenu',
		type: 'POST',
		contentType: 'application/json', 
		accept: 'application/json',
		data: JSON.stringify([structureName]),
		cache: false,
		success: function(structureDataMenu) {
			$.each(structureDataMenu, function(index, NewSDPair){
				var ColName = NewSDPair.ColName;
				var description = NewSDPair.Description;
				if (ColName && description){
					var title = NewSDPair.StructDataName + ": " + description;
				} else if (ColName=='""'){
					var title = "None: This clears circles and makes letters black.";
				} else {
					var title = "Data Description is missing.";
				}
				$("#StructDataBubbles").append($('<h3 class="dataBubble ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom" style="text-align:center;padding:0.2em">')
				.text(NewSDPair.StructDataName).attr('name',ColName).attr('title',title).data("ColName",ColName).data("OverRideColors",NewSDPair.ColorList).data("indexMode",NewSDPair.IndexMode).data("rePlaceData",NewSDPair.ExtraArg));
			})
		}
	})
	
	/*$.each(rvDataSets, function(index,rvds){
		var SDList = rvds.SpeciesEntry.StructDataMenu.split(";");
		if (SDList[0] != "") {
			for (var ii = 0; ii < SDList.length; ii++) {
				var NewSDPair = SDList[ii].split(":");
				var ColName = NewSDPair[1].match(/[^\'\\,]+/);
				var result = $.grep(rvds.DataDescriptions, function(e){ return e.ColName === ColName[0]; });
				if (ColName[0] && result[0]){
					var title = NewSDPair[0] + ": " + result[0].Description;
				} else if (ColName[0]=='""'){
					var title = "None: This clears circles and makes letters black.";
				} else {
					var title = "Data Description is missing.";
				}
				$("#StructDataBubbles").append($('<h3 class="dataBubble ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom" style="text-align:center;padding:0.2em">')
				.text(NewSDPair[0]).attr('name',NewSDPair[1]).attr('title',title));
			}
			
		}
	})*/
}

function populateAlignmentMenu() {
	// Assumes list of alignments (entropies) are identical for both structures.
	// var AlnList = rvDataSets[0].SpeciesEntry.AlnMenu.split(";");
	// $.each(rvDataSets, function(index,rvds){
		// if (AlnList[0] != "") {
			// for (var ii = 0; ii < AlnList.length; ii++) {
				// var NewAlnPair = AlnList[ii].split(":");
				// var ColName = NewAlnPair[1].match(/[^\'\\,]+/);
				// var result = $.grep(rvds.DataDescriptions, function(e){ return e.ColName === ColName[0]; });
				// if (result[0]){
					// var title = NewAlnPair[0] + ": " + result[0].Description;
				// } else {
					// var title = "Data Description is missing.";
				// }
				// $("#AlnBubbles").append($('<h3 class="dataBubble ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom" style="text-align:center;padding:0.2em">')
				// .text(NewAlnPair[0]).attr('name',NewAlnPair[1]).attr('title',title));
			// }
		// }
	// })
}

function populateProteinMenu() {
	// var pl = document.getElementById("ProtList");
	// pl.options.length = 0;
	// var title='';
	// $.each(rvDataSets, function(speciesIndex,rvds){
		// $.each(rvds.SpeciesEntry.Molecule_Names_rProtein, function(){
			// var ColName = ["All_Proteins"];
			// var result = $.grep(rvds.DataDescriptions, function(e){ return e.ColName === ColName[0]; });
			// if (result[0]){
				// title = "All_Proteins" + ": " + result[0].Description;
			// } else {
				// title = "Data Description is missing.";
			// }
		// })
		// $('#ProtList').append('<optgroup label="Proteins' + '_Struct_' + (speciesIndex + 1) + '" id="ProtList' + speciesIndex +  '" />');
		// $.each(rvds.SpeciesEntry.Molecule_Names_rProtein, function (index, value) {
			// $('#ProtList' + speciesIndex).append(new Option(rvds.SpeciesEntry.Molecule_Names_rProtein[index], rvds.SpeciesEntry.internal_protein_names[index]));
		// });
		// rvds.SpeciesEntry["SubunitProtChains"] = rvds.SpeciesEntry.SubunitProtChains;
	// })
	// $("#ProtList").multiselect("refresh");
	// return title;
}

function populateDomainHelixMenu() {
	$.each(rvDataSets, function(speciesIndex,rvds){
		// Do the Domains
		var DomainList_AN = new Array;
		var DomainList_RN = new Array;
		var DomainSelections = new Array;
		
		for (var i = 0; i < rvds.Residues.length; i++) {
			DomainList_AN[i] = rvds.Residues[i].Domain_AN;
			DomainList_RN[i] = rvds.Residues[i].Domain_RN;
			if (!DomainSelections[DomainList_AN[i]]) {
				DomainSelections[DomainList_AN[i]] = new Array;
			}
			DomainSelections[DomainList_AN[i]].push(rvds.Residues[i].uResName);

		}
		
		var DomainList_ANU = $.grep(DomainList_AN, function (v, k) {
			return $.inArray(v, DomainList_AN) === k;
		});
		var DomainList_RNU = $.grep(DomainList_RN, function (v, k) {
			return $.inArray(v, DomainList_RN) === k;
		});
		
		$('#selectByDomainHelix').append('<optgroup label="Domains' + '_Struct_' + (speciesIndex + 1) + '" id="domainsList' + speciesIndex +  '" />');
		$.each(DomainList_ANU, function (i, val) {
			if (DomainList_RNU[i] && DomainList_RNU[i].indexOf("S") >= 0) {
				$('#domainsList' + speciesIndex).append(new Option(DomainList_RNU[i], DomainSelections[val]));
			} else {
				$('#domainsList' + speciesIndex).append(new Option("Domain " + DomainList_RNU[i], DomainSelections[val]));
			}
			//console.log(DomainSelections[val]);
		});
		
		// Do the Helices
		var HelixList = new Array;
		var HelixSelections = new Array;
		
		for (var i = 0; i < rvds.Residues.length; i++) {
			HelixList[i] = rvds.Residues[i].Helix_Num;
			if (!HelixSelections[HelixList[i]]) {
				HelixSelections[HelixList[i]] = new Array;
			}
			
			HelixSelections[HelixList[i]].push(rvds.Residues[i].uResName);
		}
		
		var HelixListU = $.grep(HelixList, function (v, k) {
				return $.inArray(v, HelixList) === k;
			});
		
		$('#selectByDomainHelix').append('<optgroup label="Helicies' + '_Struct_' + (speciesIndex + 1) + '" id="heliciesList' + speciesIndex +  '" />');
		
		$.each(HelixListU, function (i, val) {
			$('#heliciesList' + speciesIndex).append(new Option("Helix " + val, HelixSelections[val]));
		});
	})
	// Refresh Menu
	$("#selectByDomainHelix").multiselect("refresh");
}


///////////////////////////////////////////////////////////////////////////////


////////////////////////////////// Canvas Functions ///////////////////////////
function watermark(usetime) {
	if (rvDataSets[0].SpeciesEntry.MapType && rvDataSets[0].SpeciesEntry.MapType != "None") {
		var output='';
		$.each(rvDataSets, function (index,rvds){
			var h = (rvds.SpeciesEntry.Orientation == "portrait") ? 779 : 612;
			var d = new Date();
			var df;
			
			var x = 10 + rvds.PageOffset[0];
			var y = h - 20;
			var Message = "A " + rvds.SpeciesEntry.MapType + " secondary structure, generated by RiboVision.";
			var df = d.toLocaleString().indexOf("G");
			var LabelLayers = rvds.getLayerByType("labels");
			var targetLayer = LabelLayers[0];
			targetLayer.CanvasContext.textAlign = 'left';
			targetLayer.CanvasContext.font = '10pt "Myriad Pro", Calibri, Arial';
			targetLayer.CanvasContext.fillStyle = "#FF5500";
			targetLayer.CanvasContext.fillText(Message, x, y);
			
			output += '<g id="g_WaterMark_' + (index + 1) + '">\n';
			output += '<text id="WaterMark" transform="matrix(1 0 0 1 ' + (parseFloat(x) - 1.262).toFixed(3) + ' ' + (parseFloat(y) + 1.145).toFixed(3) + ')" fill="#f6a828" font-family="Myriad Pro" font-size="10">' + Message + '</text>\n';
			if (usetime) {
				targetLayer.CanvasContext.fillText("Saved on " + d.toLocaleString().slice(0, df - 1), x + 75, y + 15);
				output += '<text id="Date" transform="matrix(1 0 0 1 ' + (parseFloat(x) - 1.262 + 75).toFixed(3) + ' ' + (parseFloat(y) + 15 + 1.145).toFixed(3) + ')" fill="#f6a828" font-family="Myriad Pro" font-size="8">' + "Saved on " + d.toLocaleString().slice(0, df - 1) + '</text>\n';
			}
			output += '</g>\n';
		});
		return output;
	}
}

function canvas_arrow(fromx, fromy, tox, toy) {
	
	//rvDataSets[0].Layers[rvDataSets[0].LastLayer].CanvasContext = ResidueLayer.getContext("2d");
	
	var headlen = 10; // length of head in pixels
	var angle = Math.atan2(toy - fromy, tox - fromx);
	rvDataSets[0].Layers[0].CanvasContext.moveTo(fromx, fromy);
	rvDataSets[0].Layers[0].CanvasContext.lineTo(tox, toy);
	rvDataSets[0].Layers[0].CanvasContext.lineTo(tox - headlen * Math.cos(angle - Math.PI / 6), toy - headlen * Math.sin(angle - Math.PI / 6));
	rvDataSets[0].Layers[0].CanvasContext.moveTo(tox, toy);
	rvDataSets[0].Layers[0].CanvasContext.lineTo(tox - headlen * Math.cos(angle + Math.PI / 6), toy - headlen * Math.sin(angle + Math.PI / 6));
}

function welcomeScreen() {
	var image_width=1.0*parseFloat($("#canvasDiv").css('width'));
	var image_height= image_width * 550/733;
	
	var scale_factor = parseFloat($("#canvasDiv").css('width')) / 733;
	// New Welcome Screen
	var img = new Image();
	img.onload = function() {
		if (canvas2DSupported) {
			rvDataSets[0].Layers[0].clearCanvas();
			rvDataSets[0].Layers[0].CanvasContext.drawImage(img, -1*(rvDataSets[0].HighlightLayer.Canvas.width - 612)/2,-1*(rvDataSets[0].HighlightLayer.Canvas.height - 792-242)/2,image_width,image_height);
		}
	}
	img.src = "/static/ribovision/images/RiboVisionLogoHig.png"; //

}
///////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////// Math //////////////////////////////////
function d2h(d) {
	return d.toString(16);
};
function h2d(h) {
	return parseInt(h, 16);
};
function rgb2hex(rgb) {
    rgb = rgb.match(/^rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*(\d+))?\)$/);
    function hex(x) {
        return ("0" + parseInt(x).toString(16)).slice(-2);
    }
    return "#" + hex(rgb[1]) + hex(rgb[2]) + hex(rgb[3]);
}
///////////////////////////////////////////////////////////////////////////////

function changeLineOpacity(opacity){
	document.getElementById('lineOpacity').innerHTML = "Line Opacity: " + Math.round(opacity * 100) + "%";
	$.each(rvDataSets[0].BasePairs, function (ind, item) {
		item.opacity = opacity;
	});
	rvDataSets[0].drawBasePairs("lines");
}
////////////////Nav Line ///////

function drawNavLine(){
		if (typeof d3 === 'undefined'){return;};
		if($('input[name="nl"][value=off]').is(':checked')){
			return;
		}
		$('#NavLineDiv').empty(); //clean div before draw new graph
		var data = [];
		var selectedData=[];
		var selectedDataX=[];
		var selectedDataY=[];
		var maxdata = undefined; 
		var mindata = undefined; 
		
		var targetLayer=rvDataSets[0].getSelectedLayer();
		if (targetLayer===false){
			return;
		}
		var linename = targetLayer.DataLabel;
		var	w = 1.00 * $('#NavLineDiv').innerWidth();
		var h = 0.95 * $('#NavLineDiv').innerHeight();
		var	MarginXL = 60;
		var MarginXR = 120;
		var MarginYT = 40;
		var MarginYB = 40;
		
		var maxdata2 = d3.max($.map(targetLayer.Data, function(d) { return parseFloat(d); }));
		var mindata2 = d3.min($.map(targetLayer.Data, function(d) { return parseFloat(d); }));
		
		if (maxdata2 !== undefined){
			maxdata = maxdata2;
		} 
		if (mindata2 !== undefined){
			mindata = mindata2;
		} 
		if (targetLayer.Type == "selected"){
			maxdata=1;
			mindata=0;
		}
		if (targetLayer.DataLabel === "Protein Contacts"){
			maxdata=1;
			mindata=0;
		}
		var	xScale = d3.scale.linear().domain([0, rvDataSets[0].Residues.length]).range([0 + MarginXL, w - MarginXR]);
		var	yScale = d3.scale.linear().domain([mindata, maxdata]).range([h - MarginYB,0 + MarginYT ]);
		var NavLine = d3.select("#NavLineDiv")
			.append("svg:svg")
			.attr("width", w)
			.attr("height", h)

		var g = NavLine.append("svg:g")
			.attr("width", w)
			.attr("height", h);
			//.attr("transform", "translate(0, " + 200+")");
			
			
		var line = d3.svg.line()
			.defined(function(d) { 
				return (((d!==undefined) && d!=="") ? !isNaN(d) : false) 
			})			
			.x(function(d,i) { return xScale(i); })
			.y(function(d) { return yScale(d); });	
		
		var GraphData = [];
		if (targetLayer.Type === "selected"){
			$.each(targetLayer.Data, function (index,value){
				if (value === false){
					GraphData[index]=0;
				} else {
					GraphData[index]=1;
				}
			});
			linename = "Selected Residues";
		} else if ((targetLayer.DataLabel === "empty data") || (targetLayer.DataLabel === "None")){
			$.each(targetLayer.Data, function (index,value){
				GraphData[index]=0;
			});
		} else if (targetLayer.DataLabel === "Protein Contacts"){
			$.each(targetLayer.Data, function (index,value){
					if (value === " "){
						GraphData[index]=0;
					} else {
						GraphData[index]=1;
					}
				});
			linename = "Protein Contacts";
		} else {
			GraphData = targetLayer.Data;
		}
		
		g.append("svg:path").attr("d", line(GraphData)).style("stroke", targetLayer.Color);
		//Axes
		var xAxis = d3.svg.axis()
			  .scale(xScale)
			  .orient("bottom")
			  .ticks(20);  //Set rough # of ticks
			  
		NavLine.append("g")
			.attr("class", "axis")  //Assign "axis" class
			.attr("transform", "translate(0," + (h - MarginYB) + ")")
			.call(xAxis);
			
		var yAxis = d3.svg.axis()
			  .scale(yScale)
			  .orient("left")
			  .ticks(5);
		
		NavLine.append("g")
			.attr("class", "axis")
			.attr("transform", "translate(" + MarginXL + ",0)")
			.call(yAxis);
		
		//XLabel			
		g.append("text")
		  .attr("x", (w - MarginXR-MarginXL)/2 + MarginXL)
		  .attr("y", h-MarginYB/4)
		  .attr("text-anchor", "middle")
		  .text("Nucleotide Number");	
		  
		//add legend to the navline 
		 g.append("text")
		  .attr("x", MarginXL/4)
		  .attr("y", h/2)
		  .attr("text-anchor", "middle")
		  .attr("transform", "rotate(-90 " + "," + MarginXL/4 + "," + h/2 + ")")
		  .text(linename);	
		
}

function addPopUpWindowResidue(Sele){
	
	// if (rvDataSets[Sele[1]].Residues[Sele[0]].resNum.indexOf(":") >= 0 ){
		// var ResName = rvDataSets[Sele[1]].Residues[Sele[0]].resNum;
	// } else {
		// var ResName = rvDataSets[Sele[1]].SpeciesEntry.Molecule_Names[rvDataSets[Sele[1]].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[Sele[1]].Residues[Sele[0]].ChainID)] +
		// ":" + rvDataSets[Sele[1]].Residues[Sele[0]].resNum;
	// }
	
	var ResName = rvDataSets[Sele[1]].Residues[Sele[0]].uResName;
	var dobj = $.grep(rvDataSets[Sele[1]].ConservationTable, function(e){ return e.resNum == ResName; })[0];
	//var dobj = rvDataSets[Sele[1]].ConservationTable[Sele[0]];
	if (dobj){
		//round the number to two decimal places
		var Hnum = dobj.Shannon * 1;
		var Hn = Hnum.toFixed(2);	
		drawConGraph(dobj);
	} else {
		var ConsensusSymbol = "n/a";
		var Hn = "n/a";
	}
		
	var targetLayer=rvDataSets[Sele[1]].getSelectedLayer();
	
	$('#resName').html(rvDataSets[Sele[1]].Residues[Sele[0]].resName + rvDataSets[Sele[1]].Residues[Sele[0]].resNum +
		" (" + rvDataSets[Sele[1]].Residues[Sele[0]].molName + " " + rvDataSets[Sele[1]].Residues[Sele[0]].MoleculeType +")");
	
	$('#conSeqLetter').html("Consensus: " + ConsensusSymbol);
	$('#activeData').html("Selected Data: " + targetLayer.Data[Sele[0]]);
	$("#conPercentage").html("Shannon Entropy: " + Hn);
}
function drawConGraph(dobj){
	//Width and height
	var Xoffset = 40;
	var Yoffset = 20;
	//var barHeight = 120;
	//var barWidth = 200;
		
	//var w = barWidth + Xoffset;
	//var h = barHeight + Yoffset;
	
	var w = 192;
	var h = 128;
	
	//var barPadding = 10;
	var barPaddingPer = 20;
	//var Xoffset = 40;
	//var padding = 30;
	//var Xpadding = 30;
	var barColors = ["green","blue","black","red","orange"];
	
	//round the number to two decimal places
	var Anum = dobj.A*100;
	var An = Anum.toFixed(1);
	var Cnum = dobj.C*100;
	var Cn = Cnum.toFixed(1);
	var Gnum = dobj.G*100;
	var Gn = Gnum.toFixed(1);
	var Unum = dobj.U*100;
	var Un = Unum.toFixed(1); 
	var Gpnum = dobj.Gaps*100;
	var Gpn = Gpnum.toFixed(1);	 
	var dataset = [An,Cn,Gn,Un,Gpn];
	var sLabels = ["A","C","G","U","gaps"];
	var lenDataSet = dataset.length;

	
	if (typeof d3 === 'undefined'){return;};
	//Remove old SVG
	d3.select("#ResidueTipContent svg").remove();
	//Create SVG element
	var svg = d3.select("#ResidueTipContent")
		.append("svg")
		.attr("width", w)
		.attr("height", h);
	
	var barWidth = (w - Xoffset) / (lenDataSet + ( barPaddingPer/100 * (lenDataSet -1) ) );
	//Scales
	//var xScale = d3.scale.linear();				
	
	var xScale = d3.scale.linear()
							 .domain([0, lenDataSet - 1])
							 .range([Xoffset, w - barWidth]);
							
	var yScale = d3.scale.linear()
						.domain([0, 100])
						.range([h - 2*Yoffset,0]);
		
	svg.selectAll("rect")
	   .data(dataset)
	   .enter()
	   .append("rect")
	   .attr("x", function(d, i) {
			return xScale(i);
	  })
	   .attr("y", function(d) {
			return Yoffset + yScale(d);
	   })
	   .attr("width", barWidth)
	   .attr("height", function(d) {
			return h - 2*Yoffset - yScale(d);
	   })
	   .attr("fill", function(d, i) {
		 return barColors[i];
		});

	svg.selectAll("text.number")
	   .data(dataset)		 
	   .enter()
	   .append("text")
	   .text(function(d) {
			return d;
	   })
	   .attr("text-anchor", "middle")
	   .attr("x", function(d, i) {
			return xScale(i) + barWidth/2;
	   })
	   .attr("y", function(d) {
			return Yoffset + yScale(d) - 5;
	   })
	   .attr("font-family", "sans-serif")
	   .attr("font-size", "10px")
	   .attr("fill", "black")
	   .attr('class','number')
	   .text(String);
	   
	   //Define X axis
			
		var xAxis = d3.svg.axis()
						  .scale(xScale)
						  .orient("bottom")
						  .ticks(5);

	   //Define Y axis
		var yAxis = d3.svg.axis()
				  .scale(yScale)
				  .orient("left")
				  .ticks(5);
		xAxis.tickFormat(function(d,i){
				return sLabels[i];
		});
		
		//Create X axis
		
		svg.append("g")
			.attr("class", "axis")
			.attr("transform", "translate(" + barWidth/2 + "," + (h - Yoffset) + ")")
			.call(xAxis);
		//Create Y axis
		svg.append("g")
			.attr("class", "axis")
			.attr("transform", "translate(" + 0.8* Xoffset + "," + Yoffset + ")")
			.call(yAxis);	
}	

function addPopUpWindowLine(SeleLine){
	
	var j = ActiveBasePairSet[SeleLine].residue_i;
	var k = ActiveBasePairSet[SeleLine].residue_j;
	
	var residue_j = MainResidueMap[j];
	var residue_k = MainResidueMap[k];
	
	var targetLayer = rvDataSets[0].getLayerByType("lines");
	
	/* if (residue_j.resNum.indexOf(":") >= 0 ){
		var ResName1 = residue_j.resNum;
	} else {
		// var ResName1 = rvDataSets[SeleLine[1]].SpeciesEntry.Molecule_Names[rvDataSets[SeleLine[1]].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[SeleLine[1]].Residues[j].ChainID)] +
		// ":" + rvDataSets[SeleLine[1]].Residues[j].resNum;
	}
	if (residue_k.resNum.indexOf(":") >= 0 ){
		var ResName2 = residue_j.resNum;
	} else {
		//var ResName2 = rvDataSets[SeleLine[1]].SpeciesEntry.Molecule_Names[rvDataSets[SeleLine[1]].SpeciesEntry.RNA_Chains.indexOf(rvDataSets[SeleLine[1]].Residues[k].ChainID)] +
		//":" + rvDataSets[SeleLine[1]].Residues[k].resNum;
	} */
	$('#BasePairType').html("Interaction Type: " +  targetLayer[0].DataLabel);
	if (ActiveBasePairSet[SeleLine].ProteinName){
		$('#BasePairSubType').html("Interaction Subtype: " + ActiveBasePairSet[SeleLine].ProteinName + " " + ActiveBasePairSet[SeleLine].bp_type);
	} else {
		$('#BasePairSubType').html("Interaction Subtype: " + ActiveBasePairSet[SeleLine].bp_type);
	}
	
	addPopUpWindowResidue([residue_j.index,residue_j.rvds_index]);
	$("#iResidueTipA").html($("#residuetip").find("#ResidueTipContent").html());
	addPopUpWindowResidue([residue_k.index,residue_k.rvds_index]);
	$("#iResidueTipB").html($("#residuetip").find("#ResidueTipContent").html());
}
//////////End of navline functions////

function UpdateLocalStorage(SaveStateFileName){
	var rvSaveState = {};
	if (localStorageAvailable){
		if($("input[name='LayersCheck']").attr("checked")){
			//localStorage.setItem("rvLayers",JSON.stringify(rvDataSets[0].Layers));
			rvSaveState["rvLayers"] = JSON.stringify(rvDataSets[0].Layers);
		}
		if($("input[name='SelectionsCheck']").attr("checked")){
			//localStorage.setItem("rvSelections",JSON.stringify(rvDataSets[0].Selections));
			rvSaveState["rvSelections"] = JSON.stringify(rvDataSets[0].Selections);
		}
		if($("input[name='LastSpeciesCheck']").attr("checked")){
			//localStorage.setItem("rvLastSpecies",rvDataSets[0].Name);
			rvSaveState["rvLastSpecies"] = rvDataSets[0].Name;
		}
		if($("input[name='PanelSizesCheck']").attr("checked")){
			var po = {
				PanelDivide : PanelDivide,
				TopDivide : TopDivide
			}
			//localStorage.setItem("rvPanelSizes",JSON.stringify(po));
			rvSaveState["rvPanelSizes"] = JSON.stringify(po);
		}
		if($("input[name='MouseModeCheck']").attr("checked")){
			//localStorage.setItem("rvMouseMode",onebuttonmode);	
			rvSaveState["rvMouseMode"] = onebuttonmode;
		}
		if($("input[name='CanvasOrientationCheck']").attr("checked")){
			//localStorage.setItem("rvView",JSON.stringify(rvViews[0]));
			rvSaveState["rvView"] = JSON.stringify(rvViews[0]);			
		}
		if($("input[name='JmolOrientationCheck']").attr("checked")){
			if($('input[name="3dp"][value=on]').is(':checked')){
				rvSaveState["rvJmolOrientation"] = Jmol.evaluateVar(myJmol,"script('show orientation')");	
			}			
		}
		localStorage.setItem(SaveStateFileName,JSON.stringify(rvSaveState));
		var RSL = localStorage["RV_Session_List"];
		if(RSL){
			var RSLa = JSON.parse(RSL);
			RSLa.push(SaveStateFileName);
			var RSLaU;
			RSLaU = $.grep(RSLa, function (v, k) {
				return $.inArray(v, RSLa) === k;
			});
			localStorage.setItem("RV_Session_List",JSON.stringify(RSLaU));
		} else {
			localStorage.setItem("RV_Session_List",JSON.stringify([SaveStateFileName]));
		}
		$("#SessionList").text(localStorage["RV_Session_List"].replace(/[\[\]"]/g,"").replace(/,/g,", "));
	}
}
function RestoreLocalStorage(SaveStateFileName) { 
	var DoneLoading = $.Deferred();
	var DoneLoading2 = $.Deferred();
	
	if (localStorageAvailable){
		rvSaveState = JSON.parse(localStorage[SaveStateFileName]);
		if($("input[name='LastSpeciesCheck']").attr("checked")){
			loadSpecies(rvSaveState.rvLastSpecies,[],DoneLoading,DoneLoading2);
		} else {
			DoneLoading.resolve();
		}
	}
	DoneLoading.done(function() {
		processRvState(rvSaveState);
		//RestoreLocalStorage2();
	});
	DoneLoading2.done(function() {
		$.each(rvDataSets[0].Selections, function (key, value){
			updateSelectionDiv(value.Name,0);
		});
		rvDataSets[0].drawSelection("selected");
		updateModel();
		update3Dcolors();
		if($("input[name='JmolOrientationCheck']").attr("checked")){
			if($('input[name="3dp"][value=on]').is(':checked')){
				var a = rvSaveState.rvJmolOrientation.match(/reset[^\n]+/);
				Jmol.script(myJmol, a[0]);
			}
			
		}
	});
}
function RestoreLocalStorage2(rvSaveState) {
	processRvState(rvSaveState);
}

function rvSaveManager(rvAction,rvLocation) {
	var SaveStateFileName = $("#SaveStateFileName").val();
	
	switch (rvAction) {
		case "Save":
			switch (rvLocation) {
				case "LocalStorage":
					//alert(SaveStateFileName);
					UpdateLocalStorage(SaveStateFileName);
					break;
				case "File":
					saveRvState(SaveStateFileName);
					break;		
				case "Server":
					storeRvState(SaveStateFileName);
					break;
				default:
					alert("shouldn't happen right now");
			}
			break;
		case "Restore":
			switch (rvLocation) {
				case "LocalStorage":
					//alert(SaveStateFileName);
					RestoreLocalStorage(SaveStateFileName);
					break;
				case "File":
					$("#dialog-restore-state").dialog("open");
					break;		
				case "Server":
					retrieveRvState(SaveStateFileName);
					break;
				default:
					alert("shouldn't happen right now");
			}
			break;
		default: 
			alert("shouldn't happen right now");
		}
}

function processRvState(rvSaveState) {
	if($("input[name='LayersCheck']").attr("checked")){
		var data = JSON.parse(rvSaveState.rvLayers);
		rvDataSets[0].Layers=[];
		$.each(data, function (index, value) {
			rvDataSets[0].Layers[index] = rvDataSets[0].HighlightLayer.fromJSON(value);
		});
		
		$.each(rvDataSets[0].Layers,function (index, value) {
			value.updateZIndex(index);
		});
		// Restore Selected and Linked
		var selectedLayer = rvDataSets[0].getSelectedLayer();
		var linkedLayer = rvDataSets[0].getLinkedLayer();
		resizeElements(true);
		$(".oneLayerGroup").remove();
		//remove minilayers
		$(".miniLayerName").remove();
		// Put in Layers
		$.each(rvDataSets[0].Layers, function (key, value){
			LayerMenu(value, key);
		});
		RefreshLayerMenu();
		
		$(".oneLayerGroup" + "[name=" + selectedLayer.LayerName + "]").find(".selectLayerRadioBtn").attr("checked","checked");
		rvDataSets[0].selectLayer(selectedLayer.LayerName);
		$(".oneLayerGroup" + "[name=" + linkedLayer.LayerName + "]").find(".mappingRadioBtn").attr("checked","checked");
		rvDataSets[0].linkLayer(linkedLayer.LayerName);
		
		//Refresh Linked MiniLayer
		var linkedLayer = rvDataSets[0].getLinkedLayer();
		$("#LinkSection").find(".miniLayerName").remove();
		$("#LinkSection").append($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
		.text(linkedLayer.LayerName).attr('name',linkedLayer.LayerName).droppable({
			drop: function (event,ui) {
				ProcessBubbleDrop(event,ui);
			}
		}));
		var	targetLayer=rvDataSets[0].getLayerByType("lines");
		rvDataSets[0].BasePairs=targetLayer[0].Data;
	}
	if($("input[name='SelectionsCheck']").attr("checked")){
		rvDataSets[0].Selections = JSON.parse(rvSaveState.rvSelections);
		$(".oneSelectionGroup").remove();
		// Put in Selections
		$.each(rvDataSets[0].Selections, function (key, value){
			SelectionMenu(value, key);
		});
		//Default check first selection. Come back to these to restore saved state
		$("#SelectionPanel div").first().next().find(".selectSelectionRadioBtn").attr("checked", "checked");
		RefreshSelectionMenu();
		var ret = false;
		$.each(rvDataSets[0].Selections, function (key, value) {
			if (value.Selected) {
				ret = value.Name;
				return false;
			}
		});
		$(".oneSelectionGroup[name='" + ret + "']").find(".selectSelectionRadioBtn").trigger("click");
	}
	if($("input[name='PanelSizesCheck']").attr("checked")){
		var po = JSON.parse(rvSaveState.rvPanelSizes);
		PanelDivide = po.PanelDivide;
		TopDivide = po.TopDivide;
		$( "#canvasPorportionSlider" ).slider("value",PanelDivide);
		$( "#topPorportionSlider" ).slider("value",TopDivide);					
	}
	if($("input[name='MouseModeCheck']").attr("checked")){
		$("#buttonmode").find("input[value='" + rvSaveState.rvMouseMode + "']").trigger("click");
	}
	if($("input[name='CanvasOrientationCheck']").attr("checked")){
		//localStorage.setItem("rvView",7);
		rvViews[0] = rvViews[0].fromJSON(rvSaveState.rvView);
		rvViews[0].restore();
	}
	if($("input[name='JmolOrientationCheck']").attr("checked")){
		//localStorage.setItem("rvJmolOrientation",8);
		if($('input[name="3dp"][value=on]').is(':checked')){
			var a = rvSaveState.rvJmolOrientation.match(/reset[^\n]+/);
			Jmol.script(myJmol, a[0]);
		}
	}
	
	$.each(rvDataSets, function (index, value) {
		value.refreshCanvases();
	});
	
	if(!$("input[name='LastSpeciesCheck']").attr("checked")){
		updateModel();
		update3Dcolors();
	}
	
	//InitRibovision2(true);
}
