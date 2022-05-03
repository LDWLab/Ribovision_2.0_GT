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


/////////////////////////// Classes ///////////////////////////////////////////
function RvLayer(rvds,LayerName, CanvasName, Data, Filled, ScaleFactor, Type, Color,OutLineMode) {
	//Properties
	this.LayerName = LayerName;
	this.SetNumber = rvds.SetNumber;
	this.CanvasName = CanvasName + "_" + this.SetNumber;
	if (document.getElementById(this.CanvasName) == null){
		if (this.Type === "selected"){
			$("#canvasDiv").append($('<canvas id="' + this.CanvasName + '" style="z-index:' + ( rvds.LastLayer + 1 + 800) + ';"></canvas>')); 
		} else {			
			$("#canvasDiv").append($('<canvas id="' + this.CanvasName + '" style="z-index:' + ( rvds.LastLayer + 1 ) + ';"></canvas>')); 
		}
	}
	this.Canvas = document.getElementById(this.CanvasName);
	if (this.Canvas.getContext){
		this.CanvasContext = this.Canvas.getContext("2d");
	}
	this.Data = Data;
	this.dataLayerColors = [];
	this.Filled = Filled;
	this.ScaleFactor = ScaleFactor;
	this.LinearGradients = [];
	this.Type = Type;
	this.zIndex = rvds.LastLayer + 1;
	this.Visible = true;
	this.Selected = false;
	this.Linked = false;
	this.ColorLayer = [];
	this.ColorGradientMode = "Matched";
	this.DataLabel = "None";
	if (Color) {
		this.Color = Color;
	} else {
		this.Color = "#0000FF";
	}
	if (OutLineMode==undefined) {
		this.OutLineMode=true;
	} else {
		this.OutLineMode=OutLineMode;
	}
	if (this.Type === "lines") {
		this.ColorLayer = "gray_lines";
		this.ColorGradientMode = "Matched";
	}
	
	//Methods
	this.updateZIndex = function(zIndex){
		if (this.Type === "selected"){
			this.Canvas.style.zIndex = zIndex + 800;
		} else {
			this.Canvas.style.zIndex = zIndex;
		}
		this.zIndex = zIndex;
	}
	this.toJSON = function () {
		return {
			LayerName: this.LayerName,
			CanvasName: this.CanvasName, 
			Data: this.Data, 
			dataLayerColors: this.dataLayerColors,
			Filled: this.Filled,
			ScaleFactor: this.ScaleFactor,
			LinearGradients: this.LinearGradients,
			Type: this.Type,
			zIndex: this.zIndex,
			Visible: this.Visible,
			Selected: this.Selected,
			Linked: this.Linked,
			ColorLayer: this.ColorLayer,
			ColorGradientMode: this.ColorGradientMode,
			DataLabel: this.DataLabel,
			Color: this.Color,
			ColorLayer: this.ColorLayer,
			ColorGradientMode: this.ColorGradientMode
		};
	};
	
	this.fromJSON = function (json) {
		//var data = JSON.parse(json);
		var data = json;
		var e = new RvLayer(data.LayerName, data.CanvasName, data.Data, data.Filled, data.ScaleFactor, data.Type, data.Color);
		e.dataLayerColors = data.dataLayerColors;
		e.LinearGradients = data.LinearGradients;
		e.zIndex = data.zIndex;
		e.Visible = data.Visible;
		e.Selected = data.Selected;
		e.Linked = data.Linked;
		e.ColorLayer = data.ColorLayer;
		e.ColorGradientMode = data.ColorGradientMode;
		e.DataLabel = data.DataLabel;
		e.ColorLayer = data.ColorLayer;
		e.ColorGradientMode = data.ColorGradientMode;
		return e;
	};
	this.clearCanvas = function () {
		if (canvas2DSupported){
			this.CanvasContext.setTransform(1, 0, 0, 1, 0, 0);
			this.CanvasContext.clearRect(0, 0, rvViews[0].width, rvViews[0].height);
			this.CanvasContext.setTransform(rvViews[0].scale, 0, 0, rvViews[0].scale, rvViews[0].x, rvViews[0].y);
		}
	};
	this.addLinearGradient = function (LinearGradient) {
		this.LinearGradients.push(LinearGradient);
	};
	this.deleteLayer = function () {
		$(this.Canvas).remove();
	};
	this.clearData = function () {
		this.dataLayerColors = new Array;
		this.Data = new Array;
		for (var jj = 0; jj < rvds.Residues.length; jj++) {
			this.dataLayerColors[jj] = undefined;
			this.Data[jj] = undefined;
		}
	}
	this.setVisibility = function (visibility_prop) {
		switch (visibility_prop) {
		case "hidden":
			$(this.Canvas).css("visibility", "hidden");
			this.Visible = false;
			break;
		case "visible":
			$(this.Canvas).css("visibility", "visible");
			this.Visible = true;
			break;
		default:
			alert("this shouldn't happen");
		}
	}
	this.clearAll = function (){
		switch (this.Type){
			case "circles":
				this.DataLabel = "None";
				$("[name=" + this.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(this.DataLabel);
				this.clearData();
				drawNavLine();
				rvds.clearCanvas(this.LayerName);
				break;
			case "contour":
				this.DataLabel = "None";
				$("[name=" + this.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(this.DataLabel);
				this.clearData();
				drawNavLine();
				rvds.clearCanvas(this.LayerName);
				break;	
			case "residues":
				this.DataLabel = "None";
				$("[name=" + this.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(this.DataLabel);
				this.clearData();
				clearColor(false);
				drawNavLine();
				break;
			case "lines":
				this.DataLabel = "None";
				$("[name=" + this.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(this.DataLabel);
				//$(this).find(".layerContent").find("span[name=DataLabel]").text("None"));
				//$(this).parent().parent().find(".DataDescription").text("Empty Data");
				drawNavLine();
				refreshBasePairs("clear_lines");
				break;
			default:
		}			
		//this.DataLabel = "None";
	}
}

function RvSelection(rvds,SelectionName,rvResidues,rvColor,rvResidues_rProtein){
	//Properties
	this.Name = SelectionName;
	this.zIndex = rvds.Selections.length;
	if (rvResidues) { 
		this.Residues = rvResidues;
	} else {
		this.Residues = [];
	}
	if (rvResidues_rProtein) { 
		this.Residues_rProtein = rvResidues_rProtein;
	} else {
		this.Residues_rProtein = [];
	}
	if (rvColor) {
		this.Color = rvColor;
	} else {
		this.Color = "#940B06";
	}
	this.Selected = false;
	
	//Methods
}
function rvDataSet(DataSetName,SetNumber) {
	//Properties
	this.Name = DataSetName;
	this.SetNumber = SetNumber;
	this.PageOffset = [];
	this.Layers = [];
	this.HighlightLayer = [];
	this.Residues = [];
	this.ResidueList = [];
	this.ContourLinePoints = [];
	this.ExtraContourLineSegments = [];
	this.SequenceList = "";
	this.rvTextLabels = [];
	this.rvLineLabels = [];
	//this.BasePairs = [];
	this.FullBasePairSet = [];
	this.CustomData = [];
	this.SpeciesEntry = [];
	this.Selections = [];
	this.LastLayer = -1;
	this.LayerTypes = ['circles', 'lines', 'labels', 'residues', 'contour', 'selected'];
	this.ConservationTable = [];
	this.DataDescriptions = [];
	this.ExtraPyMOLScript ='';
	this.ColorProteins = [];
	this.Font_Size_SVG = [];
	this.Font_Size_Canvas = [];
	this.Circle_Radius = [];
	//Methods
	this.toJSON = function () {
		return {
			DataSetName: this.Name,
			Layers: this.Layers,
			//this.HighlightLayer = [];
			Residues: this.Residues,
			ResidueList: this.ResidueList,
			rvTextLabels: this.rvTextLabels,
			rvLineLabels: this.rvLineLabels,
			BasePairs: this.BasePairs,
			CustomData: this.CustomData,
			SpeciesEntry: this.SpeciesEntry,
			Selections: this.Selections,
			LastLayer: this.LastLayer,
			LayerTypes: this.LayerTypes,
			ConservationTable: this.ConservationTable,
			DataDescriptions: this.DataDescriptions
		};
	};
	this.fromJSON = function (json) {
		var data = JSON.parse(json);
		var e = new rvDataSet(DataSetName);
		e.Residues = data.Residues;
		e.ResidueList = data.ResidueList;
		e.rvTextLabels = data.rvTextLabels;
		e.rvLineLabels = data.rvLineLabels;
		e.BasePairs = data.BasePairs;
		e.CustomData = data.CustomData;
		e.SpeciesEntry = data.SpeciesEntry;
		e.Selections = data.Selections;
		e.LastLayer = data.LastLayer;
		e.ConservationTable = data.ConservationTable;
		e.DataDescriptions = data.DataDescriptions;
		e.addHighlightLayer("HighlightLayer", "HighlightLayer", [], false, 1.176, 'highlight');
		$.each(data.Layers, function (index, value) {
			e.Layers[index] = e.HighlightLayer.fromJSON(value);
		});
		return e;
	};
	this.addLayers = function (rvLayers) {
		this.Layers = rvLayers;
		this.LastLayer = this.Layers.length - 1;
	};
	this.addLayer = function (LayerName, CanvasName, Data, Filled, ScaleFactor, Type, Color, OutLineMode) {
		var b = new RvLayer(this,LayerName, CanvasName, Data, Filled, ScaleFactor, Type, Color, OutLineMode);
		this.Layers[this.Layers.length] = b;
		this.LastLayer = this.Layers.length - 1;
	};
	this.addHighlightLayer = function (LayerName, CanvasName, Data, Filled, ScaleFactor, Type) {
		var b = new RvLayer(this,LayerName, CanvasName, Data, Filled, ScaleFactor, Type);
		this.HighlightLayer = b;
		this.HighlightLayer.Canvas.style.zIndex = 990 + this.SetNumber;
	};
	this.addResidues = function (rvResidues) {
		this.Residues = rvResidues;
		this.SequenceList = makeSequenceList(rvResidues);
		this.updateRNAchains();
	};
	this.addLabels = function (rvTextLabels, rvLineLabels, rvExtraLabels) {
		if (rvTextLabels !== undefined){
			this.rvTextLabels = rvTextLabels;
		}
		if (rvLineLabels !== undefined){
			this.rvLineLabels = rvLineLabels;
		}
		if (rvExtraLabels !== undefined){
			this.rvExtraLabels = rvExtraLabels;
		}
	};
	/*
	this.addBasePairs = function (BasePairs) {
		this.BasePairs = BasePairs;
	};
	
	this.addSelected = function (Selected) {
		this.Selected = Selected;
	};*/
	this.addCustomData = function (CustomData) {
		var rvds=this;
		if(rvds.SpeciesEntry.Molecule_Names == "custom"){
			// Come back, add real custom support, multiple structures, filtering, etc
			rvds.CustomData = CustomData;
		} else {
			//var molecules_names = rvds.SpeciesEntry.Molecule_Names.concat(rvds.SpeciesEntry.Molecule_Names_rProtein)
			var molecules_names = rvds.SpeciesEntry.RNA_Names
			rvds.CustomData = $.grep(CustomData, function(value,index){
				//var mol_name = value.resNum.split(':')[0];
				return $.inArray(value.resNum.split(':')[0],molecules_names) >=0;
			})
		}
		// Copy extra parameters to customdata to support multiple datasets
		if(rvds.CustomData.length > 0){
			if (rvds.CustomData[0].SwitchPoint != undefined & rvds.CustomData[0].SwitchPoint == ""){
				rvds.CustomData[0].SwitchPoint = CustomData[0].SwitchPoint;
			}
			if (rvds.CustomData[0].TwoColorMode != undefined & rvds.CustomData[0].TwoColorMode == ""){
				rvds.CustomData[0].TwoColorMode = CustomData[0].TwoColorMode;
				rvds.CustomData[1].TwoColorMode = CustomData[1].TwoColorMode;
			}
		}
	};
	this.addSpeciesEntry = function (SpeciesEntry) {
		this.SpeciesEntry = SpeciesEntry;
		// Set FontSize
		this.Font_Size_Canvas = this.SpeciesEntry.Font_Size_Canvas;
		this.Font_Size_SVG = this.SpeciesEntry.Font_Size_SVG;
		this.Circle_Radius = this.SpeciesEntry.Circle_Radius;
	};
	this.updateRNAchains = function (SpeciesEntry) {
		//this.SpeciesEntry.Molecule_Names = this.SpeciesEntry.Molecule_Names.split(";");		
		//this.SpeciesEntry.RNA_Chains = this.SpeciesEntry.RNA_Chains.split(";");
		var molName = [...new Set(this.Residues.map(item => item.molName))];
		this.SpeciesEntry.RNA_Names = molName;
		var ChainName = [...new Set(this.Residues.map(item => item.ChainName))];
		this.SpeciesEntry.RNA_Chains = ChainName;
	
	};
	this.updateProtchains = function (SpeciesEntry) {
		//this.SpeciesEntry.RNA_Chains_rProtein = this.SpeciesEntry.RNA_Chains_rProtein.split(";");
		//this.SpeciesEntry.Molecule_Names_rProtein = this.SpeciesEntry.Molecule_Names_rProtein.split(";");
		//this.SpeciesEntry.internal_protein_names = this.SpeciesEntry.internal_protein_names.split(";");
	};
	this.addSelection = function (Name, rvResidues, rvColor) {
		if (!Name) {
			Name = "Selection_" + (this.Selections.length + 1);
		}
		this.Selections.unshift(new RvSelection(this,Name,rvResidues,rvColor));
	};
	this.sort = function () {
		this.Layers.sort(function (a, b) {
			return (Number(a.zIndex) - Number(b.zIndex));
		});
		$.each(this.Layers, function (key, value) {
			this.updateZIndex(key);
		});
	};
	this.SelectionsSort = function () {
		this.Selections.sort(function (a, b) {
			return (Number(a.zIndex) - Number(b.zIndex));
		});
		$.each(this.Selections, function (key, value) {
			this.zIndex=key;
		});
	};
	this.clearCanvas = function (layer) {
		//var LayerTypes=['circles','lines','labels','residues','contour','selected'];
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					value.clearCanvas();
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					value.clearCanvas();
				}
			});
		}
	};
	this.refreshResiduesExpanded = function (layer) {
		var rvds = this;
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					refreshLayer.call(rvds,value);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					refreshLayer.call(rvds,value);
				}
			});
		}
	};
	this.drawLabels = function (layer,drawExtra) {
		var rvds = this;
		this.clearCanvas(layer);
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawLabels.call(rvds,value,drawExtra);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawLabels.call(rvds,value,drawExtra);
				}
			});
		}
	};
	this.clearData = function (layer) {
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					value.clearData();
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					value.clearData();
				}
			});
		}
	};
	this.drawResidues = function (layer, dataIndices, ColorArray, noClear) {
		//this.clearCanvas(layer);
		var ind = $.inArray(layer, this.LayerTypes);
		var rvds = this;
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawResidues.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawResidues.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		}
	};
	this.drawContourLines = function (layer, dataIndices, ColorArray, noClear) {
		//this.clearCanvas(layer);
		var rvds = this;
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawContourLine.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawContourLine.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		}
	};
	this.drawSelection = function (layer,SeleName) {
		var rvds = this;
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawSelection.call(rvds,value,SeleName);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawSelection.call(rvds,value,SeleName);
				}
			});
		}
	};
	this.drawDataCircles = function (layer, dataIndices, ColorArray, noClear) {
		var rvds = this;
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawDataCircles.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawDataCircles.call(rvds,value, dataIndices, ColorArray, noClear);
				}
			});
		}
	};
	this.drawBasePairs = function (layer, colorLayer) {
		var rvds = this;
		var ind = $.inArray(layer, this.LayerTypes);
		if (ind >= 0) {
			$.each(this.Layers, function (key, value) {
				if (value.Type === layer) {
					drawBasePairs.call(rvds,value, colorLayer);
				}
			});
		} else {
			$.each(this.Layers, function (key, value) {
				if (value.LayerName === layer) {
					drawBasePairs.call(rvds,value, colorLayer);
				}
			});
		}
	};
	this.getLayer = function (layer) {
		var ret = false;
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				ret = value;
				return false;
			}
		});
		return ret;
	};
	this.getLayerByType = function (layer) {
		if (!$.isArray(layer)){
			layer=[layer];
		}
		var ret = [];
		$.each(this.Layers, function (key, value) {
			if ($.inArray(value.Type, layer) >=0) {
				ret.push(value);
			}
		});
		if (ret.length >0){
			return ret;
		}
		return false;
	};
	this.getSelectedLayer = function () {
		var ret = false;
		$.each(this.Layers, function (key, value) {
			if (value.Selected) {
				ret = value;
				return false;
			}
		});
		return ret;
	};
	this.getLinkedLayer = function () {
		var ret = false;
		$.each(this.Layers, function (key, value) {
			if (value.Linked) {
				ret = value;
				return false;
			}
		});
		return ret;
	};
	this.linkLayer = function (layer) {
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				value.Linked = true;
			} else {
				value.Linked = false;
			}
		});
	};
	this.selectLayer = function (layer) {
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				value.Selected = true;
			} else {
				value.Selected = false;
			}
		});
	};
	this.isUniqueLayer = function (layer) {
		var ret = true;
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				ret = false;
				return false;
			}
		});
		return ret;
	}
	this.deleteLayer = function (layer) {
		var rvds = this;
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				value.deleteLayer();
				rvds.Layers.splice(key, 1);
				return false;
			}
		});
		$.each(this.Layers, function (key, value) {
			if (value.LayerName === layer) {
				value.updateZIndex(key);
			}
		});
		this.LastLayer = this.Layers.length - 1;
	};
	this.getSelection = function (rvSelectionName) {
		var ret = false;
		$.each(this.Selections, function (key, value) {
			if (value.Name === rvSelectionName) {
				ret = value;
				return false;
			}
		});
		return ret;
	};
	this.deleteSelection = function (selection) {
		var rvds = this;
		$.each(this.Selections, function (key, value) {
			if (value.Name === selection) {
				rvds.Selections.splice(key, 1);
				return false;
			}
		});
	};
	this.isUniqueSelection = function (selection) {
		var ret = true;
		$.each(this.Selections, function (key, value) {
			if (value.Name === selection) {
				ret = false;
				return false;
			}
		});
		return ret;
	}
	
	this.refreshCanvases = function () {
		this.drawResidues("residues");
		this.drawSelection("selected");
		this.refreshResiduesExpanded("circles");
		this.refreshResiduesExpanded("contour");
		this.drawLabels("labels");
		this.drawBasePairs("lines");
	}
	this.makeContourLinePoints = function (){
		var rvds = this;
		var cutoff=11.5;
		var ContourLinePoints = [];
		var ExtraContourLineSegments = [];
		var numResidues=this.Residues.length;
		
		$.each(this.Residues, function( index , value ) {
			//Special case for firstspeciesIndex
			if (index == 0) {
				var diffxy=[];
				diffxy[0]=Number(rvds.Residues[index + 1].X) - Number(rvds.Residues[index].X);
				diffxy[1]=Number(rvds.Residues[index + 1].Y) - Number(rvds.Residues[index].Y);
				var clp = {
					X1 : Number(rvds.Residues[index].X)-diffxy[0]/2 + rvds.PageOffset[0],
					X2 : Number(rvds.Residues[index].X) + rvds.PageOffset[0],
					X3 : (Number(rvds.Residues[index + 1].X) + Number(rvds.Residues[index].X))/2 + rvds.PageOffset[0],
					Y1 : Number(rvds.Residues[index].Y)-diffxy[1]/2,
					Y2 : Number(rvds.Residues[index].Y),
					Y3 : (Number(rvds.Residues[index + 1].Y) + Number(rvds.Residues[index].Y))/2,
				};
			} else if (index == numResidues - 1) {
			//Special case for last
				var diffxy=[];
				diffxy[0]=Number(rvds.Residues[index].X) - Number(rvds.Residues[index - 1].X);
				diffxy[1]=Number(rvds.Residues[index].Y) - Number(rvds.Residues[index - 1].Y);
				var clp = {
					X1 : (Number(rvds.Residues[index - 1].X) + Number(rvds.Residues[index].X))/2 + rvds.PageOffset[0],
					X2 : Number(rvds.Residues[index].X) + rvds.PageOffset[0],
					X3 : Number(rvds.Residues[index].X)-diffxy[0]/2 + rvds.PageOffset[0],
					Y1 : (Number(rvds.Residues[index - 1].Y) + Number(rvds.Residues[index].Y))/2,
					Y2 : Number(rvds.Residues[index].Y),
					Y3 : Number(rvds.Residues[index].Y)-diffxy[1]/2,
				};
			
			} else {
				//Unlikely for the first or last to need special adjustment, so do it only here. 
				//We need to adjust the line lengths if either (distance is larger than a set cutoff) 
				// OR (molecule name changes). 
				
				var clp = {
					X1 : (Number(rvds.Residues[index - 1].X) + Number(rvds.Residues[index].X))/2 + rvds.PageOffset[0],
					X2 : Number(rvds.Residues[index].X) + rvds.PageOffset[0],
					X3 : (Number(rvds.Residues[index + 1].X) + Number(rvds.Residues[index].X))/2 + rvds.PageOffset[0],
					Y1 : (Number(rvds.Residues[index - 1].Y) + Number(rvds.Residues[index].Y))/2,
					Y2 : Number(rvds.Residues[index].Y),
					Y3 : (Number(rvds.Residues[index + 1].Y) + Number(rvds.Residues[index].Y))/2,
				};
				
				// We overwrite clp accordingly
				var firstMatches=rvds.ResidueList[index - 1].match(/[^:]+/)[0] == rvds.ResidueList[index].match(/[^:]+/)[0];
				var secondMatches=rvds.ResidueList[index + 1].match(/[^:]+/)[0] == rvds.ResidueList[index].match(/[^:]+/)[0];

				var dist1 = Math.sqrt((clp.X2-clp.X1)*(clp.X2-clp.X1) + (clp.Y2-clp.Y1)*(clp.Y2-clp.Y1));
				if (dist1 > cutoff || !firstMatches){
					//console.log(dist1 + " " + index);
					clp.X1=clp.X2;
					clp.Y1=clp.Y2;
				}
				if (dist1 > cutoff && firstMatches){
					// Add Extra line segment
					var ecls = {
						X1 : ContourLinePoints[index - 1].X3,
						X2 : clp.X1,
						Y1 : ContourLinePoints[index - 1].Y3,
						Y2 : clp.Y1,
					};
					ExtraContourLineSegments.push(ecls);
				}
				
				var dist2 = Math.sqrt((clp.X3-clp.X2)*(clp.X3-clp.X2) + (clp.Y3-clp.Y2)*(clp.Y3-clp.Y2));
				if (dist2 > cutoff || !secondMatches){
					//console.log(dist2 + " " + index);
					clp.X3=clp.X2;
					clp.Y3=clp.Y2;
				}
			}
			ContourLinePoints[index]=clp;
		});
		this.ExtraContourLineSegments=ExtraContourLineSegments;
		this.ContourLinePoints=ContourLinePoints;
	}
	this.makeResidueList = function () {
		var ResidueListLocal = [],j;
		for (j = 0; j < this.Residues.length; j++) {
			ResidueListLocal[j] = this.Residues[j].uResName;
		}
		this.ResidueList=ResidueListLocal;
	}
	// Private functions, kinda
	
	function makeSequenceList(rvResidues) {
		var SequenceListLocal = "",j;
		for (j = 0; j < rvResidues.length; j++) {
			SequenceListLocal = SequenceListLocal.concat(rvResidues[j].resName);
		}
		return SequenceListLocal;
	}
	function refreshLayer(targetLayer) {
		var rvds = this;
		//var outLineMode=true;// Come back and make this optional
		if(this.Circle_Radius){
			var CircleSize = this.Circle_Radius;
		} else {
			var CircleSize = 1.7;
		}
		if (rvds.Residues !== undefined && targetLayer.Type === "circles") {
			targetLayer.clearCanvas();
			for (var i = rvds.Residues.length - 1; i >= 0; i--) {
				if (targetLayer.dataLayerColors[i] != undefined && targetLayer.dataLayerColors[i] != '#858585') {
					targetLayer.CanvasContext.beginPath();
					targetLayer.CanvasContext.arc(parseFloat(rvds.Residues[i].X) + rvds.PageOffset[0], rvds.Residues[i].Y, (targetLayer.ScaleFactor * CircleSize), 0, 2 * Math.PI, false);
					targetLayer.CanvasContext.closePath();
					targetLayer.CanvasContext.strokeStyle = targetLayer.dataLayerColors[i];
					targetLayer.CanvasContext.stroke();
					if (targetLayer.Filled) {
						targetLayer.CanvasContext.fillStyle = targetLayer.dataLayerColors[i];
						targetLayer.CanvasContext.fill();
					}
				}
			}
		}
		// I do not know where these magic numbers, 0.05 and 0.3 come from. They make the contour lines look correct. 
		if (rvds.Residues !== undefined && targetLayer.Type === "contour") {
			targetLayer.clearCanvas();
			targetLayer.CanvasContext.lineCap = 'round';
			if (rvds.Residues && rvds.Residues.length > 0) {
				if (targetLayer.OutLineMode){
					/* $.each(rvds.ExtraContourLineSegments, function (index, value) {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.lineJoin = "round";  
						targetLayer.CanvasContext.moveTo(value.X1 - .05, value.Y1 - .3);
						targetLayer.CanvasContext.lineTo(value.X2 - .05, value.Y2 - .3);
						targetLayer.CanvasContext.setLineDash([4,4]);
						targetLayer.CanvasContext.strokeStyle = '#000000';	
						targetLayer.CanvasContext.lineWidth = 1.0;					
						targetLayer.CanvasContext.stroke();
						targetLayer.CanvasContext.closePath();
					}); */
					$.each(rvds.Residues, function (index, value) {
						if (targetLayer.dataLayerColors[index] != undefined && targetLayer.dataLayerColors[index] != '#858585') {
							targetLayer.CanvasContext.beginPath();
							targetLayer.CanvasContext.lineJoin = "round";  
							targetLayer.CanvasContext.moveTo(rvds.ContourLinePoints[index].X1 - .05, rvds.ContourLinePoints[index].Y1 - .3);
							targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X2 - .05, rvds.ContourLinePoints[index].Y2 - .3);
							targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X3 - .05, rvds.ContourLinePoints[index].Y3 - .3);
							targetLayer.CanvasContext.setLineDash([]);
							targetLayer.CanvasContext.strokeStyle = '#000000';	
							targetLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.5 * 1.9 * CircleSize;
							// Scale factor is to make lines different styles. 1.5 is to make the outline. 1.9 is to convert from circlesize to contourline thickness. 		
							targetLayer.CanvasContext.stroke();
							targetLayer.CanvasContext.closePath();
						}
					});
					//alert(targetLayer.OutLineMode)
				}
				$.each(rvds.Residues, function (index, value) {
					if (targetLayer.dataLayerColors[index] != undefined && targetLayer.dataLayerColors[index] != '#858585') {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.lineJoin = "round";  
						targetLayer.CanvasContext.moveTo(rvds.ContourLinePoints[index].X1 - .05, rvds.ContourLinePoints[index].Y1 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X2 - .05, rvds.ContourLinePoints[index].Y2 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X3 - .05, rvds.ContourLinePoints[index].Y3 - .3);
						targetLayer.CanvasContext.setLineDash([]);
						targetLayer.CanvasContext.strokeStyle = targetLayer.dataLayerColors[index];	
						targetLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.0 * 1.9 * CircleSize;					
						targetLayer.CanvasContext.stroke();
						targetLayer.CanvasContext.closePath();
					}
				});
				
			}
		}
	}
	function drawLabels(targetLayer,drawExtra) {
		if (!canvas2DSupported){return};
		targetLayer.CanvasContext.textAlign = 'left';
		targetLayer.CanvasContext.lineCap = 'round';

		if (this.rvTextLabels != undefined) {
			var n = watermark(false);
			$.each(this.ExtraContourLineSegments, function (index, value) {
				targetLayer.CanvasContext.beginPath();
				targetLayer.CanvasContext.lineJoin = "round";  
				targetLayer.CanvasContext.moveTo(value.X1 - .05, value.Y1 - .3);
				targetLayer.CanvasContext.lineTo(value.X2 - .05, value.Y2 - .3);
				targetLayer.CanvasContext.setLineDash([4,4]);
				targetLayer.CanvasContext.strokeStyle = '#000000';	
				targetLayer.CanvasContext.lineWidth = 1.0;					
				targetLayer.CanvasContext.stroke();
				targetLayer.CanvasContext.closePath();
			});
	
			for (var i = 0; i < this.rvTextLabels.length; i++) {
				targetLayer.CanvasContext.font = (0.70 * this.rvTextLabels[i].FontSize) + 'pt "Myriad Pro", Calibri, Arial';
				targetLayer.CanvasContext.fillStyle = this.rvTextLabels[i].Fill;
				targetLayer.CanvasContext.fillText(this.rvTextLabels[i].LabelText, parseFloat(this.rvTextLabels[i].X) + this.PageOffset[0], this.rvTextLabels[i].Y);
			}
			
			targetLayer.CanvasContext.setLineDash([]);
			targetLayer.CanvasContext.strokeStyle = "rgba(35,31,32,32)";
			targetLayer.CanvasContext.lineWidth = .5;
			
			for (var i = 0; i < this.rvLineLabels.length; i++) {
				targetLayer.CanvasContext.beginPath();
				targetLayer.CanvasContext.moveTo(parseFloat(this.rvLineLabels[i].X1) + this.PageOffset[0], this.rvLineLabels[i].Y1);
				targetLayer.CanvasContext.lineTo(parseFloat(this.rvLineLabels[i].X2) + this.PageOffset[0], this.rvLineLabels[i].Y2);
				targetLayer.CanvasContext.closePath();
				targetLayer.CanvasContext.stroke();
			}
		}
		/*
		var data = "data:image/svg+xml," + "<svg xmlns='http://www.w3.org/2000/svg' width='612' height='792'>" + "\n";
		for (var j = 0 ;  j < this.rvExtraLabels.length ; j++){
			data +=	this.rvExtraLabels[j].SVGLine;
		}
		data +="</svg>";
		//alert(data);
		var img = new Image();
		//img.src = data;
		img.src = "js/RiboVision/SC_28S_Struct_Dash_Lines_m.svg";
		img.onload = function() { targetLayer.CanvasContext.drawImage(img, 0, 0,612,792); };
		//targetLayer.CanvasContext.drawImage(img, 0, 0);
		*/
	}
	function drawResidues(targetLayer, dataIndices, ColorArray, noClear) {
		if (this.Font_Size_Canvas){
			var FontSize=this.Font_Size_Canvas;
		} else {
			var FontSize=3.1;
		}
		//var resMod = $('input[name="ptmod"][value=on]').is(':checked');
		if (targetLayer.Type === "residues") {
			targetLayer.clearCanvas();
			
			if (!noClear) {
				//targetLayer.clearCanvas();
				//targetLayer.dataLayerColors = [];
			}
			if (this.Residues && this.Residues.length > 0) {
				if (dataIndices && ColorArray) {
					for (var i = 0; i < this.Residues.length; i++) {
						if (ColorArray[dataIndices[i]]) {
							this.Residues[i].color = ColorArray[dataIndices[i]];
							targetLayer.dataLayerColors[i] = ColorArray[dataIndices[i]];
						//} else if (dataIndices[i]!=undefined && isNaN(dataIndices[i])){ Do not remember why this was here. Try without.
						} else if (i in dataIndices){
							this.Residues[i].color = "#000000";
							targetLayer.dataLayerColors[i] = "#000000";
						}
					}
				} else {
					for (var i = 0; i < this.Residues.length; i++) {
						this.Residues[i].color = targetLayer.dataLayerColors[i];
					}
				}
				targetLayer.CanvasContext.strokeStyle = "#000000";
				targetLayer.CanvasContext.textBaseline = "middle";
				targetLayer.CanvasContext.textAlign = "center";
				for (var i = this.Residues.length - 1; i >= 0; i--) {
					targetLayer.CanvasContext.fillStyle = (targetLayer.dataLayerColors[i] || "#000000");
					targetLayer.CanvasContext.font = this.Residues[i]["font-weight"] + " " + FontSize + 'pt "Myriad Pro", Calibri, Arial';
					//if (resMod){
						//targetLayer.CanvasContext.fillText(this.Residues[i].modResName, this.Residues[i].X, this.Residues[i].Y);
					//} else {
					// Not adjusting Y with PageOffset for now, since it won't do anything. 
						targetLayer.CanvasContext.fillText(this.Residues[i].resName, parseFloat(this.Residues[i].X) + this.PageOffset[0], this.Residues[i].Y);
					//}
				}
			} else {
				welcomeScreen();
			}
		}
	}
	function drawContourLine(targetLayer, dataIndices, ColorArray, noClear) {
		var rvds = this;
		//var outLineMode=true;
		if(this.Circle_Radius){
			var CircleSize = this.Circle_Radius;
		} else {
			var CircleSize = 1.7;
		}
		
		if (targetLayer.Type === "contour") {
			//targetLayer.clearCanvas();
			if (!noClear) {
				targetLayer.clearCanvas();
				targetLayer.dataLayerColors = [];
			}
			// I do not know where these magic numbers, 0.05 and 0.3 come from. They make the contour lines look correct. 
			if (this.Residues && this.Residues.length > 0) {
				if (targetLayer.OutLineMode){
					/* $.each(rvds.ExtraContourLineSegments, function (index, value) {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.lineJoin = "round";  
						targetLayer.CanvasContext.moveTo(value.X1 - .05, value.Y1 - .3);
						targetLayer.CanvasContext.lineTo(value.X2 - .05, value.Y2 - .3);
						targetLayer.CanvasContext.setLineDash([4,4]);
						targetLayer.CanvasContext.strokeStyle = '#000000';	
						targetLayer.CanvasContext.lineWidth = 1.0;					
						targetLayer.CanvasContext.stroke();
						targetLayer.CanvasContext.closePath();
					}); */
					$.each(rvds.Residues, function (index, value) {
						if (ColorArray ) {
							targetLayer.CanvasContext.beginPath();
							targetLayer.CanvasContext.lineJoin = "round";  
							targetLayer.CanvasContext.moveTo(rvds.ContourLinePoints[index].X1 - .05, rvds.ContourLinePoints[index].Y1 - .3);
							targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X2 - .05, rvds.ContourLinePoints[index].Y2 - .3);
							targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X3 - .05, rvds.ContourLinePoints[index].Y3 - .3);
							targetLayer.CanvasContext.setLineDash([]);
							targetLayer.CanvasContext.strokeStyle = '#000000';	
							targetLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.5 * 1.9 * CircleSize				
							targetLayer.CanvasContext.stroke();
							targetLayer.CanvasContext.closePath();
						}
					});
				}
				$.each(this.Residues, function (index, value) {
					if (dataIndices && ColorArray && ColorArray[dataIndices[index]] != undefined && ColorArray[dataIndices[index]] != '#858585') {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.lineJoin = "round";  
						targetLayer.CanvasContext.moveTo(rvds.ContourLinePoints[index].X1 - .05, rvds.ContourLinePoints[index].Y1 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X2 - .05, rvds.ContourLinePoints[index].Y2 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X3 - .05, rvds.ContourLinePoints[index].Y3 - .3);
						targetLayer.CanvasContext.setLineDash([]);
						targetLayer.CanvasContext.strokeStyle = ColorArray[dataIndices[index]];	
						targetLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.0 * 1.9 * CircleSize;					
						targetLayer.CanvasContext.stroke();
						targetLayer.CanvasContext.closePath();
						targetLayer.dataLayerColors[index] = ColorArray[dataIndices[index]];
					} else if (!dataIndices && !ColorArray && targetLayer.dataLayerColors[index] && targetLayer.dataLayerColors[index] != '#858585') {	
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.lineJoin = "round";  
						targetLayer.CanvasContext.moveTo(rvds.ContourLinePoints[index].X1 - .05, rvds.ContourLinePoints[index].Y1 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X2 - .05, rvds.ContourLinePoints[index].Y2 - .3);
						targetLayer.CanvasContext.lineTo(rvds.ContourLinePoints[index].X3 - .05, rvds.ContourLinePoints[index].Y3 - .3);
						targetLayer.CanvasContext.setLineDash([]);
						targetLayer.CanvasContext.strokeStyle = targetLayer.dataLayerColors[index];	
						targetLayer.CanvasContext.lineWidth = targetLayer.ScaleFactor * 1.0 * 1.9 * CircleSize;					
						targetLayer.CanvasContext.stroke();
						targetLayer.CanvasContext.closePath();
					} else 	{
						//targetLayer.dataLayerColors[index] = undefined;
					}
				});
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	function drawSelection(targetLayer,SeleName) {
		var SelectionList =[];
		if (!SeleName){
			//SeleName = "Main";
			$('.checkBoxDIV-S').find(".visibilityCheckImg[value=visible]").parent().parent().each(function (index){SelectionList.push($(this).attr("name"))});
		} else {
			SelectionList[0]=SeleName;
		}
		if(this.Circle_Radius){
			var CircleSize = this.Circle_Radius;
		} else {
			var CircleSize = 1.7;
		}
		targetLayer.clearCanvas();
		targetLayer.Data = [];
		targetLayer.dataLayerColors = [];
		if (this.Residues && this.Residues.length > 0) {
			for (var i = this.Residues.length - 1; i >= 0; i--) {
				targetLayer.Data[i] = false;
				targetLayer.dataLayerColors[i] = "#858585";
			}
			for (var k = SelectionList.length - 1 ; k >= 0 ; k--){
				var targetSelection = this.getSelection(SelectionList[k]);
				for (var j = targetSelection.Residues.length - 1; j >= 0; j--) {
					targetLayer.CanvasContext.beginPath();
					targetLayer.CanvasContext.arc(parseFloat(targetSelection.Residues[j].X) + this.PageOffset[0], targetSelection.Residues[j].Y, (targetLayer.ScaleFactor * CircleSize), 0, 2 * Math.PI, false);
					targetLayer.CanvasContext.closePath(); 
					targetLayer.CanvasContext.strokeStyle = targetSelection.Color;
					targetLayer.CanvasContext.lineWidth = 0.5;
					targetLayer.CanvasContext.stroke();
					targetLayer.Data[targetSelection.Residues[j].map_Index - 1] = true;
					targetLayer.dataLayerColors[targetSelection.Residues[j].map_Index - 1] = targetSelection.Color;
				}
			}
			var linkedLayer = this.getLinkedLayer();
			if (linkedLayer.Type == "selected") {
				update3Dcolors();
			}
		}
		//Temp draw box
		//targetLayer.CanvasContext.beginPath();
		//targetLayer.CanvasContext.arc(parseFloat(targetSelection.Residues[j].X) + this.PageOffset[0], targetSelection.Residues[j].Y, (targetLayer.ScaleFactor * CircleSize), 0, 2 * Math.PI, false);
		//targetLayer.CanvasContext.closePath();
		//targetLayer.CanvasContext.strokeStyle = 'orange';
		//targetLayer.CanvasContext.rect(108-18,16,490+18,722);
		//targetLayer.CanvasContext.lineWidth = 4;
		//targetLayer.CanvasContext.stroke();
		//targetLayer.CanvasContext.strokeStyle = 'green';
		//targetLayer.CanvasContext.rect(108+612-18-36,16,490+18,722);
		//targetLayer.CanvasContext.lineWidth = 4;
		//targetLayer.CanvasContext.stroke();
		
		
	}
	function drawDataCircles(targetLayer, dataIndices, ColorArray, noClear) {
		if (targetLayer.Type === "circles") {
			if (!noClear) {
				targetLayer.clearCanvas();
				targetLayer.dataLayerColors = [];
			}
			if(this.Circle_Radius){
				var CircleSize = this.Circle_Radius;
			} else {
				var CircleSize = 1.7;
			}
			if (this.Residues != undefined) {
				for (var i = this.Residues.length - 1; i >= 0; i--) {
					if (dataIndices && ColorArray && ColorArray[dataIndices[i]] != undefined && ColorArray[dataIndices[i]] != '#858585') {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.arc(parseFloat(this.Residues[i].X) + this.PageOffset[0], this.Residues[i].Y, (targetLayer.ScaleFactor * CircleSize), 0, 2 * Math.PI, false);
						targetLayer.CanvasContext.closePath();
						targetLayer.CanvasContext.strokeStyle = ColorArray[dataIndices[i]];
						targetLayer.CanvasContext.stroke();
						if (targetLayer.Filled) {
							targetLayer.CanvasContext.fillStyle = ColorArray[dataIndices[i]];
							targetLayer.CanvasContext.fill();
						}
						targetLayer.dataLayerColors[i] = ColorArray[dataIndices[i]];
					} else if (!dataIndices && !ColorArray && targetLayer.dataLayerColors[i] && targetLayer.dataLayerColors[i] != '#858585') {
						targetLayer.CanvasContext.beginPath();
						targetLayer.CanvasContext.arc(parseFloat(this.Residues[i].X) + this.PageOffset[0], this.Residues[i].Y, (targetLayer.ScaleFactor * CircleSize), 0, 2 * Math.PI, false);
						targetLayer.CanvasContext.closePath();
						targetLayer.CanvasContext.strokeStyle = targetLayer.dataLayerColors[i];
						targetLayer.CanvasContext.stroke();
						if (targetLayer.Filled) {
							targetLayer.CanvasContext.fillStyle = targetLayer.dataLayerColors[i];
							targetLayer.CanvasContext.fill();
						}
					}
				}
			}
		} else {
			return false;
		}
	};

	function drawBasePairs(targetLayer, colorLayer){
		//alert("now mode");
		var rvds = this;
		var color1,color2;
		var zoomEnabled = $('input[name="za"][value=on]').is(':checked');
		targetLayer.clearCanvas();
		if (!colorLayer) {
			var colorLayer = targetLayer.ColorLayer;
		} else {
			targetLayer.ColorLayer = colorLayer;
		}
		targetLayer.Data=ActiveBasePairSet;
		
		if (targetLayer.Data != undefined || targetLayer.Data == []) {
			if (targetLayer.ColorGradientMode == "Matched") {
				var grd_order = [0, 1];
			} else if (targetLayer.ColorGradientMode == "Opposite") {
				var grd_order = [1, 0];
			} else {
				alert("how did we get here?");
			}
			$.each(targetLayer.Data, function(index,base_pair){
				var residue_i = MainResidueMap[base_pair.residue_i];
				var residue_j = MainResidueMap[base_pair.residue_j];
				//Come back and make zoom aware work correctly, with the same color and opacity as would be in other modes. 
				if (zoomEnabled) {
					var jkdist = Math.sqrt(((residue_i.X - residue_j.X) * (residue_i.X - residue_j.X) + (residue_i.Y - residue_j.Y) * (residue_i.Y - residue_j.Y)));
					
					if ((150 - rvViews[0].scale * 23) > jkdist) {
						base_pair.color = "rgba(35,31,32," + base_pair.opacity + ")";
						return true;
					}
					if (((residue_i.X * rvViews[0].scale + rvViews[0].x < 0) || (residue_i.X * rvViews[0].scale + rvViews[0].x > rvViews[0].clientWidth) || (residue_i.Y * rvViews[0].scale + rvViews[0].y < 0) || (residue_i.Y * rvViews[0].scale + rvViews[0].y > rvViews[0].clientHeight))
						 && ((residue_j.X * rvViews[0].scale + rvViews[0].x < 0) || (residue_j.X * rvViews[0].scale + rvViews[0].x > rvViews[0].clientWidth) || (residue_j.Y * rvViews[0].scale + rvViews[0].y < 0) || (residue_j.Y * rvViews[0].scale + rvViews[0].y > rvViews[0].clientHeight))) {
						base_pair.color = "rgba(35,31,32," + base_pair.color.opacity + ")";
						return true;
					}
				}
				if (residue_i && residue_j ) {
					switch (colorLayer.Type) {
						case undefined:
							if (colorLayer == "gray_lines") {
								var grd = rvds.HighlightLayer.CanvasContext.createLinearGradient(residue_i.X, residue_i.Y, residue_j.X, residue_j.Y);
								color1 = colorNameToHex("#231F20");
								color2 = colorNameToHex("#231F20");
								
								grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + base_pair.opacity + ")");
								grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + base_pair.opacity + ")");
								
								base_pair.color = grd;
								base_pair.color_hex = color1;
							} else if (colorLayer == "manual_coloring") {
								var grd = rvds.HighlightLayer.CanvasContext.createLinearGradient(residue_i.X, residue_i.Y, residue_j.X, residue_j.Y);
								color1 = base_pair.color_hex;
								color2 = base_pair.color_hex;
								
								grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + base_pair.opacity + ")");
								grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + base_pair.opacity + ")");
								
								base_pair.color = grd;
								base_pair.color_hex = color1;
							} else {
								alert("Invalid color mode");
							}
							break;
						case "residues":
							var grd = colorLayer.CanvasContext.createLinearGradient(residue_i.X, residue_i.Y, residue_j.X, residue_j.Y);
							
							if (rvDataSets[residue_i.rvds_index].Residues[residue_i.index].color && rvDataSets[residue_j.rvds_index].Residues[residue_j.index].color) {
								color1 = colorNameToHex(rvDataSets[residue_i.rvds_index].Residues[residue_i.index].color);
								color2 = colorNameToHex(rvDataSets[residue_j.rvds_index].Residues[residue_j.index].color);
								
								grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + base_pair.opacity + ")");
								grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + base_pair.opacity + ")");
							} else {
								color1='#231F20';
							}
							//colorLayer.addLinearGradient(grd);
							base_pair.color = grd;
							base_pair.color_hex = color1;
							break;
						case "selected":
						case "contour":
						case "circles":
							var grd = colorLayer.CanvasContext.createLinearGradient(residue_i.X, residue_i.Y, residue_j.X, residue_j.Y);
							//This will be broken for now, because of j,k indices
							if (rvDataSets[residue_i.rvds_index].Layers[colorLayer.zIndex].dataLayerColors[residue_i.index] && rvDataSets[residue_j.rvds_index].Layers[colorLayer.zIndex].dataLayerColors[residue_j.index]) {
								color1 = colorNameToHex(rvDataSets[residue_i.rvds_index].Layers[colorLayer.zIndex].dataLayerColors[residue_i.index]);
								color2 = colorNameToHex(rvDataSets[residue_j.rvds_index].Layers[colorLayer.zIndex].dataLayerColors[residue_j.index]);
								
								grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + base_pair.opacity + ")");
								grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + base_pair.opacity + ")");
							} else {
								color1='#231F20';
							}
							//colorLayer.addLinearGradient(grd);
							base_pair.color = grd;
							base_pair.color_hex = color1;
							break;
						default:
							alert("this shouldn't be happening right now.");
					}
					//Regular Mode
					
					targetLayer.CanvasContext.beginPath();
					targetLayer.CanvasContext.moveTo(residue_i.X, residue_i.Y);
					targetLayer.CanvasContext.lineTo(residue_j.X, residue_j.Y);
					targetLayer.CanvasContext.strokeStyle = base_pair.color;	
					targetLayer.CanvasContext.lineWidth = base_pair.color.lineWidth;					
					targetLayer.CanvasContext.stroke();
					targetLayer.CanvasContext.closePath();
					if (zoomEnabled && (rvViews[0].scale > 10)) {
						//draw the interaction type labels here
						var x1 = residue_i.X;
						var x2 = residue_j.X;
						var x12mid = x1 - ((x1 - x2) / 2);
						var xmid = residue_i.X - (residue_i.X - residue_j.X) / 2;
						var ymid = residue_i.Y - (residue_i.Y - residue_j.Y) / 2;
						targetLayer.CanvasContext.save();
						targetLayer.CanvasContext.lineWidth = 0.5;
						targetLayer.CanvasContext.fillStyle = "white";
						targetLayer.CanvasContext.fillRect(xmid - 2.4, ymid - .8, 4.7, 1.7);
						targetLayer.CanvasContext.strokeRect(xmid - 2.3, ymid - .7, 4.5, 1.5);
						targetLayer.CanvasContext.restore();
						targetLayer.CanvasContext.save();
						targetLayer.CanvasContext.font = ".5px Arial";
						targetLayer.CanvasContext.fillText(base_pair.bp_type, xmid - 2, ymid + .5);
						targetLayer.CanvasContext.restore();
						
					}
					
				}
			});
		}
		ActiveBasePairSet=targetLayer.Data;
	}
	/* function drawBasePairs(targetLayer, colorLayer,newmode) {
		if(true){
		//if(newmode){
			drawBasePairs2.call(this,targetLayer, colorLayer);
			return;
		}
		var color1,color2;
		var zoomEnabled = $('input[name="za"][value=on]').is(':checked');
		targetLayer.clearCanvas();
		if (!colorLayer) {
			var colorLayer = targetLayer.ColorLayer;
		} else {
			targetLayer.ColorLayer = colorLayer;
		}
		targetLayer.Data=this.BasePairs;
		
		if (targetLayer.Data != undefined || targetLayer.Data == []) {
			if (targetLayer.ColorGradientMode == "Matched") {
				var grd_order = [0, 1];
			} else if (targetLayer.ColorGradientMode == "Opposite") {
				var grd_order = [1, 0];
			} else {
				alert("how did we get here?");
			}
			for (var i = 0; i < targetLayer.Data.length; i++) {
				var j = targetLayer.Data[i].resIndex1;
				var k = targetLayer.Data[i].resIndex2;
				//Come back and make zoom aware work correctly, with the same color and opacity as would be in other modes. 
				if (zoomEnabled) {
					var jkdist = Math.sqrt(((this.Residues[j].X - this.Residues[k].X) * (this.Residues[j].X - this.Residues[k].X) + (this.Residues[j].Y - this.Residues[k].Y) * (this.Residues[j].Y - this.Residues[k].Y)));
					
					if ((150 - rvViews[0].scale * 23) > jkdist) {
						targetLayer.Data[i]["color"] = "rgba(35,31,32," + targetLayer.Data[i].opacity + ")";
						continue;
					}
					if (((this.Residues[j].X * rvViews[0].scale + rvViews[0].x < 0) || (this.Residues[j].X * rvViews[0].scale + rvViews[0].x > rvViews[0].clientWidth) || (this.Residues[j].Y * rvViews[0].scale + rvViews[0].y < 0) || (this.Residues[j].Y * rvViews[0].scale + rvViews[0].y > rvViews[0].clientHeight))
						 && ((this.Residues[k].X * rvViews[0].scale + rvViews[0].x < 0) || (this.Residues[k].X * rvViews[0].scale + rvViews[0].x > rvViews[0].clientWidth) || (this.Residues[k].Y * rvViews[0].scale + rvViews[0].y < 0) || (this.Residues[k].Y * rvViews[0].scale + rvViews[0].y > rvViews[0].clientHeight))) {
						targetLayer.Data[i]["color"] = "rgba(35,31,32," + targetLayer.Data[i].opacity + ")";
						continue;
					}
				}
				if (j >= 0 && k >= 0) {
					switch (colorLayer.Type) {
					case undefined:
						if (colorLayer == "gray_lines") {
							var grd = this.HighlightLayer.CanvasContext.createLinearGradient(this.Residues[j].X, this.Residues[j].Y, this.Residues[k].X, this.Residues[k].Y);
							color1 = colorNameToHex("#231F20");
							color2 = colorNameToHex("#231F20");
							
							grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							
							targetLayer.Data[i]["color"] = grd;
							targetLayer.Data[i]["color_hex"] = color1;
						} else if (colorLayer == "manual_coloring") {
							var grd = this.HighlightLayer.CanvasContext.createLinearGradient(this.Residues[j].X, this.Residues[j].Y, this.Residues[k].X, this.Residues[k].Y);
							color1 = targetLayer.Data[i]["color_hex"];
							color2 = targetLayer.Data[i]["color_hex"];
							
							grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							
							targetLayer.Data[i]["color"] = grd;
							targetLayer.Data[i]["color_hex"] = color1;
						} else {
							alert("Invalid color mode");
						}
						break;
					case "residues":
						var grd = colorLayer.CanvasContext.createLinearGradient(this.Residues[j].X, this.Residues[j].Y, this.Residues[k].X, this.Residues[k].Y);
						
						if (this.Residues[j].color && this.Residues[k].color) {
							color1 = colorNameToHex(this.Residues[j].color);
							color2 = colorNameToHex(this.Residues[k].color);
							
							grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
						} else {
							color1='#231F20';
						}
						//colorLayer.addLinearGradient(grd);
						targetLayer.Data[i]["color"] = grd;
						targetLayer.Data[i]["color_hex"] = color1;
						break;
					case "contour":
					case "circles":
						var grd = colorLayer.CanvasContext.createLinearGradient(this.Residues[j].X, this.Residues[j].Y, this.Residues[k].X, this.Residues[k].Y);
						if (colorLayer.dataLayerColors[j] && colorLayer.dataLayerColors[k]) {
							color1 = colorNameToHex(colorLayer.dataLayerColors[j]);
							color2 = colorNameToHex(colorLayer.dataLayerColors[k]);
							
							grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
						} else {
							color1='#231F20';
						}
						//colorLayer.addLinearGradient(grd);
						targetLayer.Data[i]["color"] = grd;
						targetLayer.Data[i]["color_hex"] = color1;
						break;
					case "selected":
						var grd = colorLayer.CanvasContext.createLinearGradient(this.Residues[j].X, this.Residues[j].Y, this.Residues[k].X, this.Residues[k].Y);
						if (colorLayer.Data[j] || colorLayer.Data[k]) {
							color1 = colorNameToHex(colorLayer.dataLayerColors[j]);
							color2 = colorNameToHex(colorLayer.dataLayerColors[k]);
							//color1 = colorNameToHex("#231F20");
							//color2 = colorNameToHex("#231F20");
							grd.addColorStop(grd_order[0], "rgba(" + h2d(color1.slice(1, 3)) + "," + h2d(color1.slice(3, 5)) + "," + h2d(color1.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
							grd.addColorStop(grd_order[1], "rgba(" + h2d(color2.slice(1, 3)) + "," + h2d(color2.slice(3, 5)) + "," + h2d(color2.slice(5)) + "," + targetLayer.Data[i].opacity + ")");
						} else {
							color1='#231F20';
						}
						//colorLayer.addLinearGradient(grd);
						targetLayer.Data[i]["color"] = grd;
						targetLayer.Data[i]["color_hex"] = color1;
						break;
					default:
						alert("this shouldn't be happening right now.");
					}
					//Regular Mode
					
					targetLayer.CanvasContext.beginPath();
					targetLayer.CanvasContext.moveTo(parseFloat(this.Residues[j].X) + this.PageOffset[0], this.Residues[j].Y);
					targetLayer.CanvasContext.lineTo(parseFloat(this.Residues[k].X) + this.PageOffset[0], this.Residues[k].Y);
					targetLayer.CanvasContext.strokeStyle = targetLayer.Data[i]["color"];	
					targetLayer.CanvasContext.lineWidth = targetLayer.Data[i].lineWidth;					
					targetLayer.CanvasContext.stroke();
					targetLayer.CanvasContext.closePath();
					if (zoomEnabled && (rvViews[0].scale > 10)) {
						//draw the interaction type labels here
						var x1 = this.Residues[j].X;
						var x2 = this.Residues[k].X;
						var x12mid = x1 - ((x1 - x2) / 2) + this.PageOffset[0];
						var xmid = this.Residues[j].X - (this.Residues[j].X - this.Residues[k].X) / 2 + this.PageOffset[0];
						var ymid = this.Residues[j].Y - (this.Residues[j].Y - this.Residues[k].Y) / 2;
						targetLayer.CanvasContext.save();
						targetLayer.CanvasContext.lineWidth = 0.5;
						targetLayer.CanvasContext.fillStyle = "white";
						targetLayer.CanvasContext.fillRect(xmid - 2.4, ymid - .8, 4.7, 1.7);
						targetLayer.CanvasContext.strokeRect(xmid - 2.3, ymid - .7, 4.5, 1.5);
						targetLayer.CanvasContext.restore();
						targetLayer.CanvasContext.save();
						targetLayer.CanvasContext.font = ".5px Arial";
						targetLayer.CanvasContext.fillText(targetLayer.Data[i].bp_type, xmid - 2, ymid + .5);
						targetLayer.CanvasContext.restore();
						
					}
				}
			}
		}
		this.BasePairs=targetLayer.Data;
	} */
};

function rvView(x, y, scale) {
	//Properties
	this.x = x;
	this.y = y;
	this.scale = scale;
	this.lastX = [];
	this.lastY = [];
	this.startX = [];
	this.startY = [];
	this.width = [];
	this.height = [];
	this.clientWidth = [];
	this.clientHeight = [];
	this.slideValue = 20;
	this.defaultScale = 1;
	this.defaultX = 0;
	this.defaultY = 0;
	this.windowHeight=0;
	this.xcorr=0;
	this.ycorr=0;
	
	//Methods
	this.zoom = function (event, delta) {
		if (($("#slider").slider("value") + delta) <= $("#slider").slider("option", "max") & ($("#slider").slider("value") + delta) >= $("#slider").slider("option", "min")) {
			this.zoomDelta(event.clientX - this.xcorr, event.clientY - this.ycorr, Math.pow(1.1, delta));
			$("#slider").slider("value", $("#slider").slider("value") + delta);
			this.slideValue = $("#slider").slider("value") + delta;
		}
	}
	this.zoomDelta = function (px, py, factor){
		//console.log([px,py]);
		this.scale *= factor;
		this.x = (this.x - px) * factor + px;
		this.y = (this.y - py) * factor + py;
		
		$.each(rvDataSets, function (index, value) {
			value.refreshCanvases();
		});
	}
	
	this.centerZoom = function (scale) {
		this.scale = scale;
		this.x = (rvDataSets[0].HighlightLayer.Canvas.width - this.scale*612)/2;
		this.y = (rvDataSets[0].HighlightLayer.Canvas.height - this.scale*792)/2;
		rvViews[0].width = rvDataSets[0].HighlightLayer.Canvas.width;
		rvViews[0].height = rvDataSets[0].HighlightLayer.Canvas.height;
		rvViews[0].clientWidth = rvDataSets[0].HighlightLayer.Canvas.clientWidth;
		rvViews[0].clientHeight = rvDataSets[0].HighlightLayer.Canvas.clientHeight;
	}
	this.pan = function (event) {
		rvViews[0].x += event.clientX - rvViews[0].lastX;
		rvViews[0].y += event.clientY - rvViews[0].lastY;
		$.each(rvDataSets, function (index, value) {
			value.refreshCanvases();
		});
		rvViews[0].lastX = event.clientX;
		rvViews[0].lastY = event.clientY;
	}
	this.panXY = function (dx,dy) {
		rvViews[0].x += dx;
		rvViews[0].y += dy;
		$.each(rvDataSets, function (index, value) {
			value.refreshCanvases();
		});
		rvViews[0].lastX += dx;
		rvViews[0].lastY += dy;
	}
	this.drag = function (event) {
		var targetSelection=rvDataSets[0].getSelection($('input:radio[name=selectedRadioS]').filter(':checked').parent().parent().attr('name'));
		//rvDataSets[0].HighlightLayer.CanvasContext.strokeStyle = "#ff0000";
		rvDataSets[0].HighlightLayer.CanvasContext.strokeStyle = targetSelection.Color;
		rvDataSets[0].HighlightLayer.CanvasContext.strokeRect(this.startX, this.startY, (event.clientX - this.xcorr - this.x) / this.scale - this.startX, (event.clientY - this.ycorr - this.y) / this.scale - this.startY);
	}
	this.restore = function () {
		$("#slider").slider("value",this.slideValue);
		$.each(rvDataSets, function (index, value) {
			value.refreshCanvases();
		});
	}
	
	this.toJSON = function () {
		return {
			x: this.x,
			y: this.y, 
			scale: this.scale, 
			lastX: this.lastX,
			lastY: this.lastY,
			startX: this.startX,
			startY: this.startY,
			width: this.width,
			height: this.height,
			clientWidth: this.clientWidth,
			clientHeight: this.clientHeight,
			slideValue : this.slideValue
		};
	}
	
	this.fromJSON = function (json) {
		var data = JSON.parse(json);
		//var data = json;
		var e = new rvView(data.x, data.y, data.scale);
		e.lastX = data.lastX;
		e.lastY = data.lastY;
		e.startX = data.startX;
		e.startY = data.startY;
		e.width = data.width;
		e.height = data.height;
		e.clientWidth = data.clientWidth;
		e.clientHeight = data.clientHeight;
		e.slideValue = data.slideValue;
		return e;
	}
	
};
///////////////////////////////////////////////////////////////////////////////
