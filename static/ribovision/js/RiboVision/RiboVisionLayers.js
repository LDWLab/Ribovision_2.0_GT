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

function InitLayers() {
	// Build Layer Menu

	// Put in Top Labels and ToolBar
	$("#LayerPanel").prepend($("<div id='tipBar'>").attr({
			'name' : 'TopLayerBar'
		}).html("C&nbspV&nbsp&nbsp&nbspS&nbsp&nbsp&nbsp&nbspL&nbsp&nbsp&nbsp&nbspLayerName&nbsp&nbsp")); // where to add letters
	$('[name=TopLayerBar]').append($('<button id="newLayer" class="toolBarBtn2" title="Create a new layer"></button>'));
	$('[name=TopLayerBar]').append($('<button id="clearLayer" class="toolBarBtn2" title="Clear the selected layer"></button>'));
	$('[name=TopLayerBar]').append($('<button id="deleteLayer" class="toolBarBtn2" title="Delete the selected layer"></button>'));
	$("#newLayer").button({
		text : false,
		icons : {
			primary : "ui-icon-document"
		}
	});
	$("#newLayer").click(function () {
		$("#dialog-addLayer").dialog("open");
	});

	$("#clearLayer").button({
		text : false,
		icons : {
			primary : "ui-icon-cancel"
		}
	});

	$("#clearLayer").click(function () {
		$.each(rvDataSets, function (index, rvds) {
			var targetLayer = rvds.getSelectedLayer();
			var linkedLayer = rvds.getLinkedLayer();
			targetLayer.clearAll();
			if (targetLayer == linkedLayer){
				update3Dcolors();
			}
		});
		
		drawNavLine();
		//$("#dialog-addSelection").dialog("open");
	});

	$("#deleteLayer").button({
		text : false,
		icons : {
			primary : "ui-icon-trash"
		}
	});

	$("#deleteLayer").click(function (event) {
		$("#dialog-confirm-delete p").text("The " + rvDataSets[0].getSelectedLayer().LayerName + " layer will be permanently deleted and cannot be recovered.");
		$("#dialog-confirm-delete").dialog('open');
	});

	$(".toolBarBtn2").css('height', $("#openLayerBtn").css('height'));

	$("#LayerPanel").sortable({
		update : function (event, ui) {
			$("#LayerPanel .layerContent").each(function (e, f) {
				$.each(rvDataSets, function (index, rvds) {
					var tl = rvds.getLayer($(this).parent().attr("name"));
					tl.updateZIndex(rvds.LastLayer - e);
				});
			});
			$.each(rvDataSets, function (index, rvds) {
				rvds.sort();
			});
			RefreshLayerMenu();
		},
		items : ".oneLayerGroup",
		axis : "y"

	});
	$("#LayerPanel").disableSelection();

	$("#LayerPreferenceDialog").dialog({
		autoOpen : false,
		show : {
			effect : "blind",
			duration : 300
		},
		modal : true,
		height : 600,
		width : 400,
		position : {
			my : "center",
			at : "center",
			of : $("#canvasDiv")
		},
		buttons : {
			"Save" : function () {
				var ret = changeCurrentLayerName();
				var ret2 = changeLayerColor();
				if (ret && ret2) {
					$("#LayerPreferenceDialog").dialog("close");
				}
			},
			Cancel : function () {
				$(this).dialog("close");
			}
		},
		open : function (event) {
			$("#myJmol_object").css("visibility", "hidden");
		},
		close : function () {
			$("#myJmol_object").css("visibility", "visible");
		}
	});

	$("#dialog-addLayer").dialog({
		resizable : false,
		autoOpen : false,
		height : "auto",
		width : 400,
		modal : true,
		buttons : {
			"Create New Layer" : function () {
				var namecheck = $("#newLayerName").val().match(/[A-z][\w-_:\.]*/);
				if (namecheck !== null && namecheck[0].length === $("#newLayerName").val().length && $("#newLayerName").val().length <= 16) {
					if (rvDataSets[0].isUniqueLayer($("#newLayerName").val())) {
						//$("#canvasDiv").append($('<canvas id="' + $("#newLayerName").val() + '" style="z-index:' + ( rvDataSets[0].LastLayer + 1 ) + ';"></canvas>'));
						$.each(rvDataSets, function (index, rvds) {
							rvds.addLayer($("#newLayerName").val(), $("#newLayerName").val(), [], true, 1.0, 'circles', $("#layerColor2").val());
						});
						resizeElements();
						LayerMenu(rvDataSets[0].getLayer($("#newLayerName").val()), (1000 + (rvDataSets[0].LastLayer + 1)));
						$(".oneLayerGroup[name=" + $("#newLayerName").val() + "]").find(".selectLayerRadioBtn").prop("checked", true);
						$(".oneLayerGroup[name=" + $("#newLayerName").val() + "]").find(".selectLayerRadioBtn").trigger("change");
						RefreshLayerMenu();
						var farbobj = $.farbtastic("#layerColorPicker2");
						farbobj.setColor(colorNameToHex($("#layerColor2").val()));
						$(this).dialog("close");
					} else {
						$("#dialog-unique-layer-error").dialog("open");
					}
				} else {
					$("#dialog-name-error").dialog("open");
				}
			},
			Cancel : function () {
				$(this).dialog("close");
			}
		},
		open : function () {
			$("#newLayerName").val("Layer_" + (rvDataSets[0].Layers.length + 1));
			$("#myJmol_object").css("visibility", "hidden");
		},
		close : function () {
			$("#myJmol_object").css("visibility", "visible");
		}
	});
	$("#newLayerName").button().addClass('ui-textfield').keydown(function (event) {
		if (event.keyCode == 13) {
			$("#dialog-addLayer").dialog("option", "buttons")['Create New Layer'].apply($("#dialog-addLayer"));
		}
	});
	$("#layerColor2").button().addClass('ui-textfield').keydown(function (event) {
		if (event.keyCode == 13) {
			$("#dialog-addLayer").dialog("option", "buttons")['Create New Layer'].apply($("#dialog-addLayer"));
		}
	});
	$("#layerColor").button().addClass('ui-textfield').keydown(function (event) {
		if (event.keyCode == 13) {
			//$("#LayerPreferenceDialog").dialog("option", "buttons")['Save'].apply($("#LayerPreferenceDialog"));
			var ret2 = changeLayerColor();
		}
	});

	$("#layerNameInput").button().addClass('ui-textfield').keydown(function (event) {
		if (event.keyCode == 13) {
			//$("#LayerPreferenceDialog").dialog("option", "buttons")['Save'].apply($("#LayerPreferenceDialog"));
			var ret = changeCurrentLayerName();
		}
	});

	//MiniLayer
	$("#MiniLayer").sortable({
		update : function (event, ui) {
			//Determine Layer to move
			var target = $(ui.item[0]).attr("name");
			if (target === $(".miniLayerName:first").attr("name")) {
				var Pos = $(".miniLayerName[name=" + $(ui.item[0]).attr("name") + "]").next().attr("name");
				//Move Layer
				$("#LayerPanel").find("[name=" + target + "]").insertBefore($("#LayerPanel").find("[name=" + Pos + "]"));
			} else {
				var Pos = $(".miniLayerName[name=" + $(ui.item[0]).attr("name") + "]").prev().attr("name");
				//Move Layer
				$("#LayerPanel").find("[name=" + target + "]").insertAfter($("#LayerPanel").find("[name=" + Pos + "]"));
			}

			//Update ZIndex
			$("#LayerPanel .layerContent").each(function (e, f) {
				var layer_this=this;
				$.each(rvDataSets, function (index, rvds) {
					var tl = rvds.getLayer($(layer_this).parent().attr("name"));
					tl.updateZIndex(rvds.LastLayer - e);
				});
			});
			$.each(rvDataSets, function (index, rvds) {
				rvds.sort();
			});
			RefreshLayerMenu();
		},
		start : function (event, ui) {
			$("#LinkSection").addClass("highlightBorder");
		},
		stop : function (event, ui) {
			$("#LinkSection").removeClass("highlightBorder");
		},
		items : ".miniLayerName",
		axis : "y"
	});
	$("#MiniLayer").disableSelection();

	$("#dialog-confirm-delete").dialog({
		resizable : false,
		autoOpen : false,
		height : "auto",
		width : 400,
		modal : true,
		buttons : {
			"Delete the Layer" : function (event) {
				$.each(rvDataSets, function (index, rvds) {
					var targetLayer = rvDataSets[0].getSelectedLayer();
					if (targetLayer.Type !== "residues") {
						rvds.deleteLayer(targetLayer.LayerName);
						$("[name=" + targetLayer.LayerName + "]").remove();
						rvds.sort();

						//check selected layer and linked layer
						var rLayers = rvds.getLayerByType("residues");
						var cLayers = rvds.getLayerByType("circles");
						var sLayers = rvds.getLayerByType("selected");
						if (rLayers.length >= 1) {
							$(".oneLayerGroup[name='" + rLayers[0].LayerName + "']").find(".selectLayerRadioBtn").prop("checked", true);
							$(".oneLayerGroup[name='" + rLayers[0].LayerName + "']").find(".selectLayerRadioBtn").trigger("change");
						} else if (cLayers.length >= 1) {
							$(".oneLayerGroup[name='" + cLayers[0].LayerName + "']").find(".selectLayerRadioBtn").prop("checked", true);
							$(".oneLayerGroup[name='" + cLayers[0].LayerName + "']").find(".selectLayerRadioBtn").trigger("change");
						} else if (sLayers.length >= 1) {
							$(".oneLayerGroup[name='" + sLayers[0].LayerName + "']").find(".selectLayerRadioBtn").prop("checked", true);
							$(".oneLayerGroup[name='" + sLayers[0].LayerName + "']").find(".selectLayerRadioBtn").trigger("change");
						} else {
							$(".oneLayerGroup[name='" + rvds.Layers[0].LayerName + "']").find(".selectLayerRadioBtn").prop("checked", true);
							$(".oneLayerGroup[name='" + rvds.Layers[0].LayerName + "']").find(".selectLayerRadioBtn").trigger("change");
						}
						if (targetLayer.Linked) {
							if (rLayers.length >= 1) {
								$(".oneLayerGroup[name='" + rLayers[0].LayerName + "']").find(".mappingRadioBtn").prop("checked", true);
								$(".oneLayerGroup[name='" + rLayers[0].LayerName + "']").find(".mappingRadioBtn").trigger("change");
							} else if (cLayers.length >= 1) {
								$(".oneLayerGroup[name='" + cLayers[0].LayerName + "']").find(".mappingRadioBtn").prop("checked", true);
								$(".oneLayerGroup[name='" + cLayers[0].LayerName + "']").find(".mappingRadioBtn").trigger("change");
							} else if (sLayers.length >= 1) {
								$(".oneLayerGroup[name='" + sLayers[0].LayerName + "']").find(".mappingRadioBtn").prop("checked", true);
								$(".oneLayerGroup[name='" + sLayers[0].LayerName + "']").find(".mappingRadioBtn").trigger("change");
							} else {
								if ($(".oneLayerGroup[name='" + rvds.Layers[0].LayerName + "']").find(".mappingRadioBtn").prop("disabled") == false) {
									$(".oneLayerGroup[name='" + rvds.Layers[0].LayerName + "']").find(".mappingRadioBtn").prop("checked", true);
									$(".oneLayerGroup[name='" + rvds.Layers[0].LayerName + "']").find(".mappingRadioBtn").trigger("change");
								}
							}
						}
						RefreshLayerMenu();
					} else {
						alert("You are not allowed to delete the residues layer");
					}
				});
				$(this).dialog("close");
			},
			Cancel : function () {
				$(this).dialog("close");
			}
		},
		open : function (event) {
			$("#myJmol_object").css("visibility", "hidden");
		},
		close : function () {
			$("#myJmol_object").css("visibility", "visible");
		}
	});

}

function changeLayerColor() {
	var newcolor,
	newcolorH,
	farbobj;
	var colorcheck = $("#layerColor").val().match(/#[\dABCDEFabcdef]{6,6}$/);
	if (colorcheck && (colorcheck[0].length === 7)) {
		newcolor = colorcheck[0];
	} else {
		newcolorH = colorNameToHex($("#layerColor").val());
		if (newcolorH === false) {
			$("#dialog-invalid-color-error").dialog("open");
			return false;
		} else {
			newcolor = $("#layerColor").val();
			farbobj = $.farbtastic("#layerColorPicker");
			farbobj.setColor(newcolorH);
		}
	}
	$($dblClickedLayer).parent().find(".colorBox").css("background", newcolor);
	$.each(rvDataSets, function (index, rvds) {
		var targetLayer = rvds.getLayer($dblClickedLayerName);
		targetLayer.Color = newcolor;
	});
	drawNavLine();
	return true;
}

//in "Layer Preferfence"
function changeCurrentLayerName() {
	if ($("#layerNameInput").val() != "" && ($("#layerNameInput").val() !== $($dblClickedLayer).parent().attr("name"))) {
		var namecheck = $("#layerNameInput").val().match(/[A-z][\w-_:\.]*/);
		if (namecheck !== null && namecheck[0].length === $("#layerNameInput").val().length && $("#layerNameInput").val().length <= 16) {
			if (rvDataSets[0].isUniqueLayer($("#layerNameInput").val())) {
				$("#MiniLayer").find("[name=" + $($dblClickedLayer).parent().attr("name") + "]").attr("name", $("#layerNameInput").val()).text($("#layerNameInput").val());
				$($dblClickedLayer).parent().attr("name", $("#layerNameInput").val());
				$($dblClickedLayer).html($("#layerNameInput").val()).prepend('<span class="ui-accordion-header-icon ui-icon ui-icon-triangle-1-e"></span>');
				$.each(rvDataSets, function (index, rvds) {
					var targetLayer = rvds.getLayer($dblClickedLayerName);
					targetLayer.LayerName = $("#layerNameInput").val();
					if (targetLayer.Linked) {
						$("#LinkSection").find(".miniLayerName").attr("name", $("#layerNameInput").val()).text($("#layerNameInput").val());
					}
				});
				RefreshLayerMenu();
				return true;
			} else {
				$("#dialog-unique-layer-error").dialog("open");
				return false;
			}
		} else {
			$("#dialog-name-error").dialog("open");
			return false;
		}
	} else {
		return true;
	}
}

function LayerMenu(Layer, key, RVcolor) {
	//console.log($count);
	var $currentLayerName = Layer.LayerName;
	//console.log($currentLayerName);
	$('[name=TopLayerBar]').after($("<div>").addClass("oneLayerGroup").attr({
			'name' : $currentLayerName
		}));
	//console.log($("<div>").addClass("oneLayerGroup").attr({'name' : $currentLayerName}).html());
	//$('.oneLayerGroup').append(("<div>").addClass("visibilityCheckbox"));
	//$('.visibilityCheckbox').append($("<input />").attr({type:'checkbox'}).addClass("layerVisiblity"))

	var $currentGroup = $(".oneLayerGroup[name=" + $currentLayerName + "]")[0];

	//adding color box
	$($currentGroup)
	.append($("<div>").addClass("colorBox"));
	$.each(rvDataSets, function (index, rvds) {
		var targetLayer = rvds.getLayer($currentLayerName);
		if (RVcolor) {
			$($currentGroup).find(".colorBox").css("background", RVcolor);
			targetLayer.Color = RVcolor;
		} else {
			$($currentGroup).find(".colorBox").css("background", targetLayer.Color);
		}
	});

	//targetLayer.Color = $($currentGroup).find(".colorBox").css("background");

	//hide and show icon: eye
	var $visibleImgPath = "/static/ribovision/images/visible.png";
	var $invisibleImgPath = "/static/ribovision/images/invisible.png";
	$($currentGroup)
	.append($("<div>").addClass("checkBoxDIV").css({
			'float' : 'left',
			'padding-top' : 5,
			'margin-left' : 5,
			'width' : 20
		}).append($("<img class='visibilityCheckImg' value='visible' title='visibility' src='/static/ribovision/images/visible.png'/>").css({
				'width' : '24px',
				'height' : 'auto',
				'margin-top' : 3
			}).click(function () {
				var $type = this.getAttribute('value');
				if ($type == 'visible') {
					this.setAttribute('value', 'invisible');
					this.setAttribute('src', $invisibleImgPath);
					var newthis=this;
					$.each(rvDataSets, function (index, rvds) {
						var targetLayer = rvds.getLayer(newthis.parentNode.parentNode.getAttribute("name"));
						targetLayer.setVisibility("hidden");
					});
				} else if ($type == 'invisible') {
					this.setAttribute('value', 'visible');
					this.setAttribute('src', $visibleImgPath);
					var newthis=this;
					$.each(rvDataSets, function (index, rvds) {
						var targetLayer = rvds.getLayer(newthis.parentNode.parentNode.getAttribute("name"));
						targetLayer.setVisibility("visible");
					});
				}
			})));

	//adding raido button for selection
	$($currentGroup)
	.append($("<div>").addClass("radioDIV").css({
			'float' : 'left',
			'padding-top' : 5,
			'padding-left' : 5,
			'width' : 20
		}).append($("<input />").attr({
				type : 'radio',
				name : 'selectedRadioL',
				title : 'select layer'
			}).addClass("selectLayerRadioBtn").change(function (event) {
				$.each(rvDataSets, function (index, rvds) {
					rvds.selectLayer($(event.currentTarget).parent().parent().attr("name"));
				});
				drawNavLine();
			})));

	//raido button for telling 2D-3D mapping
	$($currentGroup)
	.append($("<div>").addClass("radioDIV2").css({
			'float' : 'left',
			'padding-top' : 5,
			'padding-left' : 5,
			'width' : 20
		}).append($("<input />").attr({
				type : 'radio',
				name : 'mappingRadioL',
				title : 'select which layer to map into 3D',
				disabled : 'disabled'
			}).addClass("mappingRadioBtn").change(function (event) {
				$.each(rvDataSets, function (index, rvds) {
					rvds.linkLayer($(event.currentTarget).parent().parent().attr("name"));
					var linkedLayer = rvds.getLinkedLayer();
					$("#LinkSection").find(".miniLayerName").attr("name", linkedLayer.LayerName);
					$("#LinkSection").find(".miniLayerName").text(linkedLayer.LayerName);
					$("#LinkSection").find(".miniLayerName").attr("title", $("#MiniLayer .miniLayerName[name=" + linkedLayer.LayerName + "]").attr("title"));
				});
				update3Dcolors();
			})));

	//adding accordion
	$($currentGroup)
	.append($("<h3>").addClass("layerName").css({
			'margin-left' : 83
		}).append($currentLayerName)
		.dblclick(function () { //double click to open dialog to change layer name and set color
			//watch out! $currentLayerName != the layer you are clicking right now!
			//open a dialog window for changing layer names

			//get the name of the layer that just got double clicked
			$dblClickedLayerName = this.innerHTML.substring(this.innerHTML.lastIndexOf("</span>") + 7);
			$dblClickedLayer = this;
			document.getElementById("layerNameInput").value = "";
			document.getElementById("layerNameInput").placeholder = $dblClickedLayerName;
			var flc = $.farbtastic("#layerColor");
			flc.setColor(rgb2hex($($dblClickedLayer).parent().find(".colorBox").css("background-color")));

			$("#LayerPreferenceDialog").dialog("open");
		}))
	.append($("<div>").addClass("layerContent").css({
			'margin-left' : 83
		}));

	//$count++;
	$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div>').html("<b>Layer Type:</b> " + Layer.Type).append($("<br>")));
	switch (Layer.Type) {
	case "circles":
		//Data Label Section
		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div>').html("<b>Loaded Data: </b><span name='DataLabel'>" + Layer.DataLabel + "</span>").append($("<br>")).append($("<br>")));

		//Circle buttons
		$("#LayerPanel div").first().next().find(".layerContent").append($('<div id="' + 'pr-' + key + '">').text("Draw Circles as:").append($("<br>")));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="filled' + key + '" id="' + 'pr-' + key + '-1' + '" value="filled" checked="checked"> <label for="' + 'pr-' + key + '-1' + ' ">filled</label>'));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="filled' + key + '" id="' + 'pr-' + key + '-2' + '" value="unfilled"> <label for="' + 'pr-' + key + '-2' + '">unfilled</label>'));
		//$('#' + 'pr-' + key).buttonset();
		if (Layer.Filled) {
			$('input[name="filled' + key + '"][value=filled]').prop("checked", true);
		} else {
			$('input[name="filled' + key + '"][value=unfilled]').prop("checked", true);
		}
		$('input[name="filled' + key + '"]').change(function (event) {
			var filled=$(this).attr("value");
			$.each(rvDataSets, function (index, rvds) {
				var targetLayer = rvds.getLayer($currentLayerName);
				if (filled === "filled") {
					targetLayer.Filled = true;
				} else {
					targetLayer.Filled = false;
				}
				rvds.refreshResiduesExpanded($currentLayerName);
			});
		});

		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div id="' + 'prs-' + key + '">').text("Circle Size:").append($("<br>")));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="size' + key + '" id="' + 'prs-' + key + '-1' + '" value="regular" checked="checked"> <label for="' + 'prs-' + key + '-1' + ' ">regular</label>'));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="size' + key + '" id="' + 'prs-' + key + '-2' + '" value="large"> <label for="' + 'prs-' + key + '-2' + '">large</label>'));
		//$('#' + 'pr-' + key).buttonset();
		if (Layer.ScaleFactor <= 1.0) {
			$('input[name="size' + key + '"][value=regular]').prop("checked", true);
		} else {
			$('input[name="size' + key + '"][value=large]').prop("checked", true);
		}

		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("disabled", false);

		$('input[name="size' + key + '"]').change(function (event) {
			var scalefactor=$(this).attr("value");
			$.each(rvDataSets, function (index, rvds) {
				var targetLayer = rvds.getLayer($currentLayerName);
				if (scalefactor === "regular") {
					targetLayer.ScaleFactor = 1.0;
				} else {
					targetLayer.ScaleFactor = 1.2;
				}
				$.each(rvDataSets, function (index, rvds) {
					rvds.refreshResiduesExpanded($currentLayerName);
				});
			});
		});
		$("#MiniOpenLayerBtn").after($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
			.text($currentLayerName).attr('name', $currentLayerName).droppable({
				drop : function (event, ui) {
					ProcessBubbleDrop(event, ui);
				}
			}));

		break;
	case "lines":
		//Data Label Section
		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div>').html("<b>Loaded Data: </b><span name='DataLabel'>" + Layer.DataLabel + "</span>").append($("<br>")).append($("<br>")));
		$("#LayerPanel div").first().next().find(".layerContent").append($('<div name="llm">').text("Color lines like:").append($("<br>")));
		//$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<select id="' + 'llm-' + key + 'lineselect' + '" name="' + 'llm-' + key + 'lineselect' + '" multiple="multiple"></select>'));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<select multiple="multiple"></select>'));

		//console.log($("#LayerPanel div").first().next().html());

		$("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select").multiselect({
			minWidth : 160,
			multiple : false,
			header : "Select a layer",
			noneSelectedText : "Select an Option",
			selectedList : 1
		});
		$("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select").multiselect().multiselectfilter();

		//Fill Menu
		var llm = $("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select")[0];
		//var SDList=rvDataSets[0].SpeciesEntry.StructDataMenu.split(";");
		//llm.options.length = 0;
		llm.options[0] = new Option("All Gray", "gray_lines");
		llm.options[0].setAttribute("selected", "selected");

		var rLayers = rvDataSets[0].getLayerByType("residues");
		var cLayers = rvDataSets[0].getLayerByType("circles");
		var sLayers = rvDataSets[0].getLayerByType("selected");
		$.each(rLayers, function (key, value) {
			llm.options[llm.options.length] = new Option(value.LayerName, value.LayerName);
		});
		$.each(cLayers, function (key, value) {
			llm.options[llm.options.length] = new Option(value.LayerName, value.LayerName);
		});
		$.each(sLayers, function (key, value) {
			llm.options[llm.options.length] = new Option(value.LayerName, value.LayerName);
		});
		$("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select").multiselect("refresh");

		$("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select").change(function (event) {
			var array_of_checked_values = $(event.currentTarget).parent().find("select").multiselect("getChecked").map(function () {
					return this.value;
				});
			if (array_of_checked_values[0] === "gray_lines") {
				$.each(rvDataSets, function (index, rvds) {
					rvds.drawBasePairs("lines", "gray_lines");
				});
				$(event.currentTarget).parent().parent().find(':radio').attr('disabled', 'disabled');
			} else if (array_of_checked_values[0] === "manual_coloring") {
				$.each(rvDataSets, function (index, rvds) {
					rvds.drawBasePairs("lines", "manual_coloring");
				});
				$(event.currentTarget).parent().parent().find(':radio').attr('disabled', 'disabled');
			} else {
				$.each(rvDataSets, function (index, rvds) {
					rvds.drawBasePairs("lines", rvDataSets[0].getLayer(array_of_checked_values[0]));
				});
				$(event.currentTarget).parent().parent().find(':radio').removeAttr('disabled');
			}
		});

		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<br>'));
		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div id="' + 'llmg-' + key + '">').text("Line Gradient Direction:").append($("<br>")));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="color_lines_gradient' + key + '" id="' + 'llmg-' + key + '-1' + '" value="Matched" checked="checked"> <label for="' + 'llmg-' + key + '-1' + ' ">Matched</label>'));
		$("#LayerPanel div").first().next().find(".layerContent").first().find("div").last().append($('<input type="radio" name="color_lines_gradient' + key + '" id="' + 'llmg-' + key + '-2' + '" value="Opposite"> <label for="' + 'llmg-' + key + '-2' + '">Opposite</label>'));
		$('[name="color_lines_gradient' + key + '"]').attr('disabled', 'disabled');
		$('input[name="color_lines_gradient' + key + '"]').change(function (event) {
			$.each(rvDataSets, function (index, rvds) {
				var targetLayer = rvds.getLayer($('input[name="color_lines_gradient' + key + '"]' + ':checked').parent().parent().parent().find('h3').text());
				targetLayer.ColorGradientMode = $(event.currentTarget).parent().parent().find('input:checked').val();
				rvds.drawBasePairs("lines");
			});
		});
		break;
	case "residues":
		//Data Label Section
		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div>').html("<b>Loaded Data: </b><span name='DataLabel'>" + Layer.DataLabel + "</span>").append($("<br>")).append($("<br>")));

		$("#LayerPanel div").first().next().find(".selectLayerRadioBtn").attr("checked", "checked");
		$.each(rvDataSets, function (index, rvds) {
			rvds.selectLayer($("#LayerPanel div").first().next().attr("name"));
			rvds.linkLayer($("#LayerPanel div").first().next().attr("name"));
		});
		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("disabled", false);
		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("checked", true);
		$("#MiniOpenLayerBtn").after($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
			.text($currentLayerName).attr('name', $currentLayerName).droppable({
				drop : function (event, ui) {
					ProcessBubbleDrop(event, ui);
				}
			}));
		break;
	case "selected":
		//$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div id="selectDiv">'))
		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("disabled", false);
		//$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div id="' + 'sele-' + key + '">'))
		//$("#selectDiv").html(text)
		break;
	case "contour":	
		//Data Label Section
		$("#LayerPanel div").first().next().find(".layerContent").first().append($('<div>').html("<b>Loaded Data: </b><span name='DataLabel'>" + Layer.DataLabel + "</span>").append($("<br>")).append($("<br>")));

		$("#LayerPanel div").first().next().find(".selectLayerRadioBtn").attr("checked", "checked");
		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("disabled", false);
		$("#LayerPanel div").first().next().find(".radioDIV2").find('input').prop("checked", true);
		$("#MiniOpenLayerBtn").after($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
			.text($currentLayerName).attr('name', $currentLayerName).droppable({
				drop : function (event, ui) {
					ProcessBubbleDrop(event, ui);
				}
			}));
		break;
	default:
		break;
	}
	$("#LayerPanel div").first().next().find(".visibilityCheckImg").attr("value", "visible");
}
// Refresh Menu
function RefreshLayerMenu() {
	//Accordion that support multiple sections open
	$("#LayerPanel").multiAccordion();
	var targetLayers = rvDataSets[0].getLayerByType("lines");
	$.each(targetLayers, function (index, value) {
		var $selectbox = $(".oneLayerGroup[name=" + targetLayers[0].LayerName + "]").find(".layerContent").find("div[name=llm]").find("select")[0];
		//Fill Menu
		//var llm = $("#LayerPanel div").first().next().find(".layerContent").find("[name='llm']").find("select")[0];
		//var SDList=rvDataSets[0].SpeciesEntry.StructDataMenu.split(";");
		$selectbox.options.length = 0;
		$selectbox.options[0] = new Option("All Gray", "gray_lines");
		$selectbox.options[0].setAttribute("selected", "selected");
		$selectbox.options[1] = new Option("Manual Coloring", "manual_coloring");

		var dsLayers = rvDataSets[0].getLayerByType(["residues", "circles","contour", "selected"]);
		$.each(dsLayers.reverse(), function (key, value) {
			$selectbox.options[$selectbox.options.length] = new Option(value.LayerName, value.LayerName);
		});
		$(".oneLayerGroup[name=" + targetLayers[0].LayerName + "]").find(".layerContent").find("div[name=llm]").find("select").multiselect("refresh");
		/*
		// Refresh MiniLayers
		$(".miniLayerName").remove();
		var dLayers = rvDataSets[0].getLayerByType(["residues","circles"]);
		$.each(dLayers, function (key, value) {
		$("#MiniLayerLabel").after($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
		.text(value.LayerName).attr('name',value.LayerName).droppable({
		drop: function (event,ui) {
		ProcessBubbleDrop(event,ui);
		}
		}));
		});
		//Refresh Linked MiniLayer
		var linkedLayer = rvDataSets[0].getLinkedLayer();
		$("#LinkSection").append($('<h3 class="miniLayerName ui-helper-reset ui-corner-all ui-state-default ui-corner-bottom ">')
		.text(linkedLayer.LayerName).attr('name',linkedLayer.LayerName).droppable({
		drop: function (event,ui) {
		ProcessBubbleDrop(event,ui);
		}
		}));*/
	});
}

function ProcessBubbleDrop(event, ui) {
	if ($(ui.draggable[0]).hasClass("dataBubble")) {
		var targetLayer = rvDataSets[0].getLayer($(event.target).attr("name"));
		targetLayer.DataLabel = $(ui.draggable[0]).text();
		$("[name=" + targetLayer.LayerName + "]").find(".layerContent").find("span[name=DataLabel]").text(targetLayer.DataLabel);
		var ColName = $(ui.draggable[0]).attr("name").match(/[^\'\\,]+/);
		var result = $.grep(rvDataSets[0].DataDescriptions, function (e) {
				return e.ColName === ColName[0];
			});
		if (result[0]) {
			$("#MiniLayer .miniLayerName[name=" + targetLayer.LayerName + "]").attr("title", targetLayer.DataLabel + ": " + result[0].Description);
			$("#LinkSection .miniLayerName[name=" + targetLayer.LayerName + "]").attr("title", targetLayer.DataLabel + ": " + result[0].Description);
			$(ui.draggable[0]).parent().parent().find(".ManualLink").attr("href", "./Documentation/" + result[0].HelpLink + ".html");
		} else {
			$(".miniLayerName[name=" + targetLayer.LayerName + "]").attr("title", "Data Description is missing.");
			$(ui.draggable[0]).parent().parent().find(".ManualLink").attr("href", "./Documentation");
		}
		// For now, the LayerName and DataLabels are identical for both subunits. Only need to loop in ProcessBubble
		ProcessBubble($(ui.draggable[0]), targetLayer.LayerName);
		$(".oneLayerGroup[name='" + targetLayer.LayerName + "']").find(".selectLayerRadioBtn").prop("checked", true);
		$(".oneLayerGroup[name='" + targetLayer.LayerName + "']").find(".selectLayerRadioBtn").trigger("change");
		$(".oneLayerGroup[name='" + targetLayer.LayerName + "']").find(".mappingRadioBtn").prop("checked", true);
		$(".oneLayerGroup[name='" + targetLayer.LayerName + "']").find(".mappingRadioBtn").trigger("change");
	}
}
