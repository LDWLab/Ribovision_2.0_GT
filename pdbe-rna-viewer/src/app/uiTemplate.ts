import { UiActionsService } from './uiActions';
import { PluginOptions, ApiData, BanNameHelper } from './data';
import { CustomEvents } from './customEvents';


export class UiTemplateService {
    selectSections_RV1 = (window as any).selectSections_RV1;
    aaPropertyConstants = (window as any).aaPropertyConstants;
    mapped_aa_properties = (window as any).mapped_aa_properties;
    mapped_aa_properties3D = (window as any).mapped_aa_properties3D;
    aaColorData = (window as any).aaColorData;
    parsePVData = (window as any).parsePVData;
    getEntropyAnnotations = (window as any).getEntropyAnnotations;
    getTWCAnnotations = (window as any).getTWCAnnotations;
    getCustomAnnotations = (window as any).getCustomAnnotations;
    getAssociatedAnnotations = (window as any).getAssociatedAnnotations;
    getHelicalAnnotations = (window as any).getHelicalAnnotations;
    getPhaseAnnotations = (window as any).getPhaseAnnotations;
    getExpansionAnnotations = (window as any).getExpansionAnnotations;
    getAnnotationArray = (window as any).getAnnotationArray;
    rv3VUEcomponent = (window as any).vm;
    private containerElement: HTMLElement;
    private pluginOptions: PluginOptions;
    private uiActionsService: UiActionsService;
    private apiData: ApiData | undefined;


    private locations: Map<any, number[]> = new Map();
    menuStyle = 'position:relative;z-index:10;height:7%;line-height:7%;background-color:#696969;padding: 0 10px;font-size:16px; color: #efefef;display:block;';
    domainTypes: any[];
    selectedDomain: string;
    pathStrs: string[] = [];
    nucleotideStrs: string[] = [];
    circleStrs: string[] = [];
    banNameMap: Map<string, JSON | undefined> = new Map();
    baseStrs: Map<string, [boolean, string[]]> = new Map();
    nestedBaseStrs: Map<string, [boolean, string[]]> = new Map();
    displayBaseStrs: string;
    displayNestedBaseStrs: string;
    basePairIDs: any[] = [];
    mappingValue: string = '';
    showAllNucleotides: boolean = false;
    mouseOverMap = new Map<string, any>()
    eventMap = new Map<string, any>()
    toolTips: Map<number, string> = new Map();
    selectSectionsTest = [];
    defaultColours = {
        domainSelection: 'rgb(255,0,0)',
        mouseOver: 'rgb(105,105,105)',
        borderColor: 'rgb(0,0,0)',
        qualityGreen: 'rgb(0,182.85714285714286,0)',
        qualityRed: 'rgb(291.42857142857144,0,0)',
        qualityYellow: 'rgb(364.2857142857143,364.2857142857143,75.71428571428572)',
        qualityRiboVision: "rgb(203,203,203)",
        qualityOrange: 'rgb(291.42857142857144,121.42857142857143,0)',
        qualityBlank: 'rgb(255,255,255)'
    }
    mapped_chains = new Set<string>()
    mapped_modifications = new Set<string>()
    constructor(containerElement: HTMLElement, pluginOptions: PluginOptions, apiData: ApiData | undefined) {
        this.containerElement = containerElement;
        this.pluginOptions = pluginOptions;
        this.apiData = apiData;
        this.uiActionsService = new UiActionsService(this.pluginOptions.pdbId);
    }

    render(apiData: ApiData, FR3DData: any, FR3DNestedData: any, BanName: any, instance?: any) {
        this.containerElement.innerHTML =
            `<div class="pdb-rna-view-container pdb-rna-view-container-${this.pluginOptions.pdbId}">
            ${this.svgTemplate(apiData, FR3DData, FR3DNestedData)}
            ${this.title()}
            ${this.tooltip()}
            ${this.actionButtons()}
        </div>
        <div style="${this.menuStyle}">
            <div class="menuOptions" style="width:95%;float:center;display:inline-block;">
                <form>
                <input type="checkbox" id="nestedBP" name="nestedBP" style="width:5%; display: inline-block; float: left;>
                <label for="nestedBP" style="width:20%; display: inline-block; float: left;"> Only nested BPs</label>
                <select class="menuSelectbox" style="width:20%; display: inline-block;"><option value="">Nucleotides</option></select>
                <select class="mappingSelectbox" style="width:20%; display: inline-block;"><option value="">Mapping</option></select>
                <div class="multiselect">
                    <div class="selectBox" onclick="UiActionsService.showCheckboxes()" style="display:inline-block;float:right;width:25%">
                        <select id="basePairingSelectElement">
                            <option>Base Pairings</option>
                        </select>
                        <div class="overSelect"></div>
                    </div>
                    <div id="checkboxes">
                    </div>
                </div>
                </form>
            </div>
        </div>`;
        //this.getJSON(apiData)
        //this.fixOverlaps(apiData)
        this.createModeDropdown();
        this.createBPDropdown();
        this.addAnnotationHandlers();
        this.uiActionsService.applyButtonActions();
        this.addEvents(apiData, BanName);

        // Save button event listeners
        const svgSaveButton = document.getElementById(`rnaTopologySaveSVG-${this.pluginOptions.pdbId}`);
        const jsonSaveButton = document.getElementById(`rnaTopologySaveJSON-${this.pluginOptions.pdbId}`);

        if (svgSaveButton) {
            svgSaveButton.addEventListener("click", this.saveSVG.bind(this));
        } else {
            console.warn('SVG save button not found');
        }

        if (jsonSaveButton) {
            jsonSaveButton.addEventListener("click", this.saveJSON.bind(this));
        } else {
            console.warn('JSON save button not found');
        }

        // <any>document.querySelector(".saveSVG")!.addEventListener("click", this.saveSVG.bind(this));

        this.getAnnotationFromRibovision(this.mapped_aa_properties, this.mapped_aa_properties3D);

        var el = document.getElementById("topview");
        document.getElementById("TopologyFSCR")?.addEventListener("click", (event) => {
            document.fullscreenElement ? document.exitFullscreen() : el?.requestFullscreen();
        })

        if (instance) {
            CustomEvents.subscribeToComponentEvents(instance)
        }
        //this.rv3VUEcomponent.topology_loaded=true;
    }

    fixOverlaps(apiData: ApiData) {
        var svgEle: any = (<any>document.querySelector(`svg.rnaTopoSvg.rnaTopoSvg_${this.pluginOptions.pdbId}`))
        var xAdjust: number = 0;
        var yAdjust: number = 0;
        apiData.sequence.split('').forEach((char: string, i: number) => {
            if (i === 0 || i === (apiData.sequence.length)) { return }
            let ele = svgEle!.getElementsByClassName(`rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[i]}`)[0]
            let nextEle = svgEle!.getElementsByClassName(`rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[i + 1]}`)[0]
            nextEle.setAttribute("y", Number(nextEle.getAttribute("y")) + yAdjust)
            nextEle.setAttribute("x", Number(nextEle.getAttribute("x")) + xAdjust)
            let distances = this.getIntersections(ele.getBBox(), nextEle.getBBox())
            if (distances) {
                var lowest = 0;
                for (var j = 1; j < distances.length; j++) {
                    if (distances[j] < distances[lowest]) lowest = j;
                }
                if (lowest === 1) {
                    nextEle.setAttribute("y", Number(nextEle.getAttribute("y")) - Number(distances[1]));
                    yAdjust = yAdjust - Number(distances[1])
                } else if (lowest === 3) {
                    nextEle.setAttribute("y", Number(nextEle.getAttribute("y")) + Number(distances[3]));
                    yAdjust = yAdjust + Number(distances[3])
                } else if (lowest === 0) {
                    nextEle.setAttribute("x", Number(nextEle.getAttribute("x")) - Number(distances[0]));
                    xAdjust = xAdjust - Number(distances[0])
                } else if (lowest === 2) {
                    nextEle.setAttribute("x", Number(nextEle.getAttribute("x")) + Number(distances[2]));
                    xAdjust = xAdjust + Number(distances[2])
                }
            }
        });
    }
    drawCircle = function (pdbId: string, i: number, color: string) {
        const circle = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)!.getElementsByClassName(`circle_${pdbId}_${i}`)[0]
        //const nucleotide = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)!.getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${i}`)[0]
        //const BBox = nucleotide.getBBox()
        //const nx = `${BBox.x + BBox.width/2}`
        //const ny = `${BBox.y + BBox.height/2}`
        //circle.setAttribute("cx", nx)
        //circle.setAttribute("cy", ny)
        circle.setAttribute("stroke", `${color}`);
        circle.setAttribute("fill", `${color}`);
        circle.style.display = "block";
    }
    colorMap = () => {
        // const tempDomain = this.domainTypes.find(item => item.label === 'Shannon entropy');

        let longest = null;
        for (const domainType of this.domainTypes) {
            if (domainType.data && domainType.data.length > (longest ? longest.data.length : 0)) {
                longest = domainType;
            }
        }
        const allIndices = new Set();
        longest.data.forEach((val: any) => {
            if (val != undefined && val.start != undefined) {
                allIndices.add(val.start);
            }
        });

        // if (tempDomain) {
        // }
        // tempDomain.data.forEach((val:any, i:number) => {
        allIndices.forEach((index: any) => {
            // if(val != undefined && val.start != undefined) {
            UiActionsService.colorNucleotide(this.pluginOptions.pdbId, index, 'rgb(0,0,0)', undefined, this.mappingValue);
            let cPath = `circle_${this.pluginOptions.pdbId}_${index}`;
            let nPath = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${index}`;
            if (document.getElementsByClassName(nPath).length > 0) {
                (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(nPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(nPath)[0].getAttribute('onmouseout')!.split(';')[0]);
            }
            if (document.getElementsByClassName(cPath).length > 0) {
                (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(cPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(cPath)[0].getAttribute('onmouseout')!.split(';')[0]);
            }
            // }
        });
        const selectBoxEle: any = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        if (selectedValue) {
            const selectedDomain = this.domainTypes[selectedValue];
            if (selectedDomain.data) {
                selectedDomain.data.forEach((val: any, i: number) => {
                    if (val != undefined && val.color != undefined) {
                        if (val.color.indexOf('NaN') < 0) {
                            //this.drawCircle(this.pluginOptions.pdbId, val.start, val.color);
                            //UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, val.color);
                            UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, val.color, undefined, this.mappingValue);

                            if (['helix', 'aes', 'phase'].includes(selectedDomain.label.toLowerCase()) && val.start in this.rv3VUEcomponent.struct_to_alignment_mapping) {
                                let align_index = this.rv3VUEcomponent.struct_to_alignment_mapping[val.start]
                                let circlePath = `circle_${this.pluginOptions.pdbId}_${val.start}`;
                                let path = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${val.start}`;
                                if (align_index in this.rv3VUEcomponent.associatedDataCache) {
                                    for (let { type, value } of this.rv3VUEcomponent.associatedDataCache[align_index]) {
                                        if (type == selectedDomain.label) {
                                            let tooltip = type + ' ' + value;
                                            this.toolTips.set(val.start, tooltip);
                                            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute(
                                                'onmouseover',
                                                document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] +
                                                `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${val.color}','${val.color}')`
                                            );
                                            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute(
                                                'onmouseout',
                                                document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] +
                                                `;UiActionsService.hideTooltip('${circlePath}')`
                                            );
                                            if ((<HTMLElement>document.getElementsByClassName(path)[0])) {
                                                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute(
                                                    'onmouseover',
                                                    document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] +
                                                    `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${val.color}', '${val.color}')`
                                                );
                                                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute(
                                                    'onmouseout',
                                                    document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] +
                                                    `;UiActionsService.hideTooltip('${path}')`
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        } /* else {
            
            const tempDomain = this.domainTypes[1];
            if (tempDomain.data) {
                tempDomain.data.forEach((val:any, i:number) => {
                    if(val != undefined && val.start != undefined) {
                        UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, 'rgb(0,0,0)', undefined, this.mappingValue);
                    }
                });
            }
            this.colorMapContacts()
            this.colorMapModifications()
            
        }*/
    }
    colorMapHelper = () => {
        var mappingDropdown = (<HTMLSelectElement>this.containerElement.querySelector<HTMLElement>('.mappingSelectbox'));
        var num: number = +mappingDropdown.value;
        (<HTMLInputElement>document.getElementById('selectColorMappingProps')).value = mappingDropdown!.options[num].text;
        this.rv3VUEcomponent.selected_property = mappingDropdown!.options[num].text;
    }
    colorMapContacts = () => {

        if (this.rv3VUEcomponent.pchainid.length > 0 && this.rv3VUEcomponent.selected_property != "Select data") {
            this.rv3VUEcomponent.selected_property = "Select data"
        }
        this.mapped_chains.forEach((val) => {
            if (!this.rv3VUEcomponent.pchainid.includes(val)) {
                for (var i in this.rv3VUEcomponent.protein_contacts[val]) {
                    UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.protein_contacts[val][i], '#323232', undefined, this.mappingValue);
                    var nPath = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.protein_contacts[val][i]}`
                    var cPath = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.protein_contacts[val][i]}`
                    if (document.getElementsByClassName(nPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(nPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(nPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                        //document.getElementsByClassName(nPath)[0].removeEventListener('mouseover', this.mouseOverMap.get(nPath))
                        //document.getElementsByClassName(nPath)[0].removeEventListener('mouseout',this.eventMap.get(nPath))
                    }
                    if (document.getElementsByClassName(cPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(cPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(cPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                        //document.getElementsByClassName(cPath)[0].removeEventListener('mouseover', this.mouseOverMap.get(cPath))
                        //document.getElementsByClassName(cPath)[0].removeEventListener('mouseout',this.eventMap.get(cPath))
                    }
                    /*
                    var circle = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(cPath)[0];
                    if(circle) {
                        circle.style.display="none";
                    }*/
                }
                this.mapped_chains.delete(val)
            }
        });
        for (let val in this.rv3VUEcomponent.pchainid) {
            var chain = this.rv3VUEcomponent.pchainid[val];
            for (var i in this.rv3VUEcomponent.protein_contacts[chain]) {
                UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.protein_contacts[chain][i], this.rv3VUEcomponent.proteinColorMap.get(chain), undefined, this.mappingValue);
            }
            this.mapped_chains.add(chain)
        }
        this.addEvents(this.apiData!, BanNameHelper.getBanName(this.pluginOptions.pdbId, ''))
    }
    colorMapModifications = () => {
        if (this.rv3VUEcomponent.modifications.length > 0 && this.rv3VUEcomponent.selected_property != "Select data") {
            this.rv3VUEcomponent.selected_property = "Select data"
        }
        this.mapped_modifications.forEach((val) => {
            if (!this.rv3VUEcomponent.modifications.includes(val)) {
                for (var i in this.rv3VUEcomponent.modified_residues.get(val)) {
                    UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.modified_residues.get(val)[i], '#323232', undefined, this.mappingValue);
                    var nPath = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues.get(val)[i]}`
                    var cPath = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues.get(val)[i]}`
                    if (document.getElementsByClassName(nPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(nPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(nPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                        //document.getElementsByClassName(nPath)[0].removeEventListener('mouseover', this.mouseOverMap.get(nPath))
                        //document.getElementsByClassName(nPath)[0].removeEventListener('mouseout',this.eventMap.get(nPath))
                    }
                    if (document.getElementsByClassName(cPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(cPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(cPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                        //document.getElementsByClassName(cPath)[0].removeEventListener('mouseover', this.mouseOverMap.get(cPath))
                        //document.getElementsByClassName(cPath)[0].removeEventListener('mouseout',this.eventMap.get(cPath))
                    }
                    /*
                    var circle = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(cPath)[0];
                    if(circle) {
                        circle.style.display="none";
                    } */
                }
                this.mapped_modifications.delete(val)
            }
        });
        for (let val in this.rv3VUEcomponent.modifications) {
            var mod = this.rv3VUEcomponent.modifications[val];
            for (var i in this.rv3VUEcomponent.modified_residues.get(mod)) {
                UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.modified_residues.get(mod)[i], this.rv3VUEcomponent.modifiedColorMap.get(mod), undefined, this.mappingValue);
            }
            this.mapped_modifications.add(mod)
        }
        this.addEvents(this.apiData!)
    }
    changeBP(val: string, e: Event | undefined = undefined, flag: boolean | undefined = undefined) {
        //Figure out what flag is doing
        this.displayBaseStrs = '';
        this.displayNestedBaseStrs = '';
        const allBP = this.containerElement.querySelector<HTMLInputElement>('#Checkbox_All')!.checked
        if ((val == 'All')) {
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                this.baseStrs.set(key, [allBP, value[1]]);
                this.nestedBaseStrs.set(key, [allBP, this.nestedBaseStrs.get(key)![1]]);
                (<HTMLInputElement>document.getElementById(`Checkbox_${key}`))!.checked = allBP
                if (allBP) {
                    this.displayBaseStrs += value[1].join('');
                    this.displayNestedBaseStrs += this.nestedBaseStrs.get(key)![1].join('');
                }
            });
        } else {
            //console.log(flag)
            //if(flag ?? this.baseStrs.get(val)![0]) {
            if (this.baseStrs.get(val)![0]) {
                this.baseStrs.set(val, [false, this.baseStrs.get(val)![1]]);
                this.nestedBaseStrs.set(val, [false, this.nestedBaseStrs.get(val)![1]]);
            } else {
                this.baseStrs.set(val, [true, this.baseStrs.get(val)![1]]);
                this.nestedBaseStrs.set(val, [true, this.nestedBaseStrs.get(val)![1]]);
            }
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                if (value[0]) {
                    this.displayBaseStrs += value[1].join('');
                    this.displayNestedBaseStrs += this.nestedBaseStrs.get(key)![1].join('')
                }
            });
        }
        this.PathOrNucleotide();
    }
    createBPDropdown() {
        if (this.baseStrs.size > 0) {
            let optionList = '<table><tr><td><label for = "Checkbox_All"><input type="checkbox" id="Checkbox_All" /> All</label></td>';
            var i = 1;
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                if (i % 2 == 0) {
                    optionList = `${optionList}<tr>`
                }
                if (key == 'cWW') {
                    optionList = `${optionList}<td><label for = "Checkbox_${key}"><input type="checkbox" id="Checkbox_${key}" checked = true/> ${key}</label></td>`;
                } else {
                    optionList = `${optionList}<td><label for = "Checkbox_${key}"><input type="checkbox" id="Checkbox_${key}"/> ${key}</label></td>`;
                }
                if (i % 2 == 1) {
                    optionList = `${optionList}</tr>`
                }
                i += 1
            });
            optionList = `${optionList}</table>`
            const selectBoxEle = document.getElementById('checkboxes');
            selectBoxEle!.innerHTML = optionList;
            document.getElementById(`Checkbox_All`)?.addEventListener("change", this.changeBP.bind(this, "All"));
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                document.getElementById(`Checkbox_${key}`)?.addEventListener("change", this.changeBP.bind(this, key));
            });
            //selectBoxEle!.addEventListener("change", this.ChangeBP.bind(this));
        }
    }
    createModeDropdown() {
        let optionList = `<option value="0">Nucleotides</option><option value="1">Path</option><option value="2">Circle</option>`;
        const selectBoxEle = this.containerElement.querySelector<HTMLElement>('.menuSelectbox');
        selectBoxEle!.innerHTML = optionList;
        selectBoxEle!.addEventListener("change", this.PathOrNucleotide.bind(this));
        const nestedBP = this.containerElement.querySelector<HTMLElement>('#nestedBP');
        nestedBP!.addEventListener("change", this.PathOrNucleotide.bind(this));
    }
    createDomainDropdown() {
        this.selectedDomain = this.domainTypes[0];
        if (this.domainTypes.length > 0) {

            let optionList = '';
            this.domainTypes.forEach((opt: any, i: number) => {
                optionList = `${optionList}<option value="${i}">${opt.label}</option>`;
            });

            const selectBoxEle = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
            selectBoxEle!.innerHTML = optionList;
            selectBoxEle!.addEventListener("change", this.colorMapHelper.bind(this));

            //selectBoxEle!.addEventListener("change", this.updateProperty.bind(this));

            //const resetIconEle = this.containerElement.querySelector('.resetIcon');
            //resetIconEle.addEventListener("click", this.resetDisplay.bind(this));

        } else {
            this.containerElement!.querySelector<HTMLElement>('.menuOptions')!.style.display = 'none';
        }
    }
    create2D3DAnnotations(name: string, residueDetails: any,
        TWCrgbMap: Map<number, any>, TWCData: Map<number, string>, mapped_aa_properties: Map<string, Array<Array<number>>>, chain_start: number, chain_end: number) {
        const _this = this;
        TWCData.forEach(function (value, index) {
            if (chain_start <= index && index <= chain_end) {
                let rgb_color = TWCrgbMap.get(index);
                (window as any).selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.pluginOptions.entityId,
                    //start_residue_number: index, 
                    //end_residue_number: index,
                    residue_number: index,
                    color: rgb_color[1],
                    sideChain: false,
                });

                _this.defaultColours.qualityRiboVision = '#' + String(rgb_color[0].join(''));

                var colors = "rgb(" + String(rgb_color[0].join(',')) + ")"
                //_this.drawValidationShape(index, "circle", _this.defaultColours.qualityRiboVision);
                residueDetails.push({ //2d
                    start: index,
                    end: index,
                    color: colors,
                    tooltipMsg: Number.parseFloat(value).toPrecision(3),
                    tooltipPosition: "prefix"
                });

                //_this.drawValidationShape(index, "circle", colors);
            }
        })
        /*
         if (TWCData.size < mapped_aa_properties.get("Shannon entropy")!.length) {
             console.log("2D3D2", mapped_aa_properties);
             for(var i = TWCData.size - 1; i < this.mapped_aa_properties.get("Shannon entropy").length; i++) {
                 (window as any).selectSections_RV1.get(name).push({ //3d
                     entity_id: _this.pluginOptions.entityId,
                     //start_residue_number: i, 
                     //end_residue_number: i,
                     residue_number: i,
                     color: {r:255, g:255, b:255},
                     sideChain: false,
                 });
             }
         }*/
        /*
        if (TWCData.size < mapped_aa_properties.get("TwinCons")!.length) {
            for(var i = TWCData.size - 1; i < this.mapped_aa_properties.get("TwinCons").length; i++) {
                (window as any).selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.pluginOptions.entityId,
                    //start_residue_number: i, 
                    //end_residue_number: i,
                    residue_number: i,
                    color: {r:255, g:255, b:255},
                    sideChain: false,
                });
            }
        }*/
        return residueDetails;
    }
    getAnnotationFromRibovision(mapped_aa_properties: Map<string, Array<Array<number>>>, mapped_3D_properties?: Map<string, Array<Array<number>>>) {
        const _this = this;
        const start = this.apiData ? this.apiData.label_seq_ids[1] : 0
        //const end = this.apiData?this.apiData.label_seq_ids[this.apiData.label_seq_ids.length - 2]:0

        if (typeof _this.domainTypes == 'undefined') {
            _this.domainTypes = [{
                label: 'Select data',
                data: null
            }];
        }
        this.selectSections_RV1 = (window as any).selectSections_RV1
        this.selectSections_RV1.set("Select data", [{
            entity_id: _this.pluginOptions.entityId,
            focus: true
        } /*{entity_id: _this.pluginOptions.entityId,
                                start_residue_number: 0, 
                                end_residue_number: 1000,
                                color: {r: 255, g: 255, b: 255},
                                sideChain: false,pchainid
                            }*/])
        this.aaPropertyConstants = (window as any).aaPropertyConstants
        this.aaColorData = (window as any).aaColorData
        if (mapped_aa_properties) {
            mapped_aa_properties.forEach((value, index) => {
                let residueDetails: any = [{
                    //start: chainRange.start,
                    //end: chainRange.end,
                    //color: _this.defaultColours.qualityBlank,
                    //tooltipMsg: 'No data for '
                }];
                let name = index;
                let separatedData = value;
                this.selectSections_RV1.set(name, []);

                let min = Math.min(...this.aaPropertyConstants.get(name));
                let max = Math.max(...this.aaPropertyConstants.get(name));
                let colormapArray = this.aaColorData.get(name);
                let data3D;
                data3D = separatedData;
                if (mapped_3D_properties) {
                    if (mapped_3D_properties.get(name)) {
                        data3D = mapped_3D_properties.get(name);
                    }
                }
                if (name == "Shannon entropy") {
                    this.getEntropyAnnotations(data3D, min, max, this.pluginOptions.chainId);
                    //console.log('name_SE', name,  this.getEntropyAnnotations(separatedData, min, max, this.pluginOptions.chainId));
                };
                if (name == "TwinCons") {
                    this.getTWCAnnotations(data3D, min, max, this.pluginOptions.chainId);
                    //console.log('name_TWC', name, this.getTWCAnnotations(separatedData, min, max, this.pluginOptions.chainId));
                };

                if (name == "Custom Data") {

                    this.getCustomAnnotations(data3D, min, max, this.pluginOptions.chainId);
                    //console.log('name_CD', name,  this.getEntropyAnnotations(separatedData, min, max, this.pluginOptions.chainId));
                };
                if (name == "Associated Data1") {

                    this.getAssociatedAnnotations(data3D, min, max, this.pluginOptions.chainId);

                };
                if (name == "Helix" || name == "helix") {
                    this.getHelicalAnnotations(data3D, min, max, this.pluginOptions.chainId);

                };

                if (name == "Phase" || name == 'phase') {

                    this.getPhaseAnnotations(data3D, min, max, this.pluginOptions.chainId);

                };
                if (name == "AES" || name == 'aes') {
                    //console.log('name_AES', name, this.getExpansionAnnotations(separatedData, min, max, this.pluginOptions.chainId))
                    this.getExpansionAnnotations(data3D, min, max, this.pluginOptions.chainId);

                };
                const [TWCrgbMap, TWCData, TWCrgbMap3D, TWCData3D] = this.parsePVData(separatedData, min, max, colormapArray, null, data3D);
                // const [TWCrgbMap3D, TWCData3D] = this.parsePVData(data3D, min, max, colormapArray);
                const TWCData_keys = TWCData.keys();

                // just to get rid of build error
                // TWCrgbMap;
                //const last_item  = 0;
                let mapLastValue;
                let i;
                for (i = 0; i < TWCData.size; i += 1) {
                    mapLastValue = TWCData_keys.next().value
                }

                //console.log(mapLastValue);
                //console.log('RNAViewer_TWCData',TWCData, TWCData.size, TWCData[Symbol.iterator](),TWCData_keys);
                this.selectSections_RV1.get(name).push({ entity_id: _this.pluginOptions.entityId, focus: true });
                //const end = TWCData.size;
                const end = mapLastValue;
                if (void 0 !== TWCData3D) {
                    residueDetails = _this.create2D3DAnnotations(name, residueDetails,
                        TWCrgbMap3D, TWCData3D, mapped_aa_properties,
                        start, end);
                    if (0 < residueDetails.length) {
                        var current = _this.domainTypes.filter(order => (order.label === name))[0];
                        if (current && current != null) {
                            current.data = residueDetails;
                        } else {
                            _this.domainTypes.push({
                                label: name,
                                data: residueDetails
                            })
                        }
                    }
                }
                if (void 0 !== TWCData) {
                    residueDetails = _this.create2D3DAnnotations(name, residueDetails,
                        TWCrgbMap, TWCData, mapped_aa_properties,
                        start, end);
                    if (0 < residueDetails.length) {
                        var current = _this.domainTypes.filter(order => (order.label === name))[0];
                        if (current && current != null) {
                            current.data = residueDetails;
                        } else {
                            _this.domainTypes.push({
                                label: name,
                                data: residueDetails
                            })
                        }
                    }
                }
            });
            this.createDomainDropdown()
            this.rv3VUEcomponent.topology_loaded = true;
        }
        else {
            //catch block
        };
    }
    getIntersections(obj1: any, obj2: any) {
        var left1 = obj1.x
        var right1 = obj1.x + obj1.width
        var top1 = obj1.y
        var bottom1 = obj1.y + obj1.height
        var left2 = obj2.x
        var right2 = obj2.x + obj2.width
        var top2 = obj2.y
        var bottom2 = obj2.y + obj2.height
        if (left1 >= right2 || top1 >= bottom2 ||
            right1 <= left2 || bottom1 <= top2) {
            return false
        } else {
            return [right2 - left1, bottom2 - top1, right1 - left2, bottom1 - top2]
        }
    }

    calculateFontSize(apiData: ApiData) {
        let xVals: number[] = [];
        let yVals: number[] = [];
        let dist: number[] = [];
        const lastPathIndex = apiData.svg_paths.length - 1;
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            if (recordIndex === 0 || recordIndex === lastPathIndex) return;
            let pathStrParsed: string[] = pathStr.split('M').join(',').split(',')
            xVals[recordIndex] = (Number(pathStrParsed[1]) + Number(pathStrParsed[3])) / 2
            yVals[recordIndex] = (Number(pathStrParsed[2]) + Number(pathStrParsed[4])) / 2
            if (recordIndex > 1) {
                let xDiff = xVals[recordIndex] - xVals[recordIndex - 1]
                let yDiff = yVals[recordIndex] - yVals[recordIndex - 1]
                dist[recordIndex] = Math.pow((Math.pow(yDiff, 2) + Math.pow(xDiff, 2)), 0.5)
            }
        });
        var sortedDist: number[] = dist.sort((a, b) => { return a - b; })
        return 0.9 * sortedDist[Math.floor(sortedDist.length * 0.05)]
    }

    download(data: any, filename: string, type: string) {
        var file = new Blob([data], { type: type });
        var a = document.createElement("a"),
            url = URL.createObjectURL(file);
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(function () {
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);
        }, 0);
    }
    private calcBaseStrs(baseStr: any, baseMap: any, font_size: number) {
        let start: number = +baseStr.seq_id1
        let end: number = +baseStr.seq_id2
        if (baseStr && start && end) {
            let type: string = baseStr.bp
            let pathID: string = `rnaviewBP rnaviewBP_${this.pluginOptions.pdbId}_${this.pluginOptions.chainId} ${type}_${start}_${end}`
            let n1: string = baseStr.nt1
            let n2: string = baseStr.nt2
            if (type == 'cSH' || type == 'tSH' || type == 'cSW' || type == 'tSW' || type == 'cHW' || type == 'tHW') {
                let temp = end
                end = start
                start = temp
                let temp2 = n1
                n1 = n2
                n2 = temp2
                type = type.charAt(0) + type.slice(-2).split('').reverse().join('')
            }
            let x1 = this.locations.get(start)![0] + font_size / 2.5
            let x2 = this.locations.get(end)![0] + font_size / 2.5
            let y1 = this.locations.get(start)![1] - font_size / 2.5
            let y2 = this.locations.get(end)![1] - font_size / 2.5
            let distance = Math.pow(Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2), 0.5)
            let x1_prime = UiTemplateService.linearlyInterpolate(x1, x2, font_size / distance)
            let y1_prime = UiTemplateService.linearlyInterpolate(y1, y2, font_size / distance)
            let x2_prime = UiTemplateService.linearlyInterpolate(x1, x2, 1 - font_size / distance)
            let y2_prime = UiTemplateService.linearlyInterpolate(y1, y2, 1 - font_size / distance)
            let stroke = "#ccc"
            if (type.charAt(0) == 't') {
                var fill = "none"
            } else {
                var fill = "#ccc"
            }
            let xm = (x1_prime + x2_prime) / 2
            let ym = (y1_prime + y2_prime) / 2
            let distance2 = distance - 2 * font_size
            let height = font_size / 1.5
            const defaultAction = `<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${n1}${start} - ${n2}${end}; ${type}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"`
            var phi = 270 + Math.atan2((y1 - y2), (x1 - x2)) * 180 / Math.PI
            if (type == 'cWW') {
                if (n1 == 'G' && n2 == 'U' || n1 == 'U' && n2 == 'G') {
                    baseMap.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${n1}${start} - ${n2}${end}; ${type}', '${pathID}', '#000', '#000');"
                    onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${(x1_prime + x2_prime) / 2 - font_size / 4}, ${(y1_prime + y2_prime) / 2}
                    a ${font_size / 4},${font_size / 4} 0 1,0 ${font_size / 2},0
                    a ${font_size / 4},${font_size / 4} 0 1,0 ${-1 * font_size / 2},0
                    "
                    stroke="#000" stroke-width="${font_size / 6}" fill="#000"
                />`)
                } else {
                    baseMap.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt,  '${n1}${start} - ${n2}${end}; ${type}', '${pathID}', '#000', '#000');" onmouseout="UiActionsService.hideTooltip('${pathID}');" stroke-width="${font_size / 6}" data-stroke-color="#000" stroke="#000" d="M${x1_prime} ${y1_prime} ${x2_prime} ${y2_prime}"></path>`)
                }
            }
            else if (type == 'tWW') {
                let xm1 = UiTemplateService.linearlyInterpolate(x1_prime, xm, 1 - (font_size / 3) / (distance / 2))
                let ym1 = UiTemplateService.linearlyInterpolate(y1_prime, ym, 1 - (font_size / 3) / (distance / 2))
                let xm2 = UiTemplateService.linearlyInterpolate(xm, x2_prime, (font_size / 3) / (distance / 2))
                let ym2 = UiTemplateService.linearlyInterpolate(ym, y2_prime, (font_size / 3) / (distance / 2))
                baseMap.get(type)![1].push(defaultAction +
                    `d="
                    M ${x1_prime} ${y1_prime} ${xm1} ${ym1}
                    M ${xm - font_size / 3} ${ym}
                    a ${font_size / 3},${font_size / 3} 0 1,0 ${font_size / 1.5},0
                    a ${font_size / 3},${font_size / 3} 0 1,0 ${-1 * font_size / 1.5},0
                    M ${xm2} ${ym2} ${x2_prime} ${y2_prime}"
                    stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}"/>`
                )
            } else if (type == 'cSS' || type == 'tSS') {
                baseMap.get(type)![1].push(defaultAction + `
                d="
                M ${xm} ${ym + distance2 / 2} ${xm} ${ym + height / 2} 
                l ${height / 2} 0
                l -${height / 2} -${height} 
                l -${height / 2} ${height}
                l ${height / 2} 0
                M ${xm} ${ym - height / 2} ${xm} ${ym - distance2 / 2}
                "stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
            } else if (type == 'tHS' || type == 'cHS') {
                baseMap.get(type)![1].push(defaultAction + `
                d="
                M ${xm} ${ym + distance2 / 2} ${xm} ${ym + height + height / 4} 
                h -${height / 2}
                v -${height}
                h ${height}
                v ${height}
                h -${height / 2}
                M ${xm} ${ym + height / 4} ${xm} ${ym - height / 4}
                l ${height / 2} 0
                l -${height / 2} -${height} 
                l -${height / 2} ${height}
                l ${height / 2} 0
                M ${xm} ${ym - height - height / 4} ${xm} ${ym - distance2 / 2}
                "stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
            } else if (type == 'tWS' || type == 'cWS') {
                baseMap.get(type)![1].push(defaultAction + `
                d="
                M ${xm} ${ym + distance2 / 2} ${xm} ${ym + height + height / 4} 
                M ${xm - height / 2} ${ym + 3 * height / 4} 
                a ${height / 2},${height / 2} 0 1,0 ${height},0
                a ${height / 2},${height / 2} 0 1,0 ${-1 * height},0
                M ${xm} ${ym + height / 4} ${xm} ${ym - height / 4}
                l ${height / 2} 0
                l -${height / 2} -${height} 
                l -${height / 2} ${height}
                l ${height / 2} 0
                M ${xm} ${ym - height - height / 4} ${xm} ${ym - distance2 / 2}
                "stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
            } else if (type == 'tWH' || type == 'cWH') {
                baseMap.get(type)![1].push(defaultAction + `
                d="
                M ${xm} ${ym + distance2 / 2} ${xm} ${ym + height + height / 4} 
                M ${xm - height / 2} ${ym + 3 * height / 4} 
                a ${height / 2},${height / 2} 0 1,0 ${height},0
                a ${height / 2},${height / 2} 0 1,0 ${-1 * height},0
                M ${xm} ${ym + height / 4} ${xm} ${ym - height / 4}
                h -${height / 2}
                v -${height}
                h ${height}
                v ${height}
                h -${height / 2}
                M ${xm} ${ym - height - height / 4} ${xm} ${ym - distance2 / 2}
                "stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
            }
            else if (type == 'tHH' || type == 'cHH') {
                baseMap.get(type)![1].push(defaultAction + `
                d="
                M ${xm} ${ym + distance2 / 2} ${xm} ${ym + height / 2} 
                h -${height / 2}
                v -${height}
                h ${height}
                v ${height}
                h -${height / 2}
                M ${xm} ${ym - height / 2} ${xm} ${ym - distance2 / 2}
                "stroke="${stroke}" stroke-width="${font_size / 6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
            }
        }
    }

    private addAnnotationHandlers() {
        
        const annotateButton = document.getElementById(`rnaTopologyAnnotate-${this.pluginOptions.pdbId}`);
        const modal = document.getElementById(`annotationModal-${this.pluginOptions.pdbId}`);
        const addButton = document.getElementById(`addAnnotations-${this.pluginOptions.pdbId}`);
    
        if (annotateButton) {
            annotateButton.addEventListener('click', () => {
                if (modal) {
                    modal.style.display = 'block';
                }
            });
        }
    
        if (addButton) {
            addButton.addEventListener('click', () => {
                const startInput = document.getElementById(`startIndex-${this.pluginOptions.pdbId}`) as HTMLInputElement;
                const intervalInput = document.getElementById(`interval-${this.pluginOptions.pdbId}`) as HTMLInputElement;
                
                const startIndex = parseInt(startInput.value);
                const interval = parseInt(intervalInput.value);
    
                if (startIndex && interval) {
                    this.createAnnotations(startIndex, interval);
                    if (modal) {
                        modal.style.display = 'none';
                    }
                }
            });
        }
    }
    
    private createAnnotations(startIndex: number, interval: number) {
        const svg = document.querySelector(`svg.rnaTopoSvg`);
        if (!svg) return;
    
        // Remove any existing annotations
        const existingAnnotations = svg.querySelectorAll('.nucleotide-annotation');
        existingAnnotations.forEach(el => el.remove());
    
        // Get all nucleotides
        const nucleotides = Array.from(svg.querySelectorAll(`.rnaviewEle_${this.pluginOptions.pdbId}`));
        var font_size = this.apiData ? this.calculateFontSize(this.apiData) : 12;
        font_size += 4;
    
        // Step 1: Build a map of occupied regions (bounding boxes)
        const occupiedRegions: { x: number, y: number, width: number, height: number }[] = [];
        nucleotides.forEach((nucleotide) => {
            if (nucleotide instanceof SVGGraphicsElement) {
                occupiedRegions.push(nucleotide.getBBox());
            }
        });
    
        // Step 2: Create annotations at non-colliding positions
        for (let i = startIndex - 1; i < nucleotides.length; i += interval) {
            const nucleotide = nucleotides[i];
            if (nucleotide && nucleotide instanceof SVGGraphicsElement) {
                const bbox = nucleotide.getBBox();
                let x = bbox.x + bbox.width / 2;
                let y = bbox.y;
    
                // Step 3: Find a non-colliding position for annotation
                let offsetX = 0, offsetY = font_size + 5;
                let foundSpot = false;
    
                for (let attempt = 0; attempt < 3; attempt++) { // Try a few positions
                    let newX = x + offsetX;
                    let newY = y - offsetY;
    
                    const collides = occupiedRegions.some(region => 
                        newX >= region.x && newX <= region.x + region.width &&
                        newY >= region.y && newY <= region.y + region.height
                    );
    
                    if (!collides) {
                        foundSpot = true;
                        x = newX;
                        y = newY;
                        break;
                    }
    
                    // Try moving text sideways
                    offsetX = attempt % 2 === 0 ? font_size * 2 : -font_size * 2;
                    offsetY += font_size * 1.5;
                }
    
                if (!foundSpot) continue; // Skip annotation if no spot is found
    
                // Step 4: Create annotation line
                const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                line.setAttribute('x1', String(bbox.x + bbox.width / 2));
                line.setAttribute('y1', String(bbox.y));
                line.setAttribute('x2', String(x));
                line.setAttribute('y2', String(y));
                line.setAttribute('stroke', '#666');
                line.setAttribute('stroke-width', '1');
                line.setAttribute('class', 'nucleotide-annotation');
    
                // Step 5: Create annotation text
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', String(x));
                text.setAttribute('y', String(y - 2));
                text.setAttribute('text-anchor', 'middle');
                text.setAttribute('fill', '#666');
                text.setAttribute('font-size', String(font_size * 0.8));
                text.setAttribute('class', 'nucleotide-annotation');
                text.textContent = String(i + 1);
    
                // Step 6: Add annotation to SVG
                const svgGroup = svg.querySelector(`.rnaTopoSvg_${this.pluginOptions.pdbId}`);
                if (svgGroup) {
                    svgGroup.appendChild(line);
                    svgGroup.appendChild(text);
                }
    
                // Step 7: Update occupied regions
                occupiedRegions.push({ x: x - font_size, y: y - font_size, width: font_size * 2, height: font_size });
            }
        }
    }
    
    public static linearlyInterpolate(v0: number, v1: number, interpolationFactor: number): number {
        // See https://en.wikipedia.org/wiki/Linear_interpolation
        return (1 - interpolationFactor) * v0 + interpolationFactor * v1;
    }
    PathOrNucleotide(e?: Event) {
        const selectBoxEle: any = this.containerElement.querySelector<HTMLElement>('.menuSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        this.mappingValue = ''
        const nestedBP = this.containerElement.querySelector<HTMLInputElement>('#nestedBP')
        if (nestedBP!.checked) {
            var displayBP = this.displayNestedBaseStrs
        } else {
            var displayBP = this.displayBaseStrs
        }
        if (selectedValue == 0) {
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.nucleotideStrs.join('') + this.circleStrs.join('') + displayBP;
        } else if (selectedValue == 1) {
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.pathStrs.join('') + this.circleStrs.join('') + displayBP;
        } else if (selectedValue == 2) {
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.nucleotideStrs.join('') + this.circleStrs.join('') + displayBP;
            this.mappingValue = 'circle';
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.nucleotideStrs.join('') + this.circleStrs.join('') + displayBP;
            //this.circleStrs.join('') + displayBP;
        }
        //if(e) {
        //this.containerElement.querySelector<HTMLInputElement>('.mappingSelectbox')!.value="0";
        this.colorMap()
        //}
        //console.log("Path or nucleotide!")
        this.colorMapContacts()
        this.colorMapModifications()
    }


    private async addEvents(apiData: ApiData, BanName: any = undefined) {

        if (this.rv3VUEcomponent.protein_contacts) {
            let proteinContactsValues = Object.keys(this.rv3VUEcomponent.protein_contacts);
            for (const val of proteinContactsValues) {
                await BanNameHelper.getBanName(this.pluginOptions.pdbId, val);
            }
        }
        /*
        const lastPathIndex = apiData.svg_paths.length - 1;
        let strokeColor = this.pluginOptions.theme?.color || '#323232';
        
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            let isUnobserved = false;
            if(apiData.unobserved_label_seq_ids && apiData.unobserved_label_seq_ids.indexOf(apiData.label_seq_ids[recordIndex - 1]) > -1) {
                strokeColor = this.pluginOptions.theme?.unobservedColor || '#ccc';
                isUnobserved = true;
            }
            if(recordIndex === 0 || recordIndex === 1 || recordIndex === lastPathIndex + 1) return;
            let classPath = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}`
            let circlePath = `circle_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}`
            document.getElementsByClassName(circlePath)[0].addEventListener('click', UiActionsService.selectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1]))
            document.getElementsByClassName(circlePath)[0].addEventListener('mouseover', UiActionsService.selectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1], 'mouseover', isUnobserved, apiData.sequence[recordIndex - 2], undefined, this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined)) 
            document.getElementsByClassName(circlePath)[0].addEventListener('mouseout', UiActionsService.unSelectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1], isUnobserved, undefined, strokeColor))

            document.getElementsByClassName(classPath)[0].addEventListener('click', UiActionsService.selectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1]))
            document.getElementsByClassName(classPath)[0].addEventListener('mouseover', UiActionsService.selectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1], 'mouseover', isUnobserved, apiData.sequence[recordIndex - 2], undefined, this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined)) 
            document.getElementsByClassName(classPath)[0].addEventListener('mouseout', UiActionsService.unSelectNucleotide.bind(this, this.pluginOptions.pdbId, this.pluginOptions.entityId, apiData.label_seq_ids[recordIndex - 1], isUnobserved, undefined, strokeColor))
        });
        */
        /*if(!whichEvents) {
            for (let i in this.basePairIDs) {
                let pathInfo = this.basePairIDs[i]
                if(this.baseStrs.get(pathInfo[7])![0]) {
                    document.getElementsByClassName(pathInfo[0])[0].addEventListener('mouseover', (e) =>  {UiActionsService.showTooltip(e, `${pathInfo[3]}${pathInfo[1]} - ${pathInfo[4]}${pathInfo[2]}; ${pathInfo[7]}`, pathInfo[0], pathInfo[6], pathInfo[5])})
                    document.getElementsByClassName(pathInfo[0])[0].addEventListener('mouseout', UiActionsService.hideTooltip.bind(this, pathInfo[0]))
                }
            }
        }*/
        let contacts = new Map<number, string[]>()
        this.mapped_chains.forEach((val) => {
            //let tooltip = "Contacts: " + val;
            for (var i in this.rv3VUEcomponent.protein_contacts[val]) {
                let path = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.protein_contacts[val][i]}`
                if (document.getElementsByClassName(path)[0]) {
                    if (contacts.get(this.rv3VUEcomponent.protein_contacts[val][i])) {
                        let newList: string[] = contacts.get(this.rv3VUEcomponent.protein_contacts[val][i])!
                        newList.push(val)
                        contacts.set(this.rv3VUEcomponent.protein_contacts[val][i], newList)
                    } else {
                        contacts.set(this.rv3VUEcomponent.protein_contacts[val][i], [val])
                    }
                    //document.getElementsByClassName(path)[0].addEventListener('mouseover', (e) =>  {UiActionsService.showTooltip(e, tooltip, path, this.rv3VUEcomponent.proteinColorMap.get(val), this.rv3VUEcomponent.proteinColorMap.get(val))})
                    //document.getElementsByClassName(path)[0].addEventListener('mouseout', UiActionsService.hideTooltip.bind(this, path))
                }
            }
        });

        let chem_mods = new Map<number, string[]>()
        this.mapped_modifications.forEach((val) => {
            //let tooltip = "Contacts: " + val;
            for (var i in this.rv3VUEcomponent.modified_residues.get(val)) {
                //let path = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues[val][i]}`
                let path = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues.get(val)[i]}`;
                if (document.getElementsByClassName(path)[0]) {
                    if (chem_mods.get(this.rv3VUEcomponent.modified_residues.get(val)[i])) {
                        let newList: string[] = chem_mods.get(this.rv3VUEcomponent.modified_residues.get(val)[i])!
                        newList.push(val)
                        chem_mods.set(this.rv3VUEcomponent.modified_residues.get(val)[i], newList)

                    } else {
                        chem_mods.set(this.rv3VUEcomponent.modified_residues.get(val)[i], [val])

                    }
                    //document.getElementsByClassName(path)[0].addEventListener('mouseover', (e) =>  {UiActionsService.showTooltip(e, tooltip, path, this.rv3VUEcomponent.proteinColorMap.get(val), this.rv3VUEcomponent.proteinColorMap.get(val))})
                    //document.getElementsByClassName(path)[0].addEventListener('mouseout', UiActionsService.hideTooltip.bind(this, path))
                }

            }
        });
        const _this = this;
        async function processContacts(contacts: Map<number, string[]>) {
            await Promise.all(
                Array.from(contacts.entries()).map(async ([key, value]) => {
                    let tooltip = "Contact:";
                    let circlePath = `circle_${_this.pluginOptions.pdbId}_${key}`;
                    let path = `rnaviewEle rnaviewEle_${_this.pluginOptions.pdbId} rnaview_${_this.pluginOptions.pdbId}_${key}`;
                    let chain = '';
                    let promises = [];

                    for (let val of value) {
                        promises.push(BanNameHelper.getBanName(_this.pluginOptions.pdbId, val));
                        chain = val;
                    }

                    const banNames = await Promise.all(promises);

                    banNames.forEach((BanName: any, index) => {
                        if (BanName) {
                            tooltip += ` ${BanName[0]}; chain: ${chain}`;
                        } else {
                            tooltip = `chain: ${chain}`;
                        }

                        _this.toolTips.set(key, tooltip);
                        (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute(
                            'onmouseover',
                            document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] +
                            `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${_this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${_this.rv3VUEcomponent.proteinColorMap.get(chain)}')`
                        );
                        (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute(
                            'onmouseout',
                            document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] +
                            `;UiActionsService.hideTooltip('${circlePath}')`
                        );
                        if ((<HTMLElement>document.getElementsByClassName(path)[0])) {
                            (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute(
                                'onmouseover',
                                document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] +
                                `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${_this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${_this.rv3VUEcomponent.proteinColorMap.get(chain)}')`
                            );
                            (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute(
                                'onmouseout',
                                document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] +
                                `;UiActionsService.hideTooltip('${path}')`
                            );
                        }
                    });
                })
            );
        }
        await processContacts(contacts)
        /*
        contacts.forEach((value: string[], key: number) => {
            let tooltip = "Contact:"
            let circlePath = `circle_${this.pluginOptions.pdbId}_${key}`
            let path = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${key}`
            var chain = ''
            let promises = []
            
            //let BanProtName = ''
            for (let val in value) {

                promises.push(BanNameHelper.getBanName(this.pluginOptions.pdbId, value[val]));
                
                
                chain = value[val];

            }

            /*Promise.all(promises).then ((banNames : Array<any>) => {
                    
                banNames.forEach(BanName => {
                    if (BanName.length > 0) {
                        tooltip += " " +BanName[0]+"; chain: " + chain;
                    } else {
                        tooltip = "chain: " + chain;
                    }
              
                    });
            const banNames = await Promise.all(promises);

            banNames.forEach((BanName, index) => {
                if (BanName.length > 0) {
                tooltip += " " + BanName[0] + "; chain: " + chain;
                } else {
                tooltip = "chain: " + chain;
                }
                    
            this.toolTips.set(key, tooltip);
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseover', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}')`);
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseout', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${circlePath}')`);
            if((<HTMLElement>document.getElementsByClassName(path)[0])) {
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseover', document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}')`);
            //this.eventMap.set(path, UiActionsService.hideTooltip.bind(this, path))
            //document.getElementsByClassName(path)[0].addEventListener('mouseout', this.eventMap.get(path))
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseout', document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${path}')`);
            }
              });
           


        });*/
        chem_mods.forEach((value: string[], key: number) => {
            let tooltip = "Modified_nucleotide:"
            let circlePath = `circle_${this.pluginOptions.pdbId}_${key}`
            let path = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${key}`
            var mod = ''
            for (let val in value) {
                tooltip += " " + value[val];
                mod = value[val];
            }
            this.toolTips.set(key, tooltip);



            //this.mouseOverMap.set(circlePath, (e: Event) =>  {UiActionsService.showTooltip(e, tooltip, circlePath, this.rv3VUEcomponent.proteinColorMap.get(chain), this.rv3VUEcomponent.proteinColorMap.get(chain))})
            //document.getElementsByClassName(circlePath)[0].addEventListener('mouseover', this.mouseOverMap.get(circlePath))
            //(<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseover', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}')`);
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseover', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}')`);
            //this.eventMap.set(path, UiActionsService.hideTooltip.bind(this, circlePath))
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseout', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${circlePath}')`);
            //document.getElementsByClassName(circlePath)[0].addEventListener('mouseout', this.eventMap.get(circlePath))
            //this.mouseOverMap.set(path, (e: Event) =>  {UiActionsService.showTooltip(e, tooltip, path, this.rv3VUEcomponent.proteinColorMap.get(chain), this.rv3VUEcomponent.proteinColorMap.get(chain))})
            //document.getElementsByClassName(path)[0].addEventListener('mouseover', this.mouseOverMap.get(path))
            if ((<HTMLElement>document.getElementsByClassName(path)[0])) {
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseover', document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}')`);
                //this.eventMap.set(path, UiActionsService.hideTooltip.bind(this, path))
                //document.getElementsByClassName(path)[0].addEventListener('mouseout', this.eventMap.get(path))
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseout', document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${path}')`);
            }

        });
    }

    private removeEventHandlers(svgData: any) {
        var nodeList = svgData.childNodes[1].childNodes
        for (let i in nodeList) {
            var el = nodeList[i];
            if (el.attributes) {
                var attributes = [].slice.call(el!.attributes);

                for (let i = 0; i < attributes.length; i++) {
                    var att = attributes[i].name;
                    if (att.indexOf("on") === 0) {
                        el!.attributes.removeNamedItem(att);
                    }
                }
            }
        }
        return svgData
    }

    private svgTemplate(apiData: ApiData, FR3DData: any, FR3DNestedData: any): string {
        const font_size: number = this.calculateFontSize(apiData)
        var lastPathIndex = apiData.svg_paths.length - 1;
        var locations2: Map<any, number[]> = new Map();
        var locations3: Map<any, number[]> = new Map();
        if (this.pluginOptions.pdbId == "cust" || this.rv3VUEcomponent.structFailed) {
            lastPathIndex = lastPathIndex - 2
            apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
                if (recordIndex == 0 || recordIndex >= lastPathIndex + 1) return;
                let pathStrParsed: string[] = pathStr.split('M').join(',').split(',')
                let xVal: number = Number(pathStrParsed[3])
                let yVal: number = Number(pathStrParsed[4])
                let midX = (xVal);
                let midY = (yVal);
                locations2.set(apiData.label_seq_ids[recordIndex - 1], [midX, midY])
                if (recordIndex >= 2) {
                    let lastX = locations2.get(apiData.label_seq_ids[recordIndex - 2])![0]
                    let lastY = locations2.get(apiData.label_seq_ids[recordIndex - 2])![1]
                    locations3.set(apiData.label_seq_ids[recordIndex - 1], [(lastX + midX) / 2, (lastY + midY) / 2])
                }
                this.locations.set(apiData.label_seq_ids[recordIndex - 1], [xVal, yVal])
            });
        } else {
            apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
                if (recordIndex === 0 || recordIndex === 1 || recordIndex >= lastPathIndex + 1) return;
                let pathStrParsed: string[] = pathStr.split('M').join(',').split(',')
                let x1Val: number = Number(pathStrParsed[1])
                let y1Val: number = Number(pathStrParsed[2])
                let xVal: number = Number(pathStrParsed[3])
                let yVal: number = Number(pathStrParsed[4])
                let midX = (xVal + x1Val) / 2;
                let midY = (yVal + y1Val) / 2;
                locations2.set(apiData.label_seq_ids[recordIndex - 1], [midX, midY])
                this.locations.set(apiData.label_seq_ids[recordIndex - 1], [xVal, yVal])
            });
        }
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            //if(recordIndex === 0 || recordIndex === 1 || recordIndex === (lastPathIndex + 1)) return;
            if (recordIndex === 0 || recordIndex === 1 || recordIndex >= (lastPathIndex + 1)) return;
            const pathEleClass = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}`;
            let strokeColor = this.pluginOptions.theme?.color || '#323232';
            const strokeWide = this.pluginOptions.theme?.strokeWidth || '2';
            let isUnobserved = false;

            if (apiData.unobserved_label_seq_ids && apiData.unobserved_label_seq_ids.indexOf(apiData.label_seq_ids[recordIndex - 1]) > -1) {
                strokeColor = this.pluginOptions.theme?.unobservedColor || '#ccc';
                isUnobserved = true;
            }
            let pathStrParsed: string[] = pathStr.split('M').join(',').split(',')
            let xVal: number = Number(pathStrParsed[3])
            let yVal: number = Number(pathStrParsed[4])
            let deltaX: number = font_size / 2
            let deltaY: number = font_size / 2
            let newPathStr
            if (this.pluginOptions.pdbId == "cust") {
                if (recordIndex < (lastPathIndex)) {
                    let newX: number = locations3.get(apiData.label_seq_ids[recordIndex - 1])![0]
                    let newX2: number = locations3.get(apiData.label_seq_ids[recordIndex])![0]
                    let newY: number = locations3.get(apiData.label_seq_ids[recordIndex - 1])![1]
                    let newY2: number = locations3.get(apiData.label_seq_ids[recordIndex])![1]
                    newPathStr = `M${newX + deltaX},${newY - deltaY},${newX2 + deltaX},${newY2 - deltaY}`
                } else {
                    let newX: number = locations3.get(apiData.label_seq_ids[recordIndex - 1])![0]
                    let newY: number = locations3.get(apiData.label_seq_ids[recordIndex - 1])![1]
                    let newX2: number = 2 * xVal - newX
                    let newY2: number = 2 * yVal - newY
                    newPathStr = `M${newX + deltaX},${newY - deltaY},${newX2 + deltaX},${newY2 - deltaY}`
                }
            }
            else {
                if (recordIndex < (lastPathIndex)) {
                    let newX: number = locations2.get(apiData.label_seq_ids[recordIndex - 1])![0]
                    let newX2: number = locations2.get(apiData.label_seq_ids[recordIndex])![0]
                    let newY: number = locations2.get(apiData.label_seq_ids[recordIndex - 1])![1]
                    let newY2: number = locations2.get(apiData.label_seq_ids[recordIndex])![1]
                    newPathStr = `M${newX + deltaX},${newY - deltaY},${newX2 + deltaX},${newY2 - deltaY}`
                } else {
                    let newX: number = locations2.get(apiData.label_seq_ids[recordIndex - 1])![0]
                    let newY: number = locations2.get(apiData.label_seq_ids[recordIndex - 1])![1]
                    let newX2: number = 2 * xVal - newX
                    let newY2: number = 2 * yVal - newY
                    newPathStr = `M${newX + deltaX},${newY - deltaY},${newX2 + deltaX},${newY2 - deltaY}`
                }
            }
            /*
            newPathStr = pathStr
            */
            pathStr = newPathStr;
            this.pathStrs.push(
                `<path 
                    class="${pathEleClass}" stroke-width="${strokeWide}" stroke="${strokeColor}" d="${pathStr}" 
                    data-stroke-color="${strokeColor}" 
                    onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                    onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                    onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')">
                </path>`)
            this.circleStrs.push(
                `<circle class="circle_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}" cx="${xVal + deltaX}" cy="${yVal - deltaY}" r="${2 * font_size / 3}" display="none" alignment-baseline="middle" stroke-width="${font_size / 6}" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')"/>`)
            this.nucleotideStrs.push(
                `<text href="#${pathEleClass}" class="${pathEleClass}" x="${xVal}" y="${yVal}" font-size = "${font_size}px" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')">${apiData.sequence[recordIndex - 2]}</text>`)
        });




        let baseArray = FR3DData.annotations;
        let nestedBaseArray = FR3DNestedData.annotations;

        this.baseStrs.set('cWW', [true, []]);
        this.baseStrs.set('tWW', [false, []]);
        this.baseStrs.set('cWH', [false, []]);
        this.baseStrs.set('tWH', [false, []]);
        this.baseStrs.set('cWS', [false, []]);
        this.baseStrs.set('tWS', [false, []]);
        this.baseStrs.set('cHH', [false, []]);
        this.baseStrs.set('tHH', [false, []]);
        this.baseStrs.set('cHS', [false, []]);
        this.baseStrs.set('tHS', [false, []]);
        this.baseStrs.set('cSS', [false, []]);
        this.baseStrs.set('tSS', [false, []]);

        this.nestedBaseStrs.set('cWW', [true, []]);
        this.nestedBaseStrs.set('tWW', [false, []]);
        this.nestedBaseStrs.set('cWH', [false, []]);
        this.nestedBaseStrs.set('tWH', [false, []]);
        this.nestedBaseStrs.set('cWS', [false, []]);
        this.nestedBaseStrs.set('tWS', [false, []]);
        this.nestedBaseStrs.set('cHH', [false, []]);
        this.nestedBaseStrs.set('tHH', [false, []]);
        this.nestedBaseStrs.set('cHS', [false, []]);
        this.nestedBaseStrs.set('tHS', [false, []]);
        this.nestedBaseStrs.set('cSS', [false, []]);
        this.nestedBaseStrs.set('tSS', [false, []]);

        nestedBaseArray.forEach((baseStr: any) => {
            this.calcBaseStrs(baseStr, this.nestedBaseStrs, font_size)
        })
        baseArray.forEach((baseStr: any) => {
            this.calcBaseStrs(baseStr, this.baseStrs, font_size)
        })
        this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
            if (value[0]) {
                this.displayBaseStrs += value[1].join('');
            }
            if (this.nestedBaseStrs.get(key)![0]) {
                this.displayNestedBaseStrs += this.nestedBaseStrs.get(key)![1].join('');
            }
        });

        // for (let key in tokens) {
        //     if (tokens.hasOwnProperty(key)) {
        //         let tokenArray = tokens[key];
        //         let startPos = this.nucleotideStrs.indexOf(tokenArray[0]);
        //         let endPos   = this.nucleotideStrs.lastIndexOf(tokenArray[1]) + tokenArray[1].length;
        //         this.nucleotideStrs = this.nucleotideStrs.slice(0, startPos) + `<g class="rnaTopoSvg_${this.pluginOptions.pdbId}" id="${key}">` + this.nucleotideStrs.slice(startPos, endPos) + `</g>` + this.nucleotideStrs.slice(endPos)
        //     }
        // }

        return `
        <div style="width:100%;height:100%;z-index:0;position:absolute;">
            <svg preserveAspectRatio="xMidYMid meet" 
            viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
            style="width:100%;height:100%;position:relative;">
                <g class="rnaTopoSvgSelection rnaTopoSvgSelection_${this.pluginOptions.pdbId}"></g>
            </svg>
        </div>
        <div style="width:100%;height:100%;z-index:1;position:absolute;">
            <svg preserveAspectRatio="xMidYMid meet"
            viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
            style="width:100%;height:100%;position:relative;">
                <g class="rnaTopoSvgHighlight rnaTopoSvgHighlight_${this.pluginOptions.pdbId}"></g>
            </svg>
        </div>
        <div id="tooltip" display="none" style="position:absolute; display: none;"></div> 
            <div style="width:100%;height:100%;z-index:2;position:absolute;">
                <svg class="rnaTopoSvg" preserveAspectRatio="xMidYMid meet" 
                    viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
                    style="width:100%;height:100%;">
                        <g class="rnaTopoSvg_${this.pluginOptions.pdbId}">${this.nucleotideStrs.join('')}${this.circleStrs.join('')}${this.displayBaseStrs}</g>
                </svg>
            </div>`
            ;
    }



    private title(): string {
        return `<span class="pdb-rna-view-title">${this.pluginOptions.pdbId.toUpperCase()} Chain ${this.pluginOptions.chainId}</span>`;
    }

    private tooltip(): string {
        return `<span class="pdb-rna-view-tooltip" id="${this.pluginOptions.pdbId}-rnaTopologyTooltip"></span>`;
    }

    private actionButtons(): string {
        return `
        <div class="pdb-rna-view-btn-group">
            <span class="pdb-rna-view-btn" title="Add Annotations" id="rnaTopologyAnnotate-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M20,2H4A2,2 0 0,0 2,4V22L6,18H20A2,2 0 0,0 22,16V4A2,2 0 0,0 20,2M8,14H6V12H8V14M8,11H6V9H8V11M8,8H6V6H8V8M15,14H10V12H15V14M18,11H10V9H18V11M18,8H10V6H18V8Z" />
                </svg>
            </span>
            <span class="pdb-rna-view-btn" title="Save as SVG" id="rnaTopologySaveSVG-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"
                            d="M12 3v13M7 12l5 5 5-5M4 19h16"/>
                </svg>
            </span>
            <span class="pdb-rna-view-btn" title="Save as JSON" id="rnaTopologySaveJSON-${this.pluginOptions.pdbId}">
                <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" fill="none" viewBox="0 0 24 24" style="width: 24px; height: 24px;">
                    <path stroke="currentColor" stroke-linecap="round" stroke-width="2" d="M9 17.5H15"/>
                    <path stroke="currentColor" stroke-linecap="round" stroke-width="2" d="M9 13.5H15"/>
                    <path stroke="currentColor" stroke-linecap="round" stroke-width="2" d="M9 9.5H15"/>
                    <path stroke="currentColor" stroke-linecap="square" stroke-width="2" d="M4.5 3.5V20.5C4.5 21.3284 5.17157 22 6 22H18C18.8284 22 19.5 21.3284 19.5 20.5V7.42883C19.5 7.02985 19.341 6.64731 19.0583 6.36582L15.112 2.43699C14.8309 2.15712 14.4504 2 14.0537 2H6C5.17157 2 4.5 2.67157 4.5 3.5Z"/>
                </svg>
            </span>
            <span class="pdb-rna-view-btn" title="Zoom In" id="rnaTopologyZoomIn-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M9,2A7,7 0 0,1 16,9C16,10.57 15.5,12 14.61,13.19L15.41,14H16L22,20L20,22L14,16V15.41L13.19,14.61C12,15.5 10.57,16 9,16A7,7 0 0,1 2,9A7,7 0 0,1 9,2M8,5V8H5V10H8V13H10V10H13V8H10V5H8Z" />
                </svg>
            </span>
            
            <span class="pdb-rna-view-btn" title="Zoom Out" id="rnaTopologyZoomOut-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M9,2A7,7 0 0,1 16,9C16,10.57 15.5,12 14.61,13.19L15.41,14H16L22,20L20,22L14,16V15.41L13.19,14.61C12,15.5 10.57,16 9,16A7,7 0 0,1 2,9A7,7 0 0,1 9,2M5,8V10H13V8H5Z" />
                </svg>
            </span>
    
            <span class="pdb-rna-view-btn" title="Reset View" id="rnaTopologyReset-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M12,4C14.1,4 16.1,4.8 17.6,6.3C20.7,9.4 20.7,14.5 17.6,17.6C15.8,19.5 13.3,20.2 10.9,19.9L11.4,17.9C13.1,18.1 14.9,17.5 16.2,16.2C18.5,13.9 18.5,10.1 16.2,7.7C15.1,6.6 13.5,6 12,6V10.6L7,5.6L12,0.6V4M6.3,17.6C3.7,15 3.3,11 5.1,7.9L6.6,9.4C5.5,11.6 5.9,14.4 7.8,16.2C8.3,16.7 8.9,17.1 9.6,17.4L9,19.4C8,19 7.1,18.4 6.3,17.6Z" />
                </svg>
            </span>
    
            <span class="pdb-rna-view-btn" title="Toggle Fullscreen" id="TopologyFSCR">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M5,5H10V7H7V10H5V5M14,5H19V10H17V7H14V5M17,14H19V19H14V17H17V14M10,17V19H5V14H7V17H10Z" />
                </svg>
            </span>
        </div>
        <!-- Add modal dialog -->
        <div id="annotationModal-${this.pluginOptions.pdbId}" class="annotation-modal" style="display:none;position:absolute;background:white;border:1px solid #ccc;padding:15px;z-index:1000;top:50%;left:50%;transform:translate(-50%,-50%);box-shadow:0 2px 10px rgba(0,0,0,0.1);border-radius:4px;">
                <h3 style="margin:0 0 10px">Add Annotations</h3>
                <div style="margin-bottom:10px">
                    <label>Starting Index:</label>
                    <input type="number" id="startIndex-${this.pluginOptions.pdbId}" min="1" style="width:100px;margin-left:5px">
                </div>
                <div style="margin-bottom:10px">
                    <label>Interval:</label>
                    <input type="number" id="interval-${this.pluginOptions.pdbId}" min="1" style="width:100px;margin-left:5px">
                </div>
                <div style="text-align:right">
                    <button onclick="document.getElementById('annotationModal-${this.pluginOptions.pdbId}').style.display='none'" style="margin-right:5px">Cancel</button>
                    <button id="addAnnotations-${this.pluginOptions.pdbId}">Add</button>
                </div>
        </div>`;
    }

    renderError(type?: string) {
        let errorContent = `<div class="pdb-rna-view-error">Error! Something went wrong!</div>`;
        if (type === 'apiError') {
            errorContent = `<div class="pdb-rna-view-error">
                RNA topology data for ${this.pluginOptions.pdbId.toUpperCase()} | ${this.pluginOptions.entityId} | ${this.pluginOptions.chainId.toUpperCase()} is not available!
            </div>`;
        }

        this.containerElement.innerHTML = `<div class="pdb-rna-view-container">${errorContent}</div>`;
    }

    saveJSON() {
        const self = this;
        (window as any).debugRNA = self;

        const name = `${self.pluginOptions.pdbId}_${self.pluginOptions.entityId}_${self.pluginOptions.chainId}`;
        interface Class {
            name: string;
            "font-size"?: string;
            "font-weight"?: string;
            "font-style"?: string;
            "font-family"?: string;
            color?: string;
        }

        interface Residue {
            residueIndex: number;
            residueName: string;
            x: number;
            y: number;
            classes: string[];
        }

        // Define edge types for better type safety
        type EdgeType = 'W' | 'H' | 'S';
        type OrientationType = 'c' | 't';

        const getEdgeName = (edge: string): string => {
            switch (edge.toUpperCase()) {
                case 'W': return 'watson_crick';
                case 'H': return 'hoogsteen';
                case 'S': return 'sugar_edge';
                default: return '';
            }
        };

        const getOrientationName = (orientation: string): string => {
            switch (orientation.toLowerCase()) {
                case 'c': return 'cis';
                case 't': return 'trans';
                default: return '';
            }
        };

        const needsFlipping = (edge1: EdgeType, edge2: EdgeType): boolean => {
            const priority: { [key in EdgeType]: number } = { 'W': 3, 'H': 2, 'S': 1 };
            return priority[edge1] < priority[edge2];
        };

        const processBasePairInteraction = (
            index1: string | number,
            index2: string | number,
            interactionType: string
        ): BasePair => {
            // Convert indices to numbers
            const idx1 = parseInt(index1.toString());
            const idx2 = parseInt(index2.toString());

            // Extract orientation and edges
            const orientation = interactionType[0] as OrientationType;
            const edge1 = interactionType[1] as EdgeType;
            const edge2 = interactionType[2] as EdgeType;

            // Check if we need to flip based on edge priority
            const shouldFlip = needsFlipping(edge1, edge2);

            // Determine final indices and edges based on flip status
            const [finalIndex1, finalIndex2] = shouldFlip ? [idx2, idx1] : [idx1, idx2];
            const [finalEdge1, finalEdge2] = shouldFlip ? [edge2, edge1] : [edge1, edge2];

            // Build the base pair type string
            const basePairType = [
                getOrientationName(orientation),
                getEdgeName(finalEdge1),
                getEdgeName(finalEdge2)
            ].join('_');

            return {
                residueIndex1: finalIndex1,
                residueIndex2: finalIndex2,
                basePairType
            };
        };

        interface BasePair {
            residueIndex1: number;
            residueIndex2: number;
            basePairType: string;
        }

        interface RNAMolecule {
            name: string;
            sequence: Residue[];
            basePairs: BasePair[];
            labels: any[];
        }

        interface RNAComplex {
            name: string;
            rnaMolecules: RNAMolecule[];
        }

        interface XRNAJSON {
            rnaComplexes: RNAComplex[];
            classes: Class[];
        }

        function createXRNAJSON(): XRNAJSON {
            // Get sequence and positions from nucleotide strings
            const sequence: Residue[] = [];
            self.nucleotideStrs.forEach((nucleotide: string) => {
                // Extract values using regex
                const residueIndexMatch = nucleotide.match(/unSelectNucleotide\('[^']*', '[^']*', (\d+)/);
                const residueNameMatch = nucleotide.match(/>([AUGC])<\/text>/);
                const xMatch = nucleotide.match(/x="([\d.-]+)"/);
                const yMatch = nucleotide.match(/y="([\d.-]+)"/);

                if (residueIndexMatch && residueNameMatch && xMatch && yMatch) {
                    sequence.push({
                        residueIndex: parseInt(residueIndexMatch[1]),
                        residueName: residueNameMatch[1],
                        x: parseFloat(xMatch[1]),
                        y: parseFloat(yMatch[1]),
                        classes: ["font0", "color0"]
                    });
                }
            });

            // Get base pairs from baseStrs
            const basePairs: BasePair[] = [];
            self.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                if (value[0]) {
                    value[1].forEach((baseStr: string) => {

                        const tooltipMatch = baseStr.match(/showTooltip\(evt,\s*'([AUGC])(\d+)\s*-\s*([AUGC])(\d+);\s*([ctw][WHS]{2})'/i);
                        if (tooltipMatch) {
                            // const [basePairType, , index1, , index2] = tooltipMatch;
                            const [, , index1, , index2, interactionType] = tooltipMatch;


                            const processedBP = processBasePairInteraction(parseInt(index1), parseInt(index2), interactionType);
                            basePairs.push(processedBP);

                        }
                    });
                }
            });

            return {
                rnaComplexes: [
                    {
                        name: name,
                        rnaMolecules: [
                            {
                                name: name,
                                sequence: sequence,
                                basePairs: basePairs,
                                labels: []
                            }
                        ]
                    }
                ],
                classes: [
                    {
                        name: "font0",
                        "font-size": `${self.apiData ? self.calculateFontSize(self.apiData) : 10}px`,
                        "font-weight": "normal",
                        "font-style": "normal",
                        "font-family": "Helvetica"
                    },
                    {
                        name: "color0",
                        color: self.pluginOptions.theme?.color || '#323232'
                    }
                ]
            };
        };

        function downloadJSON(jsonData: XRNAJSON, name: string) {
            const jsonBlob = new Blob([JSON.stringify(jsonData, null, 2)], { type: 'application/json' });
            const jsonUrl = URL.createObjectURL(jsonBlob);
            const downloadLink = document.createElement("a");
            downloadLink.href = jsonUrl;
            downloadLink.download = name;
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }

        // Get current timestamp
        const dataJSON = createXRNAJSON();
        downloadJSON(dataJSON, `${name}.json`);


    }
    saveSVG() {
        const self = this;
        (window as any).debugRNA = self;
        const name = `${self.pluginOptions.pdbId}_${self.pluginOptions.entityId}_${self.pluginOptions.chainId}`;

        function getNode(n: any, v?: any) {
            n = document.createElementNS("http://www.w3.org/2000/svg", n);
            for (var p in v)
                n.setAttributeNS(null, p.replace(/[0-9]/g, 'o').replace(/\$/g, 'd').replace(/\[/g, 'b').replace(/[A-Z]/g, function (m, p, o, s) { return "-" + m.toLowerCase(); }), v[p]);
            return n
        }

        var svgData1 = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)
        let svgData_forsave = svgData1!.cloneNode(true);
        svgData_forsave = this.removeEventHandlers(svgData_forsave)
        var svg = getNode("svg");
        svg.appendChild(svgData_forsave);

        //   function groupSVG(svg: any){
        //     svg = svg.replace(`</g>`, ``);
        //     svg = svg.replace(/<g\s+class="rnaTopoSvg_[^"]*"\s*>/g, '');

        //     console.log('svg', JSON.stringify(svg));

        //     // process the neucliotides
        //     let tokens: { [key: string]: string[] } = {
        //         'Text': ['<text', '</text>'],
        //         'Colors': ['<circle', '</circle>'],
        //         'Bonds': ['<path', '</path>']
        //     };


        //     for (let key in tokens) {

        //         if (Object.prototype.hasOwnProperty.call(tokens, key)) {
        //             let tokenArray = tokens[key];
        //             let startPos = svg.indexOf(tokenArray[0]);
        //             let endPos = svg.lastIndexOf(tokenArray[1]) + tokenArray[1].length;
        //             svg = svg.slice(0, startPos) + `<g class="${key}" id="${key}">` + svg.slice(startPos, endPos) + `</g>` + svg.slice(endPos);
        //         }
        //     }
        //     return svg;

        //   }
        function groupSVG(svg: string) {
            // Remove existing </g> tags
            svg = svg.replace(/<\/g>/g, '');

            // Remove existing <g> tags with class attribute containing "rnaTopoSvg_"
            svg = svg.replace(/<g\s+class="rnaTopoSvg_[^"]*"\s*>/g, '');

            let tokens: { [key: string]: RegExp } = {
                'Text': /<text[\s\S]*<\/text>/g,
                'Colors': /<circle[\s\S]*<\/circle>/g,
                'Bonds': /<path[\s\S]*<\/path>/g
            };


            for (let key in tokens) {
                if (Object.prototype.hasOwnProperty.call(tokens, key)) {
                    let tokenRegExp = tokens[key];
                    let match;
                    while ((match = tokenRegExp.exec(svg)) !== null) {
                        let startPos = match.index;
                        let endPos = match.index + match[0].length;
                        svg = svg.slice(0, startPos) + `<g class="${key}" id="${key}">` + svg.slice(startPos, endPos) + `</g>` + svg.slice(endPos);
                    }
                }
            }
            return svg;
        }

        function downloadSVG(svgEl: any, name: any) {
            svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
            var svgData = groupSVG(svgEl.outerHTML);
            var preface = '<?xml version="1.0" standalone="no"?>\r\n';
            var svgBlob = new Blob([preface, svgData], { type: "image/svg+xml;charset=utf-8" });
            var svgUrl = URL.createObjectURL(svgBlob);
            var downloadLink = document.createElement("a");
            downloadLink.href = svgUrl;
            downloadLink.download = name;
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }

        downloadSVG(svg, `${name}.svg`);

    }

}