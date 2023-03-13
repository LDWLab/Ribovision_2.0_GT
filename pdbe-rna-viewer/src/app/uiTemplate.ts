import { UiActionsService } from './uiActions';
import { PluginOptions, ApiData, BanNameHelper } from './data';

export class UiTemplateService {
    /*
    The main class to create the RNA topology viewer image with nucleotides and base pairing
    */
    selectSections_RV1 = (window as any).selectSections_RV1;
    nPropertyConstants = (window as any).nPropertyConstants;
    mapped_n_properties = (window as any).mapped_n_properties;
    nColorData = (window as any).nColorData;
    parsePVData = (window as any).parsePVData;
    getEntropyAnnotations = (window as any).getEntropyAnnotations;
    getTWCAnnotations = (window as any).getTWCAnnotations;
    getCustomAnnotations = (window as any).getCustomAnnotations;
    getAnnotationArray = (window as any).getAnnotationArray;
    rv3VUEcomponent = (window as any).vm;
    private containerElement: HTMLElement;
    private pluginOptions: PluginOptions;
    private uiActionsService: UiActionsService;
    private apiData: ApiData | undefined;


    private locations: Map<any, number[]> = new Map();
    menuStyle = 'position:relative;z-index:10;height:7%;width:500px;line-height:7%;background-color:#696969;padding: 0 10px;font-size:16px; color: #efefef;display:inline-block;';
    domainTypes: any[];
    selectedDomain: string;
    pathStrs: string[] = [];
    nucleotideStrs: string[] = [];
    circleStrs: string[] = [];
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
        /*Constructor for UiTemplateService
        @param containerElement {HTMLElement}: the container for the RNA viewer
        @param pluginOptions {PluginOptions}: the plugin options
        @param apiData {ApiData|undefined}: the PDBe RNA structure API data if available
        */
        this.containerElement = containerElement;
        this.pluginOptions = pluginOptions;
        this.apiData = apiData;
        this.uiActionsService = new UiActionsService(this.pluginOptions.pdbId);
    }

    render(apiData: ApiData, FR3DData: any, FR3DNestedData: any) {
        /*Generate inner HTML for RNA topology viewer
        @param apiData {ApiData}: the PDBe RNA structure API data
        @param FR3DData {JSON}: the FR3D base pairing API data (un-nested)
        @param FR3DNestedData {JSON}: the FR3D base pairing API data (nested)
        */
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
                            <select>
                                <option>Base Pairings</option>
                            </select>
                            <div class="overSelect"></div>
                        </div>
                        <div id="checkboxes">
                        </div>
                    </div>
                    </form>
                </div>
        </div>
        `;
        this.createModeDropdown()
        this.createBPDropdown()
        this.uiActionsService.applyButtonActions();
        this.addEvents();
        <any>document.querySelector(".saveSVG")!.addEventListener("click", this.saveSVG.bind(this));
        this.getAnnotationFromRibovision(this.mapped_n_properties);
    }
    drawCircle = function(pdbId: string, i: number, color: string) {
        /*Make circle visible over nucleotide
        @param pdbId {string}: The PDB ID of the molecule
        @param i {number}: The index of the nucleotide
        @param color {string}: The color of the circle
        */
        const circle = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)!.getElementsByClassName(`circle_${pdbId}_${i}`)[0]
        circle.setAttribute("stroke", `${color}`);
        circle.setAttribute("fill", `${color}`);
        circle.style.display = "block";
    }
    colorMap = () => {
        /*
        Color nucleotides according to current mapping selections
        */
        const selectBoxEle: any = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        if (selectedValue) {
            const selectedDomain = this.domainTypes[selectedValue];
            if (selectedDomain.data) {
                selectedDomain.data.forEach((val: any, i: number) => {
                    if (val != undefined && val.color != undefined) {
                        if (val.color.indexOf('NaN') < 0) {
                            UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, val.color, undefined, this.mappingValue);
                        }
                    }
                });
            }
        } else {
            const tempDomain = this.domainTypes[1];
            if (tempDomain.data) {
                tempDomain.data.forEach((val: any, i: number) => {
                    if (val != undefined && val.start != undefined) {
                        UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, 'rgb(0,0,0)', undefined, this.mappingValue);
                    }
                });
            }
        }
    }
    colorMapHelper = () => {
        /*
        Update both color mapping dropdowns to have the same selection
        */
        var mappingDropdown = (<HTMLSelectElement>this.containerElement.querySelector<HTMLElement>('.mappingSelectbox'));
        var num: number = +mappingDropdown.value;
        (<HTMLInputElement>document.getElementById('selectColorMappingProps')).value = mappingDropdown!.options[num].text;
        this.rv3VUEcomponent.selected_property = mappingDropdown!.options[num].text;
    }
    colorMapContacts = () => {
        /*
        Color nucleotides by currently selected protein contacts
        */
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
                    }
                    if (document.getElementsByClassName(cPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(cPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(cPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                    }
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
        this.addEvents()
    }
    colorMapModifications = () => {
        /*
        Color nucleotides by currently selected modified residue types
        */
        this.mapped_modifications.forEach((val) => {
            if (!this.rv3VUEcomponent.modifications.includes(val)) {
                for (var i in this.rv3VUEcomponent.modified_residues.get(val)) {
                    UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.modified_residues.get(val)[i], '#323232', undefined, this.mappingValue);
                    var nPath = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues.get(val)[i]}`
                    var cPath = `circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.modified_residues.get(val)[i]}`
                    if (document.getElementsByClassName(nPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(nPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(nPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(nPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                    }
                    if (document.getElementsByClassName(cPath).length > 0) {
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseover', document.getElementsByClassName(cPath)[0].getAttribute('onmouseover')!.split(';')[0]);
                        (<HTMLElement>document.getElementsByClassName(cPath)[0]).setAttribute('onmouseout', document.getElementsByClassName(cPath)[0].getAttribute('onmouseout')!.split(';')[0]);
                    }
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
        this.addEvents()
    }
    changeBP(val: string) {
        /*
        Update image to reflect a change in base pairing selection
        @param val {string}: The base pairing selection that has been changed
        */
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
        /*
        Create base pairing selection dropdown HTML
        */
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
        }
    }
    createModeDropdown() {
        /*
        Create visualization mode dropdown HTML
        */
        let optionList = `<option value="0">Nucleotides</option><option value="1">Path</option><option value="2">Circle</option>`;
        const selectBoxEle = this.containerElement.querySelector<HTMLElement>('.menuSelectbox');
        selectBoxEle!.innerHTML = optionList;
        selectBoxEle!.addEventListener("change", this.PathOrNucleotide.bind(this));
        const nestedBP = this.containerElement.querySelector<HTMLElement>('#nestedBP');
        nestedBP!.addEventListener("change", this.PathOrNucleotide.bind(this));
    }
    createDomainDropdown() {
        /*
        Create color mapping annotation dropdown 
        */
        this.selectedDomain = this.domainTypes[0];
        if (this.domainTypes.length > 0) {

            let optionList = '';
            this.domainTypes.forEach((opt: any, i: number) => {
                optionList = `${optionList}<option value="${i}">${opt.label}</option>`;
            });

            const selectBoxEle = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
            selectBoxEle!.innerHTML = optionList;
            selectBoxEle!.addEventListener("change", this.colorMapHelper.bind(this));
        } else {
            this.containerElement!.querySelector<HTMLElement>('.menuOptions')!.style.display = 'none';
        }
    }
    create2D3DAnnotations(name: string, residueDetails: any,
        TWCrgbMap: Map<number, any>, TWCData: Map<number, string>, chain_start: number, chain_end: number) {
        /*
        Create residue color annotations for 2D and 3D mappings
        @param name {string}: the name of the mapping type (ex. "Shannon entropy")
        @param residueDetails {any}: the list of residue details
        @param TWCrgbMap {Map<number, any>}: the map of nucleotide indices to colors
        @param TWCData {Map<number, string>} the map of nucleotide indices to mapping values
        @param chain_start: chain starting index
        @param chain_end: chain ending index
        */
        const _this = this;
        TWCData.forEach(function(value, index) {
            if (chain_start <= index && index <= chain_end) {
                let rgb_color = TWCrgbMap.get(index);
                (window as any).selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.pluginOptions.entityId,
                    residue_number: index,
                    color: rgb_color[1],
                    sideChain: false,
                });

                _this.defaultColours.qualityRiboVision = '#' + String(rgb_color[0].join(''));

                var colors = "rgb(" + String(rgb_color[0].join(',')) + ")"

                residueDetails.push({ //2d
                    start: index,
                    end: index,
                    color: colors,
                    tooltipMsg: Number.parseFloat(value).toPrecision(3),
                    tooltipPosition: "prefix"
                });
            }
        })
        return residueDetails;
    }
    getAnnotationFromRibovision(mapped_n_properties: Map<string, Array<Array<number>>>) {
        /*
        Create residue color annotations for 2D and 3D mappings
        @param mapped_n_properties {Map<string, Array<Array<number>>>} map of nucleotide properties 
        to values for each nucleotide index
        */
        const _this = this;
        const start = this.apiData ? this.apiData.label_seq_ids[1] : 0

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
                                sideChain: false,
                            }*/])
        this.nPropertyConstants = (window as any).nPropertyConstants
        this.nColorData = (window as any).nColorData
        if (mapped_n_properties) {
            mapped_n_properties.forEach((value, index) => {
                let residueDetails: any = [{}];
                let name = index;
                let separatedData = value;

                this.selectSections_RV1.set(name, [])

                let min = Math.min(...this.nPropertyConstants.get(name));
                let max = Math.max(...this.nPropertyConstants.get(name));
                let colormapArray = this.nColorData.get(name);
                if (name == "Shannon entropy") {

                    this.getEntropyAnnotations(separatedData, min, max, this.pluginOptions.chainId);
                    //console.log('name_SE', name,  this.getEntropyAnnotations(separatedData, min, max, this.pluginOptions.chainId));
                };
                if (name == "TwinCons") {
                    this.getTWCAnnotations(separatedData, min, max, this.pluginOptions.chainId);
                };

                if (name == "Custom Data") {

                    this.getCustomAnnotations(separatedData, min, max, this.pluginOptions.chainId);
                };

                const [TWCrgbMap, TWCData] = this.parsePVData(separatedData, min, max, colormapArray);
                this.selectSections_RV1.get(name).push({ entity_id: _this.pluginOptions.entityId, focus: true });
                const end = TWCData.size;

                if (void 0 !== TWCData) {
                    residueDetails = _this.create2D3DAnnotations(name, residueDetails,
                        TWCrgbMap, TWCData,
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

    calculateFontSize(apiData: ApiData) {
        /* 
        Determines font size for nucleotides to reduce overlapping
        @param apiData: the PDBe RNA structure API data
        */
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

    private calcBaseStrs(baseStr: any, baseMap: any, font_size: number) {
        /*
        Generates base pairing HTML
        @param baseStr: base pairs from FR3D API
        @param baseMap: a map of base pairing types (ex. 'cWW') to visibility and values.
        This is filled in with this function.
        @param font_size: the size of nucleotides
        */
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
    public static linearlyInterpolate(v0: number, v1: number, interpolationFactor: number): number {
        // See https://en.wikipedia.org/wiki/Linear_interpolation
        return (1 - interpolationFactor) * v0 + interpolationFactor * v1;
    }
    PathOrNucleotide(e?: Event) {
        /*
        Switches between displaying paths and individual nucleotides
        Updates current view
        */
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
            this.mappingValue = 'circle';
        }
        if (e) {
            //this.containerElement.querySelector<HTMLInputElement>('.mappingSelectbox')!.value="0";
            this.colorMap()
        }
        this.colorMapContacts()
        this.colorMapModifications()
    }

    private addEvents() {
        /*
        Adds events for protein contacts and RNA modifications
        */
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
                }
            }
        });

        let chem_mods = new Map<number, string[]>()
        this.mapped_modifications.forEach((val) => {
            for (var i in this.rv3VUEcomponent.modified_residues.get(val)) {
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

        contacts.forEach((value: string[], key: number) => {
            let tooltip = "Contact:"
            let circlePath = `circle_${this.pluginOptions.pdbId}_${key}`
            let path = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${key}`
            var chain = ''
            let promises = []

            for (let val in value) {
                promises.push(BanNameHelper.getBanName(this.pluginOptions.pdbId, value[val]));
                chain = value[val];
            }

            Promise.all(promises).then((banNames: Array<any>) => {
                banNames.forEach(BanName => {
                    if (BanName.length > 0) {
                        tooltip += " " + BanName[0] + "; chain: " + chain;
                    } else {
                        tooltip = "chain: " + chain;
                    }
                });

                this.toolTips.set(key, tooltip);
                (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseover', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}')`);
                (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseout', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${circlePath}')`);
                if ((<HTMLElement>document.getElementsByClassName(path)[0])) {
                    (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseover', document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}', '${this.rv3VUEcomponent.proteinColorMap.get(chain)}')`);
                    (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseout', document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${path}')`);
                }
            });

        });
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
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseover', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${circlePath}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}')`);
            (<HTMLElement>document.getElementsByClassName(circlePath)[0]).setAttribute('onmouseout', document.getElementsByClassName(circlePath)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${circlePath}')`);
            if ((<HTMLElement>document.getElementsByClassName(path)[0])) {
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseover', document.getElementsByClassName(path)[0].getAttribute('onmouseover')!.split(';')[0] + `;UiActionsService.showTooltip(evt, '${tooltip}', '${path}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}', '${this.rv3VUEcomponent.modifiedColorMap.get(mod)}')`);
                (<HTMLElement>document.getElementsByClassName(path)[0]).setAttribute('onmouseout', document.getElementsByClassName(path)[0].getAttribute('onmouseout')!.split(';')[0] + `;UiActionsService.hideTooltip('${path}')`);
            }

        });
    }

    private removeEventHandlers(svgData: any) {
        /*
        Removes all event handlers from SVG 
        */
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
        /*
            Generates nucleotide, circle, and path HTML to show in SVG
            @param apiData {ApiData}: the PDBe RNA structure API data
            @param FR3DData {JSON}: the FR3D base pairing API data (un-nested)
            @param FR3DNestedData {JSON}: the FR3D base pairing API data (nested)
        */
        const font_size: number = this.calculateFontSize(apiData)
        const lastPathIndex = apiData.svg_paths.length - 1;
        var locations2: Map<any, number[]> = new Map();
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            if (recordIndex === 0 || recordIndex === 1 || recordIndex === lastPathIndex + 1) return;
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
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            if (recordIndex === 0 || recordIndex === 1 || recordIndex === (lastPathIndex + 1)) return;
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

            pathStr = newPathStr;

            this.pathStrs.push(
                `<path 
                    class="${pathEleClass}" stroke-width="${strokeWide}" stroke="${strokeColor}" d="${pathStr}" 
                    data-stroke-color="${strokeColor}" 
                    onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                    onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                    onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event)">
                </path>`)
            this.circleStrs.push(
                `<circle class="circle_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}" cx="${xVal + deltaX}" cy="${yVal - deltaY}" r="${2 * font_size / 3}" display="none" alignment-baseline="middle" stroke-width="${font_size / 6}" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event)"/>`)
            this.nucleotideStrs.push(
                `<text href="#${pathEleClass}" class="${pathEleClass}" x="${xVal}" y="${yVal}" font-size = "${font_size}px" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'" + this.pluginOptions.theme.highlightColor + "'" : undefined})" 
                onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event)">${apiData.sequence[recordIndex - 2]}</text>`)
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
        /*
        Makes structure title HTML
        */
        return `<span class="pdb-rna-view-title">${this.pluginOptions.pdbId.toUpperCase()} Chain ${this.pluginOptions.chainId}</span>`;
    }

    private tooltip(): string {
        /*
        Makes tooltip HTML
        */
        return `<span class="pdb-rna-view-tooltip" id="${this.pluginOptions.pdbId}-rnaTopologyTooltip"></span>`;
    }

    private actionButtons(): string {
        /*
        Generates buttons in upper right corner of viewer
        */
        return `<div class="pdb-rna-view-btn-group">
            <span class="pdb-rna-view-btn" title="Zoom-in" id="rnaTopologySaveSVG-${this.pluginOptions.pdbId}">
                <img class="saveSVG" src="static/alignments/png/Save.png" style="height:24px; width: 24px; border:0;position: relative;cursor:pointer;" title="saveSVG"/>
            </span>
            <span class="pdb-rna-view-btn" title="Zoom-in" id="rnaTopologyZoomIn-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M15.5,14L20.5,19L19,20.5L14,15.5V14.71L13.73,14.43C12.59,15.41 11.11,16 9.5,16A6.5,6.5 0 0,1 3,9.5A6.5,6.5 0 0,1 9.5,3A6.5,6.5 0 0,1 16,9.5C16,11.11 15.41,12.59 14.43,13.73L14.71,14H15.5M9.5,14C12,14 14,12 14,9.5C14,7 12,5 9.5,5C7,5 5,7 5,9.5C5,12 7,14 9.5,14M12,10H10V12H9V10H7V9H9V7H10V9H12V10Z" />
                </svg>
            </span>
            
            <span class="pdb-rna-view-btn" title="Zoom-out" id="rnaTopologyZoomOut-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M15.5,14H14.71L14.43,13.73C15.41,12.59 16,11.11 16,9.5A6.5,6.5 0 0,0 9.5,3A6.5,6.5 0 0,0 3,9.5A6.5,6.5 0 0,0 9.5,16C11.11,16 12.59,15.41 13.73,14.43L14,14.71V15.5L19,20.5L20.5,19L15.5,14M9.5,14C7,14 5,12 5,9.5C5,7 7,5 9.5,5C12,5 14,7 14,9.5C14,12 12,14 9.5,14M7,9H12V10H7V9Z" />
                </svg>
            </span>

            <span class="pdb-rna-view-btn" title="Reset" id="rnaTopologyReset-${this.pluginOptions.pdbId}">
                <svg style="width:24px;height:24px" viewBox="0 0 24 24">
                    <path fill="currentColor" d="M12,6V9L16,5L12,1V4A8,8 0 0,0 4,12C4,13.57 4.46,15.03 5.24,16.26L6.7,14.8C6.25,13.97 6,13 6,12A6,6 0 0,1 12,6M18.76,7.74L17.3,9.2C17.74,10.04 18,11 18,12A6,6 0 0,1 12,18V15L8,19L12,23V20A8,8 0 0,0 20,12C20,10.43 19.54,8.97 18.76,7.74Z" />
                </svg>
            </span>
        </div>`;
    }

    renderError(type?: string) {
        /*
        Displays error message if rendering fails
        */
        let errorContent = `<div class="pdb-rna-view-error">Error! Something went wrong!</div>`;
        if (type === 'apiError') {
            errorContent = `<div class="pdb-rna-view-error">
                RNA topology data for ${this.pluginOptions.pdbId.toUpperCase()} | ${this.pluginOptions.entityId} | ${this.pluginOptions.chainId.toUpperCase()} is not available!
            </div>`;
        }

        this.containerElement.innerHTML = `<div class="pdb-rna-view-container">${errorContent}</div>`;
    }
    saveSVG() {
        /*
        Removes events from SVG and downloads a copy
        */
        function getNode(n: any, v?: any) {
            n = document.createElementNS("http://www.w3.org/2000/svg", n);
            for (var p in v)
                n.setAttributeNS(null, p.replace(/[0-9]/g, 'o').replace(/\$/g, 'd').replace(/\[/g, 'b').replace(/[A-Z]/g, function(m, p, o, s) { return "-" + m.toLowerCase(); }), v[p]);
            return n
        }

        var svgData1 = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)
        let svgData_forsave = svgData1!.cloneNode(true);
        svgData_forsave = this.removeEventHandlers(svgData_forsave)
        var svg = getNode("svg");
        svg.appendChild(svgData_forsave);
        function saveSvg1(svgEl: any, name: any) {
            svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
            var svgData = svgEl.outerHTML;
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
        saveSvg1(svg, 'rv3Topology.svg')
    }

}