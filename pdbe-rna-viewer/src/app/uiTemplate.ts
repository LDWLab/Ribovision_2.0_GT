import { UiActionsService } from './uiActions';
import { PluginOptions, ApiData} from './data';

export class UiTemplateService {
    selectSections_RV1 = (window as any).selectSections_RV1;
    aaPropertyConstants = (window as any).aaPropertyConstants;
    mapped_aa_properties = (window as any).mapped_aa_properties;
    aaColorData = (window as any).aaColorData;
    parsePVData = (window as any).parsePVData;
    rv3VUEcomponent = (window as any).vm;
    private containerElement: HTMLElement;
    private pluginOptions: PluginOptions;
    private uiActionsService: UiActionsService;
    private apiData: ApiData | undefined;
    private locations: Map<any, number[]> = new Map();
    menuStyle = 'position:relative;z-index:10;height:38px;width:500px;line-height:38px;background-color:#696969;padding: 0 10px;font-size:16px; color: #efefef;display:inline-block;';
    domainTypes: any[];
    selectedDomain: string;
    pathStrs: string[] = [];  
    nucleotideStrs: string[] = [];
    circleStrs: string[] = [];
    baseStrs: Map <string, [boolean, string[]]> = new Map();
    displayBaseStrs: string;
    mappingValue: string = '';
    showAllNucleotides: boolean = false;
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
    constructor(containerElement: HTMLElement, pluginOptions: PluginOptions, apiData: ApiData|undefined) {
        this.containerElement = containerElement;
        this.pluginOptions = pluginOptions;
        this.apiData = apiData;
        this.uiActionsService = new UiActionsService(this.pluginOptions.pdbId);
    }

    render(apiData: ApiData, FR3DData: any) {
        this.containerElement.innerHTML = 
        `<div class="pdb-rna-view-container pdb-rna-view-container-${this.pluginOptions.pdbId}">
            ${this.svgTemplate(apiData, FR3DData)}
            ${this.title()}
            ${this.tooltip()}
            ${this.actionButtons()}
        </div>
        <div style="${this.menuStyle}">
                <div class="menuOptions" style="float:right;margin-right: 20px;display:inline-block;">
                    <form>
                    <select class="menuSelectbox" style="margin-right: 10px; display: inline-block;"><option value="">Nucleotides</option></select>
                    <select class="mappingSelectbox" style="margin-right: 10px; display: inline-block;"><option value="">Mapping</option></select>
                    <div class="multiselect">
                        <div class="selectBox" onclick="UiActionsService.showCheckboxes()" style="display:inline-block;">
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
        //this.getJSON(apiData)
        //this.fixOverlaps(apiData)
        this.createModeDropdown()
        this.createBPDropdown()
        this.uiActionsService.applyButtonActions();
        //this.getAnnotationFromRibovision(this.mapped_aa_properties)
        //this.rv3VUEcomponent.topology_loaded=true;
    }

    fixOverlaps(apiData: ApiData) {
        var svgEle: any = (<any>document.querySelector(`svg.rnaTopoSvg.rnaTopoSvg_${this.pluginOptions.pdbId}`))
        var xAdjust: number = 0;
        var yAdjust: number = 0;
        apiData.sequence.split('').forEach((char: string, i: number) => {
            if(i === 0 || i === (apiData.sequence.length)) { return }
            let ele = svgEle!.getElementsByClassName(`rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[i]}`)[0]
            let nextEle = svgEle!.getElementsByClassName(`rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[i + 1]}`)[0]
            nextEle.setAttribute("y", Number(nextEle.getAttribute("y")) + yAdjust)
            nextEle.setAttribute("x", Number(nextEle.getAttribute("x")) + xAdjust)
            let distances = this.getIntersections(ele.getBBox(), nextEle.getBBox())
            if(distances) {
                var lowest = 0;
                for (var j = 1; j < distances.length; j++) {
                    if (distances[j] < distances[lowest]) lowest = j;
                }
                if(lowest === 1) {
                    nextEle.setAttribute("y",Number(nextEle.getAttribute("y")) - Number(distances[1]));
                    yAdjust = yAdjust - Number(distances[1])
                } else if(lowest === 3) {
                    nextEle.setAttribute("y",Number(nextEle.getAttribute("y")) + Number(distances[3]));
                    yAdjust = yAdjust + Number(distances[3])
                } else if(lowest === 0) {
                    nextEle.setAttribute("x",Number(nextEle.getAttribute("x")) - Number(distances[0]));
                    xAdjust = xAdjust - Number(distances[0])
                } else if(lowest === 2) {
                    nextEle.setAttribute("x",Number(nextEle.getAttribute("x")) + Number(distances[2]));
                    xAdjust = xAdjust + Number(distances[2])
                } 
            }
        });
    }
    drawCircle = function (pdbId: string, i: number, color: string){
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
    colorMap=() => {
        const selectBoxEle:any = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        if(selectedValue) {
            const selectedDomain = this.domainTypes[selectedValue];
            if (selectedDomain.data) {
                selectedDomain.data.forEach((val:any, i:number) => {
                    if(val != undefined && val.color != undefined) {
                        if(val.color.indexOf('NaN') < 0) {
                            //this.drawCircle(this.pluginOptions.pdbId, val.start, val.color);
                            //UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, val.color);
                            UiActionsService.colorNucleotide(this.pluginOptions.pdbId, val.start, val.color, undefined, this.mappingValue);
                        }
                    }
                });
            }
        }
    }
    colorMapContacts=() => {
        this.mapped_chains.forEach((val) => {  
            if(!this.rv3VUEcomponent.pchainid.includes(val)) {
                for(var i in this.rv3VUEcomponent.protein_contacts[val]) {
                    UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.protein_contacts[val][i], '#323232', undefined, this.mappingValue);
                    var circle = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`circle_${this.pluginOptions.pdbId}_${this.rv3VUEcomponent.protein_contacts[val][i]}`)[0];
                    if(circle) {
                        circle.style.display="none";
                    } 
                }
                this.mapped_chains.delete(val)
            }           
        });  

        for (let val in this.rv3VUEcomponent.pchainid) {
            var chain = this.rv3VUEcomponent.pchainid[val];
            for(var i in this.rv3VUEcomponent.protein_contacts[chain]) {
                UiActionsService.colorNucleotide(this.pluginOptions.pdbId, this.rv3VUEcomponent.protein_contacts[chain][i], this.rv3VUEcomponent.proteinColorMap.get(chain), undefined, this.mappingValue);
            }
            this.mapped_chains.add(chain)
        }
    }
    changeBP(val: string) {
        this.displayBaseStrs = '';
        if((val == 'All' && !this.showAllNucleotides)) {
            this.showAllNucleotides = true;
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                this.baseStrs.set(key, [true,  value[1]]);
                this.displayBaseStrs += value[1].join('');
                (<HTMLInputElement>document.getElementById(`Checkbox_${key}`))!.checked = true
            });
        } else if (val == 'All') {
            this.showAllNucleotides = false;
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                this.baseStrs.set(key, [false,  value[1]]);
                (<HTMLInputElement>document.getElementById(`Checkbox_${key}`))!.checked = false
            });
        }
        else {
            if(this.baseStrs.get(val)![0]) {
                this.baseStrs.set(val, [false,  this.baseStrs.get(val)![1]]);
            } else {
                this.baseStrs.set(val, [true,  this.baseStrs.get(val)![1]]);
            }
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                if(value[0]) {
                    this.displayBaseStrs += value[1].join('');
                }
            });
        }
        this.PathOrNucleotide();
    }
    createBPDropdown() {
        if(this.baseStrs.size > 0) {
            let optionList = '<label for = "Checkbox_All"><input type="checkbox" id="Checkbox_All" />All</label>';
            this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
                if(key == 'cWW') {
                    optionList = `${optionList}<label for = "Checkbox_${key}"><input type="checkbox" id="Checkbox_${key}" checked = true/>${key}</label>`;
                } else {
                    optionList = `${optionList}<label for = "Checkbox_${key}"><input type="checkbox" id="Checkbox_${key}"/>${key}</label>`;
                }
            });
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
    }
    createDomainDropdown () {
        this.selectedDomain = this.domainTypes[0];
        if(this.domainTypes.length > 0){

            let optionList = '';
            this.domainTypes.forEach((opt:any, i:number) => {
                optionList = `${optionList}<option value="${i}">${opt.label}</option>`;
            });

            const selectBoxEle = this.containerElement.querySelector<HTMLElement>('.mappingSelectbox');
            selectBoxEle!.innerHTML = optionList;
            selectBoxEle!.addEventListener("change", this.colorMap.bind(this));

            //const resetIconEle = this.containerElement.querySelector('.resetIcon');
            //resetIconEle.addEventListener("click", this.resetDisplay.bind(this));

        }else{
            this.containerElement!.querySelector<HTMLElement>('.menuOptions')!.style.display = 'none';
        }
    }
    create2D3DAnnotations(name: string, residueDetails: any, 
        TWCrgbMap: Map<number, any>, TWCData: Map<number, string>, mapped_aa_properties: Map<string, Array<Array<number>>>,
        chain_start: number, chain_end: number) {
        const _this = this;
        TWCData.forEach(function(value, index) {
            if (chain_start <= index && index <= chain_end){
                let rgb_color = TWCrgbMap.get(index);
                (window as any).selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.pluginOptions.entityId,
                    start_residue_number: index, 
                    end_residue_number: index,
                    color: rgb_color[1],
                    sideChain: false,
                });
                _this.defaultColours.qualityRiboVision= '#'+String(rgb_color[0].join(''));
                var colors = "rgb("+String(rgb_color[0].join(','))+")"
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
        if (TWCData.size < mapped_aa_properties.get("Shannon entropy")!.length) {
            for(var i = TWCData.size - 1; i < this.mapped_aa_properties.get("Charge").length; i++) {
                (window as any).selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.pluginOptions.entityId,
                    start_residue_number: i, 
                    end_residue_number: i,
                    color: {r:255, g:255, b:255},
                    sideChain: false,
                });
            }
        }
        return residueDetails;
    }
    getAnnotationFromRibovision(mapped_aa_properties: Map<string, Array<Array<number>>>) {
        const _this = this;
        const start = this.apiData?this.apiData.label_seq_ids[1]:0
        const end = this.apiData?this.apiData.label_seq_ids[this.apiData.label_seq_ids.length - 2]:0
        if(typeof _this.domainTypes == 'undefined'){
            _this.domainTypes = [{
                label: 'Select data',
                data: null
            }];
        }
        this.selectSections_RV1 = (window as any).selectSections_RV1
        this.aaPropertyConstants = (window as any).aaPropertyConstants
        this.aaColorData = (window as any).aaColorData
        if (mapped_aa_properties) {
            mapped_aa_properties.forEach((value, index) => {    
                let residueDetails:any = [{
                    //start: chainRange.start,
                    //end: chainRange.end,
                    //color: _this.defaultColours.qualityBlank,
                    //tooltipMsg: 'No data for '
                }];
                let name = index;
                let separatedData = value;

                this.selectSections_RV1.set(name, [])

                let min = Math.min(...this.aaPropertyConstants.get(name));
                let max = Math.max(...this.aaPropertyConstants.get(name));
                let colormapArray = this.aaColorData.get(name); 
                const [TWCrgbMap, TWCData] = this.parsePVData(separatedData, min, max, colormapArray);

                this.selectSections_RV1.get(name).push({entity_id: _this.pluginOptions.entityId, focus: true});
                
                if (void 0 !== TWCData){
                    residueDetails = _this.create2D3DAnnotations(name, residueDetails, 
                                                                TWCrgbMap, TWCData, mapped_aa_properties,
                                                                start, end);
                    if(0 < residueDetails.length){
                        var current = _this.domainTypes.filter(order => (order.label === name))[0];
                        if(current && current != null) {
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
            this.rv3VUEcomponent.topology_loaded=true;
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
            right1 <= left2 || bottom1 <= top2 ) {
          return false
        } else {
          return [right2-left1, bottom2-top1, right1-left2, bottom1-top2]
        }
      }

    calculateFontSize(apiData: ApiData) {
        let xVals: number[] = [];
        let yVals: number[] = [];
        let dist: number[] = [];
        const lastPathIndex = apiData.svg_paths.length - 1;
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            if(recordIndex === 0 || recordIndex === lastPathIndex) return;
            let pathStrParsed:string[] = pathStr.split('M').join(',').split(',')
            xVals[recordIndex] = (Number(pathStrParsed[1])+Number(pathStrParsed[3]))/2
            yVals[recordIndex] = (Number(pathStrParsed[2])+Number(pathStrParsed[4]))/2
            if(recordIndex > 1) {
                let xDiff = xVals[recordIndex] - xVals[recordIndex - 1]
                let yDiff = yVals[recordIndex] - yVals[recordIndex - 1]
                dist[recordIndex] = Math.pow((Math.pow(yDiff, 2) + Math.pow(xDiff, 2)),0.5)
            }
        });
        var sortedDist: number[] = dist.sort((a,b) => {return a - b;})
        return 0.9*sortedDist[Math.floor(sortedDist.length * 0.05)]
    }
    
    download(data: any, filename: string, type: string) {
        var file = new Blob([data], {type: type});
        var a = document.createElement("a"),
                url = URL.createObjectURL(file);
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(function() {
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);  
        }, 0); 
    }
    /*
    private getJSON(apiData: ApiData){
        let JSON = 
        `{
        "pdbId": ${this.pluginOptions.pdbId},
        "entityId": ${this.pluginOptions.entityId},
        "chainId": ${this.pluginOptions.chainId},
        "nucleotides":[
        `
        const lastPathIndex = apiData.svg_paths.length - 1;
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
        if(recordIndex === 0 || recordIndex === lastPathIndex) return;
            JSON += 
        `   {"position": ${apiData.label_seq_ids[recordIndex]}, "value": ${apiData.sequence[recordIndex - 1]}, "x": ${this.locations[recordIndex - 1][0]}, "y": ${this.locations[recordIndex - 1][1]}},
        `
        });
        JSON += 
            `]
        }`
        this.download(JSON, "Test", "JSON")
    }*/
    public static linearlyInterpolate(v0 : number, v1 : number, interpolationFactor : number) : number {
        // See https://en.wikipedia.org/wiki/Linear_interpolation
        return (1 - interpolationFactor) * v0 + interpolationFactor * v1;
    }
    PathOrNucleotide(e?: Event) {
        const selectBoxEle:any = this.containerElement.querySelector<HTMLElement>('.menuSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        this.mappingValue = ''
        if(selectedValue == 0) {
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.nucleotideStrs.join('') + this.circleStrs.join('') + this.displayBaseStrs;
        } else if(selectedValue == 1) {
            (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaTopoSvg_${this.pluginOptions.pdbId}`)[0].innerHTML = this.pathStrs.join('') + this.circleStrs.join('') + this.displayBaseStrs;
        } else if(selectedValue == 2) {
            this.mappingValue = 'circle';
        } 
        if(e) {
            //this.containerElement.querySelector<HTMLInputElement>('.mappingSelectbox')!.value="0";
            this.colorMap()
        }
        this.colorMapContacts()
    }
    private svgTemplate(apiData: ApiData, FR3DData: any): string { 
        const font_size:number = this.calculateFontSize(apiData)
        const lastPathIndex = apiData.svg_paths.length - 1;
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {

            if(recordIndex === 0 || recordIndex === 1 || recordIndex === (lastPathIndex + 1)) return;
            const pathEleClass = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}`;
            let strokeColor = this.pluginOptions.theme?.color || '#323232';
            const strokeWide = this.pluginOptions.theme?.strokeWidth || '2';
            let isUnobserved = false;
            if(apiData.unobserved_label_seq_ids && apiData.unobserved_label_seq_ids.indexOf(apiData.label_seq_ids[recordIndex - 1]) > -1) {
                strokeColor = this.pluginOptions.theme?.unobservedColor || '#ccc';
                isUnobserved = true;
            }
            let pathStrParsed:string[] = pathStr.split('M').join(',').split(',')
            let x1Val: number = Number(pathStrParsed[1]) 
            let y1Val: number = Number(pathStrParsed[2]) 
            //let xVal:number = (Number(pathStrParsed[1])+Number(pathStrParsed[3]))/2 
            let xVal:number = Number(pathStrParsed[3]) 
            //let yVal:number = (Number(pathStrParsed[2])+Number(pathStrParsed[4]))/2 
            let yVal:number = Number(pathStrParsed[4])
            let deltaX: number = font_size/2
            let deltaY: number = font_size/2
            let newPathStr = `M${x1Val + deltaX},${y1Val - deltaY},${xVal + deltaX},${yVal - deltaY}` 
            pathStr = newPathStr;
            this.locations.set(apiData.label_seq_ids[recordIndex - 1], [xVal, yVal])
            /*pathStrs.push(`<text href="#${pathEleClass}" class="${pathEleClass}" x="${xVal}" y="${yVal}" font-size = "${font_size}px" onclick="UiActionsService.selectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex - 1]}, '${apiData.sequence[recordIndex - 2]}', 'click', ${isUnobserved}, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseover="UiActionsService.selectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex - 1]}, '${apiData.sequence[recordIndex - 2]}', 'mouseover', ${isUnobserved}, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseout="UiActionsService.unSelectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, '${strokeColor}')">${apiData.sequence[recordIndex - 2]}</text>`)
        });*/
        this.pathStrs.push(
            `<path 
                class="${pathEleClass}" stroke-width="${strokeWide}" stroke="${strokeColor}" d="${pathStr}" 
                data-stroke-color="${strokeColor}" 
                onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}', 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}, event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
                onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
                onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')">
            </path>`)
        this.circleStrs.push(
            `<circle class="circle_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex - 1]}" cx="${xVal}" cy="${yVal}" r="${font_size}" display="none" alignment-baseline="middle" stroke-width="${font_size/6}" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')"/>`)
        this.nucleotideStrs.push(
            `<text href="#${pathEleClass}" class="${pathEleClass}" x="${xVal}" y="${yVal}" font-size = "${font_size}px" onclick="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'click', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseover="UiActionsService.selectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, 'mouseover', ${isUnobserved}, '${apiData.sequence[recordIndex - 2]}', event, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseout="UiActionsService.unSelectNucleotide('${this.pluginOptions.pdbId}', '${this.pluginOptions.entityId}', ${apiData.label_seq_ids[recordIndex - 1]}, ${isUnobserved}, event, '${strokeColor}')">${apiData.sequence[recordIndex - 2]}</text>`)
        });

        let baseArray = FR3DData.annotations;
        this.baseStrs.set('cWW', [true, []]);
        this.baseStrs.set('tWW', [false, []]);
        this.baseStrs.set('tHH', [false, []]);
        this.baseStrs.set('cHH', [false, []]);
        this.baseStrs.set('tSS', [false, []]);
        this.baseStrs.set('cSS', [false, []]);
        this.baseStrs.set('cWH', [false, []]);
        this.baseStrs.set('tWH', [false, []]);
        this.baseStrs.set('tWS', [false, []]);
        this.baseStrs.set('cWS', [false, []]);
        this.baseStrs.set('tHS', [false, []]);
        this.baseStrs.set('cHS', [false, []]);
        
        baseArray.forEach((baseStr: any) => {
            let start:number = +baseStr.seq_id1
            //let start:number = +baseStr[`3d_id1`];
            //let end:number = +baseStr[`3d_id2`];
            let end:number = +baseStr.seq_id2
            if(baseStr && start && end) {
                let type:string = baseStr.bp
                let pathID:string = `rnaviewBP rnaviewBP_${this.pluginOptions.pdbId}_${this.pluginOptions.chainId} ${type}_${start}_${end}`
                let n1: string = baseStr.nt1
                let n2: string = baseStr.nt2
                let x1 = this.locations.get(start)![0] + font_size/2.5
                let x2 = this.locations.get(end)![0] + font_size/2.5
                let y1 = this.locations.get(start)![1] - font_size/2.5
                let y2 = this.locations.get(end)![1] - font_size/2.5
                let distance = Math.pow(Math.pow((x1-x2),2)+ Math.pow((y1-y2),2),0.5)
                let x1_prime = UiTemplateService.linearlyInterpolate(x1, x2, font_size/distance)
                let y1_prime = UiTemplateService.linearlyInterpolate(y1, y2, font_size/distance)
                let x2_prime = UiTemplateService.linearlyInterpolate(x1, x2, 1-font_size/distance)
                let y2_prime = UiTemplateService.linearlyInterpolate(y1, y2, 1-font_size/distance)
                let stroke = "#ccc"
                if (type.charAt(0) == 't') {
                    var fill = "none"
                } else {
                    var fill = "#ccc"
                }
                let xm = (x1_prime + x2_prime)/2
                let ym = (y1_prime + y2_prime)/2
                let distance2 = distance - 2 * font_size
                let height = font_size/1.5
                if(x1 - x2 != 0) {
                    var phi = 90 + Math.atan2((y1 - y2),(x1-x2)) * 180/Math.PI
                } else {
                    var phi = 0
                }
                if(type == 'cWW'){
                    if(n1 == 'G' && n2 == 'U' || n1 == 'U' && n2 == 'G') {
                        this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '#000', '#000');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                        d="
                        M ${(x1_prime + x2_prime)/2 - font_size/4}, ${(y1_prime+y2_prime)/2}
                        a ${font_size/4},${font_size/4} 0 1,0 ${font_size/2},0
                        a ${font_size/4},${font_size/4} 0 1,0 ${-1 * font_size/2},0
                        "
                        stroke="#000" stroke-width="${font_size/6} fill="${fill}"
                    />`)
                    } else{
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt,  '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '#000', '#000');" onmouseout="UiActionsService.hideTooltip('${pathID}');" stroke-width="${font_size/6}" data-stroke-color="#000" stroke="#000" d="M${x1_prime} ${y1_prime} ${x2_prime} ${y2_prime}"></path>`)
                    } 
                } else if (type == 'tWW') {
                    let xm1 = UiTemplateService.linearlyInterpolate(x1_prime, (x1_prime + x2_prime)/2, 1-(font_size/3)/(distance/2))
                    let ym1 = UiTemplateService.linearlyInterpolate(y1_prime, (y1_prime + y2_prime)/2, 1-(font_size/3)/(distance/2))
                    let xm2 = UiTemplateService.linearlyInterpolate((x1_prime + x2_prime)/2, x2_prime, (font_size/3)/(distance/2))
                    let ym2 = UiTemplateService.linearlyInterpolate((y1_prime + y2_prime)/2, y2_prime, (font_size/3)/(distance/2))
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                        d="
                        M ${x1_prime} ${y1_prime} ${xm1} ${ym1}
                        M ${(x1_prime + x2_prime)/2 - font_size/3} ${(y1_prime + y2_prime)/2}
                        a ${font_size/3},${font_size/3} 0 1,0 ${font_size/1.5},0
                        a ${font_size/3},${font_size/3} 0 1,0 ${-1 * font_size/1.5},0
                        M ${xm2} ${ym2} ${x2_prime} ${y2_prime}
                        "
                        stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}"/>`
                    )
                } else if (type == 'cSS'||type == 'tSS') {
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${xm} ${ym+distance2/2} ${xm} ${ym+height/2} 
                    l ${height/2} 0
                    l -${height/2} -${height} 
                    l -${height/2} ${height}
                    l ${height/2} 0
                    M ${xm} ${ym - height/2} ${xm} ${ym - distance2/2}
                    "stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
                } else if (type == 'tHS'|| type == 'cHS') {
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${xm} ${ym+distance2/2} ${xm} ${ym + height + height/4} 
                    h -${height/2}
                    v -${height}
                    h ${height}
                    v ${height}
                    h -${height/2}
                    M ${xm} ${ym + height/4} ${xm} ${ym - height/4}
                    l ${height/2} 0
                    l -${height/2} -${height} 
                    l -${height/2} ${height}
                    l ${height/2} 0
                    M ${xm} ${ym - height - height/4} ${xm} ${ym - distance2/2}
                    "stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
                } else if (type == 'tWS' || type == 'cWS') {
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${xm} ${ym+distance2/2} ${xm} ${ym + height + height/4} 
                    M ${xm - height/2} ${ym + 3*height/4} 
                    a ${height/2},${height/2} 0 1,0 ${height},0
                    a ${height/2},${height/2} 0 1,0 ${-1 * height},0
                    M ${xm} ${ym + height/4} ${xm} ${ym - height/4}
                    l ${height/2} 0
                    l -${height/2} -${height} 
                    l -${height/2} ${height}
                    l ${height/2} 0
                    M ${xm} ${ym - height - height/4} ${xm} ${ym - distance2/2}
                    "stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
                } else if (type == 'tWH' || type == 'cWH') {
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${xm} ${ym+distance2/2} ${xm} ${ym + height + height/4} 
                    M ${xm - height/2} ${ym + 3*height/4} 
                    a ${height/2},${height/2} 0 1,0 ${height},0
                    a ${height/2},${height/2} 0 1,0 ${-1 * height},0
                    M ${xm} ${ym + height/4} ${xm} ${ym - height/4}
                    h -${height/2}
                    v -${height}
                    h ${height}
                    v ${height}
                    h -${height/2}
                    M ${xm} ${ym - height - height/4} ${xm} ${ym - distance2/2}
                    "stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
                }
                else if (type == 'tHH' || type == 'cHH' ) {
                    this.baseStrs.get(type)![1].push(`<path class="${pathID}" onmouseover="UiActionsService.showTooltip(evt, '${type} Base Pair ${n1}${start} - ${n2}${end}', '${pathID}', '${stroke}', '${fill}');" onmouseout="UiActionsService.hideTooltip('${pathID}');"
                    d="
                    M ${xm} ${ym+distance2/2} ${xm} ${ym+height/2} 
                    h -${height/2}
                    v -${height}
                    h ${height}
                    v ${height}
                    h -${height/2}
                    M ${xm} ${ym - height/2} ${xm} ${ym - distance2/2}
                    "stroke="${stroke}" stroke-width="${font_size/6}" fill = "${fill}" transform = "rotate(${phi} ${xm} ${ym})"/>`)
                }
            }
        });
        this.baseStrs.forEach((value: [boolean, string[]], key: string) => {
            if(value[0]) {
                this.displayBaseStrs += value[1].join('');
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
        return  `<span class="pdb-rna-view-title">${this.pluginOptions.pdbId.toUpperCase()} Chain ${this.pluginOptions.chainId}</span>`;
    }

    private tooltip(): string {
        return  `<span class="pdb-rna-view-tooltip" id="${this.pluginOptions.pdbId}-rnaTopologyTooltip"></span>`;
    }

    private actionButtons(): string {
        return  `<div class="pdb-rna-view-btn-group">
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
        let errorContent = `<div class="pdb-rna-view-error">Error! Something went wrong!</div>`;
        if(type === 'apiError') {
            errorContent = `<div class="pdb-rna-view-error">
                RNA topology data for ${this.pluginOptions.pdbId.toUpperCase()} | ${this.pluginOptions.entityId} | ${this.pluginOptions.chainId.toUpperCase()} is not available!
            </div>`;
        }

        this.containerElement.innerHTML = `<div class="pdb-rna-view-container">${errorContent}</div>`;
    }
}