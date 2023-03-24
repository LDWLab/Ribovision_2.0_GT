import { CustomEvents } from './customEvents';

declare var d3: any;
export class UiActionsService {
    pdbId: string;
    static pdbevents: any = CustomEvents.create(['PDB.RNA.viewer.click', 'PDB.RNA.viewer.mouseover', 'PDB.RNA.viewer.mouseout']);
    static lastSelected: number = -1
    static selected: number = -1
    static selectedPath: string = ''
    static selectedBPColor: string = ''
    static fillBPColor: string = ''
    static selectedColor: string|null = '' 
    static expanded: boolean = false
    constructor(pdbId: string) {
        this.pdbId = pdbId;
    }

    applyButtonActions() {

        const svgEle = d3.select('.rnaTopoSvg');
        // initialize SVG zoom behaviour and remove mouse wheel zoom events
        const zoom1 = d3.zoom().on('zoom', () => {
            d3.select(`.rnaTopoSvg_${this.pdbId}`).attr("transform", d3.event.transform)
            d3.select(`.rnaTopoSvgHighlight_${this.pdbId}`).attr("transform", d3.event.transform)
            d3.select(`.rnaTopoSvgSelection_${this.pdbId}`).attr("transform", d3.event.transform)
        });
        d3.select(`.rnaTopoSvg`).call(zoom1)
        .on('dblclick.zoom', null)
        .on('wheel.zoom', null)
        .on('mousewheel.zoom', null);

        // zoom-in button behaviour
        d3.select(`#rnaTopologyZoomIn-${this.pdbId}`)
        .on("click", () =>{
            d3.event.stopPropagation();
            zoom1.scaleBy(svgEle.transition().duration(300), 1.2);
        });

        // zoom-out button behaviour
        d3.select(`#rnaTopologyZoomOut-${this.pdbId}`)
        .on("click", () => {
            d3.event.stopPropagation();
            const rnaTopoSvg = d3.select(`.rnaTopoSvg_${this.pdbId}`);
            const transformValue = rnaTopoSvg._groups[0][0].getAttribute('transform');
            if(transformValue && transformValue !== '') {
                if(transformValue === null) return;
                const transformValMatch = +transformValue.match(/.+scale\((.*)\)/)[1];
                if(transformValMatch <= 1 || (transformValMatch - 0.3 <= 1)) {
                    svgEle.transition().duration(300).call(zoom1.transform, d3.zoomIdentity);
                    return;
                }
            }
            zoom1.scaleBy(d3.select(`.rnaTopoSvg`).transition().duration(300), 0.8);
        });

        // reset button behaviour
        d3.select(`#rnaTopologyReset-${this.pdbId}`)
        .on("click", () => {
            d3.event.stopPropagation();
            svgEle.transition().duration(300).call(zoom1.transform, d3.zoomIdentity);
        });

        // selection and highlight reset on canvas click
        d3.select(`.pdb-rna-view-container-${this.pdbId}`)
        .on("click", () => {
            d3.event.stopPropagation();
            UiActionsService.clearHighlight(this.pdbId);
            UiActionsService.clearSelection(this.pdbId);
        });

    }

    static selectNucleotide = function(pdbId: string, entityId: string, label_seq_id: number, eventType: 'click' | 'mouseover', isUnobserved: boolean, residue?: string, event?: MouseEvent, color?: string): void {
        event?.stopImmediatePropagation();
        if(<SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)) {
            UiActionsService.clearHighlight(pdbId);
            UiActionsService.colorNucleotide(pdbId, label_seq_id, color, 'highlight');
            UiActionsService.selected = label_seq_id;
            UiActionsService.lastSelected = label_seq_id;
            const ttEle = document.getElementById(`${pdbId}-rnaTopologyTooltip`);
            ttEle!.style.display = 'inline';
            ttEle!.innerHTML = `<strong>${isUnobserved ? 'Unobserved ' : ''}Residue ${residue ? residue : ''} ${label_seq_id}</strong>`;

            if(!isUnobserved) {
                const eventName: any = (eventType === 'click') ? 'PDB.RNA.viewer.click' : 'PDB.RNA.viewer.mouseover';
                const evData = { pdbId, label_seq_id, entityId }
                const textElement: any = document.querySelector(`.rnaview_${pdbId}_${label_seq_id}`);
                if(textElement) {
                    CustomEvents.dispatchCustomEvent(UiActionsService.pdbevents[eventName], evData, textElement);
                }
            }
        }
    }

    static unSelectNucleotide(pdbId: string, entityId: string, label_seq_id: number, isUnobserved: boolean, event?: any, color?: string): void {
        event?.stopImmediatePropagation();
        UiActionsService.clearHighlight(pdbId);
        const ttEle = document.getElementById(`${pdbId}-rnaTopologyTooltip`);
        if (ttEle) {
            ttEle.style.display = 'none';
        }

        if(!isUnobserved) {
            if(label_seq_id == undefined) {
                label_seq_id = this.lastSelected;
            }
            if(label_seq_id > -1) {
                const evData = { pdbId, label_seq_id, entityId }
                const textElement: any = document.querySelector(`.rnaview_${pdbId}_${label_seq_id}`);
                if(textElement) {
                    CustomEvents.dispatchCustomEvent(UiActionsService.pdbevents['PDB.RNA.viewer.mouseout'], evData, textElement);
                }
            }
        }
    }

    static colorNucleotide = function(pdbId: string, label_seq_id: number, color?: string, type?: 'selection' | 'highlight', shape?: string) {
        //const textElement: any = document.querySelector(`.rnaview_${pdbId}_${label_seq_id}`);
        //const textEleStrokeWidth = parseInt(textElement.getAttribute('stroke-width'));

        let strokeColor = color || 'orange';
        //const svgClassPrefix = (type === 'highlight') ? 'rnaTopoSvgHighlight' : 'rnaTopoSvgSelection';
        //const highlightSvg = document.querySelector(`.${svgClassPrefix}_${pdbId}`);
        /*
        const newTextEle = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        newTextEle.setAttribute('stroke-width', (textEleStrokeWidth * 3).toString());
        newTextEle.setAttribute('stroke', strokeColor);
        newTextEle.setAttribute('fill', strokeColor);
        highlightSvg!.appendChild(newTextEle);
        */
        if (shape == 'circle') {    
            const circle = <SVGSVGElement>document.querySelector(`svg.rnaTopoSvg`)!.getElementsByClassName(`circle_${pdbId}_${label_seq_id}`)[0];
            if(type == 'highlight') {
                if(circle.getAttribute("stroke")) {
                    UiActionsService.selectedColor = circle.getAttribute("stroke")
                }
            }
            if(circle) {
                circle.setAttribute("stroke", `${color}`);
                circle.setAttribute("fill", `${color}`);
                circle.style.display = "block";
            }
        } 
            let nucleotide = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${label_seq_id}`)[0];
            if(nucleotide) {
                if(nucleotide.nodeName == 'path') {
                    if(type == 'highlight') {
                        if(nucleotide.getAttribute("stroke")) {
                            UiActionsService.selectedColor = nucleotide.getAttribute("stroke")
                        } else {
                            UiActionsService.selectedColor = "#323232"
                        }
                    }
                    nucleotide.setAttribute("stroke",strokeColor);
                } else if(nucleotide.nodeName == 'text'){
                    if(type == 'highlight') {
                        if(nucleotide.getAttribute("fill")) {
                            UiActionsService.selectedColor = nucleotide.getAttribute("fill")
                        } else {
                            UiActionsService.selectedColor = "#323232"
                        }
                    }
                    nucleotide.setAttribute("fill",strokeColor);
                }
            }
        
    }
    static showCheckboxes() {
        var checkboxes = document.getElementById("checkboxes");
        if (!this.expanded) {
            checkboxes!.style.display = "block";
            this.expanded = true;
        } else {
            checkboxes!.style.display = "none";
            this.expanded = false;
        }
    }
    static clearHighlight(pdbId: string) {
        if (this.selected > -1) {
            var nucleotide = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${this.selected}`)[0]
            if(this.selectedColor && nucleotide) {
                if(nucleotide.nodeName == 'text') {
                    nucleotide.setAttribute("fill",UiActionsService.selectedColor);
                } else {
                   nucleotide.setAttribute("stroke",UiActionsService.selectedColor);
                }
            } else if(nucleotide){
                if(nucleotide.nodeName == 'text') {
                    nucleotide.setAttribute("fill","#323232");
                } else {
                   nucleotide.setAttribute("stroke","#323232");
                }
            }
            this.selected = -1;
        }
        //document.querySelector(`.rnaTopoSvgHighlight_${pdbId}`)!.innerHTML = "";
    }

    static clearSelection(pdbId: string) {
        document.querySelector(`.rnaTopoSvgSelection_${pdbId}`)? document.querySelector(`.rnaTopoSvgSelection_${pdbId}`)!.innerHTML = "":null;
    }

    static showTooltip(evt: any, text: string, path: string, color: string, fill: string) {
        const tooltip = document.getElementById("tooltip");
        tooltip!.id = "tooltip"
        tooltip!.innerHTML = text;
        tooltip!.style.display = "block";
        tooltip!.style.left = evt.layerX + 'px';
        tooltip!.style.top = evt.layerY + 'px'; 
        UiActionsService.selectedBPColor = color;
        UiActionsService.fillBPColor = fill;
        let node = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(path)[0]
        if(node.nodeName != 'text') {
            node.setAttribute("stroke", "orange");
        } 
        node.setAttribute("fill", "orange");
      }
      
    static hideTooltip(path: string) {
        var tooltip = document.getElementById("tooltip");
        tooltip!.style.display = "none";
        let node = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(path)[0]
        if(node.nodeName != 'text') {
            node.setAttribute("stroke", `${UiActionsService.selectedBPColor}`);
        } 
        node.setAttribute("fill", `${UiActionsService.fillBPColor}`);
    }

    static drawCircle(index: number, color: string, pdbId: string) {
        const circle = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`circle_${pdbId}_${index}`)[0]
        const nucleotide = (<any>document.querySelector(`svg.rnaTopoSvg`))!.getElementsByClassName(`rnaviewEle rnaviewEle_${pdbId} rnaview_${pdbId}_${index}`)[0]
        const BBox = nucleotide.getBBox()
        const nx = (BBox.x + BBox.width/2)
        const ny = (BBox.y + BBox.height/2)
        circle.setAttribute("cx", nx)
        circle.setAttribute("cy", ny)
        circle.setAttribute("stroke", `${color}`);
        circle.setAttribute("fill", `${color}`);
        circle.style.display = "block";
    }
}

(window as any).UiActionsService = UiActionsService;