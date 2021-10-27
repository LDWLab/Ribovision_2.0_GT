import { UiActionsService } from './uiActions'
import { PluginOptions, ApiData } from './data';

export class UiTemplateService {
    private containerElement: HTMLElement;
    private pluginOptions: PluginOptions;
    private uiActionsService: UiActionsService;

    constructor(containerElement: HTMLElement, pluginOptions: PluginOptions) {
        this.containerElement = containerElement;
        this.pluginOptions = pluginOptions;
        this.uiActionsService = new UiActionsService(this.pluginOptions.pdbId);
    }

    render(apiData: ApiData) {
        this.containerElement.innerHTML = `<div class="pdb-rna-view-container pdb-rna-view-container-${this.pluginOptions.pdbId}">
            ${this.svgTemplate(apiData)}
            ${this.title()}
            ${this.tooltip()}
            ${this.actionButtons()}
        </div>`;
        //this.fixOverlaps(apiData)
        this.uiActionsService.applyButtonActions();
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
                    console.log("1" + apiData.label_seq_ids[j + 1])
                    nextEle.setAttribute("y",Number(nextEle.getAttribute("y")) - Number(distances[1]));
                    yAdjust = yAdjust - Number(distances[1])
                } else if(lowest === 3) {
                    console.log("3" + apiData.label_seq_ids[j + 1])
                    nextEle.setAttribute("y",Number(nextEle.getAttribute("y")) + Number(distances[3]));
                    yAdjust = yAdjust + Number(distances[3])
                } else if(lowest === 0) {
                    console.log("0" + apiData.label_seq_ids[j + 1])
                    nextEle.setAttribute("x",Number(nextEle.getAttribute("x")) - Number(distances[0]));
                    xAdjust = xAdjust - Number(distances[0])
                } else if(lowest === 2) {
                    console.log("2" + apiData.label_seq_ids[j + 1])
                    nextEle.setAttribute("x",Number(nextEle.getAttribute("x")) + Number(distances[2]));
                    xAdjust = xAdjust + Number(distances[2])
                } 
            }
        });
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

    private svgTemplate(apiData: ApiData): string {
        let pathStrs: string[] = [];  
        const font_size:number = this.calculateFontSize(apiData)
        const lastPathIndex = apiData.svg_paths.length - 1;
        apiData.svg_paths.forEach((pathStr: string, recordIndex: number) => {
            if(recordIndex === 0 || recordIndex === lastPathIndex) return;
            const pathEleClass = `rnaviewEle rnaviewEle_${this.pluginOptions.pdbId} rnaview_${this.pluginOptions.pdbId}_${apiData.label_seq_ids[recordIndex]}`;
            let strokeColor = this.pluginOptions.theme?.color || '#323232';
            let isUnobserved = false;
            if(apiData.unobserved_label_seq_ids && apiData.unobserved_label_seq_ids.indexOf(apiData.label_seq_ids[recordIndex]) > -1) {
                strokeColor = this.pluginOptions.theme?.unobservedColor || '#ccc';
                isUnobserved = true;
            }
            let pathStrParsed:string[] = pathStr.split('M').join(',').split(',')
            let xVal:number = (Number(pathStrParsed[1])+Number(pathStrParsed[3]))/2 
            let yVal:number = (Number(pathStrParsed[2])+Number(pathStrParsed[4]))/2 
            pathStrs.push(`<text href="#${pathEleClass}" class="${pathEleClass}" x="${xVal}" y="${yVal}" font-size = "${font_size}px" onclick="UiActionsService.selectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex]}, '${apiData.sequence[recordIndex - 1]}', 'click', ${isUnobserved}, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseover="UiActionsService.selectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex]}, '${apiData.sequence[recordIndex - 1]}', 'mouseover', ${isUnobserved}, ${this.pluginOptions.theme?.highlightColor ? "'"+this.pluginOptions.theme.highlightColor+"'" : undefined})" 
            onmouseout="UiActionsService.unSelectPath(event, '${this.pluginOptions.pdbId}', ${apiData.label_seq_ids[recordIndex]}, ${isUnobserved}, '${strokeColor}')">${apiData.sequence[recordIndex - 1]}</text>`)
        });

        return `
        <div style="width:100%;height:100%;z-index:0;position:absolute;">
            <svg class="rnaTopoSvgSelection rnaTopoSvgSelection_${this.pluginOptions.pdbId}" 
            preserveAspectRatio="xMidYMid meet" 
            viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
            style="width:100%;height:100%;position:relative;cursor:default;"></svg>
        </div>
        <div style="width:100%;height:100%;z-index:1;position:absolute;">
            <svg class="rnaTopoSvgHighlight rnaTopoSvgHighlight_${this.pluginOptions.pdbId}" 
            preserveAspectRatio="xMidYMid meet" 
            viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
            style="width:100%;height:100%;position:relative;"></svg>
        </div>
        <div style="width:100%;height:100%;z-index:2;position:absolute;">
            <svg class="rnaTopoSvg rnaTopoSvg_${this.pluginOptions.pdbId}" 
                preserveAspectRatio="xMidYMid meet" 
                viewBox="0 0 ${apiData.dimensions.width} ${apiData.dimensions.height}" 
                style="width:100%;height:100%;">${pathStrs.join('')}
            </svg>
        </div>`;
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