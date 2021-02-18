
const interpolateLinearly = (window as any).interpolateLinearly;
const RdPu = (window as any).RdPu;
const YlGn = (window as any).YlGn;
const mapped_aa_properties = (window as any).mapped_aa_properties;
const aaPropertyConstants = (window as any).aaPropertyConstants;
const aaColorData = (window as any).aaColorData;
const masking_range_array = (window as any).masking_range_array;
var masked_array = (window as any).masked_array;
var parsePVData = (window as any).parsePVData;
var viewerInstance = (window as any).viewerInstance;
var selectSections_RV1 = (window as any).selectSections_RV1;
var rv3VUEcomponent = (window as any).vm;
var filterRange = (window as any).filterRange?(window as any).filterRange: "-10000,10000";

class PdbTopologyViewerPlugin { 
    
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

    displayStyle = 'border:1px solid #696969;';
    errorStyle = 'border:1px solid #696969; height:54%; padding-top:46%; text-align:center; font-weight:bold;';
    menuStyle = 'position:relative;height:38px;line-height:38px;background-color:#96c9dc;padding: 0 10px;font-size:16px; color: black;';

    sequenceArr: string[];
    entityId: string;
    entryId: string;
    //filterRange: string;
    alreadyRan: boolean;
    chainId: string;
    apiData: any;
    targetEle: HTMLElement;
    pdbevents: any

    xScale: any;
    yScale: any;
    zoom: any;
    scaledPointsArr: any[];
    domainTypes: any[];
    svgWidth = 100;
    svgHeight = 100;

    svgEle: any;

    subscribeEvents = true;

    render(target: HTMLElement, options:{entityId: string, entryId: string, filterRange?: string, chainId?: string, subscribeEvents?:boolean, displayStyle?: string, errorStyle?: string, menuStyle?: string}) {
        if(options && typeof options.displayStyle != 'undefined' && options.displayStyle != null) this.displayStyle += options.displayStyle;
        if(options && typeof options.errorStyle != 'undefined' && options.errorStyle != null) this.errorStyle += options.errorStyle;
        if(options && typeof options.menuStyle != 'undefined' && options.menuStyle != null) this.menuStyle += options.menuStyle;
        //if(options && typeof options.filterRange != 'undefined' && options.filterRange != null) filterRange = options.filterRange;
        
        this.targetEle = <HTMLElement> target;
        if(this.targetEle) this.targetEle.innerHTML = '';
        if(!target || !options || !options.entryId || !options.entityId){ 
            this.displayError('param');
            return;
        }
        if(options.subscribeEvents == false) this.subscribeEvents = false;
        this.entityId = options.entityId;
        this.entryId = options.entryId.toLowerCase();
        //If chain id is not provided then get best chain id from observed residues api
        if(typeof options.chainId == 'undefined' || options.chainId == null){
            this.getObservedResidues(this.entryId).then((result) => {
                if(typeof result != 'undefined' && typeof result[this.entryId] != 'undefined' && typeof result[this.entryId][this.entityId] != 'undefined'){
                    this.chainId = result[this.entryId][this.entityId][0].chain_id;
                    this.initPainting();
                }else{
                    this.displayError();
                }
            });
        }else{
            this.chainId = options.chainId;
            this.initPainting()
        }
        
    }

    initPainting(){
        if(!this.alreadyRan){
            this.alreadyRan = true;
            this.getApiData(this.entryId, this.chainId).then(result => {
                if(result){

                    //Validate required data in the API result set (0, 2, 4)
                    if(typeof result[0] == 'undefined' || typeof result[2] == 'undefined' || typeof result[4] == 'undefined'){ 
                        this.displayError();
                        return;
                    }
                    this.apiData = result;
                    //default pdb events
		    	    this.pdbevents = this.createNewEvent(['PDB.topologyViewer.click','PDB.topologyViewer.mouseover','PDB.topologyViewer.mouseout']);
                    this.getPDBSequenceArray(this.apiData[0][this.entryId]);
                    this.drawTopologyStructures();
                    this.createDomainDropdown();

                    if(this.subscribeEvents) this.subscribeWcEvents();

                }else{

                }
            });
        }
    }

    displayError(errType?: string){
        let errtxt = "Error: Data not available!"
        if(errType == 'param') errtxt = "Error: Invalid Parameters!"
        if(this.targetEle) this.targetEle.innerHTML = `<div style="${this.errorStyle}">${errtxt}</div>`;
    }

    createNewEvent = function(eventTypeArr: string[]){
		let eventObj:any = {};
		eventTypeArr.forEach((eventType, index) => {
			let event; 
			if (typeof MouseEvent == 'function') {
				// current standard
				event = new MouseEvent(eventType, { 'view': window, 'bubbles': true, 'cancelable': true });
			
			} else if (typeof document.createEvent == 'function') {
				// older standard
				event = document.createEvent('MouseEvents');
				event.initEvent(eventType, true /*bubbles*/, true /*cancelable*/);
			
			}
			
			eventObj[eventType] = event;
		});
		
		return eventObj;
    }
    
    async getObservedResidues(pdbId: string) {
        try {
            return await (await fetch(`https://www.ebi.ac.uk/pdbe/api/pdb/entry/observed_residues_ratio/${pdbId}`)).json();
        } catch (e) {
          console.log(`Couldn't load UniProt variants`, e);
        }
    }

    async getApiData(pdbId: string, chainId: string) {
        const dataUrls = [
            `https://www.ebi.ac.uk/pdbe/api/pdb/entry/entities/${pdbId}`,
            `https://www.ebi.ac.uk/pdbe/api/mappings/${pdbId}`,
            `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdbId}`,
            `https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${pdbId}`,
            `https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/${pdbId}/chain/${chainId}`
        ]
        return Promise.all(dataUrls.map(url => fetch(url)))
        .then(resp => Promise.all( 
                resp.map((r) => { 
                    if(r.status == 200){
                        return r.json();
                    }else{
                        return undefined;
                    }
                    
                }) 
            )
        )
    }

    getPDBSequenceArray(entities: any[]) {
        const totalEntities = entities.length;
        for(let i=0; i < totalEntities; i++){
            if(entities[i].entity_id == this.entityId){
                this.sequenceArr = entities[i].sequence.split('');
            }
        }
    }

    chunkArray(arr: any[], len: number) {
			
        let chunks = [], i = 0,	n = arr.length;
        while (i < n) {
            chunks.push(arr.slice(i, i += len));
        }
        return chunks;
    }

    getDomainRange(){
        let allCordinatesArray: any[] = [];
         //filtering by a domain

        //update domain start and end positions here
        const topologyData = this.apiData[2][this.entryId][this.entityId][this.chainId];
        const istart=Number(filterRange.split(",")[0]);
        const iend=Number(filterRange.split(",")[1]);

    //if (istart == 0 && iend == 100000) {istart = parseInt (t["terms"][0]["resnum"]); iend = parseInt (t["terms"][1]["resnum"]) }

        if ( iend > 5000) {
            var tc = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"];
            //var te = this.apiData[2][this.entryId][this.entityId][this.chainId]["extents"];
            var th = this.apiData[2][this.entryId][this.entityId][this.chainId]["helices"];
            var ts = this.apiData[2][this.entryId][this.entityId][this.chainId]["strands"];
    
            //var ths = this.apiData[2][this.entryId][this.entityId][this.chainId]["helices"][0]["start"];
            //var the = this.apiData[2][this.entryId][this.entityId][this.chainId]["helices"][0]["stop"];

    for( var k = 0; k < th.length; k++){ if ((( th[k]["start"] < istart) && ( th[k]["stop"] < istart))||(( th[k]["start"] > iend) && ( th[k]["stop"] > iend))) { this.apiData[2][this.entryId][this.entityId][this.chainId]["helices"].splice(k, 1); k--; }}
    for( var k = 0; k < ts.length; k++){ if ((( ts[k]["stop"] < istart)) || (( ts[k]["start"] > iend))) { this.apiData[2][this.entryId][this.entityId][this.chainId]["strands"].splice(k, 1); k--; }}
    for( var k = 0; k < tc.length; k++){ if ((( tc[k]["stop"] < istart) && ( tc[k]["start"] !=-1 ))||(( tc[k]["start"] > iend) && ( tc[k]["stop"] > iend))) { this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"].splice(k, 1); k--; }}

    tc = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"];
    //console.log(tc);

    for( var k = 0; k < tc.length; k++){ if (( tc[0]["start"] == -1) ) { this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"].splice(0, 1); k--; }}
    for( var k = tc.length-1; k >= 0; k--){ if (( tc[tc.length-1]["start"] == -1) ) { this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"].splice(tc.length-1, 1); k--; }}
 

    tc = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"]; 
    //var tca = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["start"];
    //var tco = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["stop"];
    //var tcp = this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["path"];


    //for( var k = 0; k < 1; k++){ if ((( tc[0]["start"] < istart) && ( tc[0]["stop"] >= istart) )) { del=istart - tc[0]["start"]-1;  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["path"].splice(0, (del-1)*2);  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["start"]=this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["start"]+del;  }}
    //for( var k = tc.length-1; k > tc.length-2; k--){ if ((( tc[tc.length-1]["start"] <= iend) &&( tc[tc.length-1]["stop"] >= iend) )) { del=tc[tc.length-1]["stop"]-1-iend; path_length=tc[tc.length-1]["path"].length; console.log(del, path_length); this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["path"].splice(path_length-3, 2);  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["stop"]=this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["stop"]-del;  }}
    
    for( var k = 0; k < 1; k++){
         if ((( tc[0]["start"] < istart) && ( tc[0]["stop"] >= istart) )){ 
            var del=istart - tc[0]["start"]-1;
            var del1=del;
            if (del>2){del1=2};  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["path"].splice(0, (del1-1)*2);  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["start"]=this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["start"]+del+1;  }}
    for( var k = tc.length-1; k > tc.length-2; k--){ if ((( tc[tc.length-1]["start"] <= iend) &&( tc[tc.length-1]["stop"] >= iend) )) { del=tc[tc.length-1]["stop"]-1-iend; del1=del; if (del>2){del1=2}; var path_length=tc[tc.length-1]["path"].length; console.log(del, path_length); this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["path"].splice(path_length-(del1-1)*2-1, (del1-1)*2);  this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["stop"]=this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][tc.length-1]["stop"]-del-1;  }}
   
    //console.log(tc);
 
    //this.apiData[2][this.entryId][this.entityId][this.chainId]["coils"][0]["start"][0]=54;
    };
    //console.log(t);
        for(let secStrType in topologyData){
        
            if(topologyData[secStrType]){
                // iterating on secondary str data array to get array spliced in x,y 
                topologyData[secStrType].forEach((secStrData: any) => {
                    if(typeof secStrData.path !== 'undefined' && secStrData.path.length > 0){
                        allCordinatesArray= allCordinatesArray.concat(this.chunkArray(secStrData.path, 2));
                    }
                });
            }
                            
        };
        
        this.xScale = d3.scaleLinear()
                        .domain([d3.min(allCordinatesArray, function(d) { return d[0]; }), d3.max(allCordinatesArray, function(d) { return d[0]; })])
                        .range([1, this.svgWidth - 1]);
        
        this.yScale = d3.scaleLinear()
                        .domain([d3.min(allCordinatesArray, function(d) { return d[1]; }), d3.max(allCordinatesArray, function(d) { return d[1]; })])
                        .range([1, this.svgHeight - 1]);
        
        this.zoom = d3.zoom()
                    .on("zoom", () => this.zoomDraw())
                    //.scaleExtent([.5, 20])  // This control how much you can unzoom (x0.5) and zoom (x20)
                    // .transform(this.xScale, this.yScale)
                      
    }

    drawStrandSubpaths(startResidueNumber:number, stopResidueNumber:number, index:number) {
        const _this = this;
        const totalAaInPath = (stopResidueNumber - startResidueNumber) + 1
        const subPathHeight = (this.scaledPointsArr[7] - this.scaledPointsArr[1])/totalAaInPath;
        
        //create subsections/paths
        let dValArr = [];
        for(let subPathIndex=0; subPathIndex<totalAaInPath; subPathIndex++){
            let subPathObj:any = {type: 'strands', elementIndex: index};
            if(subPathIndex === 0){
                subPathObj['residue_number'] = startResidueNumber;
                subPathObj['pathData'] = [
                this.scaledPointsArr[4], this.scaledPointsArr[1],
                this.scaledPointsArr[4], this.scaledPointsArr[1] + subPathHeight,
                this.scaledPointsArr[8], this.scaledPointsArr[1] + subPathHeight,
                this.scaledPointsArr[8], this.scaledPointsArr[13]
                ];
            }else{
                subPathObj['residue_number'] = startResidueNumber + subPathIndex;
                subPathObj['pathData'] = [
                dValArr[subPathIndex - 1]['pathData'][2], dValArr[subPathIndex - 1]['pathData'][3],
                dValArr[subPathIndex - 1]['pathData'][2], dValArr[subPathIndex - 1]['pathData'][3] + subPathHeight,
                dValArr[subPathIndex - 1]['pathData'][4], dValArr[subPathIndex - 1]['pathData'][5] + subPathHeight,
                dValArr[subPathIndex - 1]['pathData'][4], dValArr[subPathIndex - 1]['pathData'][5]
                ];
            }
            dValArr.push(subPathObj);
        }
        
        
        this.svgEle.selectAll('.subpath-strands'+index).remove();
        
        this.svgEle.selectAll('.subpath-strands'+index)
        .data(dValArr)
        .enter()
        .append('path')  
        .attr('class', (d: any, i:number) => { return 'strandsSubPath subpath-strands'+index+' topo_res_'+d.residue_number })
        .attr('d', (d:any,i:number) => { return 'M '+d.pathData.join(' ')+' Z' })
        .attr('stroke', '#111')
        .attr('stroke-width', '0')
        .attr('fill', 'white')
        .attr('fill-opacity','0')
        .on('mouseover', function(d:any){ _this.mouseoverAction(this, d); })
        .on('mousemove', function(d:any){ _this.mouseoverAction(this, d); })
        .on('mouseout', function(d:any){ _this.mouseoutAction(this, d); })
        .on("click", function(d:any){ _this.clickAction(d); })
        
    }

    drawStrandMaskShape(index:number) {
        let maskPointsArr = this.scaledPointsArr;
        
        const adjustmentFactor = 0.3
        let adjustIndexAddArr = [7,8,10,12];
        let adjustIndexSubtractArr = [0,1,2,3,4,5,9,11,13];
        //For arrow pointing upwards
        if(maskPointsArr[0] > maskPointsArr[6]){
            adjustIndexAddArr = [0,1,2,3,4,5,9,11,13];
            adjustIndexSubtractArr = [7,8,10,12];
        }
        
        
        let addIndexLength = adjustIndexAddArr.length;
        for(let maskPtIndex = 0; maskPtIndex < addIndexLength; maskPtIndex++){
            maskPointsArr[adjustIndexAddArr[maskPtIndex]] = maskPointsArr[adjustIndexAddArr[maskPtIndex]] + adjustmentFactor;
        }
        
        let subtractIndexLength = adjustIndexSubtractArr.length;
        for(let maskPtIndex1 = 0; maskPtIndex1 < subtractIndexLength; maskPtIndex1++){
            maskPointsArr[adjustIndexSubtractArr[maskPtIndex1]] = maskPointsArr[adjustIndexSubtractArr[maskPtIndex1]] - adjustmentFactor;
        }
        
    
        //Add the outer points          
        maskPointsArr[14] = maskPointsArr[8]
        maskPointsArr[15] = maskPointsArr[13]
        maskPointsArr[16] = maskPointsArr[8]
        maskPointsArr[17] = maskPointsArr[7]
        maskPointsArr[18] = maskPointsArr[4]
        maskPointsArr[19] = maskPointsArr[7]
        maskPointsArr[20] = maskPointsArr[4]
        maskPointsArr[21] = maskPointsArr[1]
        
        this.svgEle.selectAll('.maskpath-strands'+index).remove()
        
        this.svgEle.selectAll('.maskpath-strands'+index)
        .data([maskPointsArr])
        .enter()
        .append('path')
        .attr('class', (d:any,i:number) => { return 'strandMaskPath maskpath-strands'+index })
        .attr('d', (d:any,i:number) => { return 'M'+maskPointsArr.join(' ')+'Z' })
        .attr('stroke', '#111')
        .attr('stroke-width', 0.3)
        .attr('fill', 'white')
        .attr('stroke-opacity', 0)
        
    }

    renderTooltip(elementData: any, action: string) {
        
        let toolTipEle = d3.select('.pdbTopologyTooltip') as any;
        if(toolTipEle._groups[0][0] == null){
            toolTipEle = d3.select('body').append('div').attr('class','pdbTopologyTooltip').attr('style', 'display: none;width: auto;position: absolute;background: #fff;padding: 5px;border: 1px solid #666;border-radius: 5px;box-shadow: 5px 6px 5px 0 rgba(0,0,0,.17);font-size: .9em;color: #555;z-index: 998;');
        }
        
        if(action === 'show'){
            const x = d3.event.pageX, y = d3.event.pageY;
            
            let tooltipContent = 'Residue ' + elementData.residue_number + ' (' + this.sequenceArr[elementData.residue_number - 1] + ')'; 
            
            if(typeof elementData.tooltipMsg !== 'undefined'){
                if(typeof elementData.tooltipPosition !== 'undefined' && elementData.tooltipPosition === 'postfix'){
                    tooltipContent = tooltipContent+' '+elementData.tooltipMsg;
                }else{
                    tooltipContent = elementData.tooltipMsg+' '+tooltipContent;
                }
            }

            toolTipEle
                    .html(tooltipContent)
                    .style('display','block')
                    .style('top', y + 15 +'px')
                    .style('left', x + 10 +'px')
        }else{
            toolTipEle.style('display','none')
        }
    };
    dispatchEvent(eventType:any, eventData:any, eventElement?:HTMLElement) {
        let dispatchEventElement = this.targetEle;
        if(typeof eventElement !== 'undefined'){
            dispatchEventElement = eventElement;
        }
        if(typeof eventData !== 'undefined'){
            this.pdbevents[eventType]['eventData'] = eventData;
        }
        dispatchEventElement.dispatchEvent(this.pdbevents[eventType])
    };
    clickAction(eleObj:any) {
        //Dispatch custom click event
        console.log(eleObj);
        if(masked_array[eleObj.residue_number - 1] == true) {
            this.dispatchEvent('PDB.topologyViewer.click', {
                residueNumber: eleObj.residue_number,
                type: eleObj.type,
                entryId: this.entryId,
                entityId: this.entityId,
                filterRange: filterRange,
                chainId: this.chainId,
                // structAsymId: this.bestStructAsymId
            });
        }
    }
    mouseoverAction(eleObj:any|this, eleData:any) {

        const selectedPath = d3.select(eleObj);
        //var selectedPathData = selectedPath.data();

        //Highlight Residue
        //Check if residue is not masked
        if(eleData.tooltipMsg != 'NaN') {
            //Show Tooltip
            this.renderTooltip(eleData, 'show');
            if(eleData.type === 'strands' || eleData.type === 'helices'){
                selectedPath.attr('fill', this.defaultColours.mouseOver).attr('fill-opacity','0.3')
            }if(eleData.type === 'coils'){
                selectedPath.attr('stroke', this.defaultColours.mouseOver).attr('stroke-width', 1);
            }
            //Maybe add gradient line movement here as function of eleData.tooltipMsg
            //var newLine = document.createElementNS('http://www.w3.org/2000/svg','line');
            //newLine.setAttribute('id','line2');
            //newLine.setAttribute('x1','22');
            //newLine.setAttribute('y1','12');
            //newLine.setAttribute('x2','55');
            //newLine.setAttribute('y2','474');
            //newLine.setAttribute("stroke", "black")
            //document.querySelector('#imageaef08aed83').appendChild(newLine)
        
            //Dispatch custom mouseover event
            this.dispatchEvent('PDB.topologyViewer.mouseover', {
                residueNumber: eleData.residue_number,
                type: eleData.type,
                entryId: this.entryId,
                entityId: this.entityId,
                filterRange: filterRange,
                chainId: this.chainId,
                // structAsymId: scope.bestStructAsymId
            });
        }
    }
    mouseoutAction(eleObj:any, eleData:any) {
        let mouseOverColor = 'white';
        let fillOpacity = 0;
        let strokeOpacity = 1;
        let strokeWidth = 0.3;
        const pathElement = d3.select(eleObj);
        
        //Hide Tooltip
        this.renderTooltip('', 'hide');
        
        //if path colour is changed then get the colour
        if(pathElement.classed('coloured')){
            if(eleData.type === 'coils' && (masked_array.length != 0 && masked_array[eleData.residue_number] == false)){
                mouseOverColor = this.defaultColours.borderColor;
            }else{
                mouseOverColor = pathElement.attr('data-color');
                fillOpacity = 1;
                strokeWidth = 1;
            }
        }else{
            if(eleData.type === 'coils'){
                mouseOverColor = this.defaultColours.borderColor;
            }
        }
        
        if(eleData.type === 'strands' || eleData.type === 'helices'){
            pathElement.attr('fill',mouseOverColor).attr('fill-opacity', fillOpacity)
        }if(eleData.type === 'coils'){
            pathElement.attr('stroke',mouseOverColor).attr('stroke-opacity', strokeOpacity);
            pathElement.attr('stroke',mouseOverColor).attr('stroke-width', strokeWidth);
        }
        
        //Dispatch custom mouseover event
        this.dispatchEvent('PDB.topologyViewer.mouseout', {
            filterRange: filterRange,
            entryId: this.entryId,
            entityId: this.entityId,
            chainId: this.chainId,
            // structAsymId: scope.bestStructAsymId
        });
    }

    drawHelicesSubpaths(startResidueNumber:number, stopResidueNumber:number, index:number, curveYdiff:number) { 
        const _this = this;
        curveYdiff = 0;
        const diffVal = 5;
        let curveYdiff2 = curveYdiff - diffVal;
        if(this.scaledPointsArr[3] > this.scaledPointsArr[9]) curveYdiff2 = curveYdiff + diffVal;
        const totalAaInPath = (stopResidueNumber - startResidueNumber) + 1
        if(curveYdiff === 0) curveYdiff2 = 0;
        let subPathHeight = ((this.scaledPointsArr[9] - curveYdiff2) - this.scaledPointsArr[3])/totalAaInPath;
        let startPoint = 0;
        if(curveYdiff === 0){
            let boxHeight = (this.svgEle.select('.helices'+index).node().getBBox().height) + (subPathHeight/2);
            const singleUnitHt = boxHeight/totalAaInPath;
            //boxHeight = boxHeight - singleUnitHt; //height correction
            subPathHeight = (boxHeight - singleUnitHt/2)/totalAaInPath;
            startPoint = (subPathHeight - singleUnitHt/10);
            if(this.scaledPointsArr[3] > this.scaledPointsArr[9]){
                //startPoint = -(boxHeight + singleUnitHt/3);
                startPoint = -(boxHeight + singleUnitHt);
            }
        }

        //create subsections/paths
        let dValArr2 = [];
        let subPathObj:any = {};
        if(curveYdiff === 0){
            for(let subPathIndex=0; subPathIndex<totalAaInPath; subPathIndex++){
                subPathObj = {type: 'helices'};
                if(subPathIndex === 0){
                    if(this.scaledPointsArr[3] > this.scaledPointsArr[9]){
                        subPathObj['residue_number'] = stopResidueNumber;
                    }else{
                        subPathObj['residue_number'] = startResidueNumber;
                    }
                    subPathObj['pathData'] = [
                        this.scaledPointsArr[0], this.scaledPointsArr[3] + startPoint,
                        this.scaledPointsArr[4], this.scaledPointsArr[3] + startPoint,
                        this.scaledPointsArr[4], this.scaledPointsArr[3] + startPoint + subPathHeight,
                        this.scaledPointsArr[0], this.scaledPointsArr[3] + startPoint + subPathHeight
                    ];
                }else{
                    if(this.scaledPointsArr[3] > this.scaledPointsArr[9]){
                        subPathObj['residue_number'] = stopResidueNumber - subPathIndex;
                    }else{
                        subPathObj['residue_number'] = startResidueNumber + subPathIndex;
                    }
                    subPathObj['pathData'] = [
                        dValArr2[subPathIndex - 1]['pathData'][6], dValArr2[subPathIndex - 1]['pathData'][7],
                        dValArr2[subPathIndex - 1]['pathData'][4], dValArr2[subPathIndex - 1]['pathData'][5],
                        dValArr2[subPathIndex - 1]['pathData'][4], dValArr2[subPathIndex - 1]['pathData'][5] + subPathHeight,
                        dValArr2[subPathIndex - 1]['pathData'][6], dValArr2[subPathIndex - 1]['pathData'][5] + subPathHeight
                    ];
                }
                dValArr2.push(subPathObj);
            }
        }else{
            for(let subPathIndex=0; subPathIndex<totalAaInPath; subPathIndex++){
                subPathObj = {type: 'helices', elementIndex: index};
                if(subPathIndex === 0){
                    subPathObj['residue_number'] = startResidueNumber;
                    subPathObj['pathData'] = [
                        _this.scaledPointsArr[0], _this.scaledPointsArr[3] + curveYdiff2/2,
                        _this.scaledPointsArr[4], _this.scaledPointsArr[3] + curveYdiff2/2,
                        _this.scaledPointsArr[4], _this.scaledPointsArr[3] + subPathHeight + curveYdiff2/2,
                        _this.scaledPointsArr[0], _this.scaledPointsArr[3] + subPathHeight + curveYdiff2/2
                    ];
                }else{
                    subPathObj['residue_number'] = startResidueNumber + subPathIndex;
                    subPathObj['pathData'] = [
                        dValArr2[subPathIndex - 1]['pathData'][6], dValArr2[subPathIndex - 1]['pathData'][7],
                        dValArr2[subPathIndex - 1]['pathData'][4], dValArr2[subPathIndex - 1]['pathData'][5],
                        dValArr2[subPathIndex - 1]['pathData'][4], dValArr2[subPathIndex - 1]['pathData'][5] + subPathHeight,
                        dValArr2[subPathIndex - 1]['pathData'][6], dValArr2[subPathIndex - 1]['pathData'][5] + subPathHeight
                    ];
                }
                dValArr2.push(subPathObj);
            }
        }
        
        this.svgEle.selectAll('.subpath-helices'+index).remove();
        
        this.svgEle.selectAll('.subpath-helices'+index)
        .data(dValArr2)
        .enter()
        .append('path')
        .attr('class', function(d:any){ return 'helicesSubPath subpath-helices'+index+' topo_res_'+d.residue_number })
        .attr('d', function(d:any){ return 'M'+d.pathData.join(' ')+' Z' })
        .attr('stroke', '#111')
        .attr('stroke-width', '0')
        .attr('fill', 'white')
        .attr('fill-opacity','0')
        .on('mouseover', function(d:any){ _this.mouseoverAction(this, d); })
        .on('mousemove', function(d:any){ _this.mouseoverAction(this, d); })
        .on('mouseout', function(d:any){ _this.mouseoutAction(this, d); })
        .on("click", function(d:any) { _this.clickAction(d); })
    }

    drawHelicesMaskShape(index:number) {
        const adjustmentFactor = 0.3;
        let helicesMaskArr = [
            [this.scaledPointsArr[0]-adjustmentFactor, this.scaledPointsArr[1], 
            this.scaledPointsArr[2], this.scaledPointsArr[3]-adjustmentFactor,
            this.scaledPointsArr[4]+adjustmentFactor, this.scaledPointsArr[5],
            this.scaledPointsArr[4]+adjustmentFactor, this.scaledPointsArr[3],
            this.scaledPointsArr[0]-adjustmentFactor, this.scaledPointsArr[3]
            ],
            [this.scaledPointsArr[6]+adjustmentFactor, this.scaledPointsArr[7], 
            this.scaledPointsArr[8], this.scaledPointsArr[9]+adjustmentFactor,
            this.scaledPointsArr[10]-adjustmentFactor, this.scaledPointsArr[11],
            this.scaledPointsArr[10]-adjustmentFactor, this.scaledPointsArr[9],
            this.scaledPointsArr[6]+adjustmentFactor, this.scaledPointsArr[9]
            ]
        ]
        if(this.scaledPointsArr[3] > this.scaledPointsArr[9]){
            helicesMaskArr = [
                [this.scaledPointsArr[0]-adjustmentFactor, this.scaledPointsArr[1], 
                this.scaledPointsArr[2], this.scaledPointsArr[3]+2,
                this.scaledPointsArr[4]+adjustmentFactor, this.scaledPointsArr[5],
                this.scaledPointsArr[4]+adjustmentFactor, this.scaledPointsArr[3],
                this.scaledPointsArr[0]-adjustmentFactor, this.scaledPointsArr[3]
                ],
                [this.scaledPointsArr[6]+adjustmentFactor, this.scaledPointsArr[7], 
                this.scaledPointsArr[8], this.scaledPointsArr[9]-adjustmentFactor,
                this.scaledPointsArr[10]-adjustmentFactor, this.scaledPointsArr[11],
                this.scaledPointsArr[10]-adjustmentFactor, this.scaledPointsArr[9],
                this.scaledPointsArr[6]+adjustmentFactor, this.scaledPointsArr[9]
                ]
            ]
        }
        
        //remove old maskpath
        this.svgEle.selectAll('.maskpath-helices'+index).remove();
        
        //create new resized mask path 
        this.svgEle.selectAll('.maskpath-helices'+index)
        .data(helicesMaskArr)
        .enter()
        .append('path')
        .attr('class', function(d:any){ return 'helicesMaskPath maskpath-helices'+index })
        .attr('d', function(d:any){ 
            return 'M'+d[0]+' '+d[1]+' Q'+d[2]+' '+d[3]+' '+d[4]+' '+d[5]+' L'+d[6]+' '+d[7]+' '+d[8]+' '+d[9]+' Z';
        })
        .attr('stroke', '#111')
        .attr('stroke-width', 0.3)
        .attr('fill', 'white')
        .attr('stroke-opacity', 0)
        
    }

    drawCoilsSubpaths (startResidueNumber:number, stopResidueNumber:number, index:number) {
        const _this = this;
        
        const coilEle = this.svgEle.select('.coils'+index);
        const totalAaInPath = (stopResidueNumber - startResidueNumber) + 1
        const coilLength = coilEle.node().getTotalLength();
        const subPathLength = coilLength/totalAaInPath;
        
        let subPathCordsArr = [];
        let prevPathCord = undefined;
        let prevCordArrPositon = undefined;
        //var prevSubPathCord = [];
        let newSubPathCords:any = {};
        
        if(totalAaInPath === 1){
            newSubPathCords = {
                residue_number: startResidueNumber,
                type: 'coils',
                pathData: _this.scaledPointsArr,
                elementIndex: index
            }
            subPathCordsArr.push(newSubPathCords)
        
        }else{
            for(let subPathIndex=0; subPathIndex<totalAaInPath; subPathIndex++){
                
                const segLength = subPathLength * (subPathIndex + 1);
                const subPathCord = coilEle.node().getPointAtLength(segLength);
                const cordArrPositon = coilEle.node().getPathSegAtLength(segLength);
                newSubPathCords = {
                    residue_number: startResidueNumber + subPathIndex,
                    type: 'coils',
                    elementIndex: index
                }
            
                if(cordArrPositon === 1){
                    newSubPathCords['pathData'] = _this.scaledPointsArr.slice(0, 2);
                }else{
                
                    if(prevCordArrPositon === undefined){
                        newSubPathCords['pathData'] = _this.scaledPointsArr.slice(0, cordArrPositon * 2);
                    }else{
                        newSubPathCords['pathData'] = _this.scaledPointsArr.slice(prevCordArrPositon * 2, cordArrPositon * 2);
                        newSubPathCords['pathData'].unshift(prevPathCord.x, prevPathCord.y);
                    
                    }
                    
                    prevPathCord = subPathCord;
                    prevCordArrPositon = cordArrPositon;
                
                }
                
                newSubPathCords['pathData'] = newSubPathCords['pathData'].concat([subPathCord.x, subPathCord.y]);
                subPathCordsArr.push(newSubPathCords);
            }
        
        }
        
        if(startResidueNumber !== -1 && stopResidueNumber !== -1){
            this.svgEle.selectAll('.subpath-coils'+index).remove();
            
            this.svgEle.selectAll('.subpath-coils'+index)
                .data(subPathCordsArr)
                .enter()
                .append('path')  
                .attr('class', function(d:any){ return 'coilsSubPath subpath-coils'+index+' topo_res_'+d.residue_number })
                .attr('d', function(d:any){ return 'M '+d.pathData.join(' ') })
                .attr('stroke', this.defaultColours.borderColor)
                .attr('stroke-width', 0.3)
                .attr('fill', 'none')
                .attr('stroke-opacity','1')
                .on('mouseover', function(d:any){ _this.mouseoverAction(this, d); })
                .on('mousemove', function(d:any){ _this.mouseoverAction(this, d); })
                .on('mouseout', function(d:any){ _this.mouseoutAction(this, d); })
                .on("click", function(d:any) { _this.clickAction(d); })
            
            //Hide the main coil path
            this.svgEle.selectAll('.coils'+index).attr('stroke-opacity',0);
        }
        
        const termsData = this.apiData[2][this.entryId][this.entityId][this.chainId].terms;
        const totalCoilsInStr = this.apiData[2][this.entryId][this.entityId][this.chainId].coils.length;
        if(index === 0){
            this.svgEle.selectAll('.terminal_N').remove();
            this.svgEle.selectAll('.terminal_N')
                    .data([termsData[0]])
                    .enter()
                    .append('text')
                    .attr('class', 'terminals terminal_N')
                    .attr('text-anchor','middle')
                    .text('N')
                    .attr('x', subPathCordsArr[0]['pathData'][0])
                    .attr('y', subPathCordsArr[0]['pathData'][1])
                    .attr('stroke','#0000ff')
                    .attr('stroke-width','0.3')
                    // .attr('font-size', 3 * this.zoom.scale() +'px')
                    .attr('font-size', '3px')
                    .attr('style',"-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;")
        }else if(index === totalCoilsInStr - 1){
            const pathDataLen = subPathCordsArr[totalAaInPath - 1]['pathData'].length;
            let adjustmentFactor = -2;
            if(subPathCordsArr[totalAaInPath - 1]['pathData'][pathDataLen - 1] > subPathCordsArr[totalAaInPath - 1]['pathData'][pathDataLen - 3]){
                adjustmentFactor = 2;
            }
            this.svgEle.selectAll('.terminal_C').remove();
            this.svgEle.selectAll('.terminal_C')
                    .data([termsData[1]])
                    .enter()
                    .append('text')
                    .attr('class', 'terminals terminal_N')
                    .attr('text-anchor','middle')
                    .text('C')
                    .attr('x', subPathCordsArr[totalAaInPath - 1]['pathData'][pathDataLen - 2])
                    .attr('y', subPathCordsArr[totalAaInPath - 1]['pathData'][pathDataLen - 1] + adjustmentFactor)
                    .attr('stroke','#ff0000')
                    .attr('stroke-width','0.3')
                    // .attr('font-size', 3 * this.zoom.scale() +'px')
                    .attr('font-size', '3px')
                    .attr('style',"-webkit-tap-highlight-color: rgba(0, 0, 0, 0); text-anchor: middle; font-style: normal; font-variant: normal; font-weight: normal; font-stretch: normal; line-height: normal; font-family: Arial;")
        }
        
    }
    getAdjustedStartAndStop(secStrType: any, secStrData: any) {
        if(secStrData == undefined) {
            return null;
        }
        if(secStrType != 'helices' && secStrType != 'coils' && secStrType != 'strands') {
            return [secStrData.start, secStrData.stop];
        } else if(secStrData.start === -1) {
            return [secStrData.start, secStrData.stop];
        }
        if((secStrData.stop  < Number(filterRange.split(",")[0]) || (secStrData.start > Number(filterRange.split(",")[1])))) {
            return null;
        } else if (secStrType === 'helices' || secStrType === 'strands') {
            return [secStrData.start, secStrData.stop];
        } else if (secStrData.stop > Number(filterRange.split(",")[1])) {
            return [secStrData.start, Number(filterRange.split(",")[1])];
        } else if (secStrData.start < Number(filterRange.split(",")[0])) {
            return [Number(filterRange.split(",")[0]), secStrData.stop];
        }  
        return [secStrData.start, secStrData.stop];
    }
    drawTopologyStructures() {
        //Add container elements
        this.targetEle.innerHTML = `<div style="${this.displayStyle}">
            <div class="svgSection" style="position:relative;width:100%;"></div>
            <div style="${this.menuStyle}">
                <a style="color: black;border-bottom:none; cursor:pointer;margin-left: 16px;" target="_blank" href="https://pdbe.org/${this.entryId}">${this.entryId}</a> | <span class="menuDesc">Entity ${this.entityId} | Chain ${this.chainId.toUpperCase()}</span>
                <div class="menuOptions" style="float:right;margin-right: 20px;">
                    <select class="menuSelectbox" style="margin-right: 10px;"><option value="">Select</option></select>
                    <img class="saveSVG" src="static/alignments/png/Save.png" style="height:15px; width: 15px; border:0;position: relative;margin-right: 15px;cursor:pointer;" title="saveSVG" />\n
                    <img class="resetIcon" src="static/alignments/png/refresh.png" style="height:15px; width: 15px; border:0;position: absolute;margin-top: 11px;cursor:pointer;" title="Reset view" />
                </div>
            </div>
        </div>`;

        //Get dimensions
        let targetEleWt = this.targetEle.offsetWidth;
        let targetEleHt = this.targetEle.offsetHeight;
        if(targetEleWt == 0) targetEleWt = (this.targetEle.parentElement as HTMLElement).offsetWidth;
        if(targetEleHt == 0) targetEleHt = (this.targetEle.parentElement as HTMLElement).offsetHeight;

        if(targetEleWt <= 330) (this.targetEle.querySelector('.menuDesc') as HTMLElement).innerText = `${this.entityId} | ${this.chainId.toUpperCase()}`;
        
        //Set svg section dimensions
        const svgSection:any = this.targetEle.querySelector('.svgSection');
        const svgSectionHt = targetEleHt - 40;
        const svgSectionWt = targetEleWt;
        svgSection.style.height = svgSectionHt+'px';

        //Set svg dimensions
        const svgHt = svgSectionHt - 20;
        const svgWt = svgSectionWt - 5;
        svgSection.innerHTML = `<svg class="topoSvg" preserveAspectRatio="xMidYMid meet" viewBox="0 0 100 100" style="width:${svgWt}px;height:${svgHt}px;margin:10px 0;"></svg>`;

        this.svgEle = d3.select(this.targetEle).select('.topoSvg');
       
        this.getDomainRange();
        this.scaledPointsArr = [];
        this.svgEle.call(this.zoom).on("contextmenu", function (d:any, i:number) { d3.event.preventDefault(); }); //add zoom event and block right click event
        const topologyData = this.apiData[2][this.entryId][this.entityId][this.chainId];
        for(let secStrType in topologyData){
        // angular.forEach(this.apiResult.data[_this.entryId].topology[scope.entityId][scope.bestChainId], function(secStrArr, secStrType) {
            const secStrArr =  topologyData[secStrType];
            if(!secStrArr) return;
            //iterating on secondary str data array
            
            secStrArr.forEach((secStrData:any, secStrDataIndex: number) => {
                if(typeof secStrData.path !== 'undefined' && secStrData.path.length > 0 && this.getAdjustedStartAndStop(secStrType, secStrData) != null){  

                        var istart = this.getAdjustedStartAndStop(secStrType, secStrData)[0];
                        var istop = this.getAdjustedStartAndStop(secStrType, secStrData)[1];
                        if(istart !== -1 || (secStrType === "coils" && istart=== -1 && istop === -1 && this.getAdjustedStartAndStop(secStrType, secStrArr[secStrDataIndex + 1]) != null && this.getAdjustedStartAndStop(secStrType, secStrArr[secStrDataIndex - 1]) != null)){
                        
                        if(secStrType === 'terms'){
                            //Terms
                        }else{
                            let curveYdiff = 0
                            //modify helices path data to create a capsule like structure
                            if(secStrType === 'helices'){
                                const curveCenter = secStrData.path[0] + ((secStrData.path[2] - secStrData.path[0])/2);
                                                                    
                                curveYdiff = 2 * (secStrData.minoraxis * 1.3);
                                if(secStrData.path[1] >  secStrData.path[3]){
                                    curveYdiff = -2 * (secStrData.minoraxis * 1.3);
                                }
                                
                                const newPathCords = [
                                    secStrData.path[0], secStrData.path[1],
                                    curveCenter, secStrData.path[1] - curveYdiff,
                                    secStrData.path[2], secStrData.path[1],
                                    secStrData.path[2], secStrData.path[3],
                                    curveCenter, secStrData.path[3] + curveYdiff,
                                    secStrData.path[0], secStrData.path[3]
                                ];
                                
                                secStrData.path = newPathCords;
                            }
                            
                            secStrData.secStrType = secStrType;
                            secStrData.pathIndex = secStrDataIndex;
                            const newEle = this.svgEle.selectAll('path.'+secStrType+''+secStrDataIndex)
                            .data([secStrData])
                            .enter()
                            .append('path')  
                            .attr('class', () => {
                                if( secStrData.start=== -1 && secStrData.stop === -1 && secStrType !== 'terms'){
                                    return 'dashedEle topologyEle '+secStrType+' '+secStrType+''+secStrDataIndex+' topoEleRange_'+secStrData.start+'-'+secStrData.stop;
                                } else{
                                    return 'topologyEle '+secStrType+' '+secStrType+''+secStrDataIndex+' topoEleRange_'+istart+'-'+istop;
                                } 
                            })
                            .attr('d', (d: any) => {
                                let dVal = 'M';
                                const pathLenth = secStrData.path.length;
                                let xScaleFlag = true;
                                //if(secStrData.path[1] > secStrData.path[7]) maskDiff = 1;
                                for(let i=0; i<pathLenth; i++){
                                    if(secStrType === 'helices' && (i === 2 || i === 8)) dVal += ' Q'
                                    //if(secStrType === 'coils' && secStrData.path.length < 12 && i === 2) dVal += ' C'
                                    //if(secStrType === 'coils' && secStrData.path.length < 14 && secStrData.path.length > 12 && i === 4) dVal += ' C'
                                    if((secStrType === 'helices' && i === 6) || (secStrType === 'coils' && secStrData.path.length < 12 && i === 8)) dVal += ' L'
                                    if(xScaleFlag){
                                        const xScaleValue = this.xScale(secStrData.path[i]);
                                        dVal += ' '+xScaleValue;
                                        this.scaledPointsArr.push(xScaleValue);
                                    }else{
                                        const yScaleValue = this.yScale(secStrData.path[i]);
                                        dVal += ' '+yScaleValue;
                                        this.scaledPointsArr.push(yScaleValue);
                                    }
                                    
                                    xScaleFlag = !xScaleFlag;
                                }
                                if(secStrType === 'strands' || secStrType === 'helices') dVal += ' Z'
                                return dVal;
                            })
                            .attr('fill', 'none')
                            .attr('stroke-width', 0.3)
                            .attr('stroke', this.defaultColours.borderColor)
                        
                            if(secStrData.start === -1 && secStrData.stop === -1){
                                newEle.attr('stroke-dasharray', '0.9')
                            }
                        
                            //hightlight node calculations
                            if(secStrType === 'strands'){
                                //create subsections/paths
                                this.drawStrandSubpaths(istart, istop, secStrDataIndex)

                                //Create mask to restore shape
                                this.drawStrandMaskShape(secStrDataIndex);
                                
                                //bring original/complete helices in front newEle
                                // this.svgEle.append(newEle.node());		
                                this.svgEle._groups[0][0].append(newEle.node());						
                            }
                            
                            //for helices
                            if(secStrType === 'helices'){
                                //create subsections/paths
                                this.drawHelicesSubpaths(istart, istop, secStrDataIndex, curveYdiff)

                                //Create mask to restore shape
                                this.drawHelicesMaskShape(secStrDataIndex);
                                
                                // //bring original/complete helices in front
                                // angular.element(element[0].querySelector('.topoSvg')).append(newEle.node());
                                this.svgEle._groups[0][0].append(newEle.node());
                            }
                        
                            //for coils
                            if(secStrType === 'coils'){
                                //create subsections/paths
                                this.drawCoilsSubpaths(istart, istop, secStrDataIndex);                       
                            }
                        
                            this.scaledPointsArr = []; //empty the arr for next iteration
                        }
                    }
                }
                
            });
            
        };
        
        //bring rsrz validation circles in front
        this.svgEle._groups[0][0].append(this.svgEle.selectAll('.validationResidue').node());
    };

    zoomDraw() {

        const new_xScale = d3.event.transform.rescaleX(this.xScale);
        const new_yScale = d3.event.transform.rescaleY(this.yScale);

        // return
        const _this = this;
        
        _this.scaledPointsArr = [];
        
        const pathEle = this.svgEle.selectAll('.topologyEle');
        let pathIndex = 0;
        let pathStartResidue = 0;
        let pathStopResidue = 0;
        const curveYdiff = 0;
        pathEle.each(function(d:any){
            d3.select(d3.select(this).node()).attr('d', function(d:any){
                pathIndex = d.pathIndex;
                pathStartResidue = d.start;
                pathStopResidue = d.stop;
                
                let dVal = 'M';
                const pathLenth = d.path.length;
                let xScaleFlag = true;
                // var maskDiff = -1; //value to add/minus to show the border properly
                for(let i=0; i<pathLenth; i++){
                    if(d.secStrType === 'helices' && (i === 2 || i === 8)) dVal += ' Q'
                    //if(d.secStrType === 'coils' && d.path.length < 12 && i === 2) dVal += ' C'
                    if((d.secStrType === 'helices' && i === 6) || (d.secStrType === 'coils' && d.path.length < 12 && i === 8)) dVal += ' L'
                    if(xScaleFlag){
                        const xScaleValue = new_xScale(d.path[i])
                        dVal += ' '+xScaleValue;
                        _this.scaledPointsArr.push(xScaleValue);
                    }else{
                        const yScaleValue = new_yScale(d.path[i])
                        dVal += ' '+yScaleValue;
                        _this.scaledPointsArr.push(yScaleValue);
                    }
                    
                    xScaleFlag = !xScaleFlag;
                }
                if(d.secStrType === 'strands' || d.secStrType === 'helices') dVal += ' Z'
                    return dVal;
            });
            
            //Create mask to restore shape
            if(d.secStrType === 'helices'){
                //create subsections/paths
                _this.drawHelicesSubpaths(pathStartResidue, pathStopResidue, pathIndex, curveYdiff)
                _this.drawHelicesMaskShape(pathIndex);
                
                //bring original/complete helices in front newEle
                _this.svgEle._groups[0][0].append(d3.select(this).node());
            }else if(d.secStrType === 'strands'){
                _this.drawStrandSubpaths(pathStartResidue, pathStopResidue, pathIndex)
                _this.drawStrandMaskShape(pathIndex);
                
                //bring original/complete helices in front newEle
                _this.svgEle._groups[0][0].append(d3.select(this).node());
            }//for coils
            else if(d.secStrType === 'coils'){
                //create subsections/paths
                _this.drawCoilsSubpaths(pathStartResidue, pathStopResidue, pathIndex);
            }
            
            _this.scaledPointsArr = []; //empty the arr for next iteration
        });
        
        //scale validation - rsrz circle's
        let ValResheight = 0;
        this.svgEle
            .selectAll('.validationResidue')
                .attr('transform', function(d:any){
                    //get Shape dimesions
                    const residueEle = _this.svgEle.select('.topo_res_'+d.residue_number);
                    const dimensions = residueEle.node().getBBox();
                    const residueEleData = residueEle.data();
                    let reszEleCordinates = {x:0, y:0};
                    if(residueEleData[0].type ==='strands' || residueEleData[0].type ==='helices'){
                        reszEleCordinates = {
                            x : dimensions.x + dimensions.width/2, 
                            y : dimensions.y + dimensions.height/2
                        };
                    }else{
                        const coilCenter = residueEle.node().getPointAtLength(residueEle.node().getTotalLength()/2);
                        reszEleCordinates = {
                            x : coilCenter.x, 
                            y : coilCenter.y
                        };
                    }
                    ValResheight = dimensions.height/2;
                    return "translate(" + reszEleCordinates.x + "," + reszEleCordinates.y + ")";
                    
                })
                .attr("d", d3.symbol().type(function(d:any,i:number) { return d3.symbols[0];}).size(ValResheight));
        
        //scale selection paths
        this.svgEle
            .selectAll('.residueSelection')
            .attr('d', function(d:any){
                //assign the d attribute of the corresponding sub-path
                const dataObj = d3.select(this).data() as any;
                return _this.svgEle.select('.topo_res_'+dataObj[0].residueNumber).attr('d');
            });
                        
        //shift coilssub path to top in DOM
        this.svgEle._groups[0][0].querySelectorAll('.coilsSubPath').forEach((node:any) => this.svgEle._groups[0][0].append(node));
        
        //shift dashed paths to top in DOM
        this.svgEle._groups[0][0].querySelectorAll('.dashedEle').forEach((node:any) => this.svgEle._groups[0][0].append(node));
        
        this.displayDomain('zoom');
        
        //bring rsrz validation circles in front
        this.svgEle._groups[0][0].querySelectorAll('.validationResidue').forEach((node:any) => this.svgEle._groups[0][0].append(node));
        
        //bring selection in front
        this.svgEle._groups[0][0].querySelectorAll('.residueSelection').forEach((node:any) => this.svgEle._groups[0][0].append(node));
        
    }

    clearHighlight(){
        this.svgEle.selectAll('.residueHighlight').remove();
    }

    highlight(startResidue:number, endResidue:number, color?:{r:number, g:number, b:number} | string, eventType?:string) {
        const _this = this;
        
        let fill:any = '#000000';
        let stroke:any = '#000000';
        let strokeWidth = 0.3;
        let strokeOpacity = 0;
        
        for(let residueNumber = startResidue; residueNumber <= endResidue; residueNumber++){
            //get topology residue details
            const residueEle = this.svgEle.select('.topo_res_'+residueNumber);
            if(residueEle && residueEle._groups && residueEle._groups[0][0] == null)return; //if residue element do not exist
            const residueEleNode = residueEle.node();
            const residueEleData = residueEle.data();

            if(color){
                if(typeof color == 'string'){
                    stroke = color;
                    fill = color;
                }else{
                    stroke = d3.rgb(color.r,color.g,color.b);
                    fill = d3.rgb(color.r,color.g,color.b);
                }
            }

            if(residueEleData[0].type !=='strands' && residueEleData[0].type !=='helices'){
                fill = 'none';
                strokeWidth = 2;
                strokeOpacity = 0.5;
            }else{
                stroke = 'none';
            }
            
            this.svgEle
                .append('path')
                    .data([{residueNumber: residueNumber}])
                    .attr('class', (d:any) => {
                        if(eventType == 'click'){
                            return 'residueSelection seletectedResidue_'+residueNumber;
                        }else{
                            return 'residueHighlight highlightResidue_'+residueNumber;
                        }
                    })
                    .attr('d', residueEle.attr('d'))
                    .attr('fill', fill)
                    .attr('fill-opacity', 0.5)
                    .attr('stroke', stroke)
                    .attr('stroke-opacity', strokeOpacity)
                    .attr('stroke-width', strokeWidth)
                    .on('mouseover', (d: any) => { _this.mouseoverAction(residueEleNode, residueEleData[0]); })
                    .on('mousemove', (d: any) => { _this.mouseoverAction(residueEleNode, residueEleData[0]); })
                    .on('mouseout', (d: any) => { _this.mouseoutAction(residueEleNode, residueEleData[0]); })
                    .on("click", (d: any) => { _this.clickAction(residueEleData[0]); })
        }
        
    }

    drawValidationShape(residueNumber:number, shape:string, rgbColor:string) {
        const _this = this;
        //calculate Shape dimesions
        const residueEle = _this.svgEle.select('.topo_res_'+residueNumber) as any;
        if(residueEle._groups[0][0] == null)return; //if residue element do not exist
        const dimensions = residueEle.node().getBBox();
        const residueEleData = residueEle.data();
        let reszEleCordinates = {x:0, y:0};
        if(residueEleData[0].type ==='strands' || residueEleData[0].type ==='helices'){
            reszEleCordinates = {
                x : dimensions.x + dimensions.width/2, 
                y : dimensions.y + dimensions.height/2
            };
        }else{
            const coilCenter = residueEle.node().getPointAtLength(residueEle.node().getTotalLength()/2);
            reszEleCordinates = {
                x : coilCenter.x, 
                y : coilCenter.y
            };
        }	
        const validationResData = {
            residue_number: residueNumber,
            tooltipMsg : 'Validation issue: RSRZ <br>',
            tooltipPosition: 'prefix',	
        };
        
        this.svgEle
            .append('path')
                .attr('class', 'validationResidue rsrz_'+residueNumber)
                .data([validationResData])
                .attr('fill',rgbColor)
                .attr('stroke', '#000')
                .attr('stroke-width', 0.3)
                .attr("transform", function(d:any) { return "translate(" + reszEleCordinates.x + "," + reszEleCordinates.y + ")"; })
                .attr("d", d3.symbol().type(function(d:any,i:number) { return d3.symbols[0];}).size(dimensions.height/2))
                .style('display', 'none')
                .on('mouseover', function(d:any){ _this.mouseoverAction(this, d); })
                .on('mousemove', function(d:any){ _this.mouseoverAction(this, d); })
                .on('mouseout', function(d:any){ _this.mouseoutAction(this, d); })
                .on("click", function(d:any) { _this.clickAction(d); })
        
    }

    getAnnotationFromMappings = function () {
        const mappings = this.apiData[1];
        if(typeof mappings == 'undefined') return;
        const mappingsData = this.apiData[1][this.entryId];
        const categoryArr = ['UniProt','CATH','Pfam','SCOP','RiboVision'];
        for(let catIndex=0; catIndex < 3; catIndex++){
            if(typeof mappingsData[categoryArr[catIndex]] !== 'undefined'){
                
                if(Object.entries(mappingsData[categoryArr[catIndex]]).length !== 0){
                    let residueDetails:any = [];
                    //Iterate over mappings data to get start and end residues
                    const mappingRecords = mappingsData[categoryArr[catIndex]];
                    for(let accKey in mappingRecords){

                        mappingRecords[accKey].mappings.forEach((domainMappings:any) => {
                            if(domainMappings.entity_id == this.entityId && domainMappings.chain_id == this.chainId){
                                
                                residueDetails.push({
                                    start: domainMappings.start.residue_number,
                                    end: domainMappings.end.residue_number,
                                    color: undefined
                                });

                            }
                        });
                    }
                    
                    if(residueDetails.length > 0){
                        this.domainTypes.push(
                            {
                                label: categoryArr[catIndex],
                                data: residueDetails
                            }
                        )
                    }
                                        
                }
            }
        }
        
    }

    getChainStartAndEnd() {
        //chains array from polymerCoveragePerChain api result
        if(typeof this.apiData[4] == 'undefined') return;
        const chainsData = this.apiData[4][this.entryId].molecules[0].chains;
        
        //Iterate molecule data to get chain start and end residue
        let chainRange = {start:0, end:0}
        const totalChainsInArr = chainsData.length;
        for(let chainIndex=0; chainIndex < totalChainsInArr; chainIndex++){
            if(chainsData[chainIndex].chain_id == this.chainId){
                
                //iterate over observed array
                chainsData[chainIndex].observed.forEach((observedData:any, observedDataIndex:number) => {
                    
                    if(observedDataIndex == 0){
                        chainRange.start = observedData.start.residue_number;
                        chainRange.end = observedData.end.residue_number;
                    }else{
                        if(observedData.start.residue_number < chainRange.start){
                            chainRange.start = observedData.start.residue_number;
                        }
                        if(observedData.end.residue_number > chainRange.end){
                            chainRange.end = observedData.end.residue_number;
                        }
                    }
                    
                })
                
                break;
            }
        }
        
        return chainRange;
        
    }

    create2D3DAnnotations(name: string, residueDetails: any, 
                        TWCrgbMap: Map<number, any>, TWCData: Map<number, string>,
                        chain_start: number, chain_end: number) {
        const _this = this;
        TWCData.forEach(function(value, index) {
            if (chain_start <= index && index <= chain_end){
                let rgb_color = TWCrgbMap.get(index);
                selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.entityId,
                    start_residue_number: index, 
                    end_residue_number: index,
                    color: rgb_color[1],
                    sideChain: false,
                });
                _this.defaultColours.qualityRiboVision= "rgb("+String(rgb_color[0].join(','))+")";
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
        if (TWCData.size < mapped_aa_properties.get("Charge").length) {
            for(var i = TWCData.size - 1; i < mapped_aa_properties.get("Charge").length; i++) {
                selectSections_RV1.get(name).push({ //3d
                    entity_id: _this.entityId,
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
        const chainRange:any = this.getChainStartAndEnd();
        
        if (mapped_aa_properties) {
            mapped_aa_properties.forEach(function(value, index) {    
                let residueDetails:any = [{
                    //start: chainRange.start,
                    //end: chainRange.end,
                    //color: _this.defaultColours.qualityBlank,
                    //tooltipMsg: 'No data for '
                }];
                let name = index;
                let separatedData = value;

                selectSections_RV1.set(name, [])

                let min = Math.min(...aaPropertyConstants.get(name));
                let max = Math.max(...aaPropertyConstants.get(name));
                let colormapArray = aaColorData.get(name); 
                const [TWCrgbMap, TWCData] = parsePVData(separatedData, min, max, colormapArray);

                selectSections_RV1.get(name).push({entity_id: _this.entityId, focus: true});
                
                if (void 0 !== TWCData){
                    residueDetails = _this.create2D3DAnnotations(name, residueDetails, 
                                                                TWCrgbMap, TWCData, 
                                                                chainRange.start, chainRange.end);
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
        }
        else {
            //catch block
        };
      }
      

    getAnnotationFromOutliers() {
        const _this = this;
        const chainRange:any = this.getChainStartAndEnd();
        let residueDetails:any = [{
            start: chainRange.start,
            end: chainRange.end,
            color: _this.defaultColours.qualityGreen,
            tooltipMsg: 'No validation issue reported for '
        }];
        
        //Two temporary arrays for grouping rsrz and other outliers tooltip message  
        let rsrzTempArray:any[] = [];
        let otherOutliersTempArray = [0];
        
        //Iterate Outlier data
        if(typeof this.apiData[3] == 'undefined') return;
        const outlierData = this.apiData[3][this.entryId];
        if(typeof outlierData !== 'undefined' && typeof outlierData.molecules !== 'undefined' && outlierData.molecules.length > 0){
            outlierData.molecules.forEach((qualityData:any) => {
                if(qualityData.entity_id == this.entityId){
                    
                    //Iterate chains array in outliers
                    qualityData.chains.forEach((chainDataObj:any) => {
                        if(chainDataObj.chain_id == this.chainId){
                            
                            //Iterate models array in chains array in outliers
                            chainDataObj.models.forEach((chainModelObj:any) => {
                            
                                //Iterate residues array in models array in outliers
                                chainModelObj.residues.forEach((outlierResidue:any) => {
                                    
                                    let resColor = _this.defaultColours.qualityYellow;
                                    let issueSpell = 'issue';
                                    if(outlierResidue.outlier_types.length === 1 && outlierResidue.outlier_types[0] === 'RSRZ'){
                                        resColor = _this.defaultColours.qualityRed;
                                        _this.drawValidationShape(outlierResidue.residue_number, 'circle', resColor);
                                        
                                        //add residue number in temporary rsrz array
                                        rsrzTempArray.push(outlierResidue.residue_number)
                                        
                                        //check if residue exist in other outliers
                                        const otherOutlierIndex = otherOutliersTempArray.indexOf(outlierResidue.residue_number);
                                        if(otherOutlierIndex > -1){
                                            residueDetails[otherOutlierIndex]['tooltipMsg'] = residueDetails[otherOutlierIndex]['tooltipMsg'].replace('<br>', ', RSRZ<br>');
                                        }else{
                                            
                                            //Adding this to have tooltip on subpath with only rsrz validation 
                                            residueDetails.push({
                                                start: parseInt(outlierResidue.residue_number),
                                                end: parseInt(outlierResidue.residue_number),
                                                color: _this.defaultColours.qualityGreen,
                                                tooltipMsg: 'Validation issue: RSRZ <br>',
                                                tooltipPosition: 'prefix'
                                            });
                                            //add residue number in temporary other Outliers array
                                            otherOutliersTempArray.push(outlierResidue.residue_number)
                                        }
                                        
                                        return;
                                    }else if(outlierResidue.outlier_types.length === 1){
                                        resColor = _this.defaultColours.qualityYellow;
                                    }else if(outlierResidue.outlier_types.length === 2){
                                        resColor = _this.defaultColours.qualityOrange;
                                        issueSpell = 'issues';
                                    }else{
                                        resColor = _this.defaultColours.qualityRed;
                                        issueSpell = 'issues';
                                    }
                                    
                                    //add residue number in temporary other Outliers array
                                    otherOutliersTempArray.push(outlierResidue.residue_number)
                                    
                                    //check if residue exist in other outliers and set the tooltip message
                                    let tooltipMsgText = 'Validation '+issueSpell+': '+outlierResidue.outlier_types.join(', ')+'<br>';
                                    const rsrzTempArrayIndex = rsrzTempArray.indexOf(outlierResidue.residue_number);
                                    if(rsrzTempArrayIndex > -1){
                                        tooltipMsgText = 'Validation issues: '+outlierResidue.outlier_types.join(', ')+', RSRZ<br>';
                                    }
                                    
                                    residueDetails.push({
                                        start: parseInt(outlierResidue.residue_number),
                                        end: parseInt(outlierResidue.residue_number),
                                        color: resColor,
                                        tooltipMsg: tooltipMsgText,
                                        tooltipPosition: 'prefix'
                                    });
                                    
                                
                                });
                            
                            });
                        
                        }
                    
                    
                    });
                    
                }
            
            });
            
            if(residueDetails.length > 0){
                this.domainTypes.push(
                    {
                        label: 'Quality',
                        data: residueDetails
                    }
                )
            }
            
        }
    }
    updateProperty() {
        const selectBoxEle: any = this.targetEle.querySelector('.menuSelectbox');
        const selectedIndex = parseInt(selectBoxEle.selectedIndex);
        rv3VUEcomponent.selected_property = this.domainTypes[selectedIndex].label;
        //console.log("Selected: " + rv3VUEcomponent.selected_property);
    }
    createDomainDropdown = function () {
        
        if(typeof this.domainTypes == 'undefined'){
            this.domainTypes = [{
                label: 'Calculated data',
                data: null
            }];
            //this.getAnnotationFromMappings();
            //this.getAnnotationFromOutliers();
            this.getAnnotationFromRibovision(mapped_aa_properties);
            this.selectedDomain = this.domainTypes[0];
        }
        if(this.domainTypes.length > 1){

            let optionList = '';
            this.domainTypes.forEach((opt:any, i:number) => {
                optionList = `${optionList}<option value="${i}">${opt.label}</option>`;
            });

            const selectBoxEle = this.targetEle.querySelector('.menuSelectbox');
            selectBoxEle.innerHTML = optionList;

            selectBoxEle.addEventListener("change", this.displayDomain.bind(this));
            selectBoxEle.addEventListener("change", this.updateProperty.bind(this));
            const resetIconEle = this.targetEle.querySelector('.resetIcon');
            resetIconEle.addEventListener("click", this.resetDisplay.bind(this));
            this.targetEle.querySelector(".saveSVG").addEventListener("click", this.saveSVG.bind(this));
            rv3VUEcomponent.topology_loaded=true;
            //selectBoxEle.style.display = 'none';
        }else{
            //this.targetEle.querySelector('.menuOptions').style.display = 'none';
        }
    }
    resetTheme() {
        const _this = this;
        this.svgEle.selectAll('.coloured').each(function(d:any){
            
            const element:any = d3.select(this);
            const node = element.node();
            
            //Remover tooltip content
            element.data()[0]['tooltipMsg'] = undefined;
            element.data()[0]['tooltipPosition'] = undefined;
            
            //Set coloured flag false
            const nodeEle = d3.select(node)
                .classed('coloured', false)
            
            //Change fill and border
            const nodeClassArr = nodeEle.attr('class').split(' ');
            if(nodeClassArr.indexOf('strandsSubPath') > -1 || nodeClassArr.indexOf('helicesSubPath') > -1){
                nodeEle.attr('fill', 'white').attr('fill-opacity', 0)
            }else{
                nodeEle.attr('stroke', _this.defaultColours.borderColor).attr('stroke-width', 0.3);
            }
            
        });
        
        //hide rsrz validation circles
        this.svgEle.selectAll('.validationResidue').style('display', 'none');
    }

    changeResidueColor(residueNumber: number, rgbColor: string, tooltipContent: string, tooltipPosition: string) {

        if(typeof rgbColor === 'undefined'){
            rgbColor = this.defaultColours.domainSelection;
        }
        const residueEle = this.svgEle.select('.topo_res_'+residueNumber);
        if(residueEle._groups[0][0] == null)return; //if residue element do not exist
        residueEle.data()[0]['tooltipMsg'] = tooltipContent;
        residueEle.data()[0]['tooltipPosition'] = tooltipPosition;
        residueEle
            .attr('stroke', function(d:any){ if(d.type === 'coils'){ return rgbColor; }else{ return '#111'; }})
            .attr('stroke-width', function(d:any){
                 if(d.type === 'coils' && masked_array[residueNumber] === false){
                    return .3;
                 }else if (d.type === 'coils'){
                     return 1;
                 }else{
                     return 0;
                }
             })
            .attr('fill', function(d:any){ if(d.type === 'coils'){ return 'none'; }else{ return rgbColor; }})
            .attr('fill-opacity', function(d:any){ if(d.type === 'coils'){ return 0; }else{ return 1; }})
            .classed("coloured", true)
            .attr('data-color', rgbColor)           
    }

    updateTheme(residueDetails:any) {
        const _this = this;
        residueDetails.forEach((residueDetailsObj:any) => {
            for(let i=residueDetailsObj.start; i<=residueDetailsObj.end; i++){
                _this.changeResidueColor(i, residueDetailsObj.color, residueDetailsObj.tooltipMsg, residueDetailsObj.tooltipPosition);
            }
        });
    }
    saveSVG(){
        function getNode(n: any, v?: any) {
		  n = document.createElementNS("http://www.w3.org/2000/svg", n);
		   for (var p in v) 
		   n.setAttributeNS(null, p.replace(/[0-9]/g,'o').replace(/\$/g,'d').replace(/\[/g,'b').replace(/[A-Z]/g, function(m, p, o, s) { return "-" + m.toLowerCase(); }), v[p]);
		  return n
		}

      var svgData1=this.targetEle.querySelector(".topoSvg");
      var svgParent=this.targetEle.querySelector(".svgSection");
      let svgData_forsave = svgData1!.cloneNode(true);
      var svg = getNode("svg");
      svg.appendChild(svgData_forsave);
      if(svgData1 && svgParent) {
        svgParent.appendChild(svgData1);
      }

     function saveSvg1(svgEl: any, name: any) {
          svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
          var svgData = svgEl.outerHTML;
          var preface = '<?xml version="1.0" standalone="no"?>\r\n';
          var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
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

    displayDomain(invokedFrom?: string) {
        const selectBoxEle:any = this.targetEle.querySelector('.menuSelectbox');
        const selectedValue = parseInt(selectBoxEle.value);
        if(selectedValue) {
            const selectedDomain = this.domainTypes[selectedValue];
            const rv3AnnotationLabels = Array.from(aaPropertyConstants.keys())
            if(selectedDomain.data !== null){
                this.resetTheme();
                this.updateTheme(selectedDomain.data);
                
                //Handle custom mapping data from RV3
                if(rv3AnnotationLabels.includes(selectedDomain.label) && invokedFrom !== 'zoom'){
                    if (filterRange != "-10000,10000"){
                        let filterRangeArr = filterRange.split(",")
                        var select_sections = selectSections_RV1.get(selectedDomain.label).slice(Number(filterRangeArr[0]), Number(filterRangeArr[1])+1);
                    } else {
                        var select_sections = selectSections_RV1.get(selectedDomain.label)
                    }
                    viewerInstance.visual.select({ data: select_sections, nonSelectedColor: {r:0,g:0,b:0}})
                }
                //show rsrz validation circles if Quality
                if(selectedDomain.label === 'Quality'){
                    this.svgEle.selectAll('.validationResidue').style('display', 'block');
                }
            }else{
                if(invokedFrom !== 'zoom'){
                    this.resetTheme();
                }
            }
        }         
    }
    
    resetDisplay(){
        const selectBoxEle:any = this.targetEle.querySelector('.menuSelectbox');
        selectBoxEle.value = 0;
        this.resetTheme();
        this.displayDomain();
    }

    handleSeqViewerEvents(e:any, eType:string){
        
        if(typeof e.eventData !== 'undefined'){
            //Abort if entryid and entityid do not match
            if(e.eventData.entryId.toLowerCase() != this.entryId.toLowerCase() || e.eventData.entityId != this.entityId) return;
                    
            //Abort if chain id is different
            if(e.eventData.elementData.pathData.chain_id && e.eventData.elementData.pathData.chain_id != this.chainId) return;

            //Remove previous selection / highlight
            let selectionPathClass = 'residueSelection';
            if(eType == 'mouseover'){
                selectionPathClass = 'residueHighlight';
            }
            this.svgEle.selectAll('.'+selectionPathClass).remove();

            let startResidue:number|undefined; 
            let endResidue:number|undefined;
            if(e.eventData.residueNumber){
                startResidue = e.eventData.residueNumber;
                endResidue = e.eventData.residueNumber;
            }else if(e.eventData.elementData.pathData.start.residue_number && e.eventData.elementData.pathData.end.residue_number){
                startResidue = e.eventData.elementData.pathData.start.residue_number;
                endResidue = e.eventData.elementData.pathData.end.residue_number;
            }

            if(typeof startResidue !== 'undefined' && typeof endResidue !== 'undefined'){
                let color:{r:number, g:number, b:number} | string | undefined;
                if(e.eventData.elementData.color && e.eventData.elementData.color.length == 1){
                    color = e.eventData.elementData.color[0];
                }else{
                    color = {r: e.eventData.elementData.color[0], g: e.eventData.elementData.color[1], b: e.eventData.elementData.color[2]};
                }
                this.highlight(startResidue, endResidue, color, eType);
            }
        }
        
    }

    handleProtvistaEvents(e:any, eType:string){
        if(typeof e.detail !== 'undefined'){
            let selColor = undefined;

            //Remove previous selection / highlight
            let selectionPathClass = 'residueSelection';
            if(eType == 'mouseover'){
                selectionPathClass = 'residueHighlight';
            }
            this.svgEle.selectAll('.'+selectionPathClass).remove();

            //Abort if chain id is different
            if(typeof e.detail.feature != 'undefined'){
                if(typeof e.detail.feature.accession != 'undefined'){
                    const accessionArr = e.detail.feature.accession.split(' ');
                    if(accessionArr[0] == 'Chain' && (accessionArr[1].toLowerCase() != this.chainId.toLowerCase())) return;
                }

                if(e.detail.trackIndex > -1 && e.detail.feature.locations && e.detail.feature.locations[0].fragments[e.detail.trackIndex].color) selColor = e.detail.feature.locations[0].fragments[e.detail.trackIndex].color;
                if(typeof selColor == 'undefined' && e.detail.feature.color) selColor = e.detail.feature.color;
            }
            if(typeof selColor == 'undefined' && e.detail.color) selColor = e.detail.color;

            if(typeof selColor != 'undefined'){
                const isRgb = /rgb/g;;
                if(isRgb.test(selColor)){
                    selColor = selColor.substring(4, selColor.length - 1).split(',');
                }else{
                    selColor = [selColor];
                }
            }

            let color:{r:number, g:number, b:number} | string | undefined;
            if(selColor){
                if(selColor.length == 1){
                    color = e.eventData.elementData.color[0];
                }else{
                    color = {r: e.eventData.elementData.color[0], g: e.eventData.elementData.color[1], b: e.eventData.elementData.color[2]};
                }
            }

            //Apply new selection
            this.highlight(e.detail.start, e.detail.end, color, eType);
        }
    }

    handleMolstarEvents(e:any, eType:string){
        

        if(typeof e.eventData !== 'undefined' && Object.keys(e.eventData).length > 0 && (masked_array[e.eventData.seq_id - 1] == true || masked_array[e.eventData.seq_id - 1] == undefined)){
            //Remove previous selection / highlight
            let selectionPathClass = 'residueSelection';
            if(eType == 'mouseover'){
                selectionPathClass = 'residueHighlight';
            }
            this.svgEle.selectAll('.'+selectionPathClass).remove();

            //Abort if entryid and entityid do not match or viewer type is unipdb
            if(e.eventData.entry_id.toLowerCase() != this.entryId.toLowerCase() || e.eventData.entity_id != this.entityId) return;
            
            //Abort if chain id is different
            //if(e.eventData.label_asym_id.toLowerCase() != this.chainId.toLowerCase()) return;
            if (e.eventData.entity_id != this.entityId) return;
            //Apply new selection
            this.highlight(e.eventData.seq_id, e.eventData.seq_id, undefined, eType);
        }
        
    }

    subscribeWcEvents(){

        //sequence viewer events
        document.addEventListener('PDB.seqViewer.click', (e:any) => {
            this.handleSeqViewerEvents(e, 'click');
        });
        document.addEventListener('PDB.seqViewer.mouseover', (e:any) => {
            this.handleSeqViewerEvents(e, 'mouseover');
        });
        document.addEventListener('PDB.seqViewer.mouseout', () => {
            this.svgEle.selectAll('.residueHighlight').remove();
        });

        //litemol viewer events
        document.addEventListener('PDB.litemol.click', (e:any) => {
            this.svgEle.selectAll('.residueSelection').remove();
            //Abort if entryid and entityid do not match or viewer type is unipdb
            if(e.eventData.entryId.toLowerCase() != this.entryId.toLowerCase() || e.eventData.entityId != this.entityId) return;								
            
            //Abort if chain id is different
            if(e.eventData.chainId.toLowerCase() != this.chainId.toLowerCase()) return;
            this.highlight(e.eventData.residueNumber, e.eventData.residueNumber, undefined, 'click');
        });
        document.addEventListener('PDB.litemol.mouseover', (e:any) => {
            this.svgEle.selectAll('.residueHighlight').remove();
            //Abort if entryid and entityid do not match or viewer type is unipdb
            if(e.eventData.entryId.toLowerCase() != this.entryId.toLowerCase() || e.eventData.entityId != this.entityId) return;								
            
            //Abort if chain id is different
            if(e.eventData.chainId.toLowerCase() != this.chainId.toLowerCase()) return;
            this.highlight(e.eventData.residueNumber, e.eventData.residueNumber, undefined, 'mouseover');
        });

        //protvista viewer events
        document.addEventListener('protvista-click', (e:any) => {
            this.handleProtvistaEvents(e, 'click');
        });
        document.addEventListener('protvista-mouseover', (e:any) => {
            this.handleProtvistaEvents(e, 'mouseover');
        });
        document.addEventListener('protvista-mouseout', () => {
            this.svgEle.selectAll('.residueHighlight').remove();
        });

        //molstar viewer events
        document.addEventListener('PDB.molstar.click', (e:any) => {
            this.handleMolstarEvents(e, 'click');
        });
        document.addEventListener('PDB.molstar.mouseover', (e:any) => {
            this.handleMolstarEvents(e, 'mouseover');
        });
        document.addEventListener('PDB.molstar.mouseout', () => {
            this.svgEle.selectAll('.residueHighlight').remove();
        });
        
    }
   

}

(window as any).PdbTopologyViewerPlugin = PdbTopologyViewerPlugin;