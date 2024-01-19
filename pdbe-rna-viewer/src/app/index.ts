import { DataService, PluginOptions, ApiData, BanNameHelper} from './data';
import { UiTemplateService } from './uiTemplate';
import { UiActionsService } from './uiActions'
import { CustomEvents } from './customEvents';

class PdbRnaViewerPlugin {

    // Globals
    options: PluginOptions;
    apiData: ApiData | undefined;
    FR3DData: any;
    BanName: any;
    FR3DNestedData: any;
    targetEle: HTMLElement;
    pdbevents: any
    dataService: DataService;
    uiTemplateService: UiTemplateService;

    rvAPI = false;
    


    constructor() {
        this.dataService = new DataService();
    };
    
    async render(target: HTMLElement, options: PluginOptions) {
        if(!target || !options.pdbId || !options.chainId || !options.entityId ) {
            console.log('Invalid plugin input!');         
            return;
        }
        this.options = options;
        if (this.options.pdbId != "cust") {
            this.apiData = await this.dataService.getApiData(this.options.entityId, this.options.chainId, this.options.pdbId);
            this.FR3DData = await this.dataService.getFR3DData(this.options.pdbId, this.options.chainId);
            this.FR3DNestedData = await this.dataService.getFR3DNestedData(this.options.pdbId, this.options.chainId);
        }
        this.BanName = await BanNameHelper.getBanName(this.options.pdbId, 'H');
        this.targetEle = <HTMLElement> target;


        this.uiTemplateService = new UiTemplateService(this.targetEle, this.options, this.apiData);
        //console.log("API", this.options.pdbId, this.apiData);
        if(this.apiData) {
            // draw topology
            this.uiTemplateService.render(this.apiData, this.FR3DData, this.FR3DNestedData, this.BanName);
            CustomEvents.subscribeToComponentEvents(this);

        } else {

            if(options.rvAPI == true) this.rvAPI = true;

            if (this.rvAPI){
                const vm = (window as any).vm;
                //const dataUrls = vm.URL;
               
                //this.apiData = await (await fetch(dataUrls)).json().then((r2dtjson) => r2dtjson.RNA_2D_json) as ApiData;
                //this.FR3DData = await (await fetch(dataUrls)).json().then((r2dtjson) => r2dtjson.RNA_BP_json) as any;
                //this.FR3DNestedData = await (await fetch(dataUrls)).json().then((r2dtjson) => r2dtjson.RNA_BP_json) as any;


                
                // const {
                //     apiDataJ,
                //     FR3DDataJ, 
                //     FR3DNestedDataJ
                // } = await (await fetch(dataUrls)).json().then((r2dtjson) => ({
                //     apiDataJ : r2dtjson.RNA_2D_json as ApiData,
                //     FR3DDataJ : r2dtjson.RNA_BP_json as any, 
                //     FR3DNestedDataJ : r2dtjson.RNA_BP_json as any,
                // }));
                const r2dtjson = vm.json_structures_from_r2dt;
                this.apiData = r2dtjson.RNA_2D_json as ApiData;
                this.FR3DData = r2dtjson.RNA_BP_json as any;
                this.FR3DNestedData = r2dtjson.RNA_BP_json as any;

                // console.log('dataUrls', this.apiData);
                // console.log('dataUrls', this.FR3DData);
                // draw topology
                this.uiTemplateService.render(this.apiData, this.FR3DData, this.FR3DNestedData, this.BanName);
    
        
            };
    
        }
        // Bind to other PDB Component events
        //if(this.options.subscribeEvents){
        //    CustomEvents.subscribeToComponentEvents(this);
//}
        //console.log("render")
        //document.addEventListener("PDB.molstar.mouseover", ((e: any) => {
        //    console.log(e)
        //    console.log(this.options.chainId)
        //    if(e.eventData && e.eventData.auth_seq_id && e.eventData.auth_asym_id === this.options.chainId) {
        //        this.selectResidue(e.eventData.auth_seq_id)
                //this.clearHighlight()
        //    }
       // })),
        //document.addEventListener("PDB.molstar.mouseout", ((e: any) => {
         //   this.clearSelection(e.eventData.residueNumber)
        //}))
    }

    selectResidue(label_seq_id: number, color?: string) {
        UiActionsService.selectNucleotide(this.options.pdbId, this.options.entityId, label_seq_id, 'mouseover', false, color);
    }

    clearSelection(label_seq_id: any) {
        UiActionsService.unSelectNucleotide(this.options.pdbId, this.options.entityId, label_seq_id, false);
    }

    highlightResidue(label_seq_id: number, color?: string) {
        UiActionsService.colorNucleotide(this.options.pdbId, label_seq_id, color, 'highlight');
    }

    clearHighlight() {
        UiActionsService.clearHighlight(this.options.pdbId);
    }
}

(window as any).PdbRnaViewerPlugin = PdbRnaViewerPlugin;
