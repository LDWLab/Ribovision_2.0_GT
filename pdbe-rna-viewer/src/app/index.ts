import { DataService, PluginOptions, ApiData } from './data';
import { UiTemplateService } from './uiTemplate';
import { UiActionsService } from './uiActions'
import { CustomEvents } from './customEvents';

class PdbRnaViewerPlugin {
    /*
    Define the PDB RNA viewer plugin
    */
    // Globals
    options: PluginOptions;
    apiData: ApiData | undefined;
    FR3DData: any;
    FR3DNestedData: any;
    targetEle: HTMLElement;
    pdbevents: any
    dataService: DataService;
    uiTemplateService: UiTemplateService;

    rvAPI = false;



    constructor() {
        //constructor
        this.dataService = new DataService();
    };

    async render(target: HTMLElement, options: PluginOptions) {
        /*
        Render PDB RNA viewer plugin
        */
        if (!target || !options.pdbId || !options.chainId || !options.entityId) {
            console.log('Invalid plugin input!');
            return;
        }
        this.options = options;
        if (this.options.pdbId != "cust") {
            this.apiData = await this.dataService.getApiData(this.options.entityId, this.options.chainId, this.options.pdbId);
            this.FR3DData = await this.dataService.getFR3DData(this.options.pdbId, this.options.chainId);
            this.FR3DNestedData = await this.dataService.getFR3DNestedData(this.options.pdbId, this.options.chainId);
        }
        this.targetEle = <HTMLElement>target;


        this.uiTemplateService = new UiTemplateService(this.targetEle, this.options, this.apiData);
        if (this.apiData) {
            // draw topology
            this.uiTemplateService.render(this.apiData, this.FR3DData, this.FR3DNestedData);

            // Bind to other PDB Component events
            if (this.options.subscribeEvents) {
                CustomEvents.subscribeToComponentEvents(this);
            }

        } else {

            if (options.rvAPI == true) this.rvAPI = true;

            if (this.rvAPI) {
                const dataUrls = (window as any).vm.URL;

                const {
                    apiDataJ,
                    FR3DDataJ,
                    FR3DNestedDataJ
                } = await (await fetch(dataUrls)).json().then((r2dtjson) => ({
                    apiDataJ: r2dtjson.RNA_2D_json as ApiData,
                    FR3DDataJ: r2dtjson.RNA_BP_json as any,
                    FR3DNestedDataJ: r2dtjson.RNA_BP_json as any,
                }));
                this.apiData = apiDataJ;
                this.FR3DData = FR3DDataJ;
                this.FR3DNestedData = FR3DNestedDataJ;

                // draw topology
                this.uiTemplateService.render(this.apiData, this.FR3DData, this.FR3DNestedData);


            };

        }
    }

    selectResidue(label_seq_id: number, color?: string) {
        //Color selected residue
        UiActionsService.selectNucleotide(this.options.pdbId, this.options.entityId, label_seq_id, 'mouseover', false, color);
    }

    clearSelection(label_seq_id: number) {
        //Clear selected residue
        UiActionsService.unSelectNucleotide(this.options.pdbId, this.options.entityId, label_seq_id, false);
    }

    highlightResidue(label_seq_id: number, color?: string) {
        //Highlight selected residue
        UiActionsService.colorNucleotide(this.options.pdbId, label_seq_id, color, 'highlight');
    }

    clearHighlight() {
        //Clear highlighted residue
        UiActionsService.clearHighlight(this.options.pdbId);
    }
}
//Make plugin accessible to the window
(window as any).PdbRnaViewerPlugin = PdbRnaViewerPlugin;
