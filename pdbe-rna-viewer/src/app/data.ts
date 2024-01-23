var rv3VUEcomponent = (window as any).vm;
export type ThemeParam = {
    color?: string,
    highlightColor?: string,
    unobservedColor?: string,
    strokeWidth?: number
}

export type PluginOptions = {
    pdbId: string,
    chainId: string,
    entityId: string,
    theme?: ThemeParam,
    subscribeEvents?: boolean,
    rvAPI?: boolean
}

export type ApiData = {
    svg_paths: string[],
    dimensions: {
        width: number,
        height: number
    },
    sequence: string,
    unobserved_label_seq_ids: number[],
    auth_seq_ids: number[],
    label_seq_ids: number[],
    pdb_ins_codes: string[]
}

export class DataService {
    async getApiData(entityId: string, chainId: string, pdbId: string): Promise<ApiData | undefined> {
        try {
            const apiUrl = `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId.toLowerCase()}_${entityId}_${chainId}.json`;
            return await (await fetch(apiUrl)).json() as ApiData;
        } catch (e) { 
            this.handleError(e)
            return void 0;
        };
    }
    async getFR3DData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        try {
            const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=False`
            //const csvUrl = `http://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}`
            //const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=True`
            //const csvUrl = `http://rna.bgsu.edu/rna3dhub/pdb/${pdbId.toLowerCase()}/interactions/fr3d/basepairs/csv`;
            return await (await fetch(csvUrl)).json() as JSON;
        } catch (e) { 
            this.handleFR3DError(e)
            return void 0;
        };
    }
    async getFR3DNestedData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        try {
            const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=True`
            //const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=True`
            //const csvUrl = `http://rna.bgsu.edu/rna3dhub/pdb/${pdbId.toLowerCase()}/interactions/fr3d/basepairs/csv`;
            return await (await fetch(csvUrl)).json() as JSON;
        } catch (e) { 
            this.handleFR3DError(e)
            return void 0;
        };
    }

    private handleError(e: any): void {
        rv3VUEcomponent.structFailed = true
        console.log(`RNA topology data not available!`, e);
    }
    private handleFR3DError(e: any): void {
        console.log(`FR3D mapping data not available!`, e);
    }

}
/*export class BanNameHelper {
    public static async getBanName(pdbId : string, PchainId : string): Promise<JSON | undefined> {   
           try {
               const apiUrl = `https://api.ribosome.xyz/neo4j/get_banclass_for_chain/?pdbid=${pdbId}&auth_asym_id=${PchainId}&format=json`  
               return await (await fetch(apiUrl)).json() as JSON;  
           } catch (e) {
               //console.log(`Ban naming is not available!`, e);
               return void 0;
           };
       }
    }*/

    export class BanNameHelper {
        private static banNameMap: Map<string, JSON | undefined> = new Map();
      
        public static async getBanName(pdbId: string, PchainId: string): Promise<JSON | undefined> {
          const cacheKey = `${pdbId}_${PchainId}`;
      
          // Check if the value is in the map
          if (BanNameHelper.banNameMap.has(cacheKey)) {
         
            return BanNameHelper.banNameMap.get(cacheKey);
          }
          const vm = (window as any).vm;

          if (!vm || !vm.unfilteredChains_orig) {
            console.error("vm or vm.unfilteredChains_orig is not defined");
            return void 0; // or handle the error appropriately
          }
      
          try {
           
            var matches = vm.unfilteredChains_orig.filter(function(entry:any) {
            return entry.in_chains.includes(PchainId);
            });

            for (let chain of vm.protein_chains){
                chain.banname_2D=matches[0].molecule_name[0].replace('Large ribosomal subunit', 'LSU').replace('Small ribosomal subunit', 'SSU');
            
            }
            return matches.map(function(entry : any) {

            return entry.molecule_name[0].replace('Large ribosomal subunit protein ', '').replace('Small ribosomal subunit protein ', '');
            });       
            

          } catch (e) {
            console.error(`Error fetching ban name for ${cacheKey}`, e);
            return void 0;
          }
        }
      }
      
      
