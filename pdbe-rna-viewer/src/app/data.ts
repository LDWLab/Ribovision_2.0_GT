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
            const apiUrl = `/api-proxy/pdbe/static-entry/?pdbid=${pdbId.toLowerCase()}&entity_id=${entityId}&chain_id=${chainId}`;
            return await (await fetch(apiUrl)).json() as ApiData;
        } catch (e) { 
            this.handleError(e)
            return void 0;
        };
    }       
    async getFR3DData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        const apiUrl = `/api-proxy/fr3d/data/?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=False`;
        try {
            const resp = await fetch(apiUrl);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            return await resp.json() as JSON;
        } catch (e) {
            console.warn(`getFR3DData failed`, e);
            this.handleFR3DError(e);
            return void 0;
        }
    }
    async getFR3DNestedData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        const apiUrl = `/api-proxy/fr3d/data/?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=True`;
        try {
            const resp = await fetch(apiUrl);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            return await resp.json() as JSON;
        } catch (e) {
            console.warn(`getFR3DNestedData failed`, e);
            this.handleFR3DError(e);
            return void 0;
        }
    }

    private handleError(e: any): void {
        rv3VUEcomponent.structFailed = true
        console.log(`RNA topology data not available!`, e);
    }
    private handleFR3DError(e: any): void {
        console.log(`FR3D mapping data not available!`, e);
    }

}
    export class BanNameHelper {
        private static banNameMap: Map<string, JSON | undefined> = new Map();
        private static hasWarnedAboutVm: boolean = false;
      
        public static async getBanName(pdbId: string, PchainId: string): Promise<JSON | undefined> {
          const cacheKey = `${pdbId}_${PchainId}`;
      
          // Check if the value is in the map
          if (BanNameHelper.banNameMap.has(cacheKey)) {
         
            return BanNameHelper.banNameMap.get(cacheKey);
          }
          const vm = (window as any).vm;

          if (!vm || !vm.unfilteredChains_orig) {
            // Only warn once to avoid spamming console during initialization
            if (!BanNameHelper.hasWarnedAboutVm) {
              console.warn("vm.unfilteredChains_orig not yet available, will retry when ready");
              BanNameHelper.hasWarnedAboutVm = true;
            }
            return void 0; // Return undefined, caller will handle gracefully
          }
          
          // Reset warning flag since data is now available
          BanNameHelper.hasWarnedAboutVm = false;
      
          try {
           
            var matches = vm.unfilteredChains_orig.filter(function(entry:any) {
            return entry.in_chains && entry.in_chains.includes(PchainId);
            });

            if (!matches || matches.length === 0 || !matches[0] || !matches[0].molecule_name) {
              return void 0;
            }

            if (vm.protein_chains) {
              for (let chain of vm.protein_chains){
                  if (matches[0] && matches[0].molecule_name) {
                    chain.banname_2D=matches[0].molecule_name[0].replace('Large ribosomal subunit', 'LSU').replace('Small ribosomal subunit', 'SSU');
                  }
              }
            }
            return matches.map(function(entry : any) {
                if(entry.molecule_name) {
                    return entry.molecule_name[0].replace('Large ribosomal subunit protein ', '').replace('Small ribosomal subunit protein ', '');
                } else if (entry[0] && entry[0].molecule_name) {
                    return entry[0].molecule_name[0].replace('Large ribosomal subunit protein ', '').replace('Small ribosomal subunit protein ', '');
                }
                return void 0;
            });       
            

          } catch (e) {
            // Silently handle errors during initialization
            return void 0;
          }
        }
      }
      
      
