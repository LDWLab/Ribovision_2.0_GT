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
        /*
        Get RNA structure topology data from PDBe API
        @param entityId {string}: The entity ID of the PDB molecule
        @param chainId {string}: The chain ID of the PDB molecule
        @param pdbId {string}: The PDB ID of the PDB molecule
        @return {ApiData | undefined}: The PDBe topology API data or undefined if unavailable
        */
        try {
            const apiUrl = `https://www.ebi.ac.uk/pdbe/static/entry/${pdbId.toLowerCase()}_${entityId}_${chainId}.json`;
            return await (await fetch(apiUrl)).json() as ApiData;
        } catch (e) {
            this.handleError(e)
            return void 0;
        };
    }
    async getFR3DData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        /*
        Get non-nested RNA base pairing data from FR3D API
        @param pdbId {string}: The PDB ID of the PDB molecule
        @param chainId {string}: The chain ID of the PDB molecule
        @return {JSON | undefined}: The FR3D API data or undefined if unavailable
        */
        try {
            const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=False`
            return await (await fetch(csvUrl)).json() as JSON;
        } catch (e) {
            this.handleFR3DError(e)
            return void 0;
        };
    }
    async getFR3DNestedData(pdbId: string, chainId: string): Promise<JSON | undefined> {
        /*
        Get nested RNA base pairing data from FR3D API
        @param pdbId {string}: The PDB ID of the PDB molecule
        @param chainId {string}: The chain ID of the PDB molecule
        @return {JSON | undefined}: The FR3D API data or undefined if unavailable
        */
        try {
            const csvUrl = `https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=${pdbId.toLowerCase()}&chain=${chainId}&only_nested=True`
            return await (await fetch(csvUrl)).json() as JSON;
        } catch (e) {
            this.handleFR3DError(e)
            return void 0;
        };
    }

    private handleError(e: any): void {
        /*
        Error message for PDBe API
        */
        console.log(`RNA topology data not available!`, e);
    }
    private handleFR3DError(e: any): void {
        /*
        Error message for FR3D API
        */
        console.log(`FR3D mapping data not available!`, e);
    }

}
export class BanNameHelper {
    /*
    BanNameHelper class
    */
    public static async getBanName(pdbId: string, PchainId: string): Promise<JSON | undefined> {
        /*
        Returns ban name for PDB chain if available from ribosome.xyz
        @param pdbId {string}: The PDB ID of the PDB molecule
        @param PchainId {string}: The protein chain ID
        @return {JSON | undefined}: The ribosome.xyz data or undefined if unavailable
        */
        try {
            const apiUrl = `https://api.ribosome.xyz/neo4j/get_banclass_for_chain/?pdbid=${pdbId}&auth_asym_id=${PchainId}&format=json`
            return await (await fetch(apiUrl)).json() as JSON;
        } catch (e) {
            //console.log(`Ban naming is not available!`, e);
            return void 0;
        };
    }
}
