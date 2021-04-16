import CIFTools from './CIFTools.js' //From here https://github.com/dsehnal/CIFTools.js

export function generateChainsFromLiteMol(pdbid){
    var oReq = new XMLHttpRequest();
    oReq.open("GET", `https://coords.litemol.org/${pdbid.toLowerCase()}/assembly?id=1&lowPrecisionCoords=1&encoding=BCIF`, true);
    oReq.responseType = "arraybuffer";
    oReq.onload = function(oEvent) {
        if(vm.unfilteredChains){return;}
        var parsed = CIFTools.Binary.parse(oReq.response);
        if (parsed.isError) {
            console.log(parsed.toString());
            vm.pdbid = null;
            return;
        }
        var data = parsed.result.dataBlocks[0];
        var _entity_poly = data.getCategory('_entity_poly')
        var sequenceList = _entity_poly.getColumn("pdbx_seq_one_letter_code_can").data;
        var chainList = _entity_poly.getColumn("pdbx_strand_id").data;
        var entityList = _entity_poly.getColumn("entity_id").data;
        var typeList = _entity_poly.getColumn("type").data;
        var descriptList = data.getCategory('_entity').getColumn('pdbx_description').data;
        var _atom_site = data.getCategory('_atom_site');
        var formattedChains = new Array(entityList.length);
        entityList.forEach((entity, i)=>{
            let chainID = chainList[i].split(',')[0];
            let startPos = _atom_site.getColumn('auth_asym_id').data.indexOf(chainID);
            let endPos = startPos + _atom_site.getColumn('auth_asym_id').data.filter(x => x==chainID).length - 1;
            let startIndex = _atom_site.getColumn('auth_seq_id').getInteger(startPos);
            let endIndex = _atom_site.getColumn('auth_seq_id').getInteger(endPos);
            formattedChains[i] = {
                molecule_type: typeList[i],
                sequence: sequenceList[i].replace('\n',''),
                entity_id: parseInt(entity),
                in_chains: [chainID],
                molecule_name: [descriptList[i]],
                source: [{
                    mappings: [{
                        start:{residue_number:startIndex},
                        end: {residue_number:endIndex}
                    }]
                }]
            }
        })
        vm.unfilteredChains = formattedChains;
    };
    oReq.send();
}