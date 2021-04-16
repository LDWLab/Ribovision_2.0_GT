import CIFTools from './CIFTools.js' //From here https://github.com/dsehnal/CIFTools.js

var fetchParseLiteMolData = function (url){
    return new Promise((resolve, reject) => {
        var oReq = new XMLHttpRequest();
        oReq.open("GET", url , true);
        oReq.responseType = "arraybuffer";
        oReq.onload = function(oEvent) {
            var parsed = CIFTools.Binary.parse(oReq.response);
            if (parsed.isError) {
                console.log(parsed.toString());
                vm.pdbid = null;
                reject(parsed);
            }
            resolve(parsed.result.dataBlocks[0]);
        };
        oReq.send();
    })
}

var parseCifResults = function (data, vueAssignmentName){
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
    vm[vueAssignmentName] = formattedChains;
}

export function generateChainsFromLiteMol(url, vueAssignmentName){
    fetchParseLiteMolData(url).then(data=>{
        parseCifResults(data, vueAssignmentName)
    })
}

export function generateStartEndFromLiteMol(url){
    fetchParseLiteMolData(url).then(data=>{
        var _atom_site = data.getCategory('_atom_site');
        vm.pdbStart = _atom_site.getColumn('auth_seq_id').getInteger(0);
        vm.pdbEnd = _atom_site.getColumn('auth_seq_id').getInteger(_atom_site.getColumn('auth_seq_id').data.length-1);
    })
}