export function filterAvailablePolymers(chain_list, aln_id, vueObj) {
    var nomIDs = [];
    var geneDescs = [];
    var cleanedChains = [];
    let url = `/desire-api/alignments/${aln_id}/?format=json`;
    ajax(url).then( aln_data => {
        for (let ix =0; ix < aln_data["polymers"].length; ix++){
            let desirePolymerName = aln_data["polymers"][ix]["genedescription"].trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'');
            let desireNomList = aln_data["polymers"][ix]["nomgd"].split('/');
            nomIDs.push(desireNomList[desireNomList.length - 2]);
            geneDescs.push(desirePolymerName);
        }
        for (let i = 0; i < chain_list.length; i++) {
            let chain_listI = chain_list[i]
            if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
            if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
            if (chain_listI["molecule_type"].toLowerCase() == "polyribonucleotide") {continue;}
            let pdbePolymerNames = chain_listI["molecule_name"];
            for (let nameIx =0; nameIx < pdbePolymerNames.length; nameIx++){
                let polName = pdbePolymerNames[nameIx].replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'');
                const nameIndex = cleanedChains.findIndex(chain => chain.text === polName);
                if (nameIndex === -1){
                    try{
                        cleanedChains.push({
                            text: polName,
                            value: chain_listI["in_chains"][0],
                            entityID: chain_listI["entity_id"],
                            startIndex: chain_listI.source[0].mappings[0].start.residue_number,
                            endIndex: chain_listI.source[0].mappings[0].end.residue_number,
                            sequence: chain_listI["sequence"]
                        });
                    }catch(err){console.log(err);}
                }
            }
        }
        var uniqGeneDescs = geneDescs.filter((v, i, a) => a.indexOf(v) === i);
        var uniqNomIDs = nomIDs.filter((v, i, a) => a.indexOf(v) === i);
        var filteredChains = cleanedChains.filter(e => uniqGeneDescs.indexOf(e.text) !== -1);
        if (filteredChains.length === 0) {
            searchNewNames(uniqNomIDs, cleanedChains, vueObj);
            return;
        } else if (filteredChains.length === 2){
            //Cleanup here
            vueObj.hide_chains = null;
        } else {
            vueObj.hide_chains = null;
        }
        vueObj.chains = filteredChains;
    });
}

var searchNewNames = function(nomIDS, chainList, vueObj){
    let url = `/desire-api/old-nomenclatures/?format=json&nn_fk__in=${nomIDS.join(',')}`
    ajax(url).then( nomData => {
        var matchingChains = [];
        var filteredChains = [];
        nomData.results.forEach(nomResult => {
            var oldName = nomResult.old_name.replace('L0','L').replace('S0','S');
            matchingChains.push(chainList.filter(e => 
                e.text.split(' ').some((elt) => elt == oldName)
            ));
        });
        if (matchingChains.flat().length === 0){
            var elt = document.querySelector("#onFailedChains");
            elt.innerHTML  = "Couldn't find a matching chain!<br>Try a different PDB ID."
            vueObj.pdbid = null;
            filteredChains.push({text: "Couldn't find polymers from this structure!", value: null})
        } else {
            vueObj.hide_chains = null;
            filteredChains = matchingChains.flat().filter((v,i,a) => a.indexOf(v) === i);
        }
        vueObj.chains = filteredChains;
    });
}