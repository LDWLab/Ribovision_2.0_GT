import {ajaxProper} from './ajaxProper.js'

export function parseRNAchains(chainList){
    let entities = [];
    chainList.forEach(function(entry){
        if (entry.molecule_type == "polyribonucleotide"){
            entities.push({ 
                entityID: entry["entity_id"], 
                chainID: entry["in_chains"][0]
            });
            vm.rnaChains.push({
                text: entry["molecule_name"][0],
                value: entry["in_chains"][0],
                sequence: entry["sequence"],
                entityID: entry["entity_id"],
            });
        }
    });
    postRNACIFdata (vm.pdbid, entities);
}

var postRNACIFdata = function (pdbID, entities){
    vm.postedRNAEntities = false;
    let parseURL = `custom-struc-data/${pdbID}`;
    var stringEntities = JSON.stringify(entities); 
    ajaxProper({
        url: parseURL,
        type: 'POST',
        dataType: 'json',
        postData: {"entities": stringEntities}
    }).then (parsedResponse => {
        if (parsedResponse == "Success!"){
            vm.postedRNAEntities = true;
        }
    }).catch(error => {
        console.log(error);
    });
}

var postedRNAEntities = function (successPost){
    //watch
    if (successPost){
        let tempEntities = vm.chains.filter(obj => {
            return obj["value"] == vm.chainid[0];
        });
        let entities = [];
        tempEntities.forEach(function(ent){
            entities.push({ entityID: ent["entityID"], chainID: ent["value"] })
        });
        var stringEntities = JSON.stringify(entities); 
        let parseURL = `custom-struc-data/${vm.pdbid}`;
        ajaxProper({
            url: parseURL,
            type: 'POST',
            dataType: 'json',
            postData: {"entities": stringEntities}
        }).then (parsedResponse => {
            if (parsedResponse == "Success!"){
                vm.repostedPDBentity = true;
            }
        }).catch(error => {
            console.log(error);
        });
    }
}

var reloadMolStarWithRNA = function (checkedRNA){
    //method
    if (checkedRNA){
        vm.loadingRNAProt = true;
        let entities = vm.chains.filter(obj => {
            return obj["value"] == vm.chainid[0];
        });
        var rnaString = ''
        vm.rnaChains.forEach(function(rnaChain){
            rnaString += `${vm.pdbid}-${rnaChain.entityID}-${rnaChain.value},`
        })
        var coordURL = `/custom-struc-data/${vm.pdbid}-${entities[0].entityID}-${vm.chainid[0]},${rnaString.slice(0, -1)}`;
        var options = {
                customData: { url: coordURL,
                            format: "cif", 
                            binary: false },
            hideCanvasControls: ["expand", "selection", " animation"],
            assemblyId: '1',
            hideControls: true,
            subscribeEvents: true,
            visualStyle: "gaussian-surface",
            bgColor: {r:255,g:255,b:255},
        }
        viewerInstance.visual.update(options);
        viewerInstance.events.loadComplete.subscribe(() => { 
            vm.loadingRNAProt = false;
            vm.$nextTick(function(){
                viewerInstance.visual.focus([{entity_id: String(entities[0].entityID)}]);
            })
        });
    } else {
        if (vm.pdbid&&vm.chainid&&vm.chains){
            let entities = vm.chains.filter(obj => {
                return obj["value"] == vm.chainid[0];
            });
            vm.showPDBViewer(vm.pdbid, vm.chainid[0], entities[0].entityID);
        }
        
    }
}

{/*
<p><div v-if="!(postedRNAEntities===null)&&postedRNAEntities===false">
    Parsing RNA structures <img src="static/img/loading.gif" alt="Loading RNA structure" style="height:25px;">
</div></p>
    
<div v-if="topology_loaded">
<p><div v-if="postedRNAEntities&&repostedPDBentity" class="checkbox" id="showRNAmolecules">
    <label><input type="checkbox" v-model="checkedRNA" v-on:change="reloadMolStarWithRNA(checkedRNA)">
        Show RNA in 3D viewer</label>
</div>
<div v-if="loadingRNAProt">
    Loading RNA-prot structure <img src="static/img/loading.gif" alt="Loading RNA-prot structure" style="height:25px;">
</div></p>
</div> 

*/}