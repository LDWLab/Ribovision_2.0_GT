export async function getBanName(pdbId , PchainId){   
           try {
               const apiUrl = `https://api.ribosome.xyz/neo4j/get_banclass_for_chain/?pdbid=${pdbId}&auth_asym_id=${PchainId}&format=json`  
               return await (await fetch(apiUrl)).json();  
           } catch (e) {
               //console.log(`Ban naming is not available!`, e);
               return void 0;
           };
       }