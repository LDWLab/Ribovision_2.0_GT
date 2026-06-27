export async function getBanName(pdbId , PchainId){
           try {
               const apiUrl = `/api-proxy/ban-name/?pdbid=${pdbId}&auth_asym_id=${PchainId}`
               return await (await fetch(apiUrl)).json();
           } catch (e) {
               //console.log(`Ban naming is not available!`, e);
               return void 0;
           };
       }