export async function getBanName(pdbId , PchainId){   
           try {
               const apiUrl = `/extapi/ribosome/banclass/${pdbId}/${PchainId}`  
               return await (await fetch(apiUrl)).json();  
           } catch (e) {
               //console.log(`Ban naming is not available!`, e);
               return void 0;
           };
       }