function include(file) {
      
    var script = document.createElement('script');
    script.src = file;
    script.type = 'text/javascript';
    script.defer = true;
      
    document.getElementsByTagName('head').item(0).appendChild(script);
      
    }
      
/* Include Many js files */
import $ from '/home/hmccann3/Ribovision_3/Ribovision_3.0_GT/static/alignments/external/jquery.min.js' 
export async function map_seq_aln (fasta, struc_id) {
    var largestKey = -10
    var smallestKey = -10
    var gaps = [1]
    structMappingAndData = await doAjax(fasta, struc_id)
    var struct_mapping = structMappingAndData["structureMapping"];
    largestKey = Math.max(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
    smallestKey = Math.min(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
    gaps = structMappingAndData["gapsInStruc"]
    //return [largestKey, smallestKey, gaps.length]
}

async function doAjax(fasta, struc_id) {
    return $.ajax('/mapSeqAln/', {fasta, struc_id})
}

