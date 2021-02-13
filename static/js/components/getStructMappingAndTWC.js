export function getStructMappingAndTWC (fasta, struc_id, vueObj){
    ajax('/mapSeqAln/', {fasta, struc_id}).then(struct_mapping=>{
        vueObj.structure_mapping = struct_mapping;
        if (struct_mapping['BadMappingPositions']){vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];}
        var mapped_aa_properties = mapAAProps(vueObj.aa_properties, struct_mapping);
        var topviewer = document.getElementById("PdbeTopViewer");
        if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
            if (vueObj.unmappedTWCdata) {
                mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_aa_properties);
            }
        }
        window.mapped_aa_properties = mapped_aa_properties;
        topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
        topviewer.pluginInstance.createDomainDropdown();
    }).catch(error => {
        //var topview = document.querySelector('#topview');
        console.log(error);
        //vueObj.topology_loaded = 'error';
        //topview.innerHTML = "Failed to load the viewer!<br>Try another structure."
    });
}

const build_mapped_props = function(mapped_props, twcDataUnmapped, structure_mapping){
    mapped_props.set("TwinCons", [])
    for (let i = 0; i < twcDataUnmapped.length; i++) {
        let mappedI0 = structure_mapping[twcDataUnmapped[i][0]];
        if (mappedI0) {
            mapped_props.get("TwinCons").push([mappedI0, twcDataUnmapped[i][1]]);
        }
    }
    return mapped_props;
}