$(function() {
    $('#user_input_pdb_form').submit(function() {
        find_chains_from_PDBe_api($('form').serializeObject()['pdbid'].toLowerCase())
        return false;
    });
});

$.fn.serializeObject = function(){
    var o = {};
    var a = this.serializeArray();
    $.each(a, function() {
        if (o[this.name] !== undefined) {
            if (!o[this.name].push) {
                o[this.name] = [o[this.name]];
            }
            o[this.name].push(this.value || '');
        } else {
            o[this.name] = this.value || '';
        }
    });
    return o;
};

function find_chains_from_PDBe_api(struc_id) {
    pdb = struc_id
    $.ajax({
    url: 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/'+struc_id,
    type: "GET",
    dataType: "json",
    success: function (data) {
        build_chain_selector(data, struc_id);
    },
    error: function (error) {
        console.log(`Error ${error}`);}
    });
}

function build_chain_selector(struc_data, struc_id){
    if (struc_id.length == 0) document.getElementById("chain_selector").innerHTML = "<option></option>";
    else {
        var catOptions = "";
        for (var i = 0; i < struc_data[struc_id.toLowerCase()].length; i++) {
            if (struc_data[struc_id.toLowerCase()][i]["molecule_type"] == "Bound"){
                continue;
            }
            if (struc_data[struc_id.toLowerCase()][i]["molecule_type"] == "Water"){
                continue;
            }
            catOptions += "<option value="+struc_data[struc_id.toLowerCase()][i]["in_chains"][0]+">" + struc_data[struc_id.toLowerCase()][i]["molecule_name"][0] + "</option>";
        }
        document.getElementById("chain_selector").innerHTML = catOptions;
    }
    return false;
}