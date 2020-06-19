
function ajax(url) {
	return new Promise((resolve, reject) => {
		$.ajax({
			url: url,
			type: 'GET',
			dataType: "json",
			success: function(data) {
				resolve(data)
			},
			error: function (error) {
				console.log(`Error ${error}`);
				reject(error)
			}
		})
	})
}

Vue.component('treeselect', VueTreeselect.Treeselect,)

new Vue({
	el: '#phylo_tree_dropdown',
	delimiters: ['[[',']]'],
	data: {
		value: null,
		valuef: null,
		options: null,
		alignments: null,
		pdbid: null,
		chains: null,
		chainid: null
	},
	methods: {
		limiter(e) {
			if(e.length > 2) {
				alert('You can only select two groups!')
				e.pop()
			}
		},
		loadOptions({action, callback}) {
			if (action === "LOAD_ROOT_OPTIONS"){
				ajax('/alignments/showTaxonomy').then(data =>{
					data.isDisabled = true,
					this.options = [data];
					callback();
				}).catch(error => {
					console.log(error)
				})
			}
		},
		loadData:function(value) {
			this.alignments = null;
			ajax('/alignments/fold-api/'+value)
			.then(data=>{
				//Fix up our data
				var fpa = data["Folds to polymers to alignments"]
				var fpa_viz = [];
				Object.keys(fpa).forEach(fkey => {
					Object.keys(fpa[fkey]).forEach(pkey => {
						fpa[fkey][pkey].forEach(function (akey){
							fpa_viz.push({
								text:  'Alignment '.concat(akey[1],'; fold ',fkey),
								value: fkey.concat(',',akey)
							});
						});
					});
				});
				var temp_arr = fpa_viz
				fpa_viz = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
				this.alignments = fpa_viz
			});
		},
		getPDBchains(pdbid){
			if (pdbid.length === 4){
				ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/'+pdbid.toLowerCase())
					.then(struc_data =>{
						var chain_options = [];
						for (var i = 0; i < struc_data[pdbid.toLowerCase()].length; i++) {
							if (struc_data[pdbid.toLowerCase()][i]["molecule_type"] == "Bound"){
								continue;}
							if (struc_data[pdbid.toLowerCase()][i]["molecule_type"] == "Water"){
								continue;}
							chain_options.push({
								text: struc_data[pdbid.toLowerCase()][i]["molecule_name"][0],
								value: struc_data[pdbid.toLowerCase()][i]["in_chains"][0]
							})
						}
						var temp_arr = chain_options
						chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
						this.chains = chain_options
					}).catch(error => {
						alert("No such pdb id: "+pdbid+".",error)
					})
			}
		},
		showAlignment(aln_id){
			var url = '/paralog-aln-api/'+aln_id.split(',')[1]
			ajax(url).then(fasta =>{
				var opts = {
					el: document.getElementById("alnDiv"),
					seqs: msa.io.fasta.parse(fasta),
					colorscheme: {
						scheme: "clustal2",
					},
					zoomer: {
						// general
						alignmentWidth: 1000,
						alignmentHeight: 750,
						columnWidth: 15,
						rowHeight: 15,
						labelNameLength: 300,
						autoResize: true, // only for the width
					},
					// smaller menu for JSBin
					menu: "small",
					bootstrapMenu: true
				};

				var m = new msa.msa(opts);
				m.render();
			})
		}
	}
})
