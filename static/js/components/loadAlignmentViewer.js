import {AlnViewer} from './AlignmentViewer.js'
import ReactDOM from 'react-dom';
import React from "react";
export function loadAlignmentViewer (fasta){
    if (fasta){
        var main_elmnt = document.querySelector(".alignment_section");
        var msaHeight = main_elmnt.offsetHeight * 0.8;
        if (msaHeight > 17*(vm.fastaSeqNames.length+2)){
            msaHeight = 17*(vm.fastaSeqNames.length+2);
        }
        let seqsForMSAViewer = parseFastaSeqForMSAViewer(fasta);
        var msaOptions = {
            sequences: seqsForMSAViewer,
            colorScheme: vm.colorScheme,
            height: msaHeight,
            width: main_elmnt.offsetWidth * 0.7,
            tileHeight: 17,
            tileWidth: 17,
        };
        window.msaOptions = msaOptions;
        ReactDOM.render(
            <AlnViewer ref={(PVAlnViewer) => {window.PVAlnViewer = PVAlnViewer}}/>,
            document.getElementById('alnDiv')
          );
    }
}