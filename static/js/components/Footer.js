export function addFooterImages (divID){
    const htmlTOinject = `
    <div class="white-box" style="float: left;">
        <a href="http://apollo.chemistry.gatech.edu/RiboVision/" target="_blank">
            <img 
                style="height:75px; padding:5px;"
                src="static/ribovision/images/RiboVisionLogo.png" 
                title="RiboVision" 
                alt="RiboVision Logo">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="https://ribosome.xyz/home/" target="_blank">
            <img 
                style="height:75px; padding:5px;"
                src="static/alignments/png/riboXYZ.png" 
                title="Ribosome.xyz" 
                alt="Comprehensive Resource for Ribosome Structures">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="http://ww2.chemistry.gatech.edu/~lw26/index.html" target="_blank">
            <img 
                style="height:75px; padding:10px;padding-top:5px;"
                src="static/alignments/png/Williams_lab-logo-01.png" 
                title="The Williams Lab" 
                alt="Williams Logo">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="https://cos.gatech.edu/" target="_blank">
            <img 
                style="height:75px; padding:10px;padding-top:5px;"
                src="static/alignments/png/GT-logo.png" 
                title="Georgia Institute of Technology" 
                alt="GaTech Logo">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="http://cool.gatech.edu/" target="_blank">
            <img 
                style="height:75px; padding:10px;padding-top:2px;"
                src="static/ribovision/images/cool_logo.png" 
                title="Center for the Origin Of Life (COOL)" 
                alt="COOL Logo">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="http://prebioticchem.info//" target="_blank">
            <img 
                style="height:75px; padding:10px;padding-top:7px;"
                src="static/alignments/png/PCE3_logo-dark.png" 
                title="Prebiotic Chemistry and Early Earth Environments Consortium" 
                alt="PCE3 Logo">
        </a>
    </div>
    <p style="padding:5px;"></p>
    <div class="white-box" style="float: left;">
        <a href="https://astrobiology.nasa.gov/" target="_blank">
            <img 
                style="height:75px; padding:5px;padding-right:2px;"
                src="/static/ribovision/images/NASALogo.png" 
                class="ComboLogo" 
                title="NASA Astrobiology Institute (NAI)" alt="NASA Logo">
        </a>
    </div>
    <a href="mailto:RiboZones@gmail.com?subject=ProteoVision%20question">
        <div class="white-box" style="position: absolute;right: 35px;height:78px; padding:5px;">
            <div style="padding-top:20px"><b>Contact us</b></div>
        </div>
    </a>
</div>`
    var injectionDiv = document.getElementById(divID)
    injectionDiv.innerHTML = htmlTOinject;
}