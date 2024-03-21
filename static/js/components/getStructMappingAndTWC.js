import { loadAlignmentViewer } from './loadAlignmentViewer.js'
import { ajaxProper } from './ajaxProper.js'


const typeMappings = {
  'aes': { '1a': '1', '2': '1', '10': '1', '13': '1', '7': '1', '17': '1', '22': '1', '0': '2', '7b': '2', '9': '2', '11': '2', '16': '2', '20': '2', '18': '2', '3a': '2', '3b': '3', '5a': '3', '4': '3', '15': '3', '24': '3', '6': '3', '23': '3', '27': '3', '1': '4', '4a': '4', '7a': '4', '8': '4', '14a': '4', '25': '4', '26': '4', '21': '4', '1b': '5', '3': '5', '5': '5', '12': '5', '14': '5', '19': '5', '20a': '5' },
  'AES': { '1.0': '1', '12.0': '1', '19.0': '1', '22.0': '1', '7.0': '1', '34.0': '1', '35.0': '1', '49.0': '1', '28.0': '1', '20.0': '1', '52.0': '1', '58.0': '1', '37.0': '1', '53': '2', '47': '3', '8': '2', '54': '2', '10': '2', '16': '2', '32': '2', '5': '2', '12a': '2', '18': '2', '50': '2', '4a': '3', '21': '3', '4': '3', '26': '3', '14': '3', '25': '3', '33': '3', '38': '3', '11': '3', '51': '3', '39': '3', '40': '3', '6': '4', '6a': '4', '55': '4', '45': '4', '10a': '4', '23': '4', '31': '4', '56': '4', '15': '4', '2': '4', '43': '4', '46': '4', '30': '4', '59': '4', '41': '5', '24': '5', '17': '5', '48': '5', '9': '5', '32a': '5', '57': '5', '15a': '5', '29': '5', '44': '5', '36': '5', '3': '5' },
  'helix': { '2': '1', '6': '1', '9_3a': '1', '9_3b': '1', '13': '1', '17': '1', '21_6d': '1', '24': '1', '26a': '1', '33': '1', '36': '1', '40': '1', '42': '1', '44': '1', '41_10': '1', '3': '2', '6a': '2', '10': '2', '11': '2', '14': '2', '15': '2', '21_6b': '2', '22': '2', '27': '2', '28': '2', '31': '2', '32': '2', '37': '2', '38': '2', '41': '2', '45': '2', '4': '3', '7': '3', '9': '3', '19': '3', '21': '3', '21_3c': '3', '23': '3', '26': '3', '29': '3', '34': '3', '39_9': '3', '44_12': '3', '1': '4', '5': '4', '8': '4', '12': '4', '16': '4', '18': '4', '20': '4', '21_6a': '4', '23a': '4', '25': '4', '26_7': '4', '30': '4', '35': '4', '39': '4', '43': '4' },
  'Helix': { '5S1': '1', '5S2': '3', '5S3': '2', '5S5': '2', '5S4': '4', '5': '1', '9': '1', '13': '1', '18': '1', '21': '1', '23': '1', '25': '1', '25_7b': '1', '27': '1', '30': '1', '35a': '1', '36': '1', '38a': '1', '41': '1', '43a': '1', '47': '1', '49': '1', '51': '1', '55': '1', '58': '1', '64': '1', '78': '1', '79': '1', '74': '1', '82': '1', '86': '1', '90': '1', '94': '1', '98_39b': '1', '63_27': '1', '79_31a': '1', '2': '2', '6': '2', '10': '2', '14': '2', '19': '2', '25_7a': '2', '28': '2', '33': '2', '38': '2', '40': '2', '25a': '2', '42': '2', '52': '2', '49a': '2', '54_20a': '2', '59': '2', '63': '2', '61': '2', '67': '2', '75': '2', '79_31': '2', '83': '2', '88': '2', '92': '2', '97': '2', '99': '2', '31_9': '2', '3': '3', '7': '3', '11': '3', '15': '3', '19a': '3', '22': '3', '25_7d': '3', '26': '3', '32': '3', '34': '3', '39': '3', '44': '3', '45': '3', '48': '3', '49b': '3', '53': '3', '56': '3', '63a': '3', '63_27b': '4', '66': '3', '69': '3', '70??': '3', '76': '3', '79_31c': '3', '80': '3', '84': '3', '91': '3', '93': '3', '95': '3', '98': '3', '101': '3', '9_3': '3', '4': '4', '8': '4', '10_4': '4', '12': '4', '16': '4', '20': '4', '24': '4', '25_7c': '4', '31': '4', '35': '4', '37': '4', '38_12': '4', '43': '4', '46': '4', '50': '4', '52_19': '4', '54': '4', '57': '4', '60': '4', '62': '4', '68': '4', '73': '4', '77': '4', '79_31b': '4', '81': '4', '85': '4', '87': '4', '89': '4', '96': '4', '100': '4', '98_39a': '4' }
}
async function getBanName(pdbId, PchainId) {
  try {
    const apiUrl = `https://api.ribosome.xyz/neo4j/get_banclass_for_chain/?pdbid=${pdbId}&auth_asym_id=${PchainId}&format=json`
    return await (await fetch(apiUrl)).json();
  } catch (e) {
    //console.log(`Ban naming is not available!`, e);
    return void 0;
  };
}

export function getStructMappingAndTWC(fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj, full_sequence_from_pdb = "") {
  vm.structFailed = false
  vm.sequence = ebi_sequence;
  if (vm.fasta_data) {
    let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
    vm.fasta_data = cleanFasta;
  };
  const postData = {
    fasta,
    struc_id,
    "cif_mode_flag": vm.user_uploaded_cif_flag,
    //hardcoded_structure: full_sequence_from_pdb
    hardcoded_structure: vm.customFullSequence
  };
  let structMapping3D;
  // ajax('/mapSeqAln/', postData).then(structMappingAndData=>{
  //   structMapping3D = structMappingAndData["structureMapping"];
  // });
  ajax('/mapSeqAlnOrig/', { fasta, ebi_sequence, startIndex: 1 }).then(structMappingAndData => {
    var struct_mapping = structMappingAndData["structureMapping"];
    vm.struct_to_alignment_mapping = Object.fromEntries(Object.entries(struct_mapping).map(([key, value]) => [value, key]));
    const associatedDataMappedPerType = {
      // "AES" : [[2, 1], [3, 51], [4, 101]]
    };
    console.log("struct_mapping 1", JSON.stringify(struct_mapping));
    for (let [alignmentIndexAsString, structureIndex] of Object.entries(struct_mapping)) {
      let alignmentIndex = Number.parseInt(alignmentIndexAsString);
      // associatedDataCache has correct boundaries
      let associatedDataCache = vm.associatedDataCache;
      //todo:  wrong below 
      if (alignmentIndex in associatedDataCache) {
        for (let { type, value } of associatedDataCache[alignmentIndex]) {
          if (!(type in associatedDataMappedPerType)) {
            associatedDataMappedPerType[type] = [];
          }

          if (type in typeMappings) {
            let selectedDataDict = typeMappings[type] || {};
            value = selectedDataDict[value] || 0;
          } else {
            if (value.length === 0) {
              value = "0";
            }
            value = Number.parseInt(value);
          }
          let associatedDataI = [
            structureIndex,
            value
          ];
          associatedDataMappedPerType[type].push(associatedDataI);
        }
      }
    }
    // todo: insert color patching function here.  
    // associatedDataMappedPerType.helix[0]=[6, '1']
    const fix_colors = require('./graphColorPrediction.js');
    vm.AD_headers = [];
    // console.log(associatedDataMappedPerType);
    // vm.associatedDataMappedPerType = associatedDataMappedPerType;
    vm.associatedDataMappedPerType = fix_colors(
      viewerInstanceTop.viewInstance.uiTemplateService.apiData.sequence,
      viewerInstanceTop.viewInstance.uiTemplateService.baseStrs.get('cWW')[1],
      associatedDataMappedPerType
    );
    vm.associatedDataMappedPerType_2D = vm.associatedDataMappedPerType;
    vm.associatedDataMappedPerType_3D = {"helix":[[1,"4"],[2,"4"],[3,"4"],[4,"4"],[5,"4"],[6,"4"],[7,"4"],[8,"4"],[9,"4"],[10,"4"],[11,"4"],[12,"4"],[13,"4"],[14,"1"],[15,"1"],[16,"2"],[17,"2"],[18,"2"],[19,"2"],[20,"2"],[21,"4"],[22,"4"],[23,"4"],[24,"4"],[25,"4"],[26,"4"],[27,"2"],[28,"4"],[29,"4"],[30,"4"],[31,"4"],[32,"2"],[33,"4"],[34,"4"],[35,"4"],[36,"4"],[37,"4"],[38,"3"],[39,"3"],[40,"3"],[41,"3"],[42,"3"],[43,"3"],[44,"3"],[45,"3"],[46,"3"],[47,"3"],[48,"4"],[49,"4"],[50,"4"],[51,"4"],[52,"4"],[53,"4"],[54,"4"],[55,"4"],[56,"4"],[57,"4"],[58,"4"],[59,"2"],[60,"4"],[61,"4"],[62,"1"],[63,"1"],[64,"1"],[65,"1"],[66,"1"],[67,"1"],[68,"1"],[69,"1"],[70,"1"],[71,"1"],[72,"1"],[73,"1"],[74,"1"],[75,"1"],[76,"1"],[77,"1"],[78,"1"],[79,"1"],[80,"1"],[81,"1"],[82,"1"],[83,"1"],[84,"1"],[85,"1"],[86,"1"],[87,"1"],[88,"1"],[89,"1"],[90,"1"],[91,"1"],[92,"1"],[93,"1"],[94,"1"],[95,"1"],[96,"1"],[97,"1"],[98,"1"],[99,"1"],[100,"1"],[101,"1"],[102,"1"],[103,"1"],[104,"1"],[105,"1"],[106,"1"],[107,"2"],[108,"2"],[109,"4"],[110,"2"],[111,"2"],[112,"2"],[113,"2"],[114,"2"],[115,"2"],[116,"3"],[117,"3"],[118,"3"],[119,"3"],[120,"3"],[121,"3"],[122,"3"],[123,"3"],[124,"3"],[125,"3"],[126,"3"],[127,"3"],[128,"3"],[129,"3"],[130,"3"],[131,"3"],[132,"3"],[133,"3"],[134,"3"],[135,"3"],[136,"3"],[137,"3"],[138,"3"],[139,"3"],[140,"4"],[141,"4"],[142,"4"],[143,"4"],[144,"4"],[145,"4"],[146,"4"],[147,"4"],[148,"4"],[149,"4"],[150,"4"],[151,"4"],[152,"4"],[153,"4"],[154,"4"],[155,"4"],[156,"4"],[157,"4"],[158,"4"],[159,"4"],[160,"4"],[161,"4"],[162,"4"],[163,"4"],[164,"4"],[165,"4"],[166,"4"],[167,"4"],[168,"4"],[169,"4"],[170,"4"],[171,"4"],[172,"4"],[173,"4"],[174,"4"],[175,"4"],[176,"4"],[177,"4"],[178,"2"],[179,"2"],[180,"2"],[181,"2"],[182,"3"],[183,"2"],[184,"3"],[185,"3"],[186,"3"],[187,"3"],[188,"3"],[189,"3"],[190,"3"],[191,"3"],[192,"3"],[193,"3"],[194,"3"],[195,"3"],[196,"2"],[197,"2"],[198,"2"],[199,"2"],[200,"2"],[201,"2"],[202,"2"],[203,"2"],[204,"2"],[205,"3"],[206,"3"],[207,"3"],[208,"3"],[209,"3"],[210,"3"],[211,"3"],[212,"3"],[213,"3"],[214,"3"],[215,"3"],[216,"3"],[217,"3"],[218,"3"],[219,"3"],[220,"3"],[221,"3"],[222,"3"],[223,"3"],[224,"3"],[225,"3"],[226,"3"],[227,"3"],[228,"3"],[229,"3"],[230,"3"],[231,"3"],[232,"3"],[233,"3"],[234,"3"],[235,"3"],[236,"3"],[237,"2"],[238,"2"],[239,"2"],[240,"2"],[241,"2"],[242,"2"],[243,"2"],[244,"2"],[245,"2"],[246,"2"],[247,"2"],[248,"2"],[249,"2"],[250,"2"],[251,"2"],[252,"2"],[253,"2"],[254,"2"],[255,"2"],[256,"2"],[257,"2"],[258,"2"],[259,"2"],[260,"2"],[261,"2"],[262,"2"],[263,"2"],[264,"2"],[265,"2"],[266,"2"],[267,"2"],[268,"2"],[269,"2"],[270,"2"],[271,"2"],[272,"2"],[273,"2"],[274,"2"],[275,"2"],[276,"2"],[277,"2"],[278,"2"],[279,"2"],[280,"2"],[281,"2"],[282,"2"],[283,"2"],[284,"2"],[285,"4"],[286,"4"],[287,"4"],[288,"4"],[289,"4"],[290,"4"],[291,"4"],[292,"4"],[293,"4"],[294,"4"],[295,"4"],[296,"4"],[297,"4"],[298,"4"],[299,"4"],[300,"4"],[301,"4"],[302,"4"],[303,"4"],[304,"4"],[305,"4"],[306,"4"],[307,"4"],[308,"4"],[309,"4"],[310,"4"],[311,"2"],[312,"2"],[313,"1"],[314,"1"],[315,"1"],[316,"1"],[317,"1"],[318,"1"],[319,"1"],[320,"1"],[321,"1"],[322,"1"],[323,"1"],[324,"1"],[325,"1"],[326,"1"],[327,"1"],[328,"1"],[329,"1"],[330,"1"],[331,"1"],[332,"1"],[333,"1"],[334,"1"],[335,"2"],[336,"2"],[337,"2"],[338,"2"],[339,"2"],[340,"2"],[341,"2"],[342,"2"],[343,"2"],[344,"2"],[345,"2"],[346,"2"],[347,"2"],[348,"2"],[349,"2"],[350,"2"],[351,"2"],[352,"2"],[353,"4"],[354,"4"],[355,"4"],[356,"4"],[357,"4"],[358,"4"],[359,"4"],[360,"4"],[361,"4"],[362,"4"],[363,"2"],[364,"2"],[365,"2"],[366,"2"],[367,"2"],[368,"2"],[369,"2"],[370,"2"],[371,"2"],[372,"2"],[373,"2"],[374,"2"],[375,"2"],[376,"2"],[377,"2"],[378,"2"],[379,"2"],[380,"2"],[381,"2"],[382,"2"],[383,"2"],[384,"2"],[385,"2"],[386,"2"],[387,"2"],[388,"2"],[389,"2"],[390,"2"],[391,"2"],[392,"3"],[393,"2"],[394,"3"],[395,"3"],[396,"3"],[397,"3"],[398,"3"],[399,"3"],[400,"3"],[401,"3"],[402,"4"],[403,"4"],[404,"4"],[405,"4"],[406,"4"],[407,"4"],[408,"4"],[409,"4"],[410,"4"],[411,"4"],[412,"4"],[413,"4"],[414,"4"],[415,"4"],[416,"4"],[417,"4"],[418,"4"],[419,"4"],[420,"4"],[421,"4"],[422,"4"],[423,"4"],[424,"4"],[425,"4"],[426,"4"],[427,"4"],[428,"4"],[429,"4"],[430,"4"],[431,"4"],[432,"4"],[433,"4"],[434,"4"],[435,"4"],[436,"4"],[437,"1"],[438,"1"],[439,"1"],[440,"1"],[441,"1"],[442,"1"],[443,"1"],[444,"1"],[445,"1"],[446,"1"],[447,"1"],[448,"1"],[449,"1"],[450,"1"],[451,"1"],[452,"1"],[453,"1"],[454,"1"],[455,"1"],[456,"1"],[457,"1"],[458,"1"],[459,"1"],[460,"1"],[461,"1"],[462,"1"],[463,"1"],[464,"1"],[465,"1"],[466,"1"],[467,"1"],[468,"1"],[469,"1"],[470,"1"],[471,"1"],[472,"1"],[473,"1"],[474,"1"],[475,"1"],[476,"1"],[477,"1"],[478,"1"],[479,"1"],[480,"1"],[481,"1"],[482,"1"],[483,"1"],[484,"4"],[485,"4"],[486,"4"],[487,"4"],[488,"4"],[489,"1"],[490,"1"],[491,"1"],[492,"1"],[493,"1"],[494,"1"],[495,"4"],[496,"4"],[497,"4"],[498,"4"],[499,"4"],[500,"4"],[501,"4"],[502,"4"],[503,"4"],[504,"4"],[505,"4"],[506,"4"],[507,"4"],[508,"4"],[509,"4"],[510,"4"],[511,"4"],[512,"4"],[513,"4"],[514,"4"],[515,"4"],[516,"4"],[517,"4"],[518,"4"],[519,"4"],[520,"4"],[521,"4"],[522,"4"],[523,"4"],[524,"4"],[525,"4"],[526,"4"],[527,"4"],[528,"4"],[529,"4"],[530,"4"],[531,"4"],[532,"4"],[533,"4"],[534,"4"],[535,"4"],[536,"4"],[537,"4"],[538,"4"],[539,"4"],[540,"4"],[541,"4"],[542,"4"],[543,"4"],[544,"4"],[545,"4"],[546,"4"],[547,"4"],[548,"4"],[549,"4"],[550,"4"],[551,"2"],[552,"3"],[553,"3"],[554,"2"],[555,"2"],[556,"2"],[557,"2"],[558,"2"],[559,"2"],[560,"2"],[561,"2"],[562,"2"],[563,"2"],[564,"2"],[565,"3"],[566,"3"],[567,"3"],[568,"2"],[569,"2"],[570,"3"],[571,"3"],[572,"3"],[573,"3"],[574,"3"],[575,"3"],[576,"3"],[577,"3"],[578,"3"],[579,"4"],[580,"3"],[581,"3"],[582,"3"],[583,"3"],[584,"3"],[585,"3"],[586,"3"],[587,"3"],[588,"3"],[589,"3"],[590,"4"],[591,"4"],[592,"4"],[593,"4"],[594,"4"],[595,"4"],[596,"3"],[597,"3"],[598,"3"],[599,"3"],[600,"3"],[601,"3"],[602,"3"],[603,"3"],[604,"3"],[605,"3"],[606,"3"],[607,"3"],[608,"3"],[609,"3"],[610,"3"],[611,"3"],[612,"3"],[613,"3"],[614,"3"],[615,"3"],[616,"3"],[617,"3"],[618,"3"],[619,"3"],[620,"3"],[621,"3"],[622,"3"],[623,"3"],[624,"3"],[625,"3"],[626,"3"],[627,"3"],[628,"3"],[629,"3"],[630,"3"],[631,"3"],[632,"3"],[633,"3"],[634,"3"],[635,"3"],[636,"3"],[637,"3"],[638,"3"],[639,"3"],[640,"3"],[641,"4"],[642,"3"],[643,"4"],[644,"4"],[645,"4"],[646,"4"],[647,"2"],[648,"2"],[649,"2"],[650,"4"],[651,"4"],[652,"4"],[653,"4"],[654,"4"],[655,"4"],[656,"4"],[657,"4"],[658,"2"],[659,"2"],[660,"2"],[661,"2"],[662,"3"],[663,"3"],[664,"3"],[665,"3"],[666,"3"],[667,"3"],[668,"2"],[669,"2"],[670,"2"],[671,"2"],[672,"3"],[673,"3"],[674,"3"],[675,"3"],[676,"3"],[677,"3"],[678,"3"],[679,"3"],[680,"3"],[681,"3"],[682,"3"],[683,"3"],[684,"3"],[685,"3"],[686,"3"],[687,"3"],[688,"3"],[689,"3"],[690,"3"],[691,"3"],[692,"3"],[693,"3"],[694,"3"],[695,"3"],[696,"3"],[697,"3"],[698,"3"],[699,"3"],[700,"3"],[701,"4"],[702,"3"],[703,"3"],[704,"3"],[705,"3"],[706,"3"],[707,"3"],[708,"3"],[709,"3"],[710,"3"],[711,"3"],[712,"3"],[713,"3"],[714,"3"],[715,"3"],[716,"3"],[717,"4"],[718,"4"],[719,"4"],[720,"4"],[721,"4"],[722,"4"],[723,"4"],[724,"4"],[725,"4"],[726,"4"],[727,"4"],[728,"4"],[729,"4"],[730,"4"],[731,"4"],[732,"4"],[733,"4"],[734,"4"],[735,"4"],[736,"2"],[737,"2"],[738,"2"],[739,"2"],[740,"2"],[741,"2"],[742,"2"],[743,"2"],[744,"2"],[745,"2"],[746,"2"],[747,"2"],[748,"2"],[749,"2"],[750,"2"],[751,"2"],[752,"2"],[753,"1"],[754,"1"],[755,"1"],[756,"1"],[757,"1"],[758,"1"],[759,"4"],[760,"4"],[761,"4"],[762,"1"],[763,"1"],[764,"1"],[765,"1"],[766,"1"],[767,"1"],[768,"1"],[769,"1"],[770,"1"],[771,"1"],[772,"1"],[773,"1"],[774,"1"],[775,"1"],[776,"1"],[777,"1"],[778,"1"],[779,"1"],[780,"1"],[781,"1"],[782,"1"],[783,"1"],[784,"1"],[785,"1"],[786,"1"],[787,"1"],[788,"1"],[789,"1"],[790,"1"],[791,"1"],[792,"1"],[793,"1"],[794,"1"],[795,"1"],[796,"1"],[797,"1"],[798,"1"],[799,"1"],[800,"1"],[801,"1"],[802,"1"],[803,"1"],[804,"1"],[805,"3"],[806,"3"],[807,"3"],[808,"3"],[809,"3"],[810,"3"],[811,"4"],[812,"4"],[813,"3"],[814,"3"],[815,"3"],[816,"3"],[817,"3"],[818,"3"],[819,"3"],[820,"3"],[821,"3"],[822,"3"],[823,"3"],[824,"3"],[825,"3"],[826,"3"],[827,"3"],[828,"3"],[829,"3"],[830,"3"],[831,"3"],[832,"3"],[833,"3"],[834,"3"],[835,"3"],[836,"3"],[837,"3"],[838,"3"],[839,"3"],[840,"3"],[841,"4"],[842,"3"],[843,"3"],[844,"3"],[845,"3"],[846,"3"],[847,"3"],[848,"3"],[849,"3"],[850,"3"],[851,"3"],[852,"3"],[853,"3"],[854,"3"],[855,"3"],[856,"3"],[857,"3"],[858,"3"],[859,"3"],[860,"3"],[861,"3"],[862,"1"],[863,"2"],[864,"2"],[865,"2"],[866,"1"],[867,"1"],[868,"1"],[869,"1"],[870,"1"],[871,"1"],[872,"4"],[873,"4"],[874,"4"],[875,"4"],[876,"4"],[877,"4"],[878,"4"],[879,"4"],[880,"4"],[881,"4"],[882,"4"],[883,"4"],[884,"2"],[885,"2"],[886,"2"],[887,"2"],[888,"2"],[889,"2"],[890,"2"],[891,"2"],[892,"2"],[893,"2"],[894,"2"],[895,"2"],[896,"2"],[897,"2"],[898,"2"],[899,"2"],[900,"2"],[901,"2"],[902,"2"],[903,"2"],[904,"2"],[905,"4"],[906,"4"],[907,"4"],[908,"4"],[909,"4"],[910,"4"],[911,"4"],[912,"1"],[913,"4"],[914,"4"],[915,"1"],[916,"1"],[917,"1"],[918,"1"],[919,"1"],[920,"1"],[921,"1"],[922,"2"],[923,"2"],[924,"2"],[925,"2"],[926,"2"],[927,"2"],[928,"2"],[929,"2"],[930,"3"],[931,"3"],[932,"3"],[933,"3"],[934,"2"],[935,"2"],[936,"2"],[937,"2"],[938,"2"],[939,"2"],[940,"4"],[941,"4"],[942,"4"],[943,"4"],[944,"4"],[945,"4"],[946,"4"],[947,"4"],[948,"4"],[949,"4"],[950,"4"],[951,"4"],[952,"0"],[953,"0"],[954,"0"],[955,"0"],[956,"0"],[957,"0"],[958,"0"],[959,"0"],[960,"0"],[961,"0"],[962,"3"],[963,"3"],[964,"3"],[965,"3"],[966,"3"],[967,"3"],[968,"3"],[969,"2"],[970,"2"],[971,"2"],[972,"2"],[973,"2"],[974,"2"],[975,"2"],[976,"2"],[977,"0"],[978,"0"],[979,"0"],[980,"0"],[981,"0"],[982,"1"],[983,"0"],[984,"0"],[985,"0"],[986,"1"],[987,"1"],[988,"1"],[989,"1"],[990,"1"],[991,"1"],[992,"2"],[993,"2"],[994,"2"],[995,"2"],[996,"1"],[997,"1"],[998,"1"],[999,"1"],[1000,"1"],[1001,"1"],[1002,"1"],[1003,"1"],[1004,"1"],[1005,"1"],[1006,"1"],[1007,"1"],[1008,"1"],[1009,"1"],[1010,"1"],[1011,"1"],[1012,"1"],[1013,"1"],[1014,"1"],[1015,"1"],[1016,"1"],[1017,"1"],[1018,"1"],[1019,"1"],[1020,"1"],[1021,"1"],[1022,"1"],[1023,"1"],[1024,"1"],[1025,"0"],[1026,"0"],[1027,"0"],[1028,"0"],[1029,"2"],[1030,"2"],[1031,"2"],[1032,"2"],[1033,"2"],[1034,"2"],[1035,"2"],[1036,"1"],[1037,"1"],[1038,"1"],[1039,"1"],[1040,"1"],[1041,"1"],[1042,"1"],[1043,"1"],[1044,"1"],[1045,"3"],[1046,"1"],[1047,"3"],[1048,"3"],[1049,"3"],[1050,"3"],[1051,"3"],[1052,"3"],[1053,"3"],[1054,"3"],[1055,"3"],[1056,"3"],[1057,"3"],[1058,"3"],[1059,"3"],[1060,"3"],[1061,"3"],[1062,"3"],[1063,"3"],[1064,"3"],[1065,"3"],[1066,"3"],[1067,"3"],[1068,"3"],[1069,"3"],[1070,"1"],[1071,"1"],[1072,"1"],[1073,"1"],[1074,"1"],[1075,"1"],[1076,"1"],[1077,"1"],[1078,"1"],[1079,"1"],[1080,"1"],[1081,"1"],[1082,"1"],[1083,"1"],[1084,"1"],[1085,"3"],[1086,"3"],[1087,"3"],[1088,"3"],[1089,"3"],[1090,"3"],[1091,"2"],[1092,"2"],[1093,"2"],[1094,"2"],[1095,"2"],[1096,"1"],[1097,"1"],[1098,"1"],[1099,"1"],[1100,"4"],[1101,"4"],[1102,"4"],[1103,"4"],[1104,"3"],[1105,"3"],[1106,"3"],[1107,"3"],[1108,"2"],[1109,"2"],[1110,"2"],[1111,"2"],[1112,"2"],[1113,"2"],[1114,"2"],[1115,"4"],[1116,"4"],[1117,"4"],[1118,"4"],[1119,"4"],[1120,"4"],[1121,"4"],[1122,"4"],[1123,"4"],[1124,"4"],[1125,"4"],[1126,"4"],[1127,"4"],[1128,"2"],[1129,"3"],[1130,"3"],[1131,"3"],[1132,"3"],[1133,"3"],[1134,"3"],[1135,"3"],[1136,"3"],[1137,"3"],[1138,"4"],[1139,"4"],[1140,"4"],[1141,"4"],[1142,"4"],[1143,"4"],[1144,"4"],[1145,"4"],[1146,"4"],[1147,"4"],[1148,"4"],[1149,"4"],[1150,"4"],[1151,"4"],[1152,"4"],[1153,"4"],[1154,"4"],[1155,"4"],[1156,"4"],[1157,"4"],[1158,"1"],[1159,"1"],[1160,"1"],[1161,"1"],[1162,"1"],[1163,"1"],[1164,"1"],[1165,"1"],[1166,"1"],[1167,"1"],[1168,"1"],[1169,"1"],[1170,"1"],[1171,"1"],[1172,"1"],[1173,"1"],[1174,"1"],[1175,"1"],[1176,"1"],[1177,"1"],[1178,"1"],[1179,"1"],[1180,"1"],[1181,"1"],[1182,"1"],[1183,"1"],[1184,"1"],[1185,"1"],[1186,"1"],[1187,"2"],[1188,"2"],[1189,"2"],[1190,"2"],[1191,"2"],[1192,"2"],[1193,"3"],[1194,"3"],[1195,"3"],[1196,"3"],[1197,"3"],[1198,"3"],[1199,"3"],[1200,"3"],[1201,"3"],[1202,"3"],[1203,"3"],[1204,"3"],[1205,"3"],[1206,"3"],[1207,"3"],[1208,"3"],[1209,"3"],[1210,"3"],[1211,"3"],[1212,"3"],[1213,"2"],[1214,"2"],[1215,"2"],[1216,"2"],[1217,"2"],[1218,"2"],[1219,"2"],[1220,"2"],[1221,"2"],[1222,"2"],[1223,"1"],[1224,"1"],[1225,"1"],[1226,"1"],[1227,"1"],[1228,"1"],[1229,"1"],[1230,"4"],[1231,"4"],[1232,"4"],[1233,"4"],[1234,"4"],[1235,"4"],[1236,"4"],[1237,"1"],[1238,"4"],[1239,"4"],[1240,"2"],[1241,"2"],[1242,"2"],[1243,"2"],[1244,"2"],[1245,"2"],[1246,"2"],[1247,"2"],[1248,"2"],[1249,"2"],[1250,"2"],[1251,"2"],[1252,"2"],[1253,"2"],[1254,"2"],[1255,"2"],[1256,"2"],[1257,"2"],[1258,"2"],[1259,"2"],[1260,"1"],[1261,"1"],[1262,"1"],[1263,"1"],[1264,"1"],[1265,"1"],[1266,"1"],[1267,"1"],[1268,"1"],[1269,"1"],[1270,"1"],[1271,"1"],[1272,"1"],[1273,"1"],[1274,"1"],[1275,"1"],[1276,"1"],[1277,"1"],[1278,"1"],[1279,"2"],[1280,"2"],[1281,"2"],[1282,"2"],[1283,"2"],[1284,"2"],[1285,"2"],[1286,"2"],[1287,"2"],[1288,"2"],[1289,"2"],[1290,"2"],[1291,"2"],[1292,"2"],[1293,"1"],[1294,"1"],[1295,"1"],[1296,"1"],[1297,"1"],[1298,"1"],[1299,"1"],[1300,"1"],[1301,"1"],[1302,"1"],[1303,"1"],[1304,"1"],[1305,"1"],[1306,"1"],[1307,"1"],[1308,"1"],[1309,"1"],[1310,"1"],[1311,"1"],[1312,"1"],[1313,"1"],[1314,"1"],[1315,"1"],[1316,"2"],[1317,"1"],[1318,"1"],[1319,"2"],[1320,"1"],[1321,"1"],[1322,"1"],[1323,"1"],[1324,"1"],[1325,"1"],[1326,"1"],[1327,"1"],[1328,"1"],[1329,"1"],[1330,"1"],[1331,"1"],[1332,"4"],[1333,"4"],[1334,"4"],[1335,"4"],[1336,"4"],[1337,"4"],[1338,"4"],[1339,"3"],[1340,"3"],[1341,"3"],[1342,"3"],[1343,"3"],[1344,"3"],[1345,"3"],[1346,"3"],[1347,"3"],[1348,"3"],[1349,"4"],[1350,"4"],[1351,"4"],[1352,"4"],[1353,"4"],[1354,"4"],[1355,"4"],[1356,"4"],[1357,"4"],[1358,"4"],[1359,"4"],[1360,"4"],[1361,"4"],[1362,"4"],[1363,"4"],[1364,"4"],[1365,"4"],[1366,"4"],[1367,"4"],[1368,"4"],[1369,"4"],[1370,"4"],[1371,"4"],[1372,"4"],[1373,"4"],[1374,"4"],[1375,"2"],[1376,"2"],[1377,"2"],[1378,"2"],[1379,"2"],[1380,"2"],[1381,"2"],[1382,"2"],[1383,"2"],[1384,"2"],[1385,"2"],[1386,"2"],[1387,"1"],[1388,"1"],[1389,"2"],[1390,"1"],[1391,"2"],[1392,"1"],[1393,"1"],[1394,"1"],[1395,"1"],[1396,"2"],[1397,"1"],[1398,"1"],[1399,"1"],[1400,"1"],[1401,"1"],[1402,"1"],[1403,"1"],[1404,"1"],[1405,"1"],[1406,"1"],[1407,"1"],[1408,"1"],[1409,"1"],[1410,"1"],[1411,"1"],[1412,"1"],[1413,"1"],[1414,"1"],[1415,"1"],[1416,"1"],[1417,"1"],[1418,"1"],[1419,"1"],[1420,"1"],[1421,"1"],[1422,"1"],[1423,"1"],[1424,"1"],[1425,"1"],[1426,"1"],[1427,"1"],[1428,"1"],[1429,"1"],[1430,"1"],[1431,"1"],[1432,"1"],[1433,"1"],[1434,"1"],[1435,"1"],[1436,"1"],[1437,"1"],[1438,"1"],[1439,"1"],[1440,"1"],[1441,"1"],[1442,"1"],[1443,"3"],[1444,"1"],[1445,"3"],[1446,"3"],[1447,"1"],[1448,"1"],[1449,"1"],[1450,"1"],[1451,"1"],[1452,"1"],[1453,"1"],[1454,"1"],[1455,"1"],[1456,"1"],[1457,"1"],[1458,"1"],[1459,"1"],[1460,"1"],[1461,"1"],[1462,"1"],[1463,"1"],[1464,"1"],[1465,"1"],[1466,"1"],[1467,"1"],[1468,"1"],[1469,"1"],[1470,"1"],[1471,"1"],[1472,"1"],[1473,"1"],[1474,"1"],[1475,"1"],[1476,"1"],[1477,"1"],[1478,"1"],[1479,"1"],[1480,"1"],[1481,"1"],[1482,"1"],[1483,"1"],[1484,"1"],[1485,"1"],[1486,"1"],[1487,"1"],[1488,"1"],[1489,"1"],[1490,"1"],[1491,"1"],[1492,"1"],[1493,"1"],[1494,"1"],[1495,"1"],[1496,"1"],[1497,"1"],[1498,"1"],[1499,"1"],[1500,"1"],[1501,"1"],[1502,"1"],[1503,"1"],[1504,"1"],[1505,"1"],[1506,"1"],[1507,"2"],[1508,"2"],[1509,"2"],[1510,"2"],[1511,"2"],[1512,"2"],[1513,"2"],[1514,"2"],[1515,"2"],[1516,"2"],[1517,"2"],[1518,"2"],[1519,"2"],[1520,"2"],[1521,"2"],[1522,"2"]]}; 
    
    // for (let a = 0; a < 2; a++){
    //   vm.associatedDataMappedPerType = fix_colors(
    //         viewerInstanceTop.viewInstance.uiTemplateService.apiData.sequence, 
    //         viewerInstanceTop.viewInstance.uiTemplateService.baseStrs.get('cWW')[1],
    //         vm.associatedDataMappedPerType
    //   );
    // }

    //const AD_header='Associated Data1';
    //const ADDataArray=[[1,1],[2,2],[3,3],[4,6],[5,9]];
    //vm.AD_headers.push(AD_header)
    //mapAssociatedData(ADDataArray, AD_header, topviewer );


    // console.log('ebi_sequence', ebi_sequence);
    var largestKey = Math.max(...Object.values(struct_mapping).filter(a => typeof (a) == "number"))
    var smallestKey = Math.min(...Object.values(struct_mapping).filter(a => typeof (a) == "number"))
    if ((largestKey != stopIndex || smallestKey != startIndex) && ebi_sequence) {
      ajax('/mapSeqAlnOrig/', { fasta, ebi_sequence, startIndex: 1 }).then(origStructMappingAndData => {
        var orig_struct_mapping = origStructMappingAndData["structureMapping"];
        
        if (structMappingAndData["gapsInStruc"] && structMappingAndData["gapsInStruc"].length > 0) {
          structMappingAndData["gapsInStruc"].forEach(function (gapTup) {
            let lowMiss = Number(_.invert(struct_mapping)[gapTup[0]]);
            let topMiss = Number(_.invert(struct_mapping)[gapTup[1]]);
            if (topMiss - lowMiss > 1) {
              for (let i = lowMiss + 1; i < topMiss; i++) {
                delete orig_struct_mapping[i];
              }
            }
          })
        }
        assignColorsAndStrucMappings(vueObj, origStructMappingAndData);
      })
    } else {
      assignColorsAndStrucMappings(vueObj, structMappingAndData);

    }
  }).catch(error => {
    vueObj.topology_loaded = 'error';
    console.log(error);
    var topview = document.querySelector('#topview');
    topview.innerHTML = "Failed to load the alignment-structure mapping!<br>Try another structure."
  });
  
}

var assignColorsAndStrucMappings = function (vueObj, struct_mapping) {
  vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];
  vueObj.fasta_data = struct_mapping["amendedAln"];
  vueObj.structure_mapping = struct_mapping["structureMapping"];
  loadAlignmentViewer(vueObj.fasta_data);
  var mapped_aa_properties = mapAAProps(vueObj.aa_properties, vueObj.structure_mapping);
  if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
    if (vueObj.unmappedTWCdata) {
      mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_aa_properties);
    }
  }
  window.mapped_aa_properties = mapped_aa_properties;
  //console.log(vueObj.fasta_data);
  vm.sequence = vueObj.fasta_data.split(' ')[vm.user_uploaded_cif_flag === null || vm.user_uploaded_cif_flag ? 1 : 2];
  let sequence2 = vm.sequence.replaceAll(/-|\n/g, "");
  vm.sequence3 = sequence2.substring("sequence".length, sequence2.indexOf(">"));
  vm.sequence4 = vm.sequence3;
  //tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
  delayedMapping();
  //retry(delayedMapping, 10, 1000);
}

var delayedMapping = function () {
  //console.log("delayed mapping")
  if (typeof viewerInstanceTop === 'undefined' || viewerInstanceTop === null || vm.structFailed) {
    tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
  } else {
    if (vm.type_tree != "upload") {
      if (vm.topology_loaded) {
        try {
          var topviewer = document.getElementById("PdbeTopViewer");
          
          for (const [type, _] of Object.entries(vm.associatedDataMappedPerType_2D)) {
            let associatedDataMappedPerTypeI_2D = vm.associatedDataMappedPerType_2D[type];
            let associatedDataMappedPerTypeI_3D = vm.associatedDataMappedPerType_3D[type];

            associatedDataMappedPerTypeI_2D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            associatedDataMappedPerTypeI_3D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            const AD_header = type;
            // const ADDataArray = associatedDataMappedPerTypeI;
            vm.AD_headers.push(AD_header);
            mapAssociatedData(associatedDataMappedPerTypeI_2D, associatedDataMappedPerTypeI_3D, AD_header, topviewer);
          }
        } catch (error) {
          console.log(error)
          console.log("Mapping associated data failed")
        }
      }
      //viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties);  

      else {
        viewerInstanceTop.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties);
        setTimeout(delayedMapping, 500);
      }
    }
  }
}

function retry(fn, maxAttempts = 1, delay = 0, attempts = 0) {
  return Promise.resolve()
    .then(sleeper(delay)).then(fn)
    .catch(err => {
      if (attempts < maxAttempts) {
        return retry(fn, maxAttempts, delay, attempts + 1)
      }
      var topview = document.querySelector('#topview');
      vm.topology_loaded = 'error';
      topview.innerHTML = "EBI topology diagram is taking too long!<br>Trying to generate topology from custom mode..."
      tryCustomTopology(vm.pdbid, vm.entityID, vm.chainid[0]);
      throw err
    })
}

var tryCustomTopology = function (pdbid, entityid, chainid) {
  vm.topology_loaded = false;
  //console.log("TCT_2", pdbid, vm.sequence4); 
  //vm.getR2DT(vm.sequence4);
  //vm.URL = `r2dt/${vm.sequence3}`
  //var postTopologyURL = `r2dt/${vm.sequence3}/`
  //console.log('eid1',entityid);

  //pdbid='cust';
  if (pdbid == '' || pdbid == 'cust') {
    pdbid = 'cust'
  }

  async function getRNAChain(pdbid) {
    try {
      const returnedObject = await ajax(`full-RNA-seq/${pdbid}/${chainid}`);
      //console.log('RNA_full_sequence_cif', pdbid);
      const result = returnedObject["RNAseq"];

      //console.log('RNA_full_sequence', result);
      return result;
    } catch (error) {
      console.error(error);
      return null;
    }
  }


  async function RNAseqCall(pdbid) {
    if (vm.user_uploaded_cif_flag == true) {
      const RNA_full_sequence = await getRNAChain(pdbid);
      //console.log('RNA_full_sequence2c', RNA_full_sequence);
      return RNA_full_sequence;
    }
    //console.log('cf',vm.user_uploaded_cif_flag)
    if (vm.user_uploaded_cif_flag == false) {
      const RNA_full_sequence = vm.customFullSequence;
      //console.log('RNA_full_sequence2p', vm.customFullSequence, RNA_full_sequence);
      return RNA_full_sequence;
    }

    //return RNA_full_sequence;
    //return vm.customFullSequence;
  }

  RNAseqCall(pdbid).then(seq1 => {
    //
    if (pdbid == "cust") {
      vm.sequence_for_r2dt = seq1
    }
    else
      vm.sequence_for_r2dt = vm.sequence3;
    //console.log('RNA_full_sequence3', seq1);
    vm.getR2DT(seq1);
    vm.URL = `r2dt/${entityid}`
    var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api="true"></pdb-rna-viewer>`
    document.getElementById('topview').innerHTML = topology_viewer;
    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");

    function success(parsedResponse) {
      if (parsedResponse == "Topology Success!") {
        //console.log("Topology Success!");          
      }
      // const banName = getBanName(pdbid, 'H')
      vm.json_structures_from_r2dt = parsedResponse;
      window.viewerInstanceTop.viewInstance.uiTemplateService.render(parsedResponse.RNA_2D_json, parsedResponse.RNA_BP_json, parsedResponse.RNA_BP_json, undefined, window.viewerInstanceTop.viewInstance);
      if (vm.structFailed) {
        vm.AD_headers = [];
        var topviewer = document.getElementById("PdbeTopViewer");
        try {
          
          for (const [type, _] of Object.entries(vm.associatedDataMappedPerType_2D)) {
            let associatedDataMappedPerTypeI_2D = vm.associatedDataMappedPerType_2D[type];
            let associatedDataMappedPerTypeI_3D = vm.associatedDataMappedPerType_3D[type];

            associatedDataMappedPerTypeI_2D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            associatedDataMappedPerTypeI_3D.sort(function (entry0, entry1) {
              return entry0[0] - entry1[0];
            });
            const AD_header = type;
            // const ADDataArray = associatedDataMappedPerTypeI;
            vm.AD_headers.push(AD_header);
            mapAssociatedData(associatedDataMappedPerTypeI_2D, associatedDataMappedPerTypeI_3D, AD_header, topviewer);
          }
        } catch (error) {
          console.log("Mapping associated data failed")
        }
      }
    }
    function handle_error(error) {
      var topview = document.querySelector('#topview');
      vm.topology_loaded = 'error';
      topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
      console.log(error.responseText);
    }
    call_r2dt("POST", success, handle_error);
  }).catch(error => {
    console.error(error);
  });
}

export function call_r2dt(request_method, success = function () {/* Do nothing. */ }, handle_error = function () {/* Do nothing. */ }) {
  let keys_as_string = JSON.stringify({
    cif_mode_flag: vm.user_uploaded_cif_flag,
    cif_file_path: vm.cif_file_path
  });
  let all_lines = [
    keys_as_string,
    "\n",
    ...vm.sequence_for_r2dt.split("\n")
  ];
  let sequence_file = new File(all_lines, "my_sequence.txt", {
    type: "text/plain"
  });
  const formData = new FormData();
  formData.append("custom_seq_file", sequence_file);
  $.ajax({
    url: vm.URL,
    data: formData,
    cache: false,
    contentType: false,
    processData: false,
    method: request_method,
    dataType: 'json',
    type: request_method, // For jQuery < 1.9
    success,
    handle_error
  });
}
//vm.getR2DT(vm.pdbSeq);
//vm.URL = `r2dt/${vm.pdbSeq}/`
//vm.getR2DT(seq1);
//vm.URL = `r2dt/${seq1}/`

//var topology_viewer = `<pdb-rna-viewer id="PdbeTopViewer" pdb-id="${pdbid}" entity-id="${entityid}" chain-id="${chainid}" rv-api="true" ></pdb-rna-viewer>` 
//document.getElementById('topview').innerHTML = topology_viewer; 
//window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
// ajaxProper({
//url: postTopologyURL,
//url: vm.URL,
//type: 'POST',
//dataType: 'json'
//}).then (parsedResponse => {
//if (parsedResponse == "Topology Success!"){
//console.log("Topology Success!");          
//}
//}).catch(error => {
//var topview = document.querySelector('#topview');
//vm.topology_loaded = 'error';
//topview.innerHTML = "Failed to generate topology from the structure file!<br>Try different PDB."
//console.log(error.responseText);
//});
//}

function sleeper(ms) {
  return function (x) {
    return new Promise(resolve => setTimeout(() => resolve(x), ms));
  };
}