var registerHoverResiData = function (e, tooltipObj) {
    if (vm.type_tree == 'upload') {
        //Figure out how to do the hover in this case;
        return;
    }
    const strainQuery = '&res__poldata__strain__strain=';
    var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(e.position) + 1)}&aln=${vm.alnobj.id}${strainQuery}${vm.fastaSeqNames[Number(e.i)]}`
    let index = url.indexOf('|');
    if (index !== -1) {
        url = url.substring(0, index)
    }
    //console.log('url', url);
    ajax(url).then(alnpos_data => {
        var alnViewCanvasEle = document.querySelector("#alnDiv canvas:nth-of-type(1)");
        var alnViewLabelsEle = document.querySelector("#alnViewerLabels");
        let boundLabelBox = alnViewLabelsEle.getBoundingClientRect();
        let boundingBox = absolutePosition(alnViewCanvasEle);
        let relativeBox = alnViewCanvasEle.getBoundingClientRect();
        //console.log('alnpos_data', alnpos_data.count );
        if (alnpos_data.count != 0) {
            ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
                /*
                  if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
                    let tooltipPosition = {
                      top: mousePos.y-boundingBox.top+15 +"px",
                      left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+8 +"px",
                    };
                    //console.log('AD1',resiData["Associated data"]);
                    if (resiData["Associated data"][0] !== undefined){
                        tooltipObj.setState({
                        
                        phase: resiData["Associated data"][0][1],
                        tooltipPosition,
                      });
                    }else{
                        tooltipObj.setState({
                        //fold: 'NA',
                        phase: 'NA',
                        tooltipPosition,
                      });
                    }
                  }
                  */
                window.ajaxRun = false;
            });
        } else {
            if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right) {
                /* let tooltipPosition = {
                     top: mousePos.y-boundingBox.top+15 +"px",
                     left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+5 +"px",
                 };*/
                window.ajaxRun = false;
                /*
                tooltipObj.setState({
                    fold: 'NA',
                    phase: 'NA',
                    tooltipPosition,
                });*/
            }
        }
    }).catch(error => {
        window.ajaxRun = false;
        console.log(error);
    })
    return true;
  };
  
  function isCorrectMask(mask_range){
      window.masking_range_array = null;
      var isCorrect = false;
      if (mask_range && mask_range.match(/^(\d+-\d+;)+$/)) {
          var temp_array = mask_range.split(';').join('-').split('-');
          temp_array = temp_array.slice(0, -1)
          var i = 0;
          isCorrect = true;
          while(i < temp_array.length) {
              if(i % 2 == 0) {
                  if(Number(temp_array[i]) > Number(temp_array[i + 1])) {
                      isCorrect = false;
                  }
              }
              i = i + 1;
          }
          window.masking_range_array = temp_array;
      }
      return isCorrect;
    };
  
    // masked_array[residueNumber] === true means "keep this residue coloured".
    // Focus/Highlight keep the listed ranges; Hide keeps everything else, so it
    // passes invert=true and reuses the exact same colouring path.
    function initializeMaskedArray(mask_pairs, invert) {
        var topviewer = document.getElementById("PdbeTopViewer")
        var domainTypes = topviewer.viewInstance.uiTemplateService.domainTypes
        let longest = null;
        for (const domainType of domainTypes) {
            if (domainType.data && domainType.data.length > (longest ? longest.data.length : 0)) {
                longest = domainType;
            }
        }
        const allIndices = new Set();
        if (longest) {
            longest.data.forEach((val) => {
                if (val != undefined && val.start != undefined) {
                    allIndices.add(val.start);
                }
            });
        }

      var masked_array = [];
      for (const j of allIndices) {
          var inRange = isResidueInRanges(j, mask_pairs);
          masked_array[j] = invert ? !inRange : inRange;
      }
      return masked_array;
  };

  // Converts the flat window.masking_range_array (e.g. [1,80,91,111])
  // produced by isCorrectMask() into pairs of [start,end] structure
  // residue numbers, e.g. [[1,80],[91,111]].
  function maskRangeArrayToPairs(flatArr) {
      var pairs = [];
      if (!flatArr) { return pairs; }
      for (var i = 0; i < flatArr.length; i = i + 2) {
          pairs.push([Number(flatArr[i]), Number(flatArr[i + 1])]);
      }
      pairs.sort(function(a, b) { return a[0] - b[0]; });
      return pairs.reduce(function(merged, range) {
          var previous = merged[merged.length - 1];
          if (previous && range[0] <= previous[1] + 1) {
              previous[1] = Math.max(previous[1], range[1]);
          } else {
              merged.push(range);
          }
          return merged;
      }, []);
  };

  function isMSAViewerReady() {
      return !!(window.PVAlnViewer && window.PVAlnViewer._isMounted);
  };

  // The 3D view is never masked with viewerInstance.visual.select(): that call
  // repaints every component with a single uniform colour and paints any entry
  // without an explicit `color` white, so it cannot preserve the per-residue
  // property colouring. All three modes instead go through applyMaskColoring(),
  // which greys out the masked residues in the data that both viewers read
  // (domainTypes/selectSections_RV1 for 2D, maskedAnnotationArray for the Mol*
  // colour themes).

  // Maps the kept structure-residue ranges to alignment column ranges
  // (via vm.structure_mapping: alnPos -> structure residue number).
  // Returns the mapped ranges as [start,end] pairs of 1-indexed
  // alignment positions.
  function collapsePositions(positions) {
      var collapsed = [];
      positions.sort(function(a, b) { return a - b; }).forEach(function(pos) {
          var last = collapsed[collapsed.length - 1];
          if (last && pos == last[1] + 1) {
              last[1] = pos;
          } else {
              collapsed.push([pos, pos]);
          }
      });
      return collapsed;
  };

  function isResidueInRanges(residueNumber, mask_pairs) {
      return mask_pairs.some(function(range) {
          return residueNumber >= range[0] && residueNumber <= range[1];
      });
  };

  function getStructureResidueNumbers() {
      var topviewer = document.getElementById("PdbeTopViewer");
      var container = document.getElementById('topview');
      var residueNumbers = [];
      var seenResidues = {};
      if (topviewer && topviewer.pdbId && container) {
          container.querySelectorAll('[class*="rnaview_' + topviewer.pdbId + '_"]').forEach(function(el) {
              var tokens = (el.getAttribute('class') || '').split(/\s+/);
              var prefix = 'rnaview_' + topviewer.pdbId + '_';
              var residueToken = tokens.find(function(token) { return token.indexOf(prefix) === 0; });
              var residueNumber = residueToken ? Number(residueToken.slice(prefix.length)) : NaN;
              if (!isNaN(residueNumber) && !seenResidues[residueNumber]) {
                  seenResidues[residueNumber] = true;
                  residueNumbers.push(residueNumber);
              }
          });
      }
      if (!residueNumbers.length && vm.structure_mapping) {
          for (var alnPos in vm.structure_mapping) {
              var residueNumber = Number(vm.structure_mapping[alnPos]);
              if (!isNaN(residueNumber) && !seenResidues[residueNumber]) {
                  seenResidues[residueNumber] = true;
                  residueNumbers.push(residueNumber);
              }
          }
      }
      return residueNumbers;
  };

  function mapStructureRangesToAlignment(mask_pairs) {
      if (!vm.structure_mapping || !mask_pairs || mask_pairs.length == 0) { return []; }
      var positions = [];
      for (var alnPos in vm.structure_mapping) {
          if (isResidueInRanges(Number(vm.structure_mapping[alnPos]), mask_pairs)) {
              positions.push(Number(alnPos));
          }
      }
      return collapsePositions(positions);
  };

  // Grays out the alignment columns whose mapped structure residue falls
  // outside of the kept ranges, using the MSAViewer `features` overlay.
  function applyMaskToMSA(mask_pairs) {
      if (!isMSAViewerReady() || !vm.fastaSeqNames) { return; }
      var keepPositions = {};
      mapStructureRangesToAlignment(mask_pairs).forEach(function(range) {
          for (var p = range[0]; p <= range[1]; p++) { keepPositions[p] = true; }
      });
      var sequenceLength = window.msaOptions && window.msaOptions.sequences && window.msaOptions.sequences.length
          ? window.msaOptions.sequences[0].sequence.length : 0;
      var hiddenPositions = [];
      for (var p = 1; p <= sequenceLength; p++) {
          if (!keepPositions[p]) { hiddenPositions.push(p); }
      }
      var maskFeatures = collapsePositions(hiddenPositions).map(function(range) {
          return {
              residues: {from: range[0], to: range[1]},
              sequences: {from: 0, to: vm.fastaSeqNames.length},
              fillColor: "rgb(232,232,232)",
              borderColor: "rgb(232,232,232)",
          };
      });
      window.PVAlnViewer.setState({maskFeatures: maskFeatures});
  };

  function clearMaskFromMSA() {
      if (!isMSAViewerReady()) { return; }
      window.PVAlnViewer.setState({maskFeatures: []});
  };

  function applyHighlightLabelsTo2D(mask_pairs) {
      var topviewer = document.getElementById("PdbeTopViewer");
      var container = document.getElementById('topview');
      if (!topviewer || !topviewer.pdbId || !container) { return; }
      var pdbId = topviewer.pdbId;
      container.querySelectorAll('.nucleotide-annotation').forEach(function(annotation) {
          var tokens = (annotation.getAttribute('class') || '').split(/\s+/);
          var prefix = 'rnaview_' + pdbId + '_';
          var residueToken = tokens.find(function(token) { return token.indexOf(prefix) === 0; });
          var residueNumber = residueToken ? Number(residueToken.slice(prefix.length)) : Number(annotation.textContent);
          if (isNaN(residueNumber)) { return; }
          var isHidden = !isResidueInRanges(residueNumber, mask_pairs);
          annotation.classList.toggle('rv3-range-label-hidden', isHidden);
          if (!residueToken && annotation.tagName.toLowerCase() === 'text') {
              var tick = annotation.previousElementSibling;
              if (tick && tick.classList.contains('nucleotide-annotation')) {
                  tick.classList.toggle('rv3-range-label-hidden', isHidden);
              }
          }
      });
      var style = document.getElementById('rv3RangeVisibilityStyle');
      if (!style) {
          style = document.createElement('style');
          style.id = 'rv3RangeVisibilityStyle';
          document.head.appendChild(style);
      }
      style.textContent = '#topview .rv3-range-hidden, #topview .rv3-range-label-hidden { visibility: hidden !important; pointer-events: none !important; }';
      window.currentHideMaskPairs = mask_pairs;
      window.currentHideMaskMode = 'highlight';
      startHideMaskObserver();
  };

  // ===================== FOCUS/HIDE MODES (crop alignment and structures) =====================

  // pdb-rna-viewer keys its SVG classes by pdbId (+ chainId for base pairs)
  // + structure residue number (see pdbe-rna-viewer/src/app/uiTemplate.ts:
  // rnaview_<pdbId>_<resi> for the backbone dot + letter,
  // circle_<pdbId>_<resi> / circle-text_<pdbId>_<resi> for the property
  // circle, and rnaviewBP_<pdbId>_<chainId> with an extra "<bpType>_<start>_<end>"
  // class for each base-pair interaction line). We hide/show these directly -
  // no rebuild of that library is needed.
  function applyHideTo2D(mask_pairs, mode) {
      var topviewer = document.getElementById("PdbeTopViewer");
      if (!topviewer || !topviewer.pdbId || !topviewer.viewInstance) { return; }
      var container = document.getElementById('topview');
      if (!container) { return; }
      var pdbId = topviewer.pdbId;
      var chainId = topviewer.chainId;
      container.querySelectorAll('[class*="rnaview_' + pdbId + '_"]').forEach(function(el) {
          var tokens = (el.getAttribute('class') || '').split(/\s+/);
          var residueToken = tokens.find(function(token) { return token.indexOf('rnaview_' + pdbId + '_') === 0; });
          if (!residueToken) { return; }
          var residueNumber = Number(residueToken.slice(('rnaview_' + pdbId + '_').length));
          if (isNaN(residueNumber)) { return; }
          var isSelected = isResidueInRanges(residueNumber, mask_pairs);
          var isHidden = mode === 'focus' ? !isSelected : isSelected;
          el.classList.toggle('rv3-range-hidden', isHidden);
          container.querySelectorAll('.circle_' + pdbId + '_' + residueNumber + ', .circle-text_' + pdbId + '_' + residueNumber).forEach(function(annotation) {
              annotation.classList.toggle('rv3-range-hidden', isHidden);
          });
      });
      container.querySelectorAll('text.nucleotide-annotation').forEach(function(annotation) {
          var hasResidueClass = Array.from(annotation.classList).some(function(className) {
              return className.indexOf('rnaview_' + pdbId + '_') === 0;
          });
          if (hasResidueClass) { return; }
          var residueNumber = Number(annotation.textContent);
          if (isNaN(residueNumber)) { return; }
          var isSelected = isResidueInRanges(residueNumber, mask_pairs);
          var isHidden = mode === 'focus' ? !isSelected : isSelected;
          annotation.classList.toggle('rv3-range-hidden', isHidden);
          var tick = annotation.previousElementSibling;
          if (tick && tick.classList.contains('nucleotide-annotation')) {
              tick.classList.toggle('rv3-range-hidden', isHidden);
          }
      });
      // Hide any base-pair interaction line/label touching a hidden residue.
      container.querySelectorAll('.rnaviewBP_' + pdbId + '_' + chainId).forEach(function(el) {
          var tokens = (el.getAttribute('class') || '').split(/\s+/);
          var bpToken = tokens.find(function(t) { return /^[A-Za-z]+_\d+_\d+$/.test(t); });
          if (!bpToken) { return; }
          var parts = bpToken.split('_');
          var startSelected = isResidueInRanges(Number(parts[1]), mask_pairs);
          var endSelected = isResidueInRanges(Number(parts[2]), mask_pairs);
          var isHidden = mode === 'focus' ? !startSelected || !endSelected : startSelected || endSelected;
          el.classList.toggle('rv3-range-hidden', isHidden);
      });
      var tooltip = document.getElementById(pdbId + '-rnaTopologyTooltip');
      if (tooltip) { tooltip.style.display = 'none'; }
      var style = document.getElementById('rv3RangeVisibilityStyle');
      if (!style) {
          style = document.createElement('style');
          style.id = 'rv3RangeVisibilityStyle';
          style.textContent = '#topview .rv3-range-hidden { visibility: hidden !important; pointer-events: none !important; }';
          document.head.appendChild(style);
      }
      window.currentHideMaskPairs = mask_pairs;
      window.currentHideMaskMode = mode;
      window.mask2DHideActive = true;
      startHideMaskObserver();
  };

  function clearHideFrom2D() {
      window.mask2DHideActive = false;
      window.currentHideMaskPairs = null;
      window.currentHideMaskMode = null;
      stopHideMaskObserver();
      var container = document.getElementById('topview');
      if (!container) { return; }
      container.querySelectorAll('.rv3-range-hidden, .rv3-range-label-hidden').forEach(function(el) {
          el.classList.remove('rv3-range-hidden');
          el.classList.remove('rv3-range-label-hidden');
      });
  };

  // The topology viewer fully re-renders its SVG whenever the layout/base-pair
  // filters change (Nucleotide/Helix/Circle, "Only nested BPs", etc), which
  // would otherwise silently undo our hide-mode visibility toggling. This
  // observer re-applies it whenever the SVG's node tree changes structurally.
  // It only watches childList/subtree (never `attributes`), so our own
  // visibility-class changes never re-trigger it.
  function startHideMaskObserver() {
      stopHideMaskObserver();
      var container = document.getElementById('topview');
      if (!container || typeof MutationObserver === 'undefined') { return; }
      window.hideMaskObserver = new MutationObserver(function() {
          if (window.hideMaskObserverTimer) { clearTimeout(window.hideMaskObserverTimer); }
          window.hideMaskObserverTimer = setTimeout(function() {
              if (window.appliedMaskMode === 'highlight' && window.currentHideMaskPairs) {
                  applyHighlightLabelsTo2D(window.currentHideMaskPairs);
              } else if ((window.appliedMaskMode === 'focus' || window.appliedMaskMode === 'hide') && window.currentHideMaskPairs) {
                  applyHideTo2D(window.currentHideMaskPairs, window.appliedMaskMode);
              }
          }, 150);
      });
      window.hideMaskObserver.observe(container, {childList: true, subtree: true});
  };

  function stopHideMaskObserver() {
      if (window.hideMaskObserver) {
          window.hideMaskObserver.disconnect();
          window.hideMaskObserver = null;
      }
      if (window.hideMaskObserverTimer) {
          clearTimeout(window.hideMaskObserverTimer);
          window.hideMaskObserverTimer = null;
      }
  };

  // Crops the MSA sequences down to only the alignment columns whose
  // mapped structure residue falls inside the kept ranges, concatenated
  // together (the removed columns are not just hidden - they are no
  // longer part of the rendered sequence at all). hideColumnMap lets
  // AlignmentViewer.js translate a rendered (cropped) column index back
  // to its real alignment position for tooltips/3D highlighting.
  function applySectionToMSA(mask_pairs, mode) {
      if (!isMSAViewerReady() || !window.msaOptions || !window.msaOptions.sequences) { return; }
      if (!window.msaOptions_fullSequences) {
          window.msaOptions_fullSequences = window.msaOptions.sequences;
      }
      var selectedPositions = {};
      mapStructureRangesToAlignment(mask_pairs).forEach(function(range) {
          for (var p = range[0]; p <= range[1]; p++) { selectedPositions[p] = true; }
      });
      var sequenceLength = window.msaOptions_fullSequences.length
          ? window.msaOptions_fullSequences[0].sequence.length : 0;
      var columnMap = [];
      for (var p = 1; p <= sequenceLength; p++) {
          if ((mode === 'focus' && selectedPositions[p]) || (mode === 'hide' && !selectedPositions[p])) {
              columnMap.push(p - 1);
          }
      }
      var croppedSequences = window.msaOptions_fullSequences.map(function(seqObj) {
          var cropped = columnMap.map(function(position) { return seqObj.sequence.charAt(position); }).join('');
          return Object.assign({}, seqObj, {sequence: cropped});
      });
      window.PVAlnViewer.setState({
          hideSequences: croppedSequences,
          hideColumnMap: columnMap,
          aaPos: 0,
          seqPos: 0,
      });
  };

  function clearHideFromMSA() {
      window.msaOptions_fullSequences = null;
      if (!isMSAViewerReady()) { return; }
      window.PVAlnViewer.setState({hideSequences: null, hideColumnMap: null});
  };

  // The residue ids held by the annotation arrays are "<chainId> <residueNumber>"
  // (see Utils.generateAnnotations). Chain ids may contain digits, so the number
  // has to be taken from the last whitespace-separated token rather than by
  // stripping every non-digit out of the whole id.
  function residueNumberFromAnnotationId(id) {
      var tokens = String(id).trim().split(/\s+/);
      return Number(tokens[tokens.length - 1]);
  };

  // window.maskedAnnotationArray is the annotation data the Mol* colour themes
  // read for as long as vm.checked_filter is true. It has to be rebuilt from
  // window.masked_array whenever the underlying annotations are regenerated,
  // otherwise the 3D view keeps colouring from a stale copy.
  function rebuildMaskedAnnotationArray() {
      var annotationArray = getAnnotationArray();
      var masked = {};
      for (var mapping in annotationArray) {
          masked[mapping] = annotationArray[mapping].map(function(entry) {
              return {
                  annotation: entry.annotation,
                  ids: entry.ids.filter(function(id) {
                      return !!window.masked_array[residueNumberFromAnnotationId(id)];
                  })
              };
          });
      }
      window.maskedAnnotationArray = masked;
  };

  // Greys out the masked-out residues everywhere the two viewers read their
  // colours from. This is the path the Highlight mode has always used; Focus
  // reuses it with the listed ranges and Hide with the inverted ranges.
  // Returns the promise of the recolouring so callers can await it.
  function applyMaskColoring(mask_pairs, mode) {
      var topviewer = document.getElementById("PdbeTopViewer");
      if (!topviewer || !topviewer.viewInstance) { return Promise.resolve(); }
      var uiTemplateService = topviewer.viewInstance.uiTemplateService;
      var selectBox = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
      var selectedIndex = selectBox ? selectBox.selectedIndex : -1;

      // Rebuild the annotations first: colorResidue() greys the domain data in
      // place, so a second mask would otherwise compound onto the first one.
      uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties, window.mapped_aa_properties3D);
      if (window.custom_prop) {
          uiTemplateService.getAnnotationFromRibovision(window.custom_prop, window.custom_prop_3D);
      }

      window.masked_array = initializeMaskedArray(mask_pairs, mode === 'hide');
      for (var index = 1; index < uiTemplateService.domainTypes.length; index++) {
          colorResidue(index, window.masked_array);
      }
      var selectedData = uiTemplateService.domainTypes[selectedIndex];
      rebuildMaskedAnnotationArray();

      if (mode === 'highlight') {
          // Focus/Hide manage base-pair visibility per interaction in
          // applyHideTo2D, so only Highlight blanks them wholesale.
          var checkBoxAll = document.querySelector("#Checkbox_All");
          if (checkBoxAll) {
              checkBoxAll.checked = false;
              uiTemplateService.changeBP("All", false);
          }
      }

      const mapped_highlights = new Map()
      mapped_highlights.set('highlight',[])
      window.aaPropertyConstants.set('highlight', [0, 5]);
      window.aaColorData.set('highlight', [custom_highlight])
      for (let i = 0; i < window.masked_array.length; i++) {
          mapped_highlights.get('highlight').push([i, window.masked_array[i] ? 5 : 0])
      }
      if (window.custom_prop) {
          window.custom_prop.set("highlight", mapped_highlights.get("highlight"))
      }
      uiTemplateService.getAnnotationFromRibovision(mapped_highlights)

      var keepsCurrentProperty = selectedData && selectedData.data
          && vm.selected_property != "highlight" && vm.selected_property != 'Select data'
          && vm.selected_property != 'Clear data' && vm.selected_property;
      if (keepsCurrentProperty) {
          return recolorTopStar(selectedData.label);
      }
      if (vm.selected_property != 'highlight') {
          // The selected_property watcher recolours; assigning it here avoids
          // running recolorTopStar twice for the same property.
          vm.selected_property = 'highlight';
          return Promise.resolve();
      }
      return recolorTopStar("highlight");
  };

  // ===================== FOCUS/HIDE 3D REMOVAL (Mol* transparency) =====================
  // applyMaskColoring() only greys out masked residues - the Mol* representation
  // still renders every atom, so Focus/Hide's "removed" context stayed visible
  // as solid grey in 3D while the 2D view actually hid it (applyHideTo2D). This
  // makes the masked-out residues fully transparent in place instead, using the
  // setTransparency/clearTransparency API added to PDBeMolstarPlugin's `visual`
  // object (mirrors the existing `select`/`highlight` param shape, see
  // pdbe-molstar's src/app/index.ts). No structure reload is involved, so this
  // does not race with applyMaskColoring the way an earlier reload-based
  // implementation used to (see the comment on invalidateMolstarColorCache in
  // RV3_helpers.js).

  // Complements a sorted, merged list of [start,end] pairs (as produced by
  // maskRangeArrayToPairs) within [1, maxBound]. maxBound only needs to be
  // larger than the structure's largest residue number - anything past the
  // real end of the chain simply matches no atoms.
  function invertRanges(pairs, maxBound) {
      var inverted = [];
      var cursor = 1;
      pairs.forEach(function(range) {
          if (range[0] > cursor) { inverted.push([cursor, range[0] - 1]); }
          cursor = Math.max(cursor, range[1] + 1);
      });
      if (cursor <= maxBound) { inverted.push([cursor, maxBound]); }
      return inverted;
  };

  // Hide mode removes exactly the listed ranges; Focus removes everything
  // outside of them, so it needs the complement instead.
  function buildTransparencyRanges(mask_pairs, mode) {
      return mode === 'hide' ? mask_pairs : invertRanges(mask_pairs, 100000);
  };

  function applyRangeTransparency(mask_pairs, mode) {
      if (!window.viewerInstance || !viewerInstance.plugin || !viewerInstance.visual.setTransparency) { return Promise.resolve(); }
      var topviewer = document.getElementById("PdbeTopViewer");
      if (!topviewer) { return Promise.resolve(); }
      var data = buildTransparencyRanges(mask_pairs, mode).map(function(range) {
          return {
              auth_asym_id: topviewer.chainId,
              start_auth_residue_number: range[0],
              end_auth_residue_number: range[1]
          };
      });
      return viewerInstance.visual.setTransparency({ data: data, value: 1 }).catch(function(err) { console.log(err); });
  };

  function clearRangeTransparency() {
      if (!window.viewerInstance || !viewerInstance.plugin || !viewerInstance.visual.clearTransparency) { return Promise.resolve(); }
      return viewerInstance.visual.clearTransparency().catch(function(err) { console.log(err); });
  };

  function handleMaskingRanges(mask_range){
    vm.masking_range = mask_range;
    window.masking_range_array = null;
    if (!isCorrectMask(mask_range)) {
        vm.correct_mask = false;
        return Promise.resolve();
    }
    var mask_pairs = maskRangeArrayToPairs(window.masking_range_array);
    var mode = vm.mask_mode === 'focus' || vm.mask_mode === 'hide' ? vm.mask_mode : 'highlight';
    window.appliedMaskMode = mode;
    window.appliedMaskPairs = mask_pairs;
    window.appliedMaskRange = mask_range;
    vm.correct_mask = true;

    clearHideFrom2D();
    clearHideFromMSA();
    clearMaskFromMSA();

    window.mask3DUpdateInProgress = true;
    window.mask3DUpdatePromise = Promise.resolve(applyMaskColoring(mask_pairs, mode))
        .catch(function(err) { console.log(err); })
        .then(function() {
            window.mask3DUpdateInProgress = false;
            if (window.appliedMaskMode !== mode) { return; }
            if (mode === 'highlight') {
                applyHighlightLabelsTo2D(mask_pairs);
                applyMaskToMSA(mask_pairs);
                return clearRangeTransparency();
            } else {
                applyHideTo2D(mask_pairs, mode);
                applySectionToMSA(mask_pairs, mode);
                return applyRangeTransparency(mask_pairs, mode);
            }
        });
    return window.mask3DUpdatePromise;
  };
  function handleDomainRange(domain_range) {
      //handleFilterRange(domain_range);
      domain_array = domain_range.split(';');
      if(domain_array.length == 2) {
          vm.masking_range = null;
          vm.checked_filter = false;
          handleFilterRange(domain_range);
      } else {
          var first = domain_array[0].split('-')[0];
          var last = domain_array[domain_array.length - 2].split('-')[1];
          var full_range = first + "-" + last + ";";
          vm.checked_filter = true;
          vm.handleMaskingRanges(domain_range);
          handleFilterRange(full_range);
      }
  }
  // The old CoordinateServer (coords.litemol.org) accepted a residueRange query
  // param, but its successor, the PDBe ModelServer (CONFIG.PDB_COORDINATE_URL),
  // has no equivalent GET range query. Ranges have to be requested by POSTing an
  // explicit list of atom_site selectors and turning the returned CIF/BCIF body
  // into an object URL molstar can load.
  var previousFilterRangeBlobURL = null;
  async function fetchResidueRangeModelServer(pdbId, authAsymId, filterRange, encoding) {
      encoding = encoding || 'bcif';
      const [start, end] = filterRange.split('-').map(Number);
      const atom_site = [];
      for (let seqId = start; seqId <= end; seqId++) {
          atom_site.push({ auth_asym_id: authAsymId, auth_seq_id: seqId });
      }
      const url = `${CONFIG.PDB_COORDINATE_URL}/${pdbId}/atoms?encoding=${encoding}`;
      const response = await fetch(url, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ atom_site: atom_site })
      });
      if (!response.ok) {
          const errorText = await response.text();
          throw new Error(`ModelServer error ${response.status}: ${errorText}`);
      }
      const buffer = await response.arrayBuffer();
      const blob = new Blob([buffer], { type: 'application/octet-stream' });
      if (previousFilterRangeBlobURL) {
          URL.revokeObjectURL(previousFilterRangeBlobURL);
      }
      previousFilterRangeBlobURL = URL.createObjectURL(blob);
      return previousFilterRangeBlobURL;
  }
  async function handleFilterRange(filter_range) {
      if (filter_range.match(/^\d+-\d+;/)) {
          var filter_range = filter_range.slice(0, -1);
          const temp_array = filter_range.split('-');
          if (Number(temp_array[0]) < Number(temp_array[1])){
              window.filterRange = temp_array.join(",");
              var topviewer = document.getElementById("PdbeTopViewer");
              var selectBoxOut = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
              var selectedIndexOut = indexMatchingText(selectBoxOut.options, vm.selected_property);
              topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
              var coordURL;
              try {
                  coordURL = await fetchResidueRangeModelServer(vm.pdbid.toLowerCase(), topviewer.chainId, filter_range);
              } catch (err) {
                  console.log(err);
                  return;
              }
              viewerInstance.visual.update({
                  customData: {
                      url: coordURL,
                      format: 'cif',
                      binary:true },
                  assemblyId: '1',
                  subscribeEvents: true,
                  bgColor: {r:255,g:255,b:255},
              });
              /*
              viewerInstance.events.loadComplete.subscribe(() => { 
                  if(!vm.selected_property){return;}
                  let rangeArr = window.filterRange.split(',');
                  let selectBox = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
                  let selectedIndex = indexMatchingText(selectBox.options, vm.selected_property);
                  let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
                  if(selectSections_RV1.get(selectedData.label)) {
                      var select_sections = selectSections_RV1.get(selectedData.label).filter(resi3D  => {
                          /*if (resi3D.start_residue_number >= Number(rangeArr[0]) && resi3D.start_residue_number <= Number(rangeArr[1])){
                              return resi3D;
                          }*/
                         /* if (resi3D.residue_number >= Number(rangeArr[0]) && resi3D.residue_number <= Number(rangeArr[1])){
                              return resi3D;
                          }
                      })
                      window.viewerInstance.visual.select({
                      data: select_sections,
                      nonSelectedColor: {r:255,g:255,b:255}});
                  }
                  if (selectedIndex > 0){
                      var selectedDomain = topviewer.pluginInstance.domainTypes[selectedIndex];
                      topviewer.pluginInstance.updateTheme(selectedDomain.data);
                  }
                  selectBox.selectedIndex = selectedIndex;
              });*/
              //topviewer.pluginInstance.alreadyRan = false;
              //topviewer.pluginInstance.initPainting(window.select_sections)
              //let selectedData = topviewer.pluginInstance.domainTypes[selectedIndexOut];
              topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
              if(selectedIndexOut > 0) {
                  topviewer.pluginInstance.updateTheme(selectedData.data);
              }
              //if(vm.correct_mask){
               //   handleMaskingRanges(vm.masking_range)
              //}
          }else{
              //Swapped start end
          }
      }else{
          //Incorrect syntax
      }
  };
  
  function colorResidue(index, masked_array) {
      viewerInstanceTop.viewInstance.uiTemplateService.domainTypes[index].data.forEach(function(resiEntry){
          if (!masked_array[resiEntry.start]){
              resiEntry.color = "rgb(232,232,232)";
              resiEntry.tooltipMsg = "NaN";
          } 
      })
      selectSections_RV1.get(viewerInstanceTop.viewInstance.uiTemplateService.domainTypes[index].label).forEach(function(resiEntry){
          if (!masked_array[resiEntry.residue_number]){
              resiEntry.color = {r: 232, g: 232, b: 232};
          }
      })
  };
  function clearInputFile(f){
    if(f) {
      if(f.value){
          try{
              f.value = ''; //for IE11, latest Chrome/Firefox/Opera...
          }catch(err){ }
          if(f.value){ //for IE5 ~ IE10
              var form = document.createElement('form'),
                  parentNode = f.parentNode, ref = f.nextSibling;
              form.appendChild(f);
              form.reset();
              parentNode.insertBefore(f,ref);
          }
      }
    }
  }
  
  function cleanCustomMap(checked_customMap){
      if (vm.uploadSession){return;}
      var topviewer = document.getElementById("PdbeTopViewer");
      //console.log("topviewer_RV310", topviewer.viewInstance.uiTemplateService.domainTypes);
      if (!topviewer || !topviewer.viewInstance.uiTemplateService.domainTypess){
          if (checked_customMap){return;}
          var sliceAvailProp = Array.prototype.slice.call(vm.available_properties).filter(availProp => {
              return vm.custom_headers.includes(availProp.Name)
          })
          const setSlice = new Set(sliceAvailProp.map(a=>{return a.Name}));
          const newArray = vm.available_properties.filter(obj => !setSlice.has(obj.Name));
          vm.available_properties = newArray;
          return;
      }
      //console.log("topviewer_RV320", topviewer.viewInstance.targetEle);
      //var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
      //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
      var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
      topviewer.viewInstance.uiTemplateService.domainTypes = topviewer.viewInstance.uiTemplateService.domainTypes.filter(obj => {
          return !vm.custom_headers.includes(obj.label)
      })
      
      var sliceChildren = Array.prototype.slice.call(selectBoxEle.childNodes).filter(optionsNode => {
          return vm.custom_headers.includes(optionsNode.label)
      })
      
      sliceChildren.forEach(function(){
          selectBoxEle.removeChild(selectBoxEle.childNodes[selectBoxEle.options.length-1]);
          vm.available_properties.splice(-1,1)
      })
  
      if (checked_customMap){return;}
      window.coilsOutOfCustom = null;
      window.custom_prop = null;
      vm.csv_data = null;
      vm.custom_headers = [];
  };
  function handleCustomMappingData(){
    const readFile = function (fileInput) {
        var reader = new FileReader();
        reader.onload = function () {
            vm.csv_data = reader.result.replace("\u00EF\u00BB\u00BF", '');
        };
        reader.readAsBinaryString(fileInput);
    };
    readFile(vm.$refs.custom_csv_file.files[0]);
};

var displayMappingDataByIndex = function (topviewer, selectedIndex) {
    //var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
    var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    //topviewer.pluginInstance.resetTheme();
    //topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data);
    window.viewerInstance.visual.select({
        data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label),
        nonSelectedColor: { r: 255, g: 255, b: 255 }
    });
    selectBoxEle.selectedIndex = selectedIndex;
    vm.selected_property = topviewer.pluginInstance.domainTypes[selectedIndex].label;
}

// Function to generate dynamic SVG color scale with improved annotations
function generateColorScaleSVG(colormap, minVal, maxVal, dataName) {
    const width = 200; 
    const height = 500;
    const margin = 20; 
    const scaleHeight = height - 2 * margin;
    const barWidth = 30; // Width of the color bar itself
    const offset = 5; // Offset for the text label
    const barCenterX = width / 2; // Center the bar horizontally
    const barLeftX = barCenterX - barWidth / 2; // Left edge of centered bar
    const fontSize = 12; // Font size for annotations

    // Generate color stops for the gradient
    let colorStops = '';
    if (colormap && colormap.length > 0) {
        colormap.forEach((stop, index) => {
            const offset = (index / (colormap.length - 1)) * 100;
            const color = stop[1];
            const r = Math.round(color[0] * 255);
            const g = Math.round(color[1] * 255);
            const b = Math.round(color[2] * 255);
            colorStops += `<stop offset="${offset}%" style="stop-color:rgb(${r},${g},${b});stop-opacity:1" />`;
        });
    }

    // Calculate annotation values
    const range = maxVal - minVal;
    const zeroPercent = minVal;
    const fiftyPercent = minVal + (range * 0.5);
    const hundredPercent = maxVal;

    // Create SVG with vertical gradient and external annotations
    const svg = `
        <svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
            <defs>
                <linearGradient id="gradient-${dataName.replace(/[^a-zA-Z0-9]/g, '')}" x1="0%" y1="100%" x2="0%" y2="0%">
                    ${colorStops}
                </linearGradient>
            </defs>

            <!-- Color bar -->
            <rect x="${barLeftX}" y="${margin}" width="${barWidth}" height="${scaleHeight}"
                  fill="url(#gradient-${dataName.replace(/[^a-zA-Z0-9]/g, '')})"
                  stroke="#333" stroke-width="1"/>

            <!-- Data name label -->
            <text x="${width/2}" y="${margin + scaleHeight + 15}" text-anchor="middle" font-family="Arial" font-size="${fontSize}" fill="#333" font-weight="bold">
                ${dataName}
            </text>

            <!-- External annotations -->
            <!-- 100% (Top) - Right side -->
            <line x1="${barLeftX + barWidth + 3}" y1="${margin + offset}" x2="${barLeftX + barWidth + 8}" y2="${margin + offset}"
                  stroke="#333" stroke-width="1"/>
            <text x="${barLeftX + barWidth + 12}" y="${margin + 4 + offset}" font-family="Arial" font-size="${fontSize}" fill="#333">
                <tspan x="${barLeftX + barWidth + 12}" dy="0">100%</tspan>
                <tspan x="${barLeftX + barWidth + 12}" dy="1.2em">(${Math.round(hundredPercent)})</tspan>
            </text>

            <!-- 50% (Middle) - Left side -->
            <line x1="${barLeftX - 3}" y1="${margin + scaleHeight/2}" x2="${barLeftX - 8}" y2="${margin + scaleHeight/2}"
                  stroke="#333" stroke-width="1"/>
            <text x="${barLeftX - 12}" y="${margin + scaleHeight/2 + 4}" text-anchor="end" font-family="Arial" font-size="${fontSize}" fill="#333">
                <tspan x="${barLeftX - 12}" dy="0">50%</tspan>
                <tspan x="${barLeftX - 12}" dy="1.2em">(${Math.round(fiftyPercent)})</tspan>
            </text>

            <!-- 0% (Bottom) - Right side -->
            <line x1="${barLeftX + barWidth + 3}" y1="${margin + scaleHeight - offset}" x2="${barLeftX + barWidth + 8}" y2="${margin + scaleHeight - offset}"
                  stroke="#333" stroke-width="1"/>
            <text x="${barLeftX + barWidth + 12}" y="${margin + scaleHeight - offset - 10}" font-family="Arial" font-size="${fontSize}" fill="#333">
                <tspan x="${barLeftX + barWidth + 12}" dy="0">(${Math.round(zeroPercent)})</tspan>
                <tspan x="${barLeftX + barWidth + 12}" dy="1.2em">0%</tspan>
            </text>
        </svg>
    `;

    // Convert SVG to data URL
    const svgBlob = new Blob([svg], {type: 'image/svg+xml'});
    return URL.createObjectURL(svgBlob);
}

var mapCustomMappingData = function (custom_data, custom_data_name, topviewer, selectedColormap) { //custom_data3D, 

    //var selectBoxEle = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    if(vm.cifPdbMode == null) {
        let mapping2D_3D = {};

        for (let [k,v] of Object.entries(vm.st_mapping2D)){
            if (Object.keys(vm.st_mapping2D).includes(k)){
                mapping2D_3D[v] = vm.st_mapping3D[k];
            }
        }

        let custom_data3D = [];
        for (let [k, [u, v]] of Object.entries(custom_data)){
            custom_data3D.push([mapping2D_3D[u], v]);
            
        }

        let vals = custom_data.map(function (v) { return v[1] });
        let indexes = custom_data.map(function (v) { return v[0] });
        window.aaColorData.set(custom_data_name, [selectedColormap || viridis]);
        window.aaPropertyConstants.set(custom_data_name, [Math.min(...vals), Math.max(...vals)]);
        //let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
        //window.coilsOutOfCustom = coilsOutOfCustom;
        //console.log('CD1', custom_data_name, custom_data );
        let custom_prop = new Map();
        let custom_prop3D = new Map();
        
        custom_prop.set(custom_data_name, custom_data);
        custom_prop3D.set(custom_data_name, custom_data3D);
        if (window.custom_prop) {
            window.custom_prop.set(custom_data_name, custom_data)
        } else {
            window.custom_prop = custom_prop;
        }
        if (window.custom_prop_3D) {
            window.custom_prop.set(custom_data_name, custom_data3D)
        } else {
            window.custom_prop_3D = custom_prop3D;
        }
        topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(custom_prop, custom_prop3D);
        //var custom_option = document.createElement("option");
        //custom_option.setAttribute("value", selectBoxEle.options.length);
        //custom_option.appendChild(document.createTextNode(custom_data_name));
        //selectBoxEle.appendChild(custom_option);
        if (!vm.available_properties.some(prop => prop.Name === custom_data_name)) {
            const minVal = Math.min(...vals);
            const maxVal = Math.max(...vals);
            const colorScaleURL = generateColorScaleSVG(selectedColormap || viridis, minVal, maxVal, custom_data_name);
            vm.available_properties.push({ Name: custom_data_name, url: colorScaleURL })
        }
        if (vm.correct_mask) {
            var j = topviewer.viewInstance.uiTemplateService.domainTypes.length - 1;
            colorResidue(j, window.masked_array);
            rebuildMaskedAnnotationArray();
        }
    } else {
            //var selectBoxEle = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
            //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
            //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');

            let mapping2D_3D = {};

            for (let [k,v] of Object.entries(vm.st_mapping2D)){
                if (Object.keys(vm.st_mapping2D).includes(k)){
                    mapping2D_3D[v] = vm.st_mapping3D[k];
                }
            }

            let custom_data3D = [];
            for (let [k, [u, v]] of Object.entries(custom_data)){
                custom_data3D.push([mapping2D_3D[u], v]);
                
            }
        
            let vals = custom_data.map(function (v) { return v[1] });
            let indexes = custom_data.map(function (v) { return v[0] });
            window.aaColorData.set(custom_data_name, [selectedColormap || viridis]);
            window.aaPropertyConstants.set(custom_data_name, [Math.min(...vals), Math.max(...vals)]);
            //let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
            //window.coilsOutOfCustom = coilsOutOfCustom;
            //console.log('CD1', custom_data_name, custom_data );
            let custom_prop = new Map();
            let custom_prop3D = new Map();
            
            custom_prop.set(custom_data_name, custom_data);
            custom_prop3D.set(custom_data_name, custom_data3D);
            
            custom_prop.set(custom_data_name, custom_data);
            if (window.custom_prop) {
                window.custom_prop.set(custom_data_name, custom_data)
            } else {
                window.custom_prop = custom_prop;
            }
            if (window.custom_prop_3D) {
                window.custom_prop.set(custom_data_name, custom_data3D)
            } else {
                window.custom_prop_3D = custom_prop3D;
            }
            topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(custom_prop, custom_prop3D);
            //var custom_option = document.createElement("option");
            //custom_option.setAttribute("value", selectBoxEle.options.length);
            //custom_option.appendChild(document.createTextNode(custom_data_name));
            //selectBoxEle.appendChild(custom_option);
            if (!vm.available_properties.some(prop => prop.Name === custom_data_name)) {
                const minVal = Math.min(...vals);
                const maxVal = Math.max(...vals);
                const colorScaleURL = generateColorScaleSVG(selectedColormap || viridis, minVal, maxVal, custom_data_name);
                vm.available_properties.push({ Name: custom_data_name, url: colorScaleURL })
            }
            if (vm.correct_mask) {
                var j = topviewer.viewInstance.uiTemplateService.domainTypes.length - 1;
                colorResidue(j, window.masked_array);
                rebuildMaskedAnnotationArray();
            }
    }
}

var mapAssociatedData = function (associated_data_2D, associated_data_3D, associated_data_name, topviewer) {

    //var selectBoxEle = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    // console.log("mapAssociatedData 3D data 100:", associated_data_3D[100]);
    let vals = associated_data_3D.map(function (v) { return v[1] });
    // let indexes = associated_data_3D.map(function (v) { return v[0] });
    //window.aaColorData.set(associated_data_name, [viridis]);
    window.aaColorData.set(associated_data_name, [rainbow]);
    window.aaPropertyConstants.set(associated_data_name, [Math.min(...vals), Math.max(...vals)]);
    //let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
    //window.coilsOutOfCustom = coilsOutOfCustom;
    // console.log('AD1', associated_data_name, associated_data);
    // if (!window.associated_prop) {
    //     window.associated_prop = new Map();
    // }
    // var associated_prop = window.associated_prop;
    var associated_prop_2D = new Map();
    var associated_prop_3D = new Map();

    associated_prop_2D.set(associated_data_name, associated_data_2D);
    associated_prop_3D.set(associated_data_name, associated_data_3D);
    // console.log('associated_data_2D', JSON.stringify(associated_data_2D));
    if (window.custom_prop) {
        window.custom_prop.set(associated_data_name, associated_data_2D)
    } else {
        window.custom_prop = associated_prop_2D;
    }
    if (window.custom_prop_3D) {
        window.custom_prop_3D.set(associated_data_name, associated_data_3D)
    } else {
        window.custom_prop_3D = associated_prop_3D;
    }
    // console.log('associated_data_3D', JSON.stringify(associated_data_3D));
    // if (window.custom_prop) {
    //     window.custom_prop.set(associated_data_name, associated_data_3D)
    // } else {
    //     window.custom_prop = associated_prop_3D;
    // }
    topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(associated_prop_2D, associated_prop_3D);
    //var custom_option = document.createElement("option");
    //custom_option.setAttribute("value", selectBoxEle.options.length);
    //custom_option.appendChild(document.createTextNode(custom_data_name));
    //selectBoxEle.appendChild(custom_option);
    //if (!vm.available_properties.some(prop => prop.Name === associated_data_name)){
    //    vm.available_properties.push({Name:associated_data_name, url:"static/alignments/svg/Custom.svg"})
    //}
    if (vm.correct_mask) {
        var j = topviewer.viewInstance.uiTemplateService.domainTypes.length - 1;
        colorResidue(j, window.masked_array);
        rebuildMaskedAnnotationArray();
    }
}


var mapHelixData = function (helix_data, helix_data_name, topviewer) {

    //var selectBoxEle = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.menuSelectbox');
    //var selectBoxEle = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');

    let vals = helix_data.map(function (v) { return v[1] });
    let indexes = helix_data.map(function (v) { return v[0] });
    window.aaColorData.set(helix_data_name, [viridis]);
    window.aaPropertyConstants.set(helix_data_name, [Math.min(...vals), Math.max(...vals)]);
    //let coilsOutOfCustom = vm.coil_residues.filter(value => !indexes.includes(value));
    //window.coilsOutOfCustom = coilsOutOfCustom;
    // console.log('AD1', associated_data_name, associated_data);
    // if (!window.associated_prop) {
    //     window.associated_prop = new Map();
    // }
    // var associated_prop = window.associated_prop;
    var helix_prop = new Map();
    helix_prop.set(helix_data_name, helix_data);
    if (window.helix_prop) {
        window.helix_prop.set(helix_data_name, helix_data)
    } else {
        window.helix_prop = helix_prop;
    }
    topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(helix_prop);
    //var custom_option = document.createElement("option");
    //custom_option.setAttribute("value", selectBoxEle.options.length);
    //custom_option.appendChild(document.createTextNode(custom_data_name));
    //selectBoxEle.appendChild(custom_option);
    if (!vm.available_properties.some(prop => prop.Name === helix_data_name)) {
        vm.available_properties.push({ Name: helix_data_name, url: "static/alignments/svg/Custom.svg" })
    }
    if (vm.correct_mask) {
        var j = topviewer.viewInstance.uiTemplateService.domainTypes.length - 1;
        colorResidue(j, window.masked_array);
        rebuildMaskedAnnotationArray();
    }
}
var getExampleFile = function (url, name) {
    $.ajax({
        url: url,
        type: 'GET',
        dataType: "text",
        success: function (data) {
            let anchor = document.createElement('a');
            anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(data);
            anchor.target = '_blank';
            anchor.download = name;
            anchor.click();
        },
    })
};

function cleanFilter(checked_filter, masking_range) {
    if (checked_filter) {
        // The Mol* colour themes read window.maskedAnnotationArray for as long
        // as vm.checked_filter is true, which is from the moment the box is
        // ticked - before any range has been applied. Seed it with the
        // unfiltered annotations so colouring keeps working until Apply.
        window.maskedAnnotationArray = getAnnotationArray();
        return;
    }
    if (masking_range == null) { return; }
    window.masked_array = [];
    window.maskedAnnotationArray = null;
    vm.masking_range = null;
    vm.correct_mask = null;
    window.appliedMaskMode = null;
    window.appliedMaskPairs = null;
    window.appliedMaskRange = null;
    clearMaskFromMSA();
    clearHideFrom2D();
    clearHideFromMSA();
    clearRangeTransparency();
    var topviewer = document.getElementById("PdbeTopViewer");
    topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(mapped_aa_properties, window.mapped_aa_properties3D);
    if(window.custom_prop) {
        topviewer.viewInstance.uiTemplateService.getAnnotationFromRibovision(window.custom_prop, window.custom_prop_3D);
    }
    var domainTypes = topviewer.viewInstance.uiTemplateService.domainTypes;
    var indexToRemove = domainTypes.findIndex(obj => obj.label === 'highlight');
    if (indexToRemove !== -1) {
        domainTypes.splice(indexToRemove, 1);
    }
    var selectElement = topviewer.viewInstance.targetEle.querySelector('.mappingSelectbox');
    var optionToRemove;
    if (selectElement) {
        Array.from(selectElement.options).forEach(function(option) {
            if (option.label === 'highlight') {
                optionToRemove = option;
            }
        });
    }
    if (optionToRemove) {
        selectElement.removeChild(optionToRemove);
    }
    vm.selected_property = "Clear data"
  };
  function cleanSelection(checked_selection, filter_range){
    if (checked_selection || filter_range == null || !vm.pdbid){return;}
    var selectBox = viewerInstanceTop.pluginInstance.targetEle.querySelector('.menuSelectbox');
    var newIndex = indexMatchingText(selectBox.options, vm.selected_property);
    vm.filter_range = null;
    window.filterRange = "-10000,10000";
    viewerInstanceTop.pluginInstance.alreadyRan = false;
    viewerInstanceTop.pluginInstance.initPainting();
    var coordURL = `https://coords.litemol.org/${vm.pdbid.toLowerCase()}/chains?entityId=${viewerInstanceTop.entityId}&authAsymId=${viewerInstanceTop.chainId}&encoding=bcif`;
    //var coordURL = `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/chains?entityId=${topviewer.entityId}&encoding=bcif`;
    viewerInstance.visual.update({
        customData: {
            url: coordURL,
            format: 'cif',
            binary: true
        },
        assemblyId: '1',
        subscribeEvents: true,
        bgColor: {r:255,g:255,b:255},
      }).finally(response => {
          viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
          if(window.custom_prop) {
              viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
          }
          if(newIndex > 0) {
              viewerInstanceTop.pluginInstance.updateTheme(viewerInstanceTop.pluginInstance.domainTypes[newIndex].data); 
          }
          if (response){
              window.viewerInstance.visual.select({data: selectSections_RV1.get(vm.selected_property), nonSelectedColor: {r:255,g:255,b:255}});
          }
          if(vm.correct_mask) {
              handleMaskingRanges(vm.masking_range)
          }
            handlePropensities(vm.checked_propensities);
      });
  };
  
  var populatePDBs = function (alndata){
      if (alndata != null){
          let alnPolurl = `/desire-api/polymers/?alns_of_polymer=${alndata.id}`
         
          ajax(alnPolurl).then(polymersForAln => {
              let trueNom = polymersForAln.results[0].nomgd.split('/')[5];
              var polNames = polymersForAln.results.map(entry => entry.genedescription.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,''));
              let url = `/desire-api/old-nomenclatures/?n_b_y_h_a=BAN&nn_fk=${trueNom}`;
          
              ajax(url).then(oldnomData => {
                  oldnomData.count=1;
                  if (oldnomData.count == 0){return;}
                  //let oldName = oldnomData.results[0].old_name.replace(/^(.{2})(0)/,"$1")
                  //let riboXYZurl = `https://api.ribosome.xyz/neo4j/gmo_nom_class/?banName=${oldName}&format=json`
                  //Use vm.alnobj.text for name of alignment
                  rna_class = []
                  if(vm.alnobj.text == "5S") {
                    rna_class = ['5S']
                  } else if (vm.alnobj.text == '5.8S') {
                    rna_class = ['5.8S']
                  } else if (vm.alnobj.text == 'LSUa' || vm.alnobj.text == 'LSUb') {
                    rna_class = ['23S']
                  } else if (vm.alnobj.text == '28S') {
                    //Should we include 25S for this?
                    rna_class = ['25S', '28S']
                  } else if(vm.alnobj.text == 'SSU') {
                    rna_class = ['16S']
                  }
                  rna_class.forEach(rnaClass => {
                //   let riboXYZurl = `https://api.ribosome.xyz/neo4j/get_rna_class/?rna_class=${rnaClass}rRNA&format=json`
                  let riboXYZurl = `https://api.ribosome.xyz/polymers/polynucleotide?rna_class=${rnaClass}rRNA&format=json`
                 
                  ajax(riboXYZurl).then(data => {
                      var pdb_entries = []
                 
                      data.forEach(function(entry){
                          let pdb_text = `${entry.parent_rcsb_id} ${entry.src_organism_names[0].slice(0,39)}`
                          
                          //let pdbxDescription = entry.protein.rcsb_pdbx_description.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'')
                          //if (polNames.includes(pdbxDescription)){
                          //pdb_entries.push({id: entry.parent_rcsb_id, name: `${entry.parent_rcsb_id} ${entry.src_organism_names[0].slice(0,39)}`});
                          pdb_entries.push({id: entry.parent_rcsb_id.toLowerCase(), name: pdb_text});
                              //pdb_entries.push({id: entry.parent_rcsb_id})
                              
  
                          //}
                      });
                      
                      if (pdb_entries.length == 0){return;}
                      vm.pdbs.push(...pdb_entries.sort((a, b) => (a.id > b.id) ? 1 : -1));
                      const pdbSet = new Set();
                      vm.pdbs = vm.pdbs.filter(entry => !pdbSet.has(entry.id) && pdbSet.add(entry.id));
                  }).catch(error => {
                      console.log(error);
                  })
              })}).catch(error => {
                  console.log(error);
              })
          }).catch(error => {
                  console.log(error);
          })
      }
  }
  
  
  /*var populatePDBs = function (alndata){
      if (alndata != null){
          let alnPolurl = `/desire-api/polymers/?alns_of_polymer=${alndata.id}`
          ajax(alnPolurl).then(polymersForAln => {
              let trueNom = polymersForAln.results[0].nomgd.split('/')[5];
              var polNames = polymersForAln.results.map(entry => entry.genedescription.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,''));
              let url = `/desire-api/old-nomenclatures/?n_b_y_h_a=BAN&nn_fk=${trueNom}`;
              ajax(url).then(oldnomData => {
                  if (oldnomData.count == 0){return;}
                  let oldName = oldnomData.results[0].old_name.replace(/^(.{2})(0)/,"$1")
                  let riboXYZurl = `https://api.ribosome.xyz/neo4j/gmo_nom_class/?banName=${oldName}&format=json`
                  ajax(riboXYZurl).then(data => {
                      var pdb_entries = []
                      data.forEach(function(entry){
                          let pdb_text = `${entry.parent} ${entry.orgname[0].slice(0,39)}`
                          let pdbxDescription = entry.protein.rcsb_pdbx_description.trim().replace(/-[\w]{1}$/,'').replace(/ubiquitin/ig,'')
                          if (polNames.includes(pdbxDescription)){
                              pdb_entries.push({id: entry.parent.toLowerCase(), name:pdb_text})
                          }
                      });
                      if (pdb_entries.length == 0){return;}
                      vm.pdbs.push(...pdb_entries.sort((a, b) => (a.id > b.id) ? 1 : -1));
                  }).catch(error => {
                      console.log(error);
                  })
              }).catch(error => {
                  console.log(error);
              })
          }).catch(error => {
                  console.log(error);
          })
      }
  }
  */
  var customFilter = function (object, result, key, value){
      if(object.hasOwnProperty(key) && object[key] == value)
          result.push(object);
      for(var i=0; i<Object.keys(object).length; i++){
          let nextObj = object[Object.keys(object)[i]];
          if(typeof nextObj == "object" && nextObj != null){
              customFilter(nextObj, result, key, value);
          }
      }
  }
  
  function parseConsecutiveIndices(structureTypeString, structureList, indicesList) {
      let structureIndex = 0;
      if (indicesList.length == 0) {
          return;
      }
      let previousIndex = indicesList[0];
      let currentStructure = [previousIndex];
      for (let i = 1; i < indicesList.length; i++) {
          let currentIndex = indicesList[i];
          if (currentIndex == previousIndex + 1) {
              currentStructure.push(currentIndex);
          } else {
              let structureObject = {};
              structureObject.text = structureTypeString + " #" + structureIndex;
              structureObject.value = structureIndex;
              structureObject.indices = currentStructure;
              structureList.push(structureObject);
              structureIndex++;
              currentStructure = [currentIndex];
          }
          previousIndex = currentIndex;
      }
  }
  
  function getPropensities(property) {
      let indices = null;
      if (vm.structure_mapping && property && property!=0) {
          var sequence_indices = property.indices;
          let alignment_indices = []
          let inverse_structure_mapping = _.invert(vm.structure_mapping);
          for (var sequence_index of sequence_indices) {
              if (inverse_structure_mapping[sequence_index]){
                  alignment_indices.push(inverse_structure_mapping[sequence_index])
              }
          }
          indices = alignment_indices.join(',')
      } else {
          indices = '';
      }
      vm.propensity_indices = indices
      vm.fasta_data
  }
  
  function handlePropensityIndicesOnTruncatedStructure(indices, startTrunc, endTrunc){
      var newIndices = '';
      var invertedMap = _.invert(vm.structure_mapping);
      if (!indices){
          tempIndices = [];
          vm.all_residues.forEach(function(resi){
              tempIndices.push(invertedMap[resi]);
          })
          indices = tempIndices.join(',');
      }
      indices.split(',').forEach(function(entry){
          if (invertedMap[startTrunc] <= Number(entry) &&  Number(entry) <= invertedMap[endTrunc]){
              newIndices+=`${entry},`;
          }
      })
      return newIndices.slice(0, -1);
  }
  
  function handlePropensities(checked_propensities) {
      if (checked_propensities) {
          var title = 'Amino Acid Frequencies'
          if (vm.type_tree == "orth"){
              var title = `${vm.alnobj.text} ${title}`;
          }
          if (vm.property){
              title += ` for ${vm.property.text}`
          }
          let indices = vm.propensity_indices;
          if (vm.selected_domain.length > 0){
              title += `<br>of ECOD domain ${vm.selected_domain[0].name}`;
              var domainIndices = '';
              var invertedMap = _.invert(vm.structure_mapping);
              vm.selected_domain[0].range.split(';').forEach(function(singleRange){
                  if (singleRange == ''){return;}
                  var startDomain = Number(invertedMap[Number(singleRange.split('-')[0])]);
                  var endDomain = Number(invertedMap[Number(singleRange.split('-')[1])]);
                  for (let i = startDomain; i <= endDomain; i++){
                      domainIndices += `${i},`;
                  }
              });
              if (!indices || indices == ''){
                  indices = domainIndices.slice(0, -1);
              } else {
                  let tempIndices = _.intersection(indices.split(','),domainIndices.slice(0, -1).split(','));
                  indices = tempIndices.join(',')
              }
          }
          if (vm.filter_range){
              var startRange = Number(vm.filter_range.split('-')[0]);
              var endRange = Number(vm.filter_range.split('-')[1].slice(0, -1));
              title += `<br>between positions ${startRange} and ${endRange}`;
              indices = handlePropensityIndicesOnTruncatedStructure(indices, startRange, endRange);
          }
          var full = ['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W'];
          let customFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
          if (indices) {
              ajax("/propensity-data-custom/", {indices, customFasta}).then(data => {
                  storeFrequencyData(data['amino acid'],title);
                  build_propensity_graph(data['amino acid'], full, title, 'total');
              });
          } else {
              ajax("/propensity-data-custom/", {customFasta}).then(data => {
                  storeFrequencyData(data['amino acid'],title);
                  build_propensity_graph(data['amino acid'], full, title, 'total');
              });
          }
      }
  }
  
  var listSecondaryStructures = function() {
      vm.substructures = []
      let coilObject = {
          text: "Coil residues",
          value: 0,
          indices: vm.coil_residues,
      },
      strandObject = {
          text: "Strand residues",
          value: 1,
          indices: vm.strand_residues,
      },
      helixObject = {
          text: "Helix residues",
          value: 2,
          indices: vm.helix_residues
      };
      Array.prototype.push.apply(vm.substructures, [coilObject, strandObject, helixObject])
  }
  
  var build_propensity_graph = function (data, amino_acids, title, div) {
      // reformat data into trace format
      var data2 = {};
      var hover = {};
      // initialize arrays for each amino acid in the JavaScript object
      for (aa of amino_acids) {
          data2[aa] = [];
          hover[aa] = [];
      }
      // for a species, add the propensities for each amino acid to thei
      for (species in data) {
          for (aa of amino_acids) {
              data2[aa].push(data[species][aa])
              hover[aa].push(data[species]['name'].replace(/_/g, ' '))
          }
      };
      // build traces
      traces = []
      for (aa of amino_acids) {
          // styling - can we make this a stripplot?
          var newtrace = {
              y: data2[aa],
              type: 'box',
              boxpoints: 'all',
              jitter: 0.2,
              pointpos: 0,
              name: aa,
              text: hover[aa]};
          traces.push(newtrace)
      };
      
      var layout = {
          title: title,
          xaxis: {title: 'Amino Acid'},
          yaxis: {title: 'Frequency'},
          hovermode: 'closest',
          hoveron: 'points',
      };
      Plotly.newPlot(div, traces, layout);
      var myPlot = document.getElementById(div);
      myPlot.scrollIntoView();
      myPlot.on('plotly_hover', function(data){
          data.points.map(function(d){
              let seqname = d.text.split(' ').slice(1,).join(' ')
              PVAlnViewer.highlightRegion({
                  sequences: {from: vm.fastaSeqNames.indexOf(seqname), to: vm.fastaSeqNames.indexOf(seqname)},
                  residues: {from: 0, to: vm.fasta_data.split('>')[1].split('\n')[1].length}
              })
          });
          //Plotly.Fx.hover(div,[
          //    { curveNumber:0, pointNumber:1 },
          //    { curveNumber:1, pointNumber:2 },
          //    { curveNumber:2, pointNumber:3 },
          //]);
      }).on('plotly_unhover', function(data){
          //
      });
  }
  
  var storeFrequencyData = function (data, title){
      csv_string = '';
      aaNames = [];
      aaList = [];
      once = true;
      for (var key of Object.keys(data)) {
          var entry = data[key];
          if (once){
              once = false;
              aaNames = Object.keys(entry).slice(1);
              csv_string += 'Species\\AA,'+aaNames.join(',')+'\n';
          }
          csv_string += `${entry["name"]},`
          aaNames.forEach(function(aa){
              csv_string += `${entry[aa]},`
          })
          csv_string.slice(0, -1);
          csv_string += '\n';
      }
      vm.freqCSV = `${title.replace('<br>',' ')}\n${csv_string}`;
  }
