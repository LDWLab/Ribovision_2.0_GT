const { map } = require("lodash");

// function extractBonds(st) {
//     let regex = /.*[ACGU]([0-9]+).*[ACGU]([0-9]+);/;
//     let matches = st.match(regex);
//     return matches ? [parseInt(matches[1]), parseInt(matches[2])] : null;
// }
function euclideanDistance(x1, y1, x2, y2) {
    return Math.sqrt(Math.pow(parseFloat(x2) - parseFloat(x1), 2) + Math.pow(parseFloat(y2) - parseFloat(y1), 2));
}

function extractBonds(st) {
    let basePairRegex = /.*[ACGU]([0-9]+).*[ACGU]([0-9]+);/;
    let basePairMatch = st.match(basePairRegex);
    let basePair = basePairMatch ? basePairMatch : null;
    
    let coordsRegex = /(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)/g;
    let coordsMatches = st.matchAll(coordsRegex);
    let coords = [];
    for (const match of coordsMatches) {
        coords.push([match[1], match[2], match[3], match[4]]);
    }
    
    let dist = 0;
    if (coords.length > 0) {
        dist = euclideanDistance(...coords[0]);
    }
    
    return [basePair[1], basePair[2], dist];
}

function createGraph(seqLen, edgeMap=[]) {
    let G = {};
    for (let i = 2; i <= seqLen; i++) {
        G[i] = [i - 1, i + 1];
    }
    G[1] = [2];
    G[seqLen] = [seqLen - 1];

    // console.log(edgeMap);
    for (let [k, v] of Object.entries(edgeMap)) {
        // console.log(k, v);
        G[k].push(v);
        G[v].push(k);
        
    }
    return G;
}

function findNeighbors(G, node, depth = 1) {
    let neighbors = [node];

    while (depth) {
        let tempNeighbors = [];
        for (let n of neighbors) {
            tempNeighbors.push(...G[n]);
        }
        neighbors.push(...tempNeighbors);
        neighbors = [...new Set(neighbors)]
        depth--;
    }
    neighbors.splice(neighbors.indexOf(node), 1);
    return neighbors;
}

function getBasePair(G, node){
    for (let n in G[node]){
        if (n != node + 1 && n != node - 1){
            return n;
        }

    }
    return -1;
}


function findMisColoredPairs(edgeMap, predictedColors){
    let misColoredPairs = [];
    for (let [u, v] of Object.entries(edgeMap)) {
        u_color = predictedColors[u];
        v_color = predictedColors[v];

        if (u_color != v_color) {
            misColoredPairs.push(parseInt(u));
            misColoredPairs.push(parseInt(v));
        }

    }
    return misColoredPairs;
}

function predictColors(G, colors, edgeMap, beacketSegement, depth=2, fix=false) {
    
    let predictedColors = {};
    
    for (let [k, [node, color]] of Object.entries(colors)){
        if (color != 'undefined'){
            predictedColors[node] = color;
        }

    }
    // console.log("Before: ", JSON.stringify(predictedColors));
    predictedColors = colorSegments(G, predictedColors, edgeMap, beacketSegement, 4);
    // console.log("After: ", JSON.stringify(predictedColors));

    let nodesToColor = findMisColoredPairs(edgeMap, predictedColors);
    // console.log("nodesToColor :", nodesToColor);
    let coloredNodes = Object.entries(predictedColors).map(([k, v]) => parseInt(k)).filter(k => !nodesToColor.includes(k))

    let nodeColorPair = [];

    let colorStart = Math.min(...Object.keys(predictedColors));
    let lastResidue = Math.max(...Object.keys(G));

    // force remove 0 as we are using 1 index
    let reverse = [...Array(colorStart).keys()].reverse()
    let forward = [...Array(lastResidue - colorStart + 1).keys()].map(i => i + colorStart) 
    let traversalOrder = reverse.concat(forward).filter(item => item !== 0)

    // console.log("traversalOrder", traversalOrder);

    for (let n of traversalOrder) {
        // no need to change as already colored
        if (coloredNodes.includes(n)){
            nodeColorPair.push([n, predictedColors[n]])    
            continue;
        }
        // console.log("n", n)
        let neighbors = findNeighbors(G, n, depth);
        // console.log("neighbors", neighbors);
        
        let neighborColors = neighbors.filter(node => (Object.keys(predictedColors).includes(node.toString()))).map(node => predictedColors[node]);
        // console.log("neighborColors", neighborColors);

        if (neighborColors.length === 0) continue;
        
        let nodeColor = neighborColors.reduce((acc, val) => (acc[val] = (acc[val] || 0) + 1, acc), {});
        // console.log("nodeColor", nodeColor);
        let maxColor = Object.keys(nodeColor).reduce((a, b) => nodeColor[a] > nodeColor[b] ? a : b);
        // console.log("maxColor", maxColor);
        predictedColors[n] = maxColor;
        nodeColorPair.push([n, maxColor]);
        
    }
    console.log(JSON.stringify(nodeColorPair));
    return nodeColorPair;
}

function getMaxColor(colorsList){
    let colorFreq = colorsList.reduce((acc, val) => (acc[val] = (acc[val] || 0) + 1, acc), {});
    let maxColor = Object.keys(colorFreq).reduce((a, b) => colorFreq[a] > colorFreq[b] ? a : b);

    return maxColor;

}

function getBPSegments(edgeList) {
    // let sortedEdgeList = edgeList.slice().sort((a, b) => a[0] - b[0]);
    let segments = [];
    let stack = [];
    let tolerance = 3;

    for (let [u, v] of Object.entries(edgeList)){
        edge = [parseInt(u), parseInt(v)]
        // console.log(edge, stack, segments);
        if (stack.length > 0) {
            let last = stack[stack.length - 1];
            if (Math.max(Math.abs(last[0] - edge[0]), Math.abs(last[1] - edge[1])) < tolerance) {
                stack.push(edge);
                continue;
            } else {
                let first = stack[0];
                let last = stack[stack.length - 1];
                // console.log("Range", first, last);
                let tempSeg = [ ...Array(last[0]+1).keys()].filter(x => x >= first[0]);
                tempSeg.push( ...[ ...Array(first[1]+1).keys()].filter(x => x >= last[1]));
                segments.push(tempSeg);
                stack = [];
            }
        }
        stack.push(edge);
    }

    if (stack.length > 0) {
        let first = stack[0];
        let last = stack[stack.length - 1];
        let tempSeg = [ ...Array(last[0]+1).keys()].filter(x => x >= first[0]);
        tempSeg.push( ...[ ...Array(first[1]+1).keys()].filter(x => x >= last[1]));

        segments.push(tempSeg);
    }
    console.log('BP segments :', JSON.stringify(segments));
    return segments;
}

function generateFinalSegmentation(allSegments, bpSegments) {
    function createSegmentMapping(allSegments, bpSegments) {
        let segmentMapping = {};

        for (let i = 0; i < bpSegments.length; i++) {
            let bpSegment = bpSegments[i];
            let bpKey = `bp${i}`;

            for (let j = 0; j < allSegments.length; j++) {
                let segment = allSegments[j];
                let segmentKey = `st${j}`;

                // Check if there is an intersection or subset
                if (segment.every(elem => bpSegment.includes(elem))) {
                    if (!segmentMapping[bpKey]) {
                        segmentMapping[bpKey] = [];
                    }
                    segmentMapping[bpKey].push(segmentKey);
                }
            }
        }

        return segmentMapping;
    }

    function createFinalSegmentation(allSegments, bpSegments, mapping) {
        let finalSegments = [];

        // Assign segments for base pair segments
        Object.entries(mapping).forEach(([bpKey, segments]) => {
            let mergedSegment = [];
            segments.forEach(segmentKey => {
                let segmentIndex = parseInt(segmentKey.substring(2));
                mergedSegment.push(...allSegments[segmentIndex]);
            });
            finalSegments.push(mergedSegment);
        });

        // Add segments for non-mapped segments
        allSegments.forEach((segment, index) => {
            let isMapped = Object.values(mapping).some(segments => segments.includes(`st${index}`));
            if (!isMapped) {
                finalSegments.push(segment);
            }
        });

        return finalSegments;
    }

    let mapping = createSegmentMapping(allSegments, bpSegments);
    let finalSegmentation = createFinalSegmentation(allSegments, bpSegments, mapping);
    
    return finalSegmentation;
}

function findActualSegments(dotBracket) {
    let stack = [];
    let segmentList = [];

    for (let i = 0; i < dotBracket.length; i++) {
        if (dotBracket[i] === "(") {
            stack.push(i);
        } else if (dotBracket[i] === ")") {
            let stackTop = stack.pop();
            segmentList.push([stackTop, i]);
        }
    }

    let actualSegments = [];

    for (let s of segmentList) {
        let values = new Set([...Array(s[1] - s[0] + 2).keys()].map(x => x + s[0] + 1));
        for (let i of actualSegments) {
            values = new Set([...values].filter(x => !i.has(x)));
        }
        actualSegments.push(values);
    }

    return actualSegments;
}

function calculateDotBracket(basePairsList, sequenceLength) {
    let dotBracket = Array(sequenceLength).fill(".");
    for (let basePair of basePairsList) {
        dotBracket[basePair[0] - 1] = "(";
        dotBracket[basePair[1] - 1] = ")";
    }
    
    return dotBracket;
}

function colorSegments(graph, colors, edgeMap, beacketSegement, minLength = 4) {
    let stack = [];
    let state = 0;
    let nodeColorPair = [];
    let allSegments = [];
    let bpSegments = getBPSegments(edgeMap);
    let predColors = colors;
    
    for (let [node, neighbours] of Object.entries(graph)) {
        if (neighbours.length !== state && stack.length >= minLength) {
            // console.log("Current segment:", stack);
            allSegments.push(stack);
            // let segColors = [];
            // for (let n of stack){
            //     segColors.push(colors[n]);
            // }
            // let maxColor = getMaxColor(segColors);
            // for (let n of stack){
            //     nodeColorPair.push([parseInt(n), maxColor]);
            // }
            stack = [];
        } 
        stack.push(parseInt(node));
        state = neighbours.length;
        
    }
    if (stack.length >= 1) {
        // console.log("Current segment:", stack);
        allSegments.push(stack);
        // let segColors = [];
        // for (let n of stack) {
        //     segColors.push(colors[n]);
        // }
        // let maxColor = getMaxColor(segColors);
        // for (let n of stack){
        //     nodeColorPair.push([parseInt(n), maxColor]);
        // }
    }
    
    // let finalSegments = generateFinalSegmentation(allSegments, bpSegments);
    // let finalSegments = bpSegments;
    // console.log("Stacked Segments: ", JSON.stringify(finalSegments));

    for (let finalSegments of [allSegments, bpSegments]){ // beacketSegement
        for (let i = 0; i < finalSegments.length; i++) {
            let segment = finalSegments[i];

            if (segment.length < minLength) {
                continue;
            }

            let segColors = [];
            for (let n of segment) {
                segColors.push(predColors[n]); 
            }
            let maxColor = getMaxColor(segColors);
            if (maxColor == "undefined"){
                continue;
            }
            for (let n of segment) {
                // nodeColorPair.push([parseInt(n), maxColor]);
                predColors[n] = maxColor;
            }
        }
    }

    // return nodeColorPair;
    return predColors;
}



function fix_colors(sequence, basePairsList, dataMapJson) {
    let result = {};
    let allEdges = {};
    // console.log("Starting to fill and fix colors");
    for (let bond of basePairsList) {
        let extracted = extractBonds(bond);
        // console.log("Bond length", extracted);
        if (extracted) {
            let [u, v, dist] = extracted;
            if (dist < 10) {
                allEdges[u] = v;
            }
        }
    }

    let G = createGraph(sequence.length, allEdges);
    
    // new addition, breaks things
    // let dotBracket = calculateDotBracket(allEdges, sequence.length);
    // console.log('dotBracket', dotBracket);
    let beacketSegement = [];
    // let beacketSegement = findActualSegments(dotBracket);
    // console.log('beacketSegement', beacketSegement);
    // remove above and `beacketSegement` from everywhere 


    let edgeMap = allEdges;
    // console.log("Graph :", JSON.stringify(G));
    for (let [k, colors] of Object.entries(dataMapJson)) {
        let fix = false;
        if (['helix', "Helix"].includes(k)) {
            edgeMap = allEdges;
        }
        else {
            continue;
            edgeMap = [];
        }
        // console.log(k);
        // let predColors = predictColors(graphStackedOnly, Object.assign({}, colors), allEdges, 4, fix);
        let predColors = predictColors(G, Object.assign({}, colors), edgeMap, beacketSegement, 4, fix);
        result[k] = predColors.sort((a, b) => a[0] - b[0]);
    }
    // console.log("Missing colors filled and fixed")
    
    return result;
}

// Exporting the function
module.exports = fix_colors;
