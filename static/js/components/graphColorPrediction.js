const { map } = require("lodash");


function euclideanDistance(x1, y1, x2, y2) {
    return Math.sqrt(Math.pow(parseFloat(x2) - parseFloat(x1), 2) + Math.pow(parseFloat(y2) - parseFloat(y1), 2));
}

function extractBonds(st) {
    let basePairRegex = /cWW_(\d+)_(\d+)/;
    
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

/**
 * Create a simple graph {node: list of neghbouring nodes} graph
 * @param {int} seqLen 
 * @param {map} edgeMap 
 * @returns 
 */
function createGraph(seqLen, edgeMap = []) {
    let G = {};
    for (let i = 2; i <= seqLen; i++) {
        G[i] = [i - 1, i + 1];
    }
    G[1] = [2];
    G[seqLen] = [seqLen - 1];

    for (let [k, v] of Object.entries(edgeMap)) {
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

/**
 * based on colors, find miscolored base paired nodes, 
 * and add both of them for recoloring
 * @param {*} edgeMap 
 * @param {*} predictedColors 
 * @returns 
 */
function findMisColoredPairs(edgeMap, predictedColors) {
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



/**
 * Given a list of values
 * return the value with the maximum occurance 
 * @param {*} colorsList 
 * @returns 
 */
function getMaxColor(colorsList) {
    colorsList = colorsList.filter(value => value != undefined);
    let maxColor = '0';
    
    if(colorsList.length > 0){
        let colorFreq = colorsList.reduce((acc, val) => (acc[val] = (acc[val] || 0) + 1, acc), {});
        maxColor = Object.keys(colorFreq).reduce((a, b) => colorFreq[a] > colorFreq[b] ? a : b);
    }

    return maxColor;

}

/**
 * NOT IN USE / REPLACED BY segmentBracketElements
 * Extracts segments from the given edge list representing base pairs.
 * @param {Object} edgeList - The edge list representing base pairs.
 * @returns {Array} An array containing segments of base pairs.
 */
function getBPSegments(edgeList) {
    // Initialize an array to store segments
    let segments = [];
    // Initialize a stack to track consecutive base pairs forming segments
    let stack = [];
    // Set tolerance for considering base pairs as part of the same segment
    let tolerance = 3;

    // Iterate through each entry in the edge list
    for (let [u, v] of Object.entries(edgeList)) {
        // Convert edge indices to integers
        edge = [parseInt(u), parseInt(v)];

        // Check if the stack is not empty
        if (stack.length > 0) {
            // Get the last element in the stack
            let last = stack[stack.length - 1];
            // Check if the difference between current edge and last edge is within tolerance
            if (Math.max(Math.abs(last[0] - edge[0]), Math.abs(last[1] - edge[1])) < tolerance) {
                // If within tolerance, push current edge to the stack and continue
                stack.push(edge);
                continue;
            } else {
                // If outside tolerance, form a segment from the stack elements
                let first = stack[0];
                let last = stack[stack.length - 1];
                let tempSeg = [...Array(last[0] + 1).keys()].filter(x => x >= first[0]);
                tempSeg.push(...[...Array(first[1] + 1).keys()].filter(x => x >= last[1]));
                segments.push(tempSeg);
                stack = []; // Clear the stack for the next segment
            }
        }
        stack.push(edge); // Push current edge to the stack
    }

    // If there are remaining edges in the stack, form a segment from them
    if (stack.length > 0) {
        let first = stack[0];
        let last = stack[stack.length - 1];
        let tempSeg = [...Array(last[0] + 1).keys()].filter(x => x >= first[0]);
        tempSeg.push(...[...Array(first[1] + 1).keys()].filter(x => x >= last[1]));
        segments.push(tempSeg);
    }

    // Return the extracted segments of base pairs
    return segments;
}


/**
 * Segments elements enclosed within parentheses in the given dot-bracket notation.
 * @param {string} dotBracket - The dot-bracket notation representing RNA secondary structure.
 * @returns {Array} An array containing segments of elements enclosed within parentheses.
 */
function segmentBracketElements(dotBracket) {
    // Initialize a stack to track opening parentheses indices
    let stack = [];
    // Initialize an array to store segments of elements enclosed within parentheses
    let segments = [];

    // Iterate through each character in the dot-bracket notation
    for (let i = 0; i < dotBracket.length; i++) {
        // If the character is an opening parenthesis, push its index to the stack
        if (dotBracket[i] === "(") {
            stack.push(i);
        } 
        // If the character is a closing parenthesis
        else if (dotBracket[i] === ")") {
            // Pop the index of the corresponding opening parenthesis from the stack
            let boundaries = [stack.pop() + 1, i + 1];
            // Check if there are existing segments
            if (segments.length !== 0) {
                let last = segments[segments.length - 1];
                // Check if the new segment is close enough to the last segment
                if (Math.abs(last[0] - boundaries[0]) < 3 && Math.abs(last[1] - boundaries[1]) < 3) {
                    // If close enough, remove the last segment
                    segments.pop();
                }
            }
            // Push the boundaries of the current segment to the segments array
            segments.push(boundaries);
        }
    }

    // Initialize an array to store the segments of elements enclosed within parentheses
    let elementsSegments = [];

    // Iterate through each segment of parentheses
    for (let s of segments) {
        // Create a set containing indices of elements within the segment
        s[0] = isNaN(s[0]) ? 1 : s[0];
        let values = new Set([...Array(s[1] - s[0] + 1).keys()].map(x => x + s[0]));

        // Iterate through existing segments
        for (let i of elementsSegments) {
            // Filter out indices already present in the existing segments
            values = new Set([...values].filter(x => !i.includes(x)));
        }

        // Convert the set to an array
        values = [...values];

        // Find the locations where there are jumps in indices
        let jumpLocations = [];
        for (let i = 0; i < values.length - 1; i++) {
            if (values[i + 1] - values[i] !== 1) {
                jumpLocations.push(i);
            }
        }

        // If there are multiple jumps, remove the corresponding elements
        if (jumpLocations.length > 1) {
            let start = jumpLocations[0] + 1;
            let end = jumpLocations[jumpLocations.length - 1] + 1;
            values.splice(start, end - start);
        }

        // Push the final segments of elements enclosed within parentheses to the elementsSegments array
        elementsSegments.push(values);
    }

    // Return the segments of elements enclosed within parentheses
    return elementsSegments;
}

/**
 * Segments stacked elements in the given graph based on a minimum length threshold.
 * @param {Object} graph - The graph representing nodes and their neighbors.
 * @param {number} minLength - The minimum length of a segment to be considered.
 * @returns {Array} An array containing segments of stacked elements.
 */
function segmentStackedElements(graph, minLength = 4) {
    // Initialize an array to store segmented stacked elements
    let stackedSegments = [];
    // Initialize a stack to track consecutive elements forming segments
    let stack = [];
    // Initialize a variable to track the state of neighbors' length
    let state = 0;

    // Iterate through each node and its neighbors in the graph
    for (let [node, neighbours] of Object.entries(graph)) {
        // Check if the number of neighbors has changed and the stack length is greater than or equal to the minimum length
        if (neighbours.length !== state && stack.length >= minLength) {
            // Push the current stack (segment) to the stackedSegments array
            stackedSegments.push(stack);
            // Reset the stack for the next segment
            stack = [];
        }
        // Push the current node to the stack
        stack.push(parseInt(node));
        // Update the state to the current number of neighbors
        state = neighbours.length;
    }

    // If there are remaining elements in the stack, push them as the last segment
    if (stack.length >= 1) {
        stackedSegments.push(stack);
    }

    // Return the segmented stacked elements
    return stackedSegments;
}


function colorSegments(colors, listOfSegements) {
    let predColors = colors;

    for (let finalSegments of listOfSegements) { // bpSegments, beacketSegement, allSegments
        for (let i = 0; i < finalSegments.length; i++) {
            let segment = finalSegments[i];

            let segColors = [];
            for (let n of segment) {
                segColors.push(predColors[n]);
            }
            
            let maxColor = getMaxColor(segColors);
            // console.log("Segments", segment, segColors, maxColor);
            if (maxColor == "undefined") {
                // continue;
                maxColor = '-1';
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


function predictColors(G, predictedColors, edgeMap, depth=2, isHelix=false) {
    let nodeColorPair = [];
    let coloredNodes = [];
    if (isHelix == true){
        let nodesToColor = findMisColoredPairs(edgeMap, predictedColors);
        coloredNodes = Object.entries(predictedColors).map(([k, v]) => parseInt(k)).filter(k => !nodesToColor.includes(k))
    }
    let colorStart = Math.min(...Object.keys(predictedColors));
    let lastResidue = Math.max(...Object.keys(G));

    // force remove 0 as we are using 1 index
    let reverse = [...Array(colorStart).keys()].reverse()
    let forward = [...Array(lastResidue - colorStart + 1).keys()].map(i => i + colorStart) 
    let traversalOrder = reverse.concat(forward).filter(item => item !== 0)


    for (let n of traversalOrder) {
        // no need to change as already colored
        
        if (coloredNodes.includes(n)){ 
            nodeColorPair.push([n, predictedColors[n]])    
            continue;
        }
        
        let neighbors = findNeighbors(G, n, depth);
        // 
        
        let neighborColors = neighbors.filter(node => (Object.keys(predictedColors).includes(node.toString()))).map(node => predictedColors[node]);
        
        if (neighborColors.length === 0) continue;

        let maxColor = getMaxColor(neighborColors);
        
        predictedColors[n] = maxColor;
        nodeColorPair.push([n, maxColor]);
        
    }
    return nodeColorPair;
}

function fix_colors(sequence, basePairsList, dataMapJson) {
    let result = {};
    let allEdges = {};
    
    let dotBracket_array = Array(sequence.length).fill(".");
    for (let bond of basePairsList) {
        let extracted = extractBonds(bond);
        if (extracted) {
            let [u, v, dist] = extracted;
            if (dist < 10) {
                allEdges[u] = v;
                dotBracket_array[parseInt(u) - 1] = "(";
                dotBracket_array[parseInt(v) - 1] = ")";
            }
        }
    }

    let G = createGraph(sequence.length, allEdges);
    
    let dotBracket = dotBracket_array.join('');
    let beacketSegment = segmentBracketElements(dotBracket);
    // let stackedSegment = segmentStackedElements(G, minLength = 4);

    
    let tempColors;
    let isHelix = false;
    for (let [k, colors] of Object.entries(dataMapJson)) {
        // remove undefined colors 
        let predictedColors = {};
    
        // remove undefined colors from map
        for (let [_, [node, color]] of Object.entries(colors)){
            if (color != 'undefined'){
                
                // if (parseInt(color) < 5){
                    predictedColors[node] = color;
                // }
            }
    
        }
        // console.log("predictedColors", predictedColors);
        
        if (['helix', "Helix"].includes(k)) {
            tempColors = colorSegments(predictedColors, [beacketSegment]);
            isHelix = true;
        }
        else {
            // tempColors = colorSegments(predictedColors, [stackedSegment]);
            tempColors = predictedColors;
            isHelix = false;
        }
        // console.log(k, tempColors);
        let predColors = predictColors(G, tempColors, allEdges, 4, isHelix);
        result[k] = predColors.sort((a, b) => a[0] - b[0]);
    }

    return result;
}

// Exporting the function
module.exports = fix_colors;
