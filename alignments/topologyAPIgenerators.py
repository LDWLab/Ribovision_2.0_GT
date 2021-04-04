import re

def generateTopologyJSONfromSVG(svgContents, pdbID, chainID, entityID):
    '''Generates topology JSON from ProOrigami svg output'''
    helices = []
    coils = []
    strands = []
    terms = []
    extents = []
    for helixMatch in re.findall(r"<rect\s+[^>]*dunnart:type\s*=\s*\"bioHelix\"[^>]*>", svgContents):
        helixXMatch = re.search(r"\s+x\s*=\s*\"(-?[\d.]*)\"", helixMatch)
        if helixXMatch is None:
            raise Exception("SVG x variable not found within dunnart bioHelix object.")
        x = float(helixXMatch.group(1))

        helixYMatch = re.search(r"\s+y\s*=\s*\"(-?[\d.]*)\"", helixMatch)
        if helixYMatch is None:
            raise Exception("SVG y variable not found within dunnart bioHelix object.")
        y = float(helixYMatch.group(1))

        helixWidthMatch = re.search(r"\s+width\s*=\s*\"([\d.]*)\"", helixMatch)
        if helixWidthMatch is None:
            raise Exception("SVG width variable not found within dunnart bioHelix object.")
        width = float(helixWidthMatch.group(1))

        helixHeightMatch = re.search(r"\s+height\s*=\s*\"([\d.]*)\"", helixMatch)
        if helixHeightMatch is None:
            raise Exception("SVG height variable not found within dunnart bioHelix object.")
        height = float(helixHeightMatch.group(1))

        helixRxMatch = re.search(r"\s+rx\s*=\s*\"([\d.]*)\"", helixMatch)
        if helixRxMatch is None:
            raise Exception("SVG rx variable not found within dunnart bioHelix object.")
        rx = float(helixRxMatch.group(1))

        helixRyMatch = re.search(r"\s+ry\s*=\s*\"([\d.])*\"", helixMatch)
        if helixRyMatch is None:
            raise Exception("SVG ry variable not found within dunnart bioHelix object.")
        ry = float(helixRyMatch.group(1))

        residueSequenceNumbers = re.search(r"\s+proorigami:residueSeqNums\s*=\s*\"[\s\d.]*\"", helixMatch)
        if residueSequenceNumbers is None:
            print ("match: " + helixMatch)
            raise Exception("Residue sequence numbers not found within dunnart bioHelix object.")
        residueSequenceNumbers = residueSequenceNumbers.group(0)
        firstQuoteIndex = residueSequenceNumbers.index('\"')
        secondQuoteIndex = residueSequenceNumbers.index('\"', firstQuoteIndex + 1)
        residueSequenceNumbers = residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex].split()
        start = int(residueSequenceNumbers[0])
        stop = int(residueSequenceNumbers[-1])
        helices.append({
            "start" : start,
            "stop" : stop,
            "path" : [x, y, x + width, y + height],
            "minoraxis" : rx,
            "majoraxis" : ry
        })
    for strandMatch in re.findall(r"<path\s+[^>]*dunnart:type\s*=\s*\"bioStrand\"[^>]*>", svgContents):
        strandPath = re.search(r"\s+d\s*=\s*\"[MLCZz.,\-\d\s]*\"", strandMatch)
        if strandPath is None:
            raise Exception("SVG path not found within dunnart bioStrand object.")
        strandPath = strandPath.group(0)
        firstQuoteIndex = strandPath.index('\"')
        secondQuoteIndex = strandPath.index('\"', firstQuoteIndex + 1)
        strandPath = strandPath[firstQuoteIndex + 1:secondQuoteIndex]
        lineCoordinates = []
        for pathElementMatch in re.finditer(r"(M|L)\s+([\d\s.,-]+)", strandPath):
            x, y = pathElementMatch.group(2).split(',')
            x = float(x)
            y = float(y)
            lineCoordinates.append(x)
            lineCoordinates.append(y)
        lineCoordinates = lineCoordinates[:-2]
        lineCoordinates = lineCoordinates[8:] + lineCoordinates[:8]
        residueSequenceNumbers = re.search(r"\s+proorigami:residueSeqNums\s*=\s*\"([\d\s])*\"", strandMatch)
        if residueSequenceNumbers is None:
            raise Exception("Residue sequence numbers not found within dunnart bioStrand object.")
        residueSequenceNumbers = residueSequenceNumbers.group(0)
        firstQuoteIndex = residueSequenceNumbers.index('\"')
        secondQuoteIndex = residueSequenceNumbers.index('\"', firstQuoteIndex + 1)
        residueSequenceNumbers = residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex].split()
        start = int(residueSequenceNumbers[0])
        stop = int(residueSequenceNumbers[-1])
        strands.append({
            "start" : start,
            "stop" : stop,
            "path" : lineCoordinates
        })
    for coilMatch in re.findall(r"<path\s+[^>]*dunnart:type\s*=\s*\"connAvoidPoly\"[^>]*>", svgContents):
        lineCoordinates = []
        # print (coilMatch)
        strandPath = re.search(r"\s+d\s*=\s*\"[MLCZz.,\-\d\s]*\"", coilMatch)
        if strandPath is None:
            raise Exception("SVG path not found within dunnart bioStrand object.")
        strandPath = strandPath.group(0)
        # print ("strandPath: ", strandPath)
        for pathElementMatch in re.finditer(r"(M|L|C)\s+([\d\s.,-]+)", strandPath):
            # print ("\tpathElementMatch: " + pathElementMatch.group(0) + " | " + pathElementMatch.group(1) + " | " + pathElementMatch.group(2))
            group1 = pathElementMatch.group(1)
            if group1 == "C":
                group2 = pathElementMatch.group(2)
                split = re.sub(r'\s+', ',', group2.strip()).split(',')
                # print ("group2: ", group2, " | group2.split(,): ", split)
                x0, y0, x1, y1, x2, y2 = split
                x0 = float(x0)
                y0 = float(y0)
                x1 = float(x1)
                y1 = float(y1)
                x2 = float(x2)
                y2 = float(y2)
                lineCoordinates.append(x0)
                lineCoordinates.append(y0)
                lineCoordinates.append(x1)
                lineCoordinates.append(y1)
                lineCoordinates.append(x2)
                lineCoordinates.append(y2)
            else:
                x, y = pathElementMatch.group(2).split(',')
                x = float(x)
                y = float(y)
                lineCoordinates.append(x)
                lineCoordinates.append(y)

        residueSequenceNumbers = re.search(r"\s+proorigami:residueSeqNums\s*=\s*\"([\d\s])*\"", coilMatch)
        if residueSequenceNumbers is None:
            raise Exception("Residue sequence numbers not found within dunnart bioStrand object.")
        residueSequenceNumbers = residueSequenceNumbers.group(0)
        firstQuoteIndex = residueSequenceNumbers.index('\"')
        secondQuoteIndex = residueSequenceNumbers.index('\"', firstQuoteIndex + 1)
        residueSequenceNumbers = residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex].split()
        # print ("__len(residueSequenceNumbers): ", len(residueSequenceNumbers), "__")
        if len(residueSequenceNumbers) > 0:
            start = int(residueSequenceNumbers[0])
            stop = int(residueSequenceNumbers[-1])
        else:
            start = -1
            stop = -1
        coils.append({
            "start" : start,
            "stop" : stop,
            "path" : lineCoordinates
        })
    tree = {
        pdbID : {
            entityID : {
                chainID : {
                    "helices" : helices,
                    "coils" : coils,
                    "strands" : strands,
                    "terms" : terms,
                    "extents" : extents
                }
            }
        }
    }
    return tree

def generateEntityJSON (pdbID, entityID, sequenceStr, start, end):
    return {pdbID:[{
                "entity_id":entityID,
                "sequence":sequenceStr,
                "source":[{
                    "mappings":[{
                        "start":{"residue_number":start},
                        "end":{"residue_number":end}
                        }]
                    }],
                "length":len(sequenceStr)
            }]}

def generatePolCoverageJSON(pdbID, chainID, entityID, startAuth, startNum, endAuth, endNum):
    return {pdbID:{
                "molecules":[{
                    "entity_id":entityID,
                    "chains":[{
                        "observed":[{
                            "start":{"author_residue_number":startAuth,"residue_number":startNum},
                            "end":{"author_residue_number":endAuth,"residue_number":endNum}
                        }],
                        "chain_id": chainID
                    }]
                }]
            }}