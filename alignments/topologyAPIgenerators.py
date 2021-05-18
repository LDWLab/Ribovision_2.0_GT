import re
from shapely.geometry import Point, Polygon, LineString

def generateTopologyJSONfromSVG(svgContents, pdbID, chainID, entityID):
    '''Generates topology JSON from ProOrigami svg output'''
    helices = []
    coils = []
    strands = []
    terms = []
    extents = []
    non_printable_data_map = {}
    maximum_stop_found = -1
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
            raise Exception("Residue sequence numbers not found within dunnart bioHelix object.")
        residueSequenceNumbers = residueSequenceNumbers.group(0)
        firstQuoteIndex = residueSequenceNumbers.index('\"')
        secondQuoteIndex = residueSequenceNumbers.index('\"', firstQuoteIndex + 1)
        residueSequenceNumbers = parseSequentialResidueSequenceNumbers(residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex])
        start = int(residueSequenceNumbers[0])
        stop = int(residueSequenceNumbers[-1])
        if stop > maximum_stop_found:
            maximum_stop_found = stop

        dunnart_reversed = re.search(r"\s+dunnart:reversed\s*=\s*\"([01]{1})\"", helixMatch)
        if dunnart_reversed is None:
            raise Exception("Dunnart:reversed flag not found within dunnart bioHelix object.")
        dunnart_reversed = dunnart_reversed.group(1)
        
        # print ("Dunnart_reversed:", dunnart_reversed, (dunnart_reversed == "0"))
        helix_reversed_flag = dunnart_reversed == "0"
        path = [x, y, x + width, y + height]
        new_helix = {
            "start" : start,
            "stop" : stop,
            "path" : path
        }
        non_printable_data_map[(start, stop)] = {
            "reversed" : helix_reversed_flag,
            "rx" : rx,
            "ry" : ry
        }
        helices.append(new_helix)
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
        # print("Before: ", lineCoordinates)
        lineCoordinates = lineCoordinates[8:] + lineCoordinates[:8]
        # print("After: ", lineCoordinates)
        residueSequenceNumbers = re.search(r"\s+proorigami:residueSeqNums\s*=\s*\"([\d\s])*\"", strandMatch)
        if residueSequenceNumbers is None:
            raise Exception("Residue sequence numbers not found within dunnart bioStrand object.")
        residueSequenceNumbers = residueSequenceNumbers.group(0)
        firstQuoteIndex = residueSequenceNumbers.index('\"')
        secondQuoteIndex = residueSequenceNumbers.index('\"', firstQuoteIndex + 1)
        residueSequenceNumbers = parseSequentialResidueSequenceNumbers(residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex])
        start = int(residueSequenceNumbers[0])
        stop = int(residueSequenceNumbers[-1])
        if stop > maximum_stop_found:
            maximum_stop_found = stop
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
        residueSequenceNumbers = parseSequentialResidueSequenceNumbers(residueSequenceNumbers[firstQuoteIndex + 1:secondQuoteIndex])
        # print ("__len(residueSequenceNumbers): ", len(residueSequenceNumbers), "__")
        if len(residueSequenceNumbers) > 0:
            start = int(residueSequenceNumbers[0])
            stop = int(residueSequenceNumbers[-1])
        else:
            start = -1
            stop = -1

        if stop > maximum_stop_found:
            maximum_stop_found = stop

        coils.append({
            "start" : start,
            "stop" : stop,
            "path" : lineCoordinates
        })
    
    n_terminus_match = re.search(r"<rect\s+[^>]*dunnart:label\s*=\s*\"N\"[^>]*>", svgContents)
    if n_terminus_match is None:
        raise Exception("N terminus not found SVG file.")
    n_terminus_match = n_terminus_match.group(0)

    xMatch = re.search(r"\s+x\s*=\s*\"([\d\-.]+)\"", n_terminus_match)
    if xMatch is None:
        raise Exception("x variable not found within n terminus")
    x = float(xMatch.group(1))

    yMatch = re.search(r"\s+y\s*=\s*\"([\d\-.]+)\"", n_terminus_match)
    if yMatch is None:
        raise Exception("y variable not found within n terminus")
    y = float(yMatch.group(1))
    
    widthMatch = re.search(r"\s+width\s*=\s*\"([\d.]+)\"", n_terminus_match)
    if widthMatch is None:
        raise Exception("width variable not found within n terminus")
    width = float(widthMatch.group(1))

    heightMatch = re.search(r"\s+height\s*=\s*\"([\d.]+)\"", n_terminus_match)
    if heightMatch is None:
        raise Exception("height variable not found within n terminus")
    height = float(heightMatch.group(1))

    # print ("N: " + n_terminus_match + "\n\tx: " + str(x) + "\n\ty: " + str(y) + "\n\twidth: " + str(width) + "\n\theight: " + str(height))
    terms.append({
        "resnum" : "1",
        "type" : "N",
        "start" : -1,
        "stop" : -1,
        "path" : [
            x, y,
            x + width, y,
            x + width, y + height,
            x, y + height
        ]
    })
    
    c_terminus_match = re.search(r"<rect\s+[^>]*dunnart:label\s*=\s*\"C\"[^>]*>", svgContents)
    if c_terminus_match is None:
        raise Exception("C terminus not found SVG file.")
    c_terminus_match = c_terminus_match.group(0)

    xMatch = re.search(r"\s+x\s*=\s*\"([\d\-.]+)\"", c_terminus_match)
    if xMatch is None:
        raise Exception("x variable not found within n terminus")
    x = float(xMatch.group(1))

    yMatch = re.search(r"\s+y\s*=\s*\"([\d\-.]+)\"", c_terminus_match)
    # print ("yMatch: " + yMatch.group(0) + " " + yMatch.group(1))
    if yMatch is None:
        raise Exception("y variable not found within n terminus")
    y = float(yMatch.group(1))
    
    widthMatch = re.search(r"\s+width\s*=\s*\"([\d.]+)\"", c_terminus_match)
    if widthMatch is None:
        raise Exception("width variable not found within n terminus")
    width = float(widthMatch.group(1))

    heightMatch = re.search(r"\s+height\s*=\s*\"([\d.]+)\"", c_terminus_match)
    if heightMatch is None:
        raise Exception("height variable not found within n terminus")
    height = float(heightMatch.group(1))

    # print ("C: " + c_terminus_match + "\n\tx: " + str(x) + "\n\ty: " + str(y) + "\n\twidth: " + str(width) + "\n\theight: " + str(height))
    terms.append({
        "resnum" : str(maximum_stop_found),
        "type" : "C",
        "start" : -1,
        "stop" : -1,
        "path" : [
            x, y,
            x + width, y,
            x + width, y + height,
            x, y + height
        ]
    })

    for coil in coils:
        coil_path = coil["path"]
        coil_start_clipped_flag = False
        coil_end_clipped_flag = False
        for strand in strands:
            last_coil_vertex_index_in_strand = -1
            first_coil_vertex_index_in_strand = -1
            strand_path = strand["path"]
            for vertex_index in range(0, len(coil_path), 2):
                x = coil_path[vertex_index]
                y = coil_path[vertex_index + 1]
                vertex_in_strand = vertex_in_polygon(x, y, strand_path)
                if not vertex_in_strand:
                    break
                last_coil_vertex_index_in_strand = vertex_index
            if last_coil_vertex_index_in_strand >= 0:
                # print ("Coil-Strand intersection found.")
                previous_x = strand_path[-2]
                previous_y = strand_path[-1]
                coil_line = LineString([(coil_path[last_coil_vertex_index_in_strand], coil_path[last_coil_vertex_index_in_strand + 1]), (coil_path[last_coil_vertex_index_in_strand + 2], coil_path[last_coil_vertex_index_in_strand + 3])])
                for vertex_index in range(0, len(strand_path), 2):
                    current_x = strand_path[vertex_index]
                    current_y = strand_path[vertex_index + 1]
                    strand_line = LineString([(previous_x, previous_y), (current_x, current_y)])
                    if coil_line.intersects(strand_line):
                        intersection = coil_line.intersection(strand_line)
                        # print ("Intersection: " + str(intersection))
                        edited_path = [intersection.x, intersection.y] + coil_path[last_coil_vertex_index_in_strand + 2:]
                        coil_path = edited_path
                        coil["path"] = edited_path
                        break
                    previous_x = current_x
                    previous_y = current_y

            for vertex_index in range(len(coil_path) - 2, 0, -2):
                x = coil_path[vertex_index]
                y = coil_path[vertex_index + 1]
                vertex_in_strand = vertex_in_polygon(x, y, strand_path)
                if not vertex_in_strand:
                    break
                first_coil_vertex_index_in_strand = vertex_index
            if first_coil_vertex_index_in_strand >= 0:
                coil_line = LineString([(coil_path[first_coil_vertex_index_in_strand - 2], coil_path[first_coil_vertex_index_in_strand - 1]), (coil_path[first_coil_vertex_index_in_strand], coil_path[first_coil_vertex_index_in_strand + 1])])
                # print ("Coil Line: " + str(coil_line))
                strand_line = LineString([(strand_path[0], strand_path[1]), (strand_path[12], strand_path[13])])
                # print ("Strand Line: " + str(strand_line))
                if coil_line.intersects(strand_line):
                    intersection = coil_line.intersection(strand_line)
                    edited_path = coil_path[:first_coil_vertex_index_in_strand] + [intersection.x, intersection.y]
                    coil_path = edited_path
                    coil["path"] = edited_path
        for helix in helices:
            last_coil_vertex_index_in_helix = -1
            first_coil_vertex_index_in_helix = -1
            min_x, min_y, max_x, max_y = helix["path"]
            helix_path = [min_x, min_y, min_x, max_y, max_x, max_y, max_x, min_y]
            # print ("Helix Path: " + str(helix_path))

            for vertex_index in range(0, len(coil_path), 2):
                x = coil_path[vertex_index]
                y = coil_path[vertex_index + 1]
                vertex_in_helix = vertex_in_polygon(x, y, helix_path)
                if not vertex_in_helix:
                    break
                last_coil_vertex_index_in_helix = vertex_index
            if last_coil_vertex_index_in_helix >= 0:
                coil_line = LineString([(coil_path[last_coil_vertex_index_in_helix], coil_path[last_coil_vertex_index_in_helix + 1]), (coil_path[last_coil_vertex_index_in_helix + 2], coil_path[last_coil_vertex_index_in_helix + 3])])
                previous_x = helix_path[-2]
                previous_y = helix_path[-1]
                for vertex_index in range(0, len(helix_path), 2):
                    current_x = helix_path[vertex_index]
                    current_y = helix_path[vertex_index + 1]
                    helix_line = LineString([(previous_x, previous_y), (current_x, current_y)])
                    if (coil_line.intersects(helix_line)):
                        intersection = coil_line.intersection(helix_line)
                        edited_path = [intersection.x, intersection.y] + coil_path[last_coil_vertex_index_in_helix + 2:]
                        # print ("Start Clipped.\n\tCoil Path: " + str(coil_path) + "\n\tEdited Path: " + str(edited_path))
                        coil_path = edited_path
                        coil["path"] = edited_path
                        break
                    previous_x = current_x
                    previous_y = current_y
            for vertex_index in range(len(coil_path) - 2, 0, -2):
                x = coil_path[vertex_index]
                y = coil_path[vertex_index + 1]
                vertex_in_helix = vertex_in_polygon(x, y, helix_path)
                if not vertex_in_helix:
                    break
                first_coil_vertex_index_in_helix = vertex_index
            if first_coil_vertex_index_in_helix >= 0:
                # print ("Vertex within helix.")
                coil_line = LineString([(coil_path[first_coil_vertex_index_in_helix - 2], coil_path[first_coil_vertex_index_in_helix - 1]), (coil_path[first_coil_vertex_index_in_helix], coil_path[first_coil_vertex_index_in_helix + 1])])
                previous_x = helix_path[-2]
                previous_y = helix_path[-1]
                for vertex_index in range(0, len(helix_path), 2):
                    current_x = helix_path[vertex_index]
                    current_y = helix_path[vertex_index + 1]
                    helix_line = LineString([(previous_x, previous_y), (current_x, current_y)])
                    # print ("\tHelix Line: " + str(helix_line))
                    if coil_line.intersects(helix_line):
                        intersection = coil_line.intersection(helix_line)
                        edited_path = coil_path[:first_coil_vertex_index_in_helix] + [intersection.x, intersection.y]
                        # print ("End Clipped.\n\tCoil Path: " + str(coil_path) + "\n\tEdited Path: " + str(edited_path))
                        coil_path = edited_path
                        coil["path"] = edited_path
                        break
                    previous_x = current_x
                    previous_y = current_y
        
        start = coil["start"]
        stop = coil["stop"]
        minimum_num_vertices_per_coil = stop - start + 3
        deficit = minimum_num_vertices_per_coil - len(coil_path) // 2
        # print ("Start: " + str(start) + "\tStop: " + str(stop) + "\tMinimum #Verts: " + str(minimum_num_vertices_per_coil) + "\tDeficit: " + str(deficit))
        if deficit > 0:
            x0 = coil_path[0]
            y0 = coil_path[1]
            x1 = coil_path[2]
            y1 = coil_path[3]
            dummy_points = []
            t = 0
            dt = 1 / (deficit + 1)
            for i in range(deficit):
                t += dt
                # print ("T: " + str(t))
                interpolated_x, interpolated_y = interpolate(x0, y0, x1, y1, t)
                dummy_points.append(round(interpolated_x, 3))
                dummy_points.append(round(interpolated_y, 3))
            coil_path = coil_path[0:2] + dummy_points + coil_path[2:]
            coil["path"] = coil_path

    for helix in helices:
        non_printable_data_map_entry = non_printable_data_map[(helix["start"], helix["stop"])]
        min_x, min_y, max_x, max_y = helix["path"]
        width = max_x - min_x
        height = max_y - min_y
        if width > height:
            smaller_dimension = height
            # helix["minoraxis"]
            dx = non_printable_data_map_entry["rx"]
            min_x += dx
            max_x -= dx
        else:
            smaller_dimension = width
            # helix["majoraxis"]
            dy = non_printable_data_map_entry["ry"]
            min_y += dy
            max_y -= dy
        helix["majoraxis"] = smaller_dimension / 2
        helix["minoraxis"] = smaller_dimension / 4
        if non_printable_data_map_entry["reversed"]:
            helix["path"] = [max_x, max_y, min_x, min_y]
        else:
            helix["path"] = [min_x, min_y, max_x, max_y]
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

def vertex_in_polygon(point_x, point_y, polygon_coordinates):
    formatted_coordinates = []
    for vertex_index in range(0, len(polygon_coordinates), 2):
        x = polygon_coordinates[vertex_index]
        y = polygon_coordinates[vertex_index + 1]
        formatted_coordinates.append((x, y))
    point = Point(point_x, point_y)
    poly = Polygon(formatted_coordinates)
    return point.within(poly)

def interpolate_1d(x0, x1, t):
    return (1 - t) * x0 + t * x1

def interpolate(x0, y0, x1, y1, t):
    return interpolate_1d(x0, x1, t), interpolate_1d(y0, y1, t)

def parseSequentialResidueSequenceNumbers(residueSequenceNumbersString):
    residueSequenceNumberStrings = residueSequenceNumbersString.split()
    residueSequenceNumbers = []
    for residueSequenceNumberString in residueSequenceNumberStrings:
        residueSequenceNumbers.append(int(residueSequenceNumberString))
    if (len(residueSequenceNumbers) > 0):
        nextResidueSequenceNumber = residueSequenceNumbers[-1]
        sequentialResidueSequenceNumbers = [nextResidueSequenceNumber]
        for i in range(len(residueSequenceNumbers) - 2, -1, -1):
            # print (i)
            residueSequenceNumber = residueSequenceNumbers[i]
            if residueSequenceNumber != nextResidueSequenceNumber - 1:
                break
            sequentialResidueSequenceNumbers = [residueSequenceNumber] + sequentialResidueSequenceNumbers
            nextResidueSequenceNumber = residueSequenceNumber
        # print ()
        return sequentialResidueSequenceNumbers
    else:
        return []

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