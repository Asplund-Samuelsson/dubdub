#!/usr/bin/env python3

import sys

# Read input and output files
infile = sys.argv[1]
outfile = sys.argv[2]

# Create list to store features in
features = []

# Iterate over lines in input Genbank file
for line in open(infile).readlines():
    # Store current sequence name
    if line.startswith("LOCUS"):
        sequence = line.split()[1]
    # Store current feature properties
    if not line.startswith(" "*6) and line.startswith(" "*5):
        # Remove leading and trailing whitespace and split
        line = line.strip().split()
        # Determine feature type and initialize feature in features list
        features.append({"sequence":sequence, "feature":line[0]})
        # Determine feature strand
        if "complement" in line[1]:
            features[-1]["strand"] = "-"
        else:
            features[-1]["strand"] = "+"
        # Remove junk from range
        line[1] = line[1].replace("join", "").replace("complement", "")
        line[1] = line[1].replace("(", "").replace(")", "")
        # Determine feature range
        range_values = line[1].replace(",", "..").split("..")
        from_to = [range_values[0], range_values[-1]]
        # Fix for "join" ranges
        if len(range_values) == 4:
            if range_values[0] < range_values[3]:
                from_to = [range_values[2], range_values[1]]
        # Store initial feature attributes
        features[-1].update({
            "start":from_to[0], "end":from_to[1],
            "pseudo":False, "product":""
        })
        # Skip features with "order"
        order = "order" in line[1]
    # Determine attributes of interest
    elif line.startswith(" "*21):
        # Skip features with "order"
        if order:
            continue
        # Remove leading and trailing whitespace
        line = line.strip()
        # Check for current attribute
        if line.startswith("/"):
            line = line.lstrip("/").split("=", maxsplit=1)
            attribute = line[0]
            # Store attribute value
            if len(line) > 1:
                features[-1][attribute] = line[1].strip('"')
            else:
                features[-1][attribute] = True
        else:
            # Continue adding to value from following rows
            if not attribute == "translation":
                features[-1][attribute] += " "
            features[-1][attribute] += line.split('=', maxsplit=1)[0].strip('"')

# Write features to GFF

outfile = open(outfile, 'w')

for i in range(len(features)):
    # Only consider coding sequences
    if not features[i]['feature'] == "CDS":
        continue
    # Write column 1: Sequence
    output = features[i]['sequence'] + "\t"
    # Write column 2: Source
    output += 'Custom' + "\t"
    # Write column 3: Type
    output += 'gene' + "\t"
    # Write column 4: Start
    output += features[i]['start'] + "\t"
    # Write column 5: End
    output += features[i]['end'] + "\t"
    # Write column 6: Score
    output += '.' + "\t"
    # Write column 7: Strand
    output += features[i]['strand'] + "\t"
    # Write column 8: Frame
    output += '0' + "\t"
    # Write column 9: Attributes
    try:
        locus_tag = features[i]['old_locus_tag']
    except KeyError:
        try:
            locus_tag = features[i]['locus_tag']
        except KeyError:
            locus_tag = features[i]['gene']
    try:
        ID = features[i]['locus_tag']
    except KeyError:
        ID = features[i]['gene']
    locus_tag = "locus_tag=" + locus_tag
    ID = "ID=" + ID
    product = "product=" + features[i]['product'].replace(";", "_")
    output += "\t" + ";".join([product, locus_tag, ID]) + "\n"
    junk = outfile.write(output)

outfile.close()
