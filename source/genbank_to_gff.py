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
            "start":from_to[0].replace("<","").replace(">",""),
            "end":from_to[1].replace("<","").replace(">",""),
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

# Count all old_locus_tag, locus_tag, and gene to find non-unique tags
tag_counts = {"old_locus_tag":{}, "locus_tag":{}, "gene":{}}

for i in range(len(features)):
    # Only consider coding sequences
    if not features[i]['feature'] == "CDS":
        continue
    # Count old_locus_tag
    try:
        tag_counts['old_locus_tag'][features[i]['old_locus_tag']] += 1
    except KeyError:
        try:
            tag_counts['old_locus_tag'][features[i]['old_locus_tag']] = 1
        except KeyError:
            pass
    # Count locus_tag
    try:
        tag_counts['locus_tag'][features[i]['locus_tag']] += 1
    except KeyError:
        try:
            tag_counts['locus_tag'][features[i]['locus_tag']] = 1
        except KeyError:
            pass
    # Count gene
    try:
        tag_counts['gene'][features[i]['gene']] += 1
    except KeyError:
        try:
            tag_counts['gene'][features[i]['gene']] = 1
        except KeyError:
            pass

# Identify all non-unique old_locus_tag, locus_tag, and gene tags
non_uniq = {
    'old_locus_tag':set(filter(
        lambda x: tag_counts['old_locus_tag'][x] > 1,
        tag_counts['old_locus_tag']
    )),
    'locus_tag':set(filter(
        lambda x: tag_counts['locus_tag'][x] > 1,
        tag_counts['locus_tag']
    )),
    'gene':set(filter(
        lambda x: tag_counts['gene'][x] > 1,
        tag_counts['gene']
    ))
}

# Rename all features that are non-unique
def rename_feature(feature, feature_type):
    try:
        if feature[feature_type] in non_uniq[feature_type]:
            feature[feature_type] = "_".join([
                feature[feature_type], feature['sequence'],
                feature['start'], feature['end'], feature['strand']
            ])
    except KeyError:
        pass
    return feature

# Write features to GFF

outfile = open(outfile, 'w')

for i in range(len(features)):
    # Select feature
    feature = features[i]
    # Only consider coding sequences
    if not feature['feature'] == "CDS":
        continue
    # Rename non-unique tags
    for feature_type in ['old_locus_tag', 'locus_tag', 'gene']:
        feature = rename_feature(feature, feature_type)
    # Write column 1: Sequence
    output = feature['sequence'] + "\t"
    # Write column 2: Source
    output += 'Custom' + "\t"
    # Write column 3: Type
    output += 'gene' + "\t"
    # Write column 4: Start
    output += feature['start'] + "\t"
    # Write column 5: End
    output += feature['end'] + "\t"
    # Write column 6: Score
    output += '.' + "\t"
    # Write column 7: Strand
    output += feature['strand'] + "\t"
    # Write column 8: Frame
    output += '0' + "\t"
    # Write column 9: Attributes
    try:
        locus_tag = feature['old_locus_tag']
    except KeyError:
        try:
            locus_tag = feature['locus_tag']
        except KeyError:
            locus_tag = feature['gene']
    try:
        ID = feature['locus_tag']
    except KeyError:
        ID = feature['gene']
    locus_tag = "locus_tag=" + locus_tag
    ID = "ID=" + ID
    product = "product=" + feature['product'].replace(";", "_")
    output += ";".join([product, locus_tag, ID]) + "\n"
    junk = outfile.write(output)

outfile.close()
