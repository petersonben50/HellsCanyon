
import os
import sys
import pandas as pd

gffName = sys.argv[1]
gene = sys.argv[2]
output = sys.argv[3]


listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
df = pd.read_csv(gffName, sep = '\t', names = listOfNames)

outputDF = pd.DataFrame(columns = ['geneID', 'scaffoldID', 'start', 'end', 'strand', 'focalGene', 'color', 'geneName', 'kofam'])

def get_element(my_list, position):
    return my_list[position]

outputDF['geneID'] = df.attributes.str.split(';').apply(get_element, position=0).str.split('=').apply(get_element, position=1)
outputDF['scaffoldID'] = df['sequence']
outputDF['start'] = df['start']
outputDF['end'] = df['end']
outputDF['strand'] = df['strand']
outputDF['focalGene'] = 'FALSE'
outputDF['color'] = 'NA'
#outputDF['locusTag'] = df.attributes.str.split(';').apply(get_element, position=1).str.split('=').apply(get_element, position=1)
outputDF['geneName'] = df.attributes.str.split(';').apply(get_element, position=13).str.split('=').apply(get_element, position=1).str.rsplit('-', 1).apply(get_element, position=0)
outputDF['kofam'] = df.attributes.str.split(';').apply(get_element, position=13).str.split('=').apply(get_element, position=1).str.rsplit('-', 1).apply(get_element, position=1)

# Identify focal gene and make it blue
genekey = gene.split("_", 1)[1].lstrip('0')
outputDF.loc[(outputDF['geneID'] == genekey), 'color'] = 'blue'
outputDF.loc[(outputDF['geneID'] == genekey), 'focalGene'] = 'True'

# Add colors for other genes
outputDF.loc[(outputDF['kofam'] == "K00520"), 'color'] = 'red'

# Transport proteins are green
outputDF.loc[(outputDF['kofam'] == "K08363"), 'color'] = 'green'
outputDF.loc[(outputDF['kofam'] == "K19058"), 'color'] = 'green'
outputDF.loc[(outputDF['kofam'] == "K19059"), 'color'] = 'green'

# Binding protein is yellow
outputDF.loc[(outputDF['kofam'] == "K08364"), 'color'] = 'yellow'

# MerR is black
outputDF.loc[(outputDF['kofam'] == "K08365"), 'color'] = 'black'
outputDF.loc[(outputDF['kofam'] == "K19057"), 'color'] = 'black'

# Output to csv
outputDF.to_csv(output, index = False)
