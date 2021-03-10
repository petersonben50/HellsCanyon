import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none' # Editable SVG text
CSV = sys.argv[1]
PDF = sys.argv[2]

def add_polygon(coordinates):
    # Plots a gene as an arrow with offsets from the end of translation
    # Reference: https://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
    x1,x2,y,h,color,strand = coordinates
    edge_color='black'
    edge_width=1
    alpha=1
    # Polygon for ORF
    points = [[x1, y-h/2],   # left bottom
              [x1, y],       # left center
              [x1, y + h/2], # left top
              [x2, y + h/2], # right top
              [x2, y],       # right center
              [x2, y - h/2]] # right bottom
    # Short genes require different offset length
    x_offset = h/3
    if x_offset > abs(x1-x2):
        x_offset = abs(x1-x2)
    # "Point" the polygon to indicate direction
    if strand == '+':
        # Arrow on right
        points[3][0] = x2 - x_offset
        points[5][0] = x2 - x_offset
    else:
        # Arrow on left
        points[0][0] = x1 + x_offset
        points[2][0] = x1 + x_offset
    # Plot parameters
    polygon = plt.Polygon(points, fc=color, edgecolor=edge_color, linewidth=edge_width, alpha=alpha)
    return plt.gca().add_patch(polygon)

def set_relative_coordinates(df, central_gene):
    # Calculate proximity to reference gene for other genes on scaffold
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    startcoord_ref = int(df.at[central_gene,'start'])
    endcoord_ref = int(df.at[central_gene,'end'])
    # Define baseline values for hit
    df.at[central_gene,'Rel Start Coord'] = 0
    df.at[central_gene,'Rel strand'] = '+'
    # Depending on strand (top=+, bottom=-)
    sign = df.at[central_gene,'strand']
    if sign == "+": # Subtract all coords from initial start coord
        vector = +1.0
        strand_key = {"+" : "+", "-" : "-"}
        df['Rel Start Coord'] = df['start'] - startcoord_ref
        df['Rel End Coord'] = (df['end'] - startcoord_ref)
        df['Rel strand'] = df['strand'].map(strand_key)
    elif sign == "-": # -1 * Subtract all coords from initial end coord (true start)
        vector = -1.0
        strand_key = {"-" : "+", "+" : "-"} # Flip strands
        df['Rel Start Coord'] = vector*(df['end'] - endcoord_ref)
        df['Rel End Coord'] = vector*(df['start'] - endcoord_ref)
        df['Rel strand'] = df['strand'].map(strand_key)
    # Assign Rel # below or above gene
    df = df.sort_values(by='Rel Start Coord')
    df.loc[(df['Rel Start Coord'] > 0), 'Rel #'] = range(1, 1 + len(df['Rel Start Coord'][(df['Rel Start Coord'] > 0)].tolist()), 1)
    df.loc[(df['Rel Start Coord'] < 0), 'Rel #'] = list(reversed(range(-1, -1 - len(df['Rel Start Coord'][(df['Rel Start Coord'] < 0)].tolist()), -1)))
    df.at[central_gene,'Rel #'] = 0
    return df

def plot_gene_clusters(df):
    # Spacing of plot
    y = 0 # initial y-coordinate
    h = 500 # height of genes
    spacing_vertical = 5 # aspect ratio between scaffolds
    # Plot genes centered on these genes
    focal_genes = df[(df['focalGene'] == True)].sort_values(by='scaffoldID').index.tolist()
    for FOCAL_GENE in focal_genes:
        # Assign new y-coordinate
        y = y - h * spacing_vertical
        xs = []
        scaffold = df.at[FOCAL_GENE,'scaffoldID']
        # Set relative coordinates for subset
        scaf = df[df['scaffoldID'] == scaffold].copy()
        scaf = set_relative_coordinates(scaf, FOCAL_GENE)
        # Plot ea. gene in scaffold
        genes_in_scaffold = scaf.index.tolist()
        for GENE in genes_in_scaffold:
            # Gene
            x1 = scaf.at[GENE,'Rel Start Coord']
            x2 = scaf.at[GENE,'Rel End Coord']
            strand = scaf.at[GENE,'Rel strand']
            color = scaf.at[GENE,'color']
            if x1 >= -10000 and x2 <= 10000:
                coordinates = [x1, x2, y, h, color, strand]
                xs.append(x2)
                # Label
                label_x = scaf.at[GENE,'Rel Start Coord'] + 0.1 * abs(scaf.at[GENE,'Rel Start Coord'] - scaf.at[GENE,'Rel End Coord'])
                label_y = y + h
                label_text = scaf.at[GENE,'geneName']
                #plt.text(label_x, label_y, label_text, size=8, rotation=45, verticalalignment='bottom')
                add_polygon(coordinates)
        # Plot genome text
        genome = scaf.at[FOCAL_GENE,'scaffoldID']
        plt.text(10000, y, genome, size=10, verticalalignment='center')
    plt.axis('scaled')
    plt.axis('off')
    return plt.savefig(PDF)

# Parse command line arguments
#parser = argparse.ArgumentParser(description='Plot gene clusters from formatted CSV file')
#parser.add_argument('--csv', dest='csv', help='CSV containing gene information (see example)')
#parser.add_argument('--output', dest='svg', help='File name for SVG output')
#args = parser.parse_args()

#CSV = args.csv
#SVG = args.svg

# Upload and format CSV
df = pd.read_csv(CSV)
df = df.drop_duplicates('geneID') # Remove duplicates
df = df.set_index('geneID') # Set index
df['color'] = df['color'].fillna('white') # Provide color if none provided

# Plot
plt.figure(figsize=(10,10))
plot_gene_clusters(df)
