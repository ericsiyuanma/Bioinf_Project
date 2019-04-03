import re
import string
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class Genome:
    def __init__(self, name, genes, sequence):
        self.name = name
        self.genes = genes
        self.sequence = sequence


class Gene:

    def __init__(self, name, gene_coords, cds_coords):
        self.name = name
        self.gene_coords = gene_coords
        self.cds_coords = cds_coords
        self.gene_sequence = ""
        self.cds_sequence_list = []
        self.size = 0

        self.gc_content_overall = 0
        self.gc_content_coding = 0
        self.gc_content_non_coding = 0

    def __str__(self):
        return "name is %s, gene coordinates are %s, cds coordinates is %s \n " \
               "sequence is %s" % (self.name, self.gene_coords, self.cds_coords, self.gene_sequence)

    def calculate_gc_contents(self):
        self.gc_content_overall = (self.gene_sequence.count('G') +
                                   self.gene_sequence.count('C')) / len(self.gene_sequence)
        self.size = int(len(self.gene_sequence))


def find_genes(file):
    genelist = []
    f = open(file, 'r+')
    readstuff = f.read()
    p = re.compile("gene\\s+\\d+\.\.\\d+\\s+.*\\s+\/locus_tag=\".+\"")
    result = p.findall(readstuff)
    for match in result:
        coords = re.findall('\d+', match)
        start = coords[0]
        end = coords[1]
        name = match.split("\"")[3].split("\"")[0]
        new_gene = {"name": name, "gene_coords_list": [int(start), int(end)]}
        genelist.append(new_gene)
    return genelist


def find_cds(file):
    cdslist = []
    f = open(file, 'r+')
    readstuff = f.read()
    p = re.compile("CDS\\s+\\d+\.\.\\d+\\s+.*\\s+\/locus_tag=\".+\"")
    p1 = re.compile("CDS\\s+join\(.+\)\\s+.*\\s+\/locus_tag=\".+\"")
    result = p.findall(readstuff)
    result1 = p1.findall(readstuff)
    for match in result:
        coords = re.findall('\d+', match)
        start = coords[0]
        end = coords[1]
        name = match.split("\"")[3].split("\"")[0]
        cdsobject = {"name": name, "cds_coords_list": [[int(start), int(end)]]}
        cdslist.append(cdsobject)

    for match1 in result1:
        coords = re.findall('\d+\.\.\d+', match1)
        pairlist = []
        for pair in coords:
            pairtemp = pair.split("..")
            pairobject = [int(pairtemp[0]), int(pairtemp[1])]
            pairlist.append(pairobject)

        name = match1.split("\"")[3].split("\"")[0]
        cdsobject = {"name": name, "cds_coords_list": pairlist}
        cdslist.append(cdsobject)
    return cdslist


def read_fasta(file):
    f = open(file, 'r+')
    readstuff = f.read()
    genome = readstuff.split("complete genome")[1].translate(str.maketrans('', '', string.whitespace))
    return genome


def map_gene_object(genelist, cdslist):
    gene_object_list = []
    for gene in genelist:

        key = gene['name']
        for cds in cdslist:
            if cds["name"] == key:
                geneobject = Gene(name=key, gene_coords=gene["gene_coords_list"], cds_coords=cds['cds_coords_list'])
                gene_object_list.append(geneobject)

    return gene_object_list


def draw_plot_something(diamonds):
    plt.clf()
    plt.scatter(x='overall_gc', y='size', data=diamonds, alpha=0.5)
    plt.show()


def main():
    genelist = find_genes("hiv1.gb")

    cdslist = find_cds("hiv1.gb")
    sequence = read_fasta("hiv1.fasta")

    gene_object_list = map_gene_object(genelist, cdslist)

    dflisttemp = []
    for genex in gene_object_list:
        genex.gene_sequence = sequence[genex.gene_coords[0] - 1:genex.gene_coords[1]]
        cds_seq_list = []
        for cds_coords in genex.cds_coords:
            cds_seq = sequence[cds_coords[0] - 1:cds_coords[1]]
            cds_seq_list.append(cds_seq)
        genex.cds_sequence_list = cds_seq_list
        genex.calculate_gc_contents()
        dflisttemp.append([genex.name, int(genex.size), genex.gc_content_overall])

    hiv1 = Genome("hiv1", gene_object_list, sequence=sequence)
    hiv1df = pd.DataFrame(np.array(dflisttemp), columns=['locus_tag', 'size', 'overall_gc'])
    hiv1dfsorted = hiv1df.sort_values(by=['size'])
    hiv1dfsorted['size'] = pd.to_numeric(hiv1dfsorted['size'])
    hiv1dfsorted['overall_gc'] = pd.to_numeric(hiv1dfsorted['overall_gc'])
    draw_plot_something(hiv1dfsorted)
    hiv1dfsorted.to_csv('hiv1.csv', encoding='utf-8', index=False)


if __name__ == "__main__":
    main()
