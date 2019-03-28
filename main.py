import re
import string


class Gene:

    def __init__(self, start, end, name):
        self.start = start
        self.end = end
        self.name = name
        self.sequence = ""
        self.cds = []

    def __str__(self):
        return "name is %s, start is %s, end is %s" % (self.name, self.start, self.end)


def find_genes(file):
    genelist = []
    f = open(file, 'r+')
    readstuff = f.read()
    p = re.compile("gene\\s+\\d+..\\d+\\s*\/gene=\".+\"")
    result = p.findall(readstuff)
    for match in result:
        coords =  re.findall('\d+', match)
        start = coords[0]
        end = coords[1]
        name = match.split("\"")[1].split("\"")[0]
        new_gene = Gene(start,end,name)
        genelist.append(new_gene)
    return genelist


def find_cds(file):
    genelist = []
    f = open(file, 'r+')
    readstuff = f.read()
    p = re.compile("CDS\\s+\\d+..\\d+\\s*\/gene=\".+\"")
    p1 = re.compile("CDS\\s+\\d+..\\d+\\s*\/gene=\".+\"")
    result = p.findall(readstuff)
    print(result)
    for match in result:
        coords =  re.findall('\d+', match)
        start = coords[0]
        end = coords[1]
        name = match.split("\"")[1].split("\"")[0]
        new_gene = Gene(start,end,name)
        genelist.append(new_gene)
    return genelist


def read_fasta(file):
    f = open(file, 'r+')
    readstuff = f.read()
    genome = readstuff.split("complete genome")[1].translate(str.maketrans('', '', string.whitespace))
    return genome


def main():
    sequence = read_fasta("hiv1.fasta")
    genelist = find_genes("hiv1.gb")
    for genex in genelist:
        genex.sequence = sequence[int(genex.start)-1:int(genex.end)]

    cdslist = find_cds("hiv1.gb")
    for genex in cdslist:
        print(genex)


if __name__ == "__main__":
    main()


