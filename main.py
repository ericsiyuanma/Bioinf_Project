import re


def read_parse_sequence():
    f = open('sequence.txt', 'r+')
    readstuff = f.read()
    result = ''.join([i for i in readstuff if i.isalpha()]).lower()
    return result


def read_parse_mrna():
    f1 = open("mRNA.txt", "r+")
    readmrna = f1.read()
    result1 = re.findall('\d+', readmrna)
    return result1


def read_parse_cds():
    f2 = open("cds.txt", "r+")
    readcds = f2.read()
    result2 = re.findall('\d+', readcds)
    return result2


def make_intron(result, result1):
    count = 0
    # go through all exon sequences
    for x in range(1, int(len(result1) / 2)+1):
        # exon sequences make uppercase
        for y in range(int(result1[count]) -1, int(result1[count+1])):
            result = result[:y] + result[y].swapcase() + result[y+1:]
        count += 2

    finalfile = open("finished.txt", "w")
    finalfile.write(result)
    return(result)


def find_three_utr(gene):
    """
    These are hardcoded for now

    :param gene:
    :return:
    """
    three = gene[79806:81188]
    finalfile1 = open("three.txt", "w")
    finalfile1.write(three)
    print(three)


def find_five_utr(gene):
    """
    Hardcoded for now
    :param gene:
    :return:
    """
    five = gene[:213] + gene[1368:1386]
    finalfile2 = open("five.txt", "w")
    finalfile2.write(five)
    print(five)



sequence = read_parse_sequence()
mrna = read_parse_mrna()
cds = read_parse_cds()
gene = make_intron(sequence, mrna)
find_three_utr(gene)
find_five_utr(gene)







