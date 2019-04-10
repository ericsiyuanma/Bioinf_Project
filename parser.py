from lib_parser import parse_gene;


def main():
    parse_gene("hiv1.gb", "hiv1.fasta", "hiv1")
    parse_gene("herpes.gb", "herpes.fasta", "herpes")
    parse_gene("pandora.gb", "pandora.fasta", "pandora")
    parse_gene("phage.gb", "phage.fasta", "phage")


if __name__ == "__main__":
    main()
