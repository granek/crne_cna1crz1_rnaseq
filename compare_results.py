import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Determine FASTQ file encoding")
    # parser.add_argument("OUTFILE", help="Name to save to.  If ends in '.gz', output will be gzipped")
    parser.add_argument("FILE1", type=file,
                        help="?????")
    parser.add_argument("FILE2", type=file,
                        help="?????")
    args = parser.parse_args()


    file1_reader = csv.reader(args.FILE1)
    file1_genes =  [row[0] for row in file1_reader if row[0]]

    file2_reader = csv.reader(args.FILE2)
    file2_genes =  [row[0] for row in file2_reader if row[0]]
    print file1_genes
    print file2_genes

if __name__ == "__main__":
    main()
