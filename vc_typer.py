import argparse
import glob
import sys
from os import path
from Bio.Blast.Applications import NcbiblastnCommandline
import datetime
from time import sleep as sl
from Bio import SeqIO

def run_main(set_wd, args):
    ##set blast_path, database list
    blast_db_path = set_blast_path(set_wd)
    databases = databases_list(args)

    ##output during analysis
    print("Blasting database alleles against genomes.")
    all = len([x for x in glob.iglob(args.strains_directory + '/*.f*')])
    genome_count = 1

    ##blast results
    blast_result_dict = {}
    for strain in glob.iglob(args.strains_directory + '/*.f*'):
        accession = strain.split('/')[-1].split('.')[0]

        # get blast results
        blast_result = blast_input_against_db(strain, blast_db_path, args)
        blast_result = blast_result_filtering(blast_result, databases, set_wd, args)
        blast_result_dict[accession] = blast_result

        print_list = []
        for i in blast_result:
            if len(blast_result[i]) == 0:
                print_list.append('')

            else:
                print_list.append(blast_result[i].split(',')[1])
        print(accession, *print_list, sep=',')

    ##writing out
    write_out(blast_result_dict, args)

# blast searching
def set_blast_path(set_wd):

    print("Analaying strains.")
    current_blast_db_path = set_wd + "/blast_db/alleles"

    return current_blast_db_path
def blast_input_against_db(strain, blast_db_path, args):

    #blast query database against genome
    blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=blast_db_path, outfmt=6, num_threads=args.threads)
    out, err = blast_string()
    out = out.splitlines()
    return out
def blast_result_filtering(blast_results_dict, databases, set_wd, args):

    # read in file of allele lengths
    gene_size_dict = input_allele_lengths(set_wd, args)

    #seperate blast results based on database
    sep_blast_results = seperate_blast_results(blast_results_dict, databases)

    #select the hits from each databaes
    result_dict = selecting_hits(sep_blast_results, gene_size_dict, args)

    return result_dict

#handling alleles
def seperate_blast_results(blast_results_dict, databases):
    result_dict = {}
    for item in databases:
        result_dict[item] = []

    for key in result_dict:
        for line in blast_results_dict:
            col = line.split()
            allele = col[1]
            if key in allele:
                result_dict[key].append(line)

    return result_dict
def selecting_hits(sep_blast_results, gene_size_dict, args):

    result_dict = {}
    for key, value in sep_blast_results.items():
        result_dict[key] = []
        format_blast_list = [x.split('\t') for x in value]
        format_blast_list.sort(key=lambda x: (float(x[2]),int(x[3])), reverse=True)

        for element in format_blast_list:
            allele_length = gene_size_dict[element[1]]
            min_return_allele_length = int(args.length / 100 * allele_length)
            min_return_coverage = int(args.coverage / 100 * allele_length)

            if (min_return_allele_length <= int(element[3])) and (min_return_coverage <= int(element[3])):
                top_hit = format_blast_list[0]
                keep = ','.join(top_hit)
                result_dict[key] = keep
                break

    return result_dict

##writing out
def write_out(results_dict, args):
    print("Writing Out Results")

    if args.overwrite:
        # create output file and fill
        with open(args.output_folder + '/blast_results.csv', 'w') as out:
            write_out_iterator(results_dict, out, args)
        with open(args.output_folder + '/simple_blast_results.csv', 'w') as out:
            simple_write_out_iterator(results_dict, out, args)

    else:
        #create new output file and fill
        time_stamp = '_'.join('_'.join(str(datetime.datetime.now()).split('.')[0].split(' ')).split(':'))
        with open(args.output_folder + '/blast_results_' + time_stamp + '.csv','w') as out:
            write_out_iterator(results_dict, out, args)
        with open(args.output_folder + '/simple_blast_results.csv', 'w') as out:
            simple_write_out_iterator(results_dict, out, args)
def write_out_iterator(results_dict, out, args):
    # column names
    headers_list = ["Accession", "query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                    "qend", "sstart", "send", "evalue", "bitscore"]

    #writing out all outputs blast results
    for key in results_dict:
        #adding column headers
        out.write(key + '\n')
        for col_head in headers_list:
            out.write(col_head + ',')
        out.write('\n')

        for strain in results_dict[key]:
            out.write(strain + ',')
            for item in results_dict[key][strain]:
                out.write(item)
            out.write('\n')

        out.write('\n')
def simple_write_out_iterator(results_dict, out, args):
    # writing out all outputs blast results
    for strain in results_dict:
        out.write(strain + ',')
        for gene in results_dict[strain]:
            if len(results_dict[strain][gene]) == 0:
                out.write(',')
            else:
                out.write(results_dict[strain][gene].split(',')[1] + ',')
        out.write('\n')

#smaller functions
def check_databases_input(args):
    database_okay = True
    accepted_databases = ['All','sero','ctxB','tcpA','rstR','sxt', 'ICE_setD', 'seventh', 'species']
    in_split = args.databases.split(',')
    for i in in_split:
        if i not in accepted_databases:
            database_okay = False
    return database_okay
def input_allele_lengths(set_wd, args):

    gene_size_dict = {}
    for record in SeqIO.parse(set_wd + "/blast_db/alleles.fasta", "fasta"):
        gene = record.id
        length = len(str(record.seq))
        gene_size_dict[gene] = length

    return gene_size_dict
def databases_list(args):

    if args.databases == 'All':
        database_list = ['sero', 'ctxB', 'tcpA', 'rstR', 'sxt', 'ICE_setD', 'seventh',
                              'species']
        return database_list

    else:
        input_database_list = args.databases.split(',')
        database_list = []

        for el in input_database_list:
            database_list.append(el)

        return database_list
def query_length_gen(col):

    if int(col[9]) > int(col[8]):
        query_length = int(col[9]) - int(col[8]) + 1
    else:
        query_length = int(col[8]) - int(col[9]) + 1

    return query_length
def parseargs():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-db", "--databases", default='All',
                        help="Databases to process. Options ctxB, tcpA, rstR, sxt, \
                             ICE_check, seventh_check, species_check, \
                             For combinations use comma seperated. Eg. ctxB,tcpA \
                             Default is all databases.")
    parser.add_argument("-dir", "--strains_directory", required=True,
                        help="A directory of strains to analyse.")
    parser.add_argument("-o", "--output_folder",
                        help="Output folder to save if not specified.")
    parser.add_argument("-f", "--overwrite", action='store_true',
                        help="Overwrite previous results output file.")
    parser.add_argument("-cov", "--coverage", default=90, type=int,
                        help="Minimum percentage of blast query length.")
    parser.add_argument("-len", "--length", default=90, type=int,
                        help="Minmum percetange gene can be missing nucleotides.")
    parser.add_argument("-t", "--threads", default=1,
                        help="Threads for computing.")



    args = parser.parse_args()

    return args

def main():

    set_wd = path.dirname(path.abspath(__file__))
    args = parseargs()
    database_input = check_databases_input(args)
    if database_input:
        run_main(set_wd, args)
    else:
        print('One of following database(s) is incorrect : ' + str(args.databases))
        sys.exit()

if __name__ == '__main__':
    main()