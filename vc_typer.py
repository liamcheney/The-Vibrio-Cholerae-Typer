import argparse
import glob
import sys
from os import path
from Bio.Blast.Applications import NcbiblastnCommandline
import multiprocessing
import datetime
import progressbar
from operator import itemgetter
from time import sleep as sl
from Bio import SeqIO

#smaller functions
def check_databases_input(args):
    database_okay = True
    accepted_databases = ['sero','ctxB','tcpA','rstR','bio','sxt', 'ICE_check', 'seventh_check', 'species_check']
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
        database_list = ['sero', 'ctxB', 'tcpA', 'rstR', 'bio', 'sxt', 'ICE_check', 'seventh_check',
                              'species_check']
        return database_list

    else:
        input_database_list = args.databases.split(',')
        database_list = []

        for el in input_database_list:
            database_list.append(el)

        return database_list
def progess_use(args):
    number_of_genomes = len(list(glob.iglob(args.strains_directory + '/*.f*'))) + 1
    bar = progressbar.ProgressBar(maxval=number_of_genomes,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

    return bar
def query_length_gen(col):

    if int(col[9]) > int(col[8]):
        query_length = int(col[9]) - int(col[8]) + 1
    else:
        query_length = int(col[8]) - int(col[9]) + 1

    return query_length
def parseargs(set_wd):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-db", "--databases", default='All',
                        help="Databases to process. Options ctxB, tcpA, rstR, bio, sxt, \
                             ICE_check, seventh_check, species_check, \
                             For combinations use comma seperated. Eg. ctxB,tcpA \
                             Default is all databases.")
    parser.add_argument("-dir", "--strains_directory", required=True,
                        help="A directory of strains to analyse.")
    parser.add_argument("-r", "--return_all_blast", action='store_true',
                        help="Return all blast hits")
    parser.add_argument("-o", "--output_folder", default=set_wd + "/output/",
                        help="Output folder to save if not specified.")
    parser.add_argument("-f", "--overwrite", action='store_true',
                        help="Overwrite previous results output file.")
    parser.add_argument("-cov", "--coverage", default=95,
                        help="Minimum percentage of blast query length.")
    parser.add_argument("-len", "--length", default=95,
                        help="Minmum percetange gene can be missing nucleotides.")
    parser.add_argument("-t", "--threads", default=1,
                        help="Threads for computing.")



    args = parser.parse_args()

    return args

#blast searching
def blast_search(set_wd, args):

    print("Analaying strains.")
    current_blast_db_path = set_wd + "/blast_db/alleles"
    results_dict = run_genomes(current_blast_db_path, args)

    return results_dict
def run_genomes(current_blast_db_path, args):
    print("Blasting database alleles against genomes.")
    bar = progess_use(args)
    bar.start()
    genome_count = 0
    blast_results_dict = {}
    for strain in glob.iglob(args.strains_directory + '/*.f*'):
        bar.update(genome_count)
        genome_count = genome_count + 1
        accession = strain.split('/')[-1].split('.')[0]
        blast_result = blast_input_against_db(strain, current_blast_db_path, args)
        blast_results_dict[accession] = blast_result
    bar.finish()

    results_dict = blast_result_filtering(blast_results_dict, args)

    return results_dict
def blast_input_against_db(strain, current_blast_db_path, args):
    #blast query database against genome

    blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=current_blast_db_path, outfmt=6, num_threads=args.threads)
    out, err = blast_string()
    return out
def blast_result_filtering(blast_results_dict, args):

    print("Processing blast outputs.")
    # read in file of allele lengths
    set_wd = path.dirname(path.abspath(__file__))
    gene_size_dict = input_allele_lengths(set_wd, args)

    #store the blast results and create a dict for all results
    returning_results = {}

    #get list of databases
    databases = databases_list(args)

    #go over each selected databases
    for item in databases:

        if item == 'bio':
            returning_results[item] = biotype_filter_results(blast_results_dict, gene_size_dict)
        else:
            returning_results[item] = remaining_filter_results(blast_results_dict, item, gene_size_dict,args)

    return returning_results

#handling biotyping
def biotype_filter_results(blast_results_dict, gene_size_dict):

    return_dict = {}

    #go over each strain result for the database
    for strain in blast_results_dict.keys():
        #remove non relevant results and organise list
        format_blast_list = biotype_format_blast_output(blast_results_dict[strain])

        # process ctxB, rstR and tcpA blast results
        result_list = biotype_blast_filter(format_blast_list, gene_size_dict)

        #determine the biotype
        return_dict[strain] = biotype_selector(result_list)

    return return_dict
def biotype_format_blast_output(results_list):

    want_database_list = ['ctxB', 'tcpA', 'rstR']

    #split the blast result
    blast_hits = results_list.split('\n')[:-1] #remove last new line char

    keep_list = []
    for number in range(0, len(blast_hits), 1):
        for item in want_database_list:
            if item in blast_hits[number].split('\t')[1].split('_')[0]:
                keep_list.append(blast_hits[number])

    #organise blast hits by precen of id
    keep_list = sorted(keep_list, key=itemgetter(2,11))

    return keep_list
def biotype_blast_filter(format_blast_list, gene_size_dict):

    result_list = []
    allele_type = []
    allele_asigned = []

    for element in format_blast_list:
        col = element.split('\t')
        allele_length = gene_size_dict[col[1]]
        query_length = int(col[3]) - int(col[4])
        if col[1].split('_')[0] not in allele_asigned:

            # ##check if hit is exact match
            if query_length == allele_length and float(col[2]) > 99.0:
                result_list.append(element)
                allele_type.append(col[1])
                allele_asigned.append(col[1].split('_')[0])

            elif query_length >= (0.99 * allele_length) and float(col[2]) > 90.0 and 'ctxB' not in col[1]:
                result_list.append(element)
                allele_type.append(col[1])
                allele_asigned.append(col[1].split('_')[0])

    result_list.insert(0, '/'.join(allele_type))
    final_list = '\t'.join(result_list)
    return final_list
def biotype_selector(result_list):

        # create function to check what biotype the strain is
        biotype_results = []

        ##possible combiniations of alleles for each biotype
        bio_classical = ['ctxB1', 'rstR_cla', 'tcpA_cla_WT']
        bio_eltor = ['ctxB3', 'rstR_el', 'tcpA_el_WT', 'tcpA_el_A226']
        bio_atypical = ['ctxB1', 'ctxB7', 'rstR_cla', 'rstR_el', 'tcpA_el_A226', 'tcpA_el_WT']
        bio_mozambique = ['ctxB1', 'tcpA_el_WT', 'rstR_cla']

        col = result_list.split('\t')
        alleles_list = col[0].split('/')

        ##check strain biotype by comparing alleles against bio_lists
        if set(alleles_list).issubset(bio_classical):
            biotype_results = 'classical' + '\t' + result_list

        elif set(alleles_list).issubset(bio_eltor):
            biotype_results = 'eltor' + '\t' + result_list

        elif set(alleles_list).issubset(bio_mozambique):
            biotype_results = 'mozambique' + '\t' + result_list

        elif set(alleles_list).issubset(bio_atypical):
            biotype_results = 'atypical' + '\t' + result_list

        return biotype_results

#handling all alleles
def remaining_filter_results(blast_results_dict, item, gene_size_dict,args):

    return_dict = {}
    #go over each strain result for the database
    for strain in blast_results_dict.keys():
        #remove non relevant results and organise list
        format_blast_list = remaining_format_blast_output(blast_results_dict[strain], item)

        #check the blast results
        return_dict[strain] = remaining_blast_filter(format_blast_list, item, gene_size_dict,args)

    return return_dict
def remaining_format_blast_output(results_list, item):

    #split the blast result
    blast_hits = results_list.split('\n')[:-1] #remove last new line char

    keep_list = []
    for number in range(0, len(blast_hits), 1):
        if item in blast_hits[number].split('\t')[1]:
            keep_list.append(blast_hits[number])

    return keep_list
def remaining_blast_filter(format_blast_list, item, gene_size_dict, args):

    result_list = []
    format_blast_list = [x.split('\t') for x in format_blast_list]
    format_blast_list.sort(key=lambda x: (float(x[2]),int(x[3])), reverse=True)

    for element in format_blast_list:
        allele_length = gene_size_dict[element[1]]
        min_return_allele_length = int(args.length / 100 * allele_length)
        min_return_coverage = int(args.coverage / 100 * allele_length)

        if args.return_all_blast:
            if (min_return_allele_length <= int(element[3])) and (min_return_coverage <= int(element[3])):
                keep = ','.join(element)
                result_list.append(keep)

        if not args.return_all_blast:
            if (min_return_allele_length <= int(element[3])) and (min_return_coverage <= int(element[3])):
                top_hit = format_blast_list[0]
                keep = ','.join(top_hit)
                result_list.append(keep)
                break

    return result_list

##writing out
def write_out(results_dict, args):
    print("Writing Out Results")

    if args.overwrite:
        # create output file and fill
        with open(args.output_folder + '/blast_results.csv', 'w') as out:
            write_out_iterator(results_dict, out, args)

    else:
        #create new output file and fill
        time_stamp = '_'.join('_'.join(str(datetime.datetime.now()).split('.')[0].split(' ')).split(':'))
        with open(args.output_folder + '/blast_results_' + time_stamp + '.csv','w') as out:
            write_out_iterator(results_dict, out, args)
def write_out_iterator(results_dict, out, args):
    # column names
    headers_list = ["Accession", "query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                    "qend", "sstart", "send", "evalue", "bitscore"]
    biotype_headers_list = ["Accession", 'biotype'] + headers_list[1:]

    for key in results_dict:

        #adding column headers
        out.write(key + '\n')
        for col_head in headers_list:
            out.write(col_head + ',')
        out.write('\n')
        print(results_dict[key])

def main():

    set_wd = path.dirname(path.abspath(__file__))
    args = parseargs(set_wd)
    database_input = check_databases_input(args)
    if database_input:
        results_dict = blast_search(set_wd, args)
        write_out(results_dict, args)
    else:
        print('One of following database(s) is incorrect : ' + str(args.databases))
        sys.exit()

if __name__ == '__main__':
    main()

#TODO make multiple databases from other tool
#TODO add abricate for antibiotic resistance
