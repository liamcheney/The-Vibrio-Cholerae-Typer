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

#smaller functions
def check_databases_input(args):
    database_okay = True
    accepted_databases = ['sero','ctxB','tcpA','rstR','bio','sxt', 'ICE_check']
    in_split = args.databases.split(',')
    for i in in_split:
        if i not in accepted_databases:
            database_okay = False
    return database_okay
def input_allele_lengths(set_wd, args):

    gene_size_dict = {}
    for line in open(set_wd + "/allele_lengths.txt",'r').read().splitlines():
        col = line.split('\t')
        gene_size_dict[col[0]] = int(col[1])
    return gene_size_dict
def databases_list(args):

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
    parser.add_argument("-db", "--databases", required=True,
                        help="Databases to process. Options ctxB, tcpA, rstR, sero, bio, sxt, ICE_check. \
                             For combinations use comma seperated. Eg. ctxB,tcpA")
    parser.add_argument("-dir", "--strains_directory", required=True,
                        help="A directory of strains to analyse.")
    parser.add_argument("-r", "--blast_return", default=1, type=int,
                        help="1: return only top blast hits, 2: return all blast hits")
    parser.add_argument("-o", "--output_folder", default=set_wd + "/output/",
                        help="Output folder to save if not specified.")
    parser.add_argument("-f", "--overwrite", default=False,
                        help="Overwride previous results output file.")


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

    cpus = multiprocessing.cpu_count()

    #blast query database against genome
    blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=current_blast_db_path, outfmt=6, perc_identity=75, num_threads=cpus)
    out, err = blast_string()

    return out
def blast_result_filtering(blast_results_dict, args):

    print("Processing blast outputs.")
    # read in file of allele lengths
    set_wd = path.dirname(path.abspath(__file__))
    gene_size_dict = input_allele_lengths(set_wd, args)

    #handle the blast results and create a dict for all results
    returning_results = {}

    #get list of databases
    databases = databases_list(args)

    #go over each selected databases
    for item in databases:

        if item == 'bio':
            returning_results[item] = biotype_filter_results(blast_results_dict, gene_size_dict)

        elif item == 'sero':
            returning_results[item] = serogroup_filter_results(blast_results_dict,item, gene_size_dict)

        elif item == 'sxt':
            returning_results[item] = sxt_filter_results(blast_results_dict, item, gene_size_dict)

        elif item == 'ICE_check':
            returning_results[item] = ICE_check_filter_results(blast_results_dict, item, gene_size_dict)

        else:
            returning_results[item] = remaining_filter_results(blast_results_dict, item, gene_size_dict)

    return returning_results


#handling SXT
def sxt_filter_results(blast_results_dict, item, gene_size_dict):
    return_dict = {}

    # go over each strain result for the database
    for strain in blast_results_dict.keys():
        # remove non relevant results and organise list
        format_blast_list = sxt_format_blast_output(blast_results_dict[strain], item)

        # check the blast results
        return_dict[strain] = sxt_blast_filter(format_blast_list, item, gene_size_dict)

    return return_dict

    return
def sxt_format_blast_output(results_list, item):

    #split the blast result
    blast_hits = results_list.split('\n')[:-1] #remove last new line char

    keep_list = []
    for number in range(0, len(blast_hits), 1):
        if item in blast_hits[number].split('\t')[1]:
            keep_list.append(blast_hits[number])

    #organise blast hits by precen of id
    keep_list = sorted(keep_list, key=itemgetter(2,11))

    return keep_list
def sxt_blast_filter(format_blast_list, item, gene_size_dict):

    result_list = []
    allele_type = []
    allele_asigned = []

    sxt_group_dict = {'sxt_ICEVchInd6':'group 4','sxt_MO10':'first SXT', 'sxt_ICEVchBan5':'group 1', 'sxt_ICEVchBan9':'group 2', 'sxt_ICEVchInd4':'group 3'}

    for element in format_blast_list:
        col = element.split('\t')
        allele_length = gene_size_dict[col[1]]
        query_length = int(col[3]) - int(col[4])
        if col[1].split('_')[0] not in allele_asigned:

            # ##check if hit is exact match
            if query_length == allele_length and float(col[2]) > 99.0:
                result_list.append(element)
                allele_type.append(sxt_group_dict['sxt_' + col[1].split('_')[1]])
                allele_asigned.append(col[1])

    result_list.insert(0, '/'.join(allele_type))
    final_list = '\t'.join(result_list)
    return final_list

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

#checking ICEs presence
def ICE_check_filter_results(blast_results_dict, item, gene_size_dict):
    return_dict = {}

    # go over each strain result for the database
    for strain in blast_results_dict.keys():
        # remove non relevant results and organise list
        format_blast_list = ICE_check_format_blast_output(blast_results_dict[strain])

        # check the blast results
        return_dict[strain] = ICE_check_blast_filter(format_blast_list, item, gene_size_dict)

    return return_dict

def ICE_check_format_blast_output(results_list):

    # split the blast result
    blast_hits = results_list.split('\n')[:-1]  # remove last new line char

    keep_list = []
    for number in range(0, len(blast_hits), 1):
        if 'setC' in blast_hits[number].split('\t')[1] or 'setD' in blast_hits[number].split('\t')[1]:
            keep_list.append(blast_hits[number])

    # organise blast hits by precen of id
    keep_list = sorted(keep_list, key=itemgetter(2, 11))

    return keep_list
def ICE_check_blast_filter(format_blast_list, item, gene_size_dict):
    result_list = []
    allele_type = []
    allele_asigned = []

    for element in format_blast_list:
        col = element.split('\t')
        allele_length = gene_size_dict[col[1]]
        query_length = int(col[3]) - int(col[4])
        allele = col[1].split('_')[0]
        if allele not in allele_asigned:

            # ##check if hit is exact match
            if query_length == allele_length and float(col[2]) > 99.0:
                result_list.append(element)
                allele_type.append(col[1].split('_')[-1])
                allele_asigned.append(allele)

    result_list.insert(0, '/'.join(allele_type))
    final_list = '\t'.join(result_list)
    return final_list

#handling serogroup
def serogroup_filter_results(blast_results_dict, args):

            # elif key == 'sero' and query_length > (0.95 * allele_length) and float(col[2]) > 95.0:
            #     result_list.append(element)
            #     allele_type.append(col[1])

    #TODO add when have multiple 100% fragments
    return
#not finished didnt need at the time when updating program

#handling all ctxB, tcpA and rstR
def remaining_filter_results(blast_results_dict, item, gene_size_dict):

    return_dict = {}
    #go over each strain result for the database
    for strain in blast_results_dict.keys():
        #remove non relevant results and organise list
        format_blast_list = remaining_format_blast_output(blast_results_dict[strain], item)

        #check the blast results
        return_dict[strain] = remaining_blast_filter(format_blast_list, item, gene_size_dict)

    return return_dict
def remaining_format_blast_output(results_list, item):

    #split the blast result
    blast_hits = results_list.split('\n')[:-1] #remove last new line char

    keep_list = []
    for number in range(0, len(blast_hits), 1):
        if item in blast_hits[number].split('\t')[1]:
            keep_list.append(blast_hits[number])

    #organise blast hits by precen of id
    keep_list = sorted(keep_list, key=itemgetter(2,11))

    return keep_list
def remaining_blast_filter(format_blast_list, item, gene_size_dict):

    result_list = []
    allele_type = []
    allele_asigned = []

    for element in format_blast_list:
        col = element.split('\t')
        allele_length = gene_size_dict[col[1]]
        query_length = int(col[3]) - int(col[4])
        if item in col[1] and col[1].split('_')[0] not in allele_asigned:

            # ##check if hit is exact match
            if query_length == allele_length and float(col[2]) > 99.0:
                result_list.append(element)
                allele_type.append(col[1])
                allele_asigned.append(item)

            elif query_length >= (0.99 * allele_length) and float(col[2]) > 90.0 and 'ctxB' not in col[1]:
                result_list.append(element)
                allele_type.append(col[1])
                allele_asigned.append(col[1].split('_')[0])

    result_list.insert(0, '/'.join(allele_type))
    final_list = '\t'.join(result_list)
    return final_list

##writing out
def write_out(results_dict, args):
    print("Writing Out Results")

    if args.overwrite:
        # create output file and fill
        with open(args.output_folder + '/blast_results_.csv', 'w') as out:
            write_out_iterator(results_dict, out, args)

    else:
        #create new output file and fill
        time_stamp = '_'.join('_'.join(str(datetime.datetime.now()).split('.')[0].split(' ')).split(':'))
        with open(args.output_folder + '/blast_results_' + time_stamp + '.csv','w') as out:
            write_out_iterator(results_dict, out, args)
def write_out_iterator(results_dict, out, args):
    # column names
    headers_list = ["Accession", "allele(s)", "query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                    "qend", "sstart", "send", "evalue", "bitscore"]
    biotype_headers_list = ["Accession", 'biotype'] + headers_list[1:]

    for key in results_dict:
        out.write(key + '\n')
        #biotype has different headers
        if key == 'bio':
            for col_head in biotype_headers_list:
                out.write(col_head + ',')
            out.write('\n')

        else:
            for col_head in headers_list:
                out.write(col_head + ',')
            out.write('\n')

        for strain in results_dict[key]:
            out.write(strain + ',')
            col = results_dict[key][strain].split('\t')
            for cell in col:
                out.write(cell + ',')
            out.write('\n')

        out.write('\n')


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

#TODO must be able to type atypical strains and then within