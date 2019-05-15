import argparse
import glob
import sys
from os import path
from Bio.Blast.Applications import NcbiblastnCommandline
import multiprocessing
import datetime
import progressbar
from operator import itemgetter

#smaller functions
def check_databases_input(args):
    database_okay = True
    accepted_databases = ['sero','ctxB','tcpA','rstR','bio']
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
    database_dict = {}
    for el in input_database_list:
        if el == 'bio':
            for i in ['ctxB', 'tcpA', 'rstR']:
                database_dict[i] = {}
            return database_dict

        else:
            database_dict[el] = {}

    return database_dict
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
                        help="Databases to process. Options ctxB, tcpA, rstR, sero, bio"
                             "For combinations use comma seperated. Eg. ctxB,tcpA")
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
    print("Processing genomes.")
    bar = progess_use(args)
    bar.start()
    genome_count = 0
    results_dict = {}
    for strain in glob.iglob(args.strains_directory + '/*.f*'):
        bar.update(genome_count)
        genome_count = genome_count + 1
        accession = strain.split('/')[-1].split('.')[0]
        results_dict[accession] = blast_input_against_db(strain, current_blast_db_path, args)
    bar.finish()

    if args.databases == 'bio':
        biotype_results_dict = biotype_filter(results_dict)
        results_dict = {**results_dict, **biotype_results_dict}

    return results_dict
def blast_input_against_db(strain, current_blast_db_path, args):

    cpus = multiprocessing.cpu_count()

    #blast query database against genome
    blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=current_blast_db_path, outfmt=6, perc_identity=90, num_threads=cpus)
    out, err = blast_string()

    #if no blast alignment found
    if out == '':
        blast_result_dict = {1:'No alignment found. Modify "-r 2" to see low quality hits.'}
        return blast_result_dict

    # format blast result to have dict with all blast hit
    else:
        blast_result_dict = format_blast_output(out, args)
        return blast_result_dict
def format_blast_output(out, args):

    # create dict of databases to fetch results from
    databases_dict = databases_list(args)

    #split the blast result
    blast_hits = out.split('\n')[:-1] #remove last new line char

    #organise blast hits by precen of id
    blast_hits = sorted(blast_hits, key=itemgetter(2,11))

    databases_dict = blast_filter(databases_dict, blast_hits, args)

    return databases_dict
def blast_filter(databases_dict, blast_hits, args):

    set_wd = path.dirname(path.abspath(__file__))

    #read in file of allele lengths
    gene_size_dict = input_allele_lengths(set_wd, args)

    for key in databases_dict:
        result_list = []
        allele_type = []

        for element in blast_hits:
            col = element.split('\t')
            allele_length = gene_size_dict[col[1]]
            query_length = int(col[3]) - int(col[4])
            if key in col[1]:
                # print(element)

                ##check if hit is exact match
                if query_length == allele_length and float(col[2]) == 100.0:
                    result_list.append(element)
                    allele_type.append(col[1])

                elif key == 'sero' and query_length > (0.95 * allele_length) and float(col[2]) > 95.0:
                    result_list.append(element)
                    allele_type.append(col[1])

        #TODO add when have multiple 100% fragments


        result_list.insert(0, '/'.join(allele_type))
        final_list = '\t'.join(result_list)
        databases_dict[key] = final_list

    return databases_dict
def biotype_filter(results_dict):

    #create function to check what biotype the strain is
    biotype_results = {}

    ##possible combiniations of alleles for each biotype
    bio_classical = ['ctxB1', 'rstR_cla', 'tcpA_cla_WT']
    bio_eltor = ['ctxB3', 'rstR_el', 'tcpA_el_WT','tcpA_el_A226']
    bio_mozam = ['ctxB1', 'rstR_cla', 'tcpA_el_A226', 'tcpA_el_WT']
    bio_atypical = ['ctxB1', 'rstR_cla', 'rstR_el', 'tcpA_el_A226', 'tcpA_el_WT']

    #for each strain with blast results
    for key in results_dict.keys():
        biotype_results[key] = {}
        #get the alleles for all blast results
        alleles_list = []
        for value in results_dict[key]:
            col = results_dict[key][value].split('\t')
            #if only a single allele present, add to list
            if len(col[0].split('/')) == 1:
                alleles_list.append(col[0])

            #if more than one allele present
            #add the first allele, check if next are different
            if len(col[0].split('/')) > 1:
                alleles_list.append(col[0].split('/')[0])
                for frag in col[0].split('/'):
                    if frag not in alleles_list:
                        alleles_list.append(frag)

        #check strain biotype by comparing alleles against bio_lists
        if set(alleles_list).issubset(bio_classical):
            biotype_results[key]['classical'] = alleles_list

        if set(alleles_list).issubset(bio_eltor):
            biotype_results[key]['eltor'] = alleles_list

        if set(alleles_list).issubset(bio_atypical):
            biotype_results[key]['atypical'] = alleles_list

    return biotype_results


#writing out
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
            write_out_iterator(results_dict, out)
def write_out_iterator(results_dict, out, args):
    # column names
    headers_list = ["Accession", "allele", "query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                    "qend", "sstart", "send", "evalue", "bitscore"]
    biotype_headers_list = ["Accession", 'biotype', 'ctxB', 'tcpA', 'rstR']

    # # write out column headers
    if args.databases == 'bio':
        for col_head in biotype_headers_list:
            out.write(col_head + ',')
        out.write('\n')
    else:
        for col_head in headers_list:
            out.write(col_head + ',')
        out.write('\n')

    for key in results_dict:

        # # write out blast hits for each genome
        for item in results_dict[key]:
            if bool(results_dict[key][item]):
                if args.databases == 'bio':
                    out.write(key + ',' + item + ',')
                    for i in results_dict[key][item]:
                        out.write(i + ',')
                    out.write('\n')

                else:
                    out.write(key + ',')
                    col = results_dict[key][item].split('\t')
                    for i in col:
                        out.write(i + ',')
                    out.write('\n')
def main():

    set_wd = path.dirname(path.abspath(__file__))
    args = parseargs(set_wd)
    database_input = check_databases_input(args)
    if database_input:
        results_dict = blast_search(set_wd, args)
        write_out(results_dict, args)
    else:
        print('One of following database(s) is incorrect :' + str(args.databases))
        sys.exit()

if __name__ == '__main__':
    main()

#TODO must be able to type atypical strains and then within