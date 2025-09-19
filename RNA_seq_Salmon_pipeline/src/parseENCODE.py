#!/usr/bin/env python

import requests, csv, argparse

def filter_bam_files(bamfiles, bamfiles_filtered, linkfiles2experiments):
    
 
    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}
    list_of_files_cleared = []

    # filter filenames of experiments
    file = open(bamfiles, "r")
    new_file = open(bamfiles_filtered, "w")

    lines = file.readlines()
    for number, line in enumerate(lines):
        file_name = line.split('/')[4]

        url = 'https://www.encodeproject.org/search/?searchTerm='+str(file_name)+'&frame=object&format=json'
    
        # GET the search result
        try:
            response = requests.get(url, headers=headers)
        except TimeoutError:
            pass
    
        # Extract the JSON response as a python dictionary
        search_results = response.json()
    
        # filtering cond //assembly, reference, file format
        file_type = search_results['@graph'][0]['file_type'] #bam
        if (file_type == 'bam'):
            output_type = search_results['@graph'][0]['output_type'] #transcriptome alignments
            assembly = search_results['@graph'][0]['assembly'] #GRCh38
	    
            try:
                annotation = search_results['@graph'][0]['genome_annotation']
            except KeyError:
                pass

            if (output_type == 'transcriptome alignments') & (assembly == 'GRCh38') & (annotation == 'V29'):
        
                replicates_num = search_results['@graph'][0]['technical_replicates'][0]
                proj_num = search_results['@graph'][0]['dataset'].split('/')[2]
            
                #save meta info: original line number, file name, accession number, and replicate number
                list_of_files_cleared.append([number, file_name, proj_num, replicates_num])
            
                #save lines that met conditions
                new_file.write(line)
            
            total = len(lines)
            bar_update = "Processed {} lines out of {}".format(number, total)
            print(bar_update)

    new_file.close()

    #save files-to-experiments table
    
    with open(linkfiles2experiments, "w") as f:
        writer = csv.writer(f)
        writer.writerows(list_of_files_cleared)   
        print('DONE')
        return



def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", "--bfiles_ori", required=True, help="original list with bam files")
    parser.add_argument("-f2", "--bfiles_new", required=True, help="new name for the filtered list with bam files")
    parser.add_argument("-f3", "--files_exp", required=True, help="name for files_ids2experiment_ids table")
    args = parser.parse_args()
    
    filter_bam_files(args.bfiles_ori, args.bfiles_new, args.files_exp)

if __name__ == "__main__":
	main()
