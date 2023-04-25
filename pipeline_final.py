import os
import glob
from Bio import SeqIO
import argparse
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET

'''
ADD ARGPARSE:
--help: gives help explaining each flag
--test: runs with test data provided in github
-i: path to SearchSRA folders
-o: path to where you want output
-c: coverage percentage minimum (defaults to 0.1)
-q: path to query fasta sequence that you initally put into SearchSRA
-e: NCBI email used to retrieve SRA files
#this cannot do multifasta file!!! You can only do single fasta files OR you can automatically concatenate the contigs
of the fna file and put NNNNNNN between the contigs
'''

def msg(name=None):
    return '''python3 pipeline_final.py -f phage_reference_file -i input_path -o output_path -e ncbi_email -c [optional] coverage_threshold'''

#This adds the argparse input for the phage reference file, searchsra output, entrz email, coverage threshold, and path to where the output will go. 
parser=argparse.ArgumentParser(usage=msg())
parser.add_argument('-o', '--output_path', action="store", metavar='<directory>', help='Directory to store resulting files (required)')
parser.add_argument('-i', '--input_path', action="store", metavar='<directory>', help='Directory to SearchSRA folders of results (required)')
parser.add_argument('-t', '--num_threads', action="store", metavar='<int>', type=int, default = 4, help='Number of processors to use (default=4)')
parser.add_argument('-c', '--coverage_threshold', action="store", metavar='<int>', default= 0.5, help='Threshold for phage coverages (default=0.5)')
parser.add_argument('-f', '--phage_reference_file', action="store", metavar='<filename>', help='Reference file of phage sequences (required)')
parser.add_argument('-e', '--entrez_email', action="store", metavar='email', help='Email to run entrez efetch on (required)')
group=parser.add_mutually_exclusive_group()
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args=parser.parse_args()

#If a certain argument is provided, provide an error statement.
if args.phage_reference_file is None:
    parser.error('Phage reference file must be provided for analysis.')
if args.output_path is None:
    parser.error('Output path must be provided.')
if args.input_path is None:
    parser.error('SearchSRA results name must be provided.')
if args.entrez_email is None:
    parser.error('Email for entrez esearch must be provided.')

#This creates the directory where the output of the pipeline will go. 
base_path = args.output_path+'/sra_quality_pipeline'
os.makedirs(base_path)
os.chdir(base_path)

#This is the path to the search sra output.
file_folder = args.input_path 

#This is the path to the phage query used in the searchsra.
query_path = args.phage_reference_file

#This creates the log file.
log_file = open(base_path+'/PipelineProject.log', 'w+')

#1: CONVERT TO PILEUP FILE
#This function will not only convert the bam files to pileup files, but also 
#it will exclude the pileup files with 0 bytes in it.
def bam_to_pileup(cwd):
    
    #This list holds the list of total bam files from the searchsra output.
    bam_files = []    
    for filename in glob.iglob(file_folder+'/**', recursive=True):
        if os.path.isfile(filename): # filter dirs
            if filename.endswith('.bam'):
                #Since the .bai files that are 24 bytes are empty, this filters them out 
                #so now only the bam files with actual information are converted into pileup files.
                if os.stat(filename+'.bai').st_size > 24:
                    bam_files.append(filename)
    
    #This folder will contain the pileup files.              
    pileup_path = cwd+'/pileup_files'
    os.makedirs(pileup_path)   
    
    #For each bam file that had actual info in it,
    for b in bam_files:  
        
        #Convert the file into a pileup file.
        filename = b[-14:-3]+'pileup'
        command = 'samtools mpileup -o '+pileup_path+'/'+filename+' '+b
        os.system(command)

#2: FILTER PILEUP FILES
def parse_pileup_files(cwd):
    
    #This holds the converted pileup files.
    pileup_files = []
    pileup_path = cwd+'/pileup_files'

    #This gets the paths to all the pileup files and puts them in a list.
    for filename in glob.iglob(pileup_path+'/**', recursive=True):
            if os.path.isfile(filename): # filter dirs
                if filename.endswith('.pileup'):
                    #If the pileup file is more than 0 bytes, append it to the list.
                    if os.stat(filename).st_size > 0:
                        pileup_files.append(filename)
    
    #This counts the number of bases in the input fasta file.         
    for record in SeqIO.parse(query_path, 'fasta'):
        query_length = len(record.seq)
        
    #This list will hold the coverage of each pileup file to the query.
    coverage = []
    
    #For each pileup file,
    for f in pileup_files:
        
        #Read in the file and count the number of rows.
        file = open(f, 'r')
        beep = file.readlines()
        num_rows = len(beep)
        
        #Calculate the query coverage by dividing #rows/#bases in query
        #and make a list of lists of the query coverage, filename.
        coverage.append([num_rows/query_length, f])
    
    #This filters the query coverage by the provided or default coverage threshold.
    filtered_coverage = [] 
    for i, j in coverage:
        if float(i) >= float(args.coverage_threshold):
            filtered_coverage.append([i,j])
    
    #This sorts the filtered query coverage in ascending order.
    filtered_coverage = sorted(filtered_coverage, key=lambda x: x[0])
    
    #This writes out the headers of the log file.
    log_file.write('SRA Number'+'\t'+'Query Coverage %'+'\t'+'Average Read Length'+'\t'+'Study Title'+'\n')
    
    #In the reverse order of the files above the provided threshold.
    for i,j in reversed(filtered_coverage):
        sra_number = j[-17:-7]
        result = sra_number.startswith('/')
        print(result)
        
        if result == True:
            sra_number = j[-16:-7]
            print(j[-16:-7])
        
        print(sra_number)
        #Convert and round the decimal to a percentage.
        percentage = round(i*100,2)
        
        #Read in the file name path as a dataframe, 
        #and calculate the average of the 4th column (read count) and write out in log file.
        pileup_df = pd.read_csv(j,header=0,delimiter="\t",names = ['no1','no2', 'no3','yes','no4', 'no5'])
        avg = pileup_df['yes'].mean()
        avg = round(avg,2)
        
        #this fetches the sra db
        Entrez.email = args.entrez_email
        
        try:
            handle = Entrez.efetch(db="sra", id=sra_number)

            #Parse XML file with ElementTree.
            tree = ET.parse(handle)
            root = tree.getroot()

            title_list = []
            
            # Gets all the titles from the sra study
            for title in root.iter('TITLE'):
                title_list.append(title.text)

            metadata = ';'.join(title_list)
            
        #Exception for if efetch doens't work or if it can't find the metadata.
        except:
            metadata = ('There is no metadata available')
        
        #This writes out the query coverage %, the SRA record, average read count in the pileup file, and the metadata.
        log_file.write(sra_number+'\t'+str(percentage)+'\t'+str(avg)+'\t'+metadata+'\n')
         

bam_to_pileup(base_path)
parse_pileup_files(base_path)
log_file.close()
