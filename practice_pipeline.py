import pysam
import os
import glob

foundation = os.getcwd()
base_path = os.getcwd()

base_path = base_path+'/sra_quality_pipeline'
os.makedirs(base_path)

file_folder = '/home/nataliestegman/COMP_483/group_project/test/bam_files_v2/1'

log_file = open(base_path+'/PipelineProject.log', 'w+')

#1: CONVERT TO PILEUP FILE
#this function will not only convert the bam files to pileup files, but also 
#it will exclude the pileup files with 0 bytes in it 

def bam_to_pileup(cwd):
    
    bam_files = []
    #pileup_files = []
    
    
    for filename in glob.iglob(file_folder+'/**', recursive=True):
        if os.path.isfile(filename): # filter dirs
            if filename.endswith('.bam'):
                #OR I could do IF .bam file is 24 bytes, next iteration and do not append to list:
                if os.stat(filename+'.bai').st_size > 24:
                    bam_files.append(filename)
                    #print(filename)
                    log_file.write(filename+'\n')
    log_file.write('\n')
    
    pileup_path = cwd+'/pileup_files'
    os.makedirs(pileup_path)       
    
    for b in bam_files:
        #convert to pileup file
        
        filename = b[64:-3]+'pileup'
        print(filename)
        command = 'samtools mpileup -o '+pileup_path+'/'+filename+' '+b
        print(command)
        #pileup_files.append(pileup_path+'/'+filename)
        os.system(command)
    #for p in pileup_files:
        #if it is empty file, then remove it from pileup_files list?
    
def parse_pileup_files(cwd):
    
    pileup_files = []
    
    pileup_path = cwd+'/pileup_files'
    
    for filename in glob.iglob(pileup_path+'/**', recursive=True):
            if os.path.isfile(filename): # filter dirs
                if filename.endswith('.pileup'):
                    pileup_files.append(filename)
                    #print(filename)
                    log_file.write(filename+'\n')
                    
    #return all the reads that are above 10.
    good_reads = []
    
    for f in pileup_files:
        for line in open(f):     # Variable 'i' named more intuitively
            split_line = line.split('\t') 
            read_number = split_line[3]
            #if read_number > 10:
                #something
            
        
    
bam_to_pileup(base_path)
parse_pileup_files(base_path)

log_file.close()

#2: FILTER TO DISCLUDE EMPTY PILEUP FILES

#3: COUNT READS (3RD OR 4TH COLUMN)

#4: WRITE OUT INTO CSV THE SRR AND THE AVG NUMBER OF READS?? IDK




'''
samfile = pysam.AlignmentFile('/home/nataliestegman/COMP_483/group_project/test/bam_files/1/SRR1051814.bam', "rb")


iter = samfile.fetch("seq1", 10, 20)
for x in iter:
    print(str(x))

#samtools view -c??
#finding the number of mapped reads
samtools flagstat file.sorted.bam
number of mapped locations:
    samtools view -F 0x04 -c file.sorted.bam

Example:
samtools sort -o '/home/nataliestegman/COMP_483/group_project/test/SRR1051814_sorted.bam' '/home/nataliestegman/COMP_483/group_project/test/bam_files/1/SRR1051814.bam' 
samtools view -F 0x04 -c '/home/nataliestegman/COMP_483/group_project/test/SRR1051814_sorted.bam'

'''

