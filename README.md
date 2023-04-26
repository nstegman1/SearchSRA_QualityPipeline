<h1 align="center">SearchSRA and Quality Pipeline</h1>

<p align="center">This is a project for COMP 383/483. This filters through the data from Search SRA and returns the SRA number that have multiple hits to the query fasta file. Please read the Github Wiki page for more in-depth information.</p>

### Before you run this, make sure

* You have ran SearchSRA with your desired fasta file (https://www.searchsra.org/).
* You have provided the fasta file that is one continuous sequence instead of a file with multiple contigs that you used to run SearchSRA with.
* You have an NCBI account.

### For Testing Purposes for Comp Bio class:

* There is sample SearchSRA output for the T3 bacteriophage under the compbio.cs.luc.edu server
* path to SearchSRA output:
```
/home/nstegman/search_sra_data/bam_files_v2
```
* This SearchSRA output is intended to be used with the T3 fasta file in this repo

### Dependencies

* os
* python
* biopython (SeqIO)
* glob
* pandas
* argparse
* searchSRA
* samtools
* ElementTree (xml.etree.ElementTree)

### Flags Description

* -f: [REQUIRED] Phage reference file that was the query file in the Search SRA run.
* -i: [REQUIRED] The output data from SearchSRA, which should contain a folder of 25 folders of .bam and .bam.bai files.
* -o: [REQUIRED] This is where you specify where you want the output of the pipeline to be.
* -e: [REQUIRED] NCBI email used to retrieve SRA files.
* -c: This is the coverage threshold, which is how well you want the SRA files to map to the query file. This can be between 0 and 1 but defaults to 0.5.

### Executing Program

* This will download the pipeline script, and sample phage genome fasta file:
```
git clone https://github.com/nstegman1/SearchSRA_QualityPipeline.git
```

* Then, you must be in the same directory where the .py pipeline file is:
```
cd insert_base_path_here/sra_quality_pipeline
```

* Now, once you know the path to your SearchSRA ouput and query fasta file, then you are ready to execute the pipeline.

```
python3 pipeline_final.py -f phage_reference_file -i input_path -o output_path -e ncbi_email -c [optional] coverage_threshold
```

### The Files Produced
* A folder containing the filtered pileup files of the SearchSRA output
* A .csv file containing the SRA record, query coverage, average read count, and study names











