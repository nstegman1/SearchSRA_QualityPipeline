<h1 align="center">SearchSRA and Quality Pipeline</h1>

<p align="center">This is a project for COMP 383/483. This filters through the data from Search SRA and returns the SRA number that have multiple hits to the query fasta file. Please read the Github Wiki page for more in-depth information.</p>

### Before you run this, make sure

* you have ran SearchSRA with your desired fasta file.
* you have provided the fasta file that is one continuous sequence instead of a file with multiple contigs that you used to run SearchSRA with. 

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

* -f: [REQUIRED] Phage reference file that was the query file in the Search SRA run
* -i: [REQUIRED] The output data from SearchSRA, which should contain a folder of 25 folders of .bam and .bam.bai files.
* -o: [REQUIRED] This is where you specify where you want the output of the pipeline to be.
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
python3 pipeline_v2.py -f phage_reference_file -i input_path -o output_path -c [optional] coverage_threshold
```

### The Files Produced
* A folder containing the filtered pileup files of the SearchSRA output
* A .csv file containing the coverage, SRA record, and average read length











