```
Usage: pb-mod2bed.py [OPTIONS] BAM

  This script converts pacbio modified bam files to expanded bed files in the
  following format: 
  
  read_name start end methylation_status pass_tag

Options:
  -n, --read_names PATH  filter analysis based on file containing reamv d names
                         per line
  -u, --can_prob FLOAT   probability threshold for canonical bases
  -m, --mod_prob FLOAT   probability threshold for modified bases
  -o, --out FILENAME     output path
  --help                 Show this message and exit.
```