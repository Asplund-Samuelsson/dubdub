![alt text](dubdub.png "DubDub")

# DubDub: Automated DubSeq data analysis

Automated [DubSeq](https://github.com/psnovichkov/DubSeq) data analysis.


## About

`./dubdub.sh -d <results directory> -f <fastq.gz directory> -l <barseq layout table> -L <library> [-U <pre-sequence> -P <pre-sequence position> -D <post-sequence>]` to run pipeline in the _results directory_ using the specified _fastq.gz directory_ barcode sequencing files and _barseq layout table_ on the specified DubSeq _library_. Parameters in brackets are optional and are by default obtained from `config.sh`.

Analysis based on [Mutalik _et al._ (2019), DOI 10.1038/s41467-018-08177-8](https://www.nature.com/articles/s41467-018-08177-8) and [DubSeq version f452baa](https://github.com/psnovichkov/DubSeq/commit/f452baab7d9c9e56150803962dc755a5a39b78fd).


## System requirements

Linux operating system (Tested on Ubuntu 20.04 LTS)

Python 3.7.6 (Tested with 3.7.6)

R ≥ 3.6.3 (Tested with 3.6.3)

GNU parallel 20161222 (Tested with 20161222)

Python libraries: ...

R libraries: tidyverse


## Installation guide

```
git clone https://github.com/Asplund-Samuelsson/dubdub
```

Installation of the pipeline takes approximately X seconds.


## Demonstration

### Commands

`command`

The demonstration example takes X to run on a Y.

### Expected output

...

### Further instructions

...

## Authors

Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)
