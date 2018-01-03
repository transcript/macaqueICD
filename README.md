# macaqueICD
Python, R, and bash scripts used for analysis of rhesus macaque gut microbiomes for ICD.  Except when described otherwise, tools used are all components of the SAMSA2 pipeline [https://github.com/transcript/samsa2](https://github.com/transcript/samsa2) .

###Base pipeline and analysis

*Script used: bash\_scripts/macaque\_master\_script.bash*

**Explanation:** This is the main pipeline script of SAMSA2, set up to run on a cluster for the macaque metatranscriptome files.  The reference databases used are NCBI's RefSeq Bacterial non-redundant database (version 102, released 22 December 2015), and TheSEED Subsystems database (retrieved 18 January 2017).

###Host read analysis

*Script used: bash\_scripts/macaque\_host\_pipeline\_script.bash*

**Explanation:** Macaque host reads were screened out of metatranscriptome later than analysis of bacterial sequences, and thus a simplified pipeline script was created for analysis.  The reference database used is NCBI's RefSeq macaque proteins database (retrieved 4 October 2017).

###Figure generation R scripts

