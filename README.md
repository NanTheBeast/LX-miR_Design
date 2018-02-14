# LX-miR Design

This repository contains all the code used to create the miRNA-targeting LX-miR CRISPR library. To keep the relationship between primary miRNAs, mature miRNAs, and sgRNAs clear, a relational (MySQL or MSSQL) database is created.

## Dependencies


* Python 2.7
* [Jupyter](http://jupyter.org) - notebooks for organizing code 
* [data_processing](https://github.com/jkurata/data_processing) - helper package with database connect and editing functions
* [MySQL.Connector](https://www.mysql.com/products/connector/) - connects to mySQL database
* [Pandas](http://pandas.pydata.org/) - dataframes 
* [Paramiko](http://www.paramiko.org/) - SSH connection
* [PyODBC](https://github.com/mkleehammer/pyodbc) - connects to MSSQL database
* [RNAfold](http://www.tbi.univie.ac.at/RNA/) - RNA secondary structure prediction
* [SciPy](https://www.scipy.org/) - scientific computing functions
* [SSHTunnel](https://pypi.python.org/pypi/sshtunnel) - SSH tunnel

## Database Summary

The code in this repository will create a database with 6 tables. Below is a summary of this database.

The <b>PrimaryMicroRNA</b> table contains information about all 1881 human primary microRNAs. The PriID is the accession number for the miRNA and is the <a href="http://www.tutorialspoint.com/sql/sql-primary-key.htm">primary key</a> for this table. The chromosome on which the primary miRNA is located is stored as 1-22 or X or Y. The StemloopSeq contains the DNA sequence of the stemloop as anotated in <a href="http://www.mirbase.org/">miRBase</a>. The LongSeq is the stemloop sequence plus an additional 20 bp on either side, which allows for all sgRNA with cutting sites in the stemloop to be identified. The RNAfold dot bracket structure of the stemloop sequence is also provided, which helps determine where in the stemloop the sgRNA will cleave. The miRNA family name is provided if the primary miRNA is known to be part of a miRNA family. The HighConfidence variable is a 'T' if the primary miRNA is considered high confidence (see <a href="http://www.mirbase.org/blog/2014/07/high-confidence-mirna-set-available-for-mirbase-21/">this blog post</a> and <a href="http://nar.oxfordjournals.org/content/42/D1/D68.full">this paper</a> for more about the high confidence anotation).

The <b>MatureMicroRNA</b> table contains information about the mature microRNAs. The MatID is the accestion number for the mature miRNA and is the <a href="http://www.tutorialspoint.com/sql/sql-primary-key.htm">primary key</a> for this table. The PriID corresponds to the PriID in the PrimaryMiRNA table and has a <a href="http://www.w3schools.com/sql/sql_foreignkey.asp">foreign key constraint</a>. 

The <b>SingleGuideRNA</b> table contains information about each unique sgRNA. This table does not contain any information about where the sgRNA targets due to the fact some sgRNAs target multiple miRNAs. The SgID is an automatically generated number will no meaning. The exclude field contains a single character describing why the sgRNA should be excluded from the pool:
* 'R': Greater than 10 targets in the hg19 genome with no mismatches
* 'T': 4 or greater T's in a row, which can lead to <a href="http://www.sciencedirect.com/science/article/pii/0092867481905225">RNA Pol III termination</a>
* 'A': All targets of the sgRNA have an NAG PAM, which is <a href="http://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html">much less effective than an NGG PAM</a>
* 'C': The cleavage site is outside of the stemloop
* 'D': Less than 10, but at least one off-target site which exactly matches the sgRNA
* 'L': Low max Zhang score (<0.2)

If multiple possible exclussion reasons apply, the reason most likely to affect the sgRNA read out was used ('R'>'T'>'A'>'C'>'D'>'L'). The presence of the exact sgRNA sequence in the <a href="http://genome-engineering.org/gecko/?page_id=15">Zhang lentiCRISPRv2 pool</a> was recorded as a 'Y' in the ZhangLibrary column. The scores for each sgRNA using various scoring methods are also included in this table. The <a href="http://crispr.mit.edu/">Zhang lab sgRNA scoring website</a> was scripted against (see <a href='#Zhang_Scoring'>Zhang Scoring</a>) and the scores were recorded in the ZhangScore column. The <a href="https://www.bioconductor.org/packages/3.3/bioc/html/CRISPRseek.html">CRISPRseek</a> scores for the first 1,000 sgRNAs are included, though running this R package was too resource intensive to perform for all sgRNAs. On-target scores for each sgRNA were provided by sgRNAScorer (see <a href="http://www.nature.com/nmeth/journal/v12/n9/full/nmeth.3473.html">this paper</a>), the orgional Doench scoring algorthm(see <a href="http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html">this paper</a>) and an updated version of Doech, <a href="http://research.microsoft.com/en-us/projects/azimuth/">Azimuth</a> (described in <a href="http://biorxiv.org/content/early/2015/06/26/021568">this bioRXiv paper</a> and <a href="http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3437.html">this Nature Biotechnology paper</a>). The number of exact matches to the sgRNA sequence in the hg19 genome, as found by alignment using <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">Bowtie2</a>, are also recorded. The sequence of the oligo to be ordered (85 nts long), which includes primer binding sites, is included in the table. The final element in the table is a boolean value stating if the give sgRNA is or is not included in the ordered miR-ome oligos. 

The <b>SgRNATargetInformation</b> table contains information about each miRNA site targeted by the sgRNAs. This table contains any data specific to a single site, not to the sgRNA which can target multiple sgRNAs. This includes the LongSg sequence: NNNN[sgRNA target 20]NGGNNN, which is used for Doench/Azimuth scoring (see above). The PriID of the targeted miRNA corresponds to the PriID in the PrimaryMiRNA table and has a <a href="http://www.w3schools.com/sql/sql_foreignkey.asp">foreign key constraint</a>. If the cleavage site falls within a mature miRNA site, the MatID also corresponds the the MatID in the MatureMiRNA table and has a foreign key constraint. The SgID also has a foreign key constraint with the SgIDs in the SingleGuideRNA table. 

The <b>OverlappingSgRNAs</b> table contains information about which sgRNAs have identical seeds. These sgRNAs which target the same miRNA are likely to have similar off-target sites and may create false positives. Therefore, only one sgRNA with any given seed is used when there are enough sgRNA to allow for the top four to be choosen. The SelectedSgID has a higher ZhangScore than the OverlappingSgID, which is discarded. 

The <b>ControlSgRNA</b> table contains the name and the sequence of the control sgRNAs used along with the ID originally used in the paper from whtich the control sgRNA was taken.

The <b>InPool</b> table contains information the sgRNAs which are included in the miRNA-ome pool. Includes the id of the sgRNA or the name of the sgRNA (for controls) and the sequence of the ssDNA oligo submitted for plate synthesis. 

### PrimaryMicroRNA

| Variable | Type |
|-|-|
| PriID | varchar(50) |
| PriMiRName | varchar(50) |
| Chr | varchar(2) |
| ChrStrand | char(1) |
| GenomeStart | int |
| GenomeEnd | int |
| StemLoopSeq | varchar(200) |
| LongSeq | varchar(250) |
| RNAfold | varchar(200) |
| MiRFamily | varchar(50) |
| HighConfidence | nchar(1) |


### MatureMicroRNA

| Variable | Type |
|-|-|
| MatID | varchar(50) |
| MatMiRName | varchar(50) |
| PriID | varchar(50) |
| MatStart | int |
| MatEnd | int |
| MatSeq | varchar(50) |

### SingleGuideRNA

| Variable | Type |
|-|-|
| SgID | int |
| SgRNA | char(20) |
| Exclude | char(1) |
| ZhangLibrary | varchar(1) |
| ZhangScore | float |
| CRISPRseek | float |
| SgRNAScorer | float |
| MaxDoenchScore | float |
| MaxAzimuthScore | float |
| NumExactMatch | int |
| hg38NumExactMatch | int |

### SgRNATargetInformation

| Variable | Type |
|-|-|
| SgID | int |
| LongSg | char(30) |
| SgStart | int |
| SgEnd | int |
| SgChr | varchar(2) |
| SgStrand | nchar(1) |
| CleavageSite | varchar(50) |
| CleaveStart | int |
| CleaveEnd | int |
| PriID | varchar(50) |
| MatID | varchar(50) |
| PAM | char(3) |
| DoenchScore | float |
| AzimuthScore | float |
| AzimuthScorev2 | float |

### OverlappingSgRNAs

| Variable | Type |
|-|-|
| SelectedSgID | int |
| OverlappingSgID | int|

### ControlSgRNA


| Variable | Type |
| -------- | ---- |
| SgRNAName | varchar(50) |
| SgRNA    | char(20) |
| LiteratureSgRNA | varchar(200) |

### InPool

| Variable | Type |
|-|-|
| SgID | int |
| SgName | varchar(50) |
| OligoSeq | char(85) |