# iMapSplice
__iMapSplice__ is an individualized RNA-seq read alignment software. It makes use of personal genomic information and performs an unbiased alignment towards a genome index carrying both the reference and any alternative bases. It is light-weight and does not require an index to be rebuilt for each individual. Importantly, it breaks the computational emphasis or dependency on the position of canonical splice site dinucleotide motifs in the reference genome and enables iMapSplice to discover personal splice junctions created through splice site mutations.

## Installation

### System Requirements
iMapSplice has been tested on Linux platforms with the following system settings.
  * Red Hat 4.4.5-6
  * gcc version 4.9.1

### Installation
1. Clone the repository:

    ```
    git clone https://github.com/xa6xa6/iMapSplice.git
    ```

2. Build and install:

    ```
    cd code
    make all
    ```

    After successfully build the code, you can find __iMapSplice__ toolchain at ``bin/``.

    ```
    ls bin/
    ```
    
### Manual
1. Building reference genome index (MapSplice index, only need to run once, shared by all the individuals)

    Note: before building index, you need to put all the sequence files of reference genome into a directory (like /PATH/hg19/). And all the   sequence files are required to be in the following format:

    (1) In "FASTA" format, with ".fa" extension.
    (2) One chromosome per sequence file.
    (3) Chromosome name in the header line (">" not included) is the same as the sequence file base name, and does not contain any blank space. E.g. If the header line is ">chr1", then the sequence file name should be "chr1.fa".
    (5) No other files in the same folder.

    
    ```
    ./buildWholeGenome <input_chromosomes_folder_path> <output_globalIndex_folder_path>
    ./build2ndLevelIndex <input_globalIndex_folder_path> <output_localIndex_folder_path>
    ```
    
2. Generating SNP-mers
   
    Generating SNP-mers with phased SNPs
    
    ```
    ./getSNPmer-phased <input_globalIndex_folder_path> <input_GAF_path> <input_SNPlist_path(phased)> <output_SNPmer_folder_path> <SNPmer_length>
    ```    
    
    __or__
    
    Generating SNP-mers with unphased SNPs
    
    ```
    ./getSNPmer-unphased <input_globalIndex_folder_path> <input_GAF_path> <input_SNPlist_path(unphased)> <output_SNPmer_folder_path> <SNPmer_length>
    ```    
       
    Note: SNP list (tab delimited, four columns):
    
    ```
    ...
    chr1 1000 A T
    chr2 2000 C G
    chrX 3000 G T 
    ...
    ```
    
    Each row records a SNP:
    Column 1: chromosome name
    Column 2: chromosome pos
    Column 3: Haplotype 1 base (phased) or reference base (unphased)
    Column 4: Haplotype 2 base (phased) or alternate base (unphased)
    
 3. Building SNP-mer index
    
    ```
    ./buildSNPmerIndex <input_SNPmer_folder_path/SNPinAnn.fa> <output_SNPmerIndex_folder_path> 
    ```
  
 4. Mapping
   
    Mapping with phased SNPs
    
    ```
    ./iMapSplice-phased -P <input_SNPlist_path(phased)> -Q <input_SNPmerIndex_folder_path> -G <input_globalIndex_folder_path> -L <input_localIndex_folder_path> -1 <read_end1> -2 <read_end2> -T <threads_num> -O <output_folder>
    ```
   
    __or__
   
    Mapping with unphased SNPs
   
    ```
    ./iMapSplice-unphased -P <input_SNPlist_path(unphased)> -Q <input_SNPmerIndex_folder_path> -G <input_globalIndex_folder_path> -L <input_localIndex_folder_path> -1 <read_end1> -2 <read_end2> -T <threads_num> -O <output_folder>
    ```
    
### License
Please refer to LICENSE.txt
