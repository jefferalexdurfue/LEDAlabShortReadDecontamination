# LEDAlabShortReadDecontamination

# 01-Quality Control
  
requires: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/, https://multiqc.info/
   
    mkdir -p project-tutorial
    cd /project-tutorial
    mkdir -p 00-raw-fastq
    cd 00-raw-fastq

load paired short reads

    project-tutorial
    └──00-raw-fastq
        ├──raw-read_R1.fastq.gz
        └──raw-read_R2.fastq.gz

create a directory for the 01-Quality-control folder and enter
    
    mkdir -p ../01-quality-control
    cd /project-tutorial/01-quality-control

    #run the script
   
    RAW-READ1="/project-tutorial/00-raw-fastq/raw-read_R1.fastq.gz"
    RAW-READ2="/project-tutorial/00-raw-fastq/raw-read_R1.fastq.gz"
    ID=species
    
    seqkit stats ${RAW-READ*} 1> ${ID}_seqkit.txt 2>&1
    
    mkdir fastQC
    fastqc -t 48 ${RAW-READ*} -o fastQC/
    
    mkdir -p MultiQC
    multiqc fastQC/. -o MultiQC/

check the files: "species_seqkit.txt" and "multiqc_report.html"

    project-tutorial
    └──00-raw-fastq
        ├──raw-read_R1.fastq.gz
        ├──raw-read_R2.fastq.gz
        └──01-Quality_Control
            ├──species_seqkit.txt
            ├──fastQC
            |   ├──raw-read_R1.fastq.gz.html
            |   ├──raw-read_R1.fastq.gz.zip
            |   ├──raw-read_R2.fastq.gz.html
            |   └──raw-read_R2.fastq.gz.zip
            └──MultiQC
                ├──multiqc_report.html
                └──multiqc_data

----
# 02-TRIMMING
requires: https://github.com/OpenGene, http://www.usadellab.org/cms/?page=trimmomatic, https://github.com/dutilh/CAT, https://github.com/shenwei356/seqkit

    #create a directory for the 02-trimming folder and enter

    mkdir -p 02-trimming
    cd 02-trimming

    #run the script

    RAW-READ1="/project-tutorial/00-raw-fastq/raw-read_R1.fastq.gz"
    RAW-READ2="/project-tutorial/00-raw-fastq/raw-read_R1.fastq.gz"
    ID=species
    adapter="~/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
    unpaired1=${ID}_R1_SE_trimmed.fastq.gz
    unpaired2=${ID}_R2_SE_trimmed.fastq.gz
    cpus=48
    
    fastp -g -l 30 -f 3 -F 3 -t 5 -T 5 \
      --thread ${cpus}\
      --cut_tail \
      --cut_tail_window_size 5 \
      --cut_tail_mean_quality 15 \
      -i ${RAW-READ1} \
      -I ${RAW-READ2} \
      -o ${ID}_R1_PE_trimmed.fastq.gz \
      -O ${ID}_R2_PE_trimmed.fastq.gz \
      --unpaired1 ${ID}_R1_SE_trimmed.fastq.gz \
      --unpaired2 ${ID}_R2_SE_trimmed.fastq.gz \
      --adapter_fasta ${adapter} \
      --json fastp.json \
      --html fastp.html

    #to merge R1_SE+R2_SE

    cat ${unpaired1} ${unpaired2} > ${ID}_SE_trimmed.fastq.gz
    seqkit stats \
      ${ID}_R*_PE_trimmed.fastq.gz \
      ${ID}_SE_trimmed.fastq.gz 1> ${ID}_TRIM_seqkit.txt 2>&1

    #to run the quality control of the trimmed reads

    mkdir -p fastQC
    fastqc -t ${cpus} \
      ${ID}_R1_SE_trimmed.fastq.gz \
      ${ID}_R2_PE_trimmed.fastq.gz \
      ${ID}_SE_trimmed.fastq.gz \
      -o fastQC/

    mkdir -p MultiQC
    multiqc fastQC/. -o MultiQC/

to check the files: "fastp.html", "species_seqkit.txt" and "multiqc_report.html"

    project-tutorial
    ├──00-raw-fastq
    ├──01-Quality_Control
    └──02-trimming
        ├──fastp.html
        ├──fastp.json
        ├──species_R1_PE_trimmed.fastq.gz
        ├──species_R1_R2_SE_trimmed.fastq.gz
        ├──species_R1_SE_trimmed.fastq.gz
        ├──species_R2_PE_trimmed.fastq.gz
        ├──species_R2_SE_trimmed.fastq.gz
        ├──species_TRIM_seqkit.txt
        ├──fastQC
        |   ├──species_R1_PE_trimmed_fastqc.html
        |   ├──species_R1_R2_SE_trimmed_fastqc.zip
        |   ├──species_R1_PE_trimmed_fastqc.zip
        |   ├──species_R2_PE_trimmed_fastqc.html
        |   ├──species_R1_R2_SE_trimmed_fastqc.html
        |   └──species_R2_PE_trimmed_fastqc.zip
        └──MultiQC
            ├──multiqc_report.html
            └──multiqc_data

----
# 03-ECR [ERROR CORRECTION READS][ALLPATHS-LG]

Whole genome shotgun assembler that can generate high quality assemblies from short reads. 
requires: ftp://ftp.broadinstitute.org/pub/crd/ALLPATHS/Release-LG/latest_source_code/LATEST_VERSION.tar.gz

    cd project-tutorial
    mkdir -p 03-ECR
    cd 03-ECR

    #run the script

    READ1_PE_TRIMMED="/project-tutorial/02-trimming/species_R1_PE_trimmed.fastq.gz"
    READ2_PE_TRIMMED="/project-tutorial/02-trimming/species_R2_PE_trimmed.fastq.gz"
    READ_SE_TRIMMED="/project-tutorial/02-trimming/species_R1_R2_SE_trimmed.fastq.gz"
    OUTPUT_ECR="/project-tutorial/03-ECR"
    cpus=48
    ID=species
    
    perl ErrorCorrectReads.pl \
        PAIRED_READS_A_IN=${READ1_PE_TRIMMED} \
        PAIRED_READS_B_IN=${READ2_PE_TRIMMED} \
        UNPAIRED_READS_IN=${READ_SE_TRIMMED} \
        PAIRED_SEP=100 \
        THREADS=${cpus} \
        PHRED_ENCODING=33 \
        READS_OUT=${OUTPUT_ECR} 1> report_ECR.txt 2>&1
    
    seqkit \
      stats ${OUTPUT_ECR}/${ID}.paired.*.fastq ${OUTPUT_ECR}/${ID}.unpaired.fastq 
      1> ${OUTPUT_ECR}/report_ECR_seqkit.txt 2>&1
  
  to check the files: "report_ECR.txt" and "report_ECR_seqkit.txt"
  
    project-tutorial
    ├──00-raw-fastq
    ├──01-Quality_Control
    ├──02-trimming
    └──03-ECR
        ├──species.fastq
        ├──species.fastq.ids
        ├──species.paired.A.fastq
        ├──species.paired.B.fastq
        ├──species.unpaired.fastq
        ├──report_ECR.txt
        └──report_ECR_seqkit.txt
----
# 04-DNAmit
requares: http://github.com/Kinggerm/GetOrganelle2

    cd project-tutorial
    mkdir -p 04-DNA-mit

    #include a reference mitogenome or a gene (e.g., COI) from a closely related species
    #run the script

    READ1_PE_CORRECTED="/project-tutorial/03-ECR/species.paired.A.fastq"
    READ2_PE_CORRECTED="/project-tutorial/03-ECR/species.paired.B.fastq"
    SED=species_mitogenome/gene.fasta
    OUTPUT_DNA-mit="/project-tutorial/04-DNA-mit"
    cpus=48
    ID=species
    
    get_organelle_from_reads.py \
      -1 ${READ1_CORRECTED} \
      -2 ${READ2_CORRECTED} \
      -t ${cpus} -R 10 -k 21,45,65,85,105 \
      -F animal_mt \
      -s ${SED} \
      -o ${OUTPUT_DNA-mit}/${ID}_mtDNA

"contigs.fasta" or "scaffolds.fasta" file will be used for mitogenome annotation or removal in the reads (see below)

    project-tutorial
    ├──00-raw-fastq
    ├──01-Quality_Control
    ├──02-trimming
    ├──03-ECR
    └──04-DNA-mit
        ├──species_mitogenome/gene.fasta
        └──species_mtDNA
            ├──extended_1_paired.fq
            ├──extended_1_unpaired.fq
            ├──extended_2_paired.fq
            ├──extended_2_unpaired.fq
            ├──get_org.log.txt
            ├──seed
            └──extended_spades
                ├──assembly_graph.fastg
                ├──assembly_graph_after_simplification.gfa
                ├──assembly_graph_with_scaffolds.gfa
                ├──before_rr.fasta
                ├──contigs.fasta
                ├──contigs.paths
                ├──dataset.info
                ├──input_dataset.yaml
                ├──params.txt
                ├──run_spades.sh
                ├──run_spades.yaml
                ├──scaffolds.fasta
                ├──scaffolds.paths
                ├──spades.log
                ├──warnings.log
                ├──corrected
                ├──K105
                ├──K21
                ├──K45
                ├──K65
                ├──K85
                ├──misc
                ├──pipeline_state
                └──tmp

----
# 05-DECON
## 01-KRAKEN2_DB
requires: https://github.com/DerrickWood/kraken2, 
create the DB libraries
REMEMBER: (i) --use-ftp command would fail! (ii) you will need a x amount of gigabases to build your own DB to identify (ID) your original NGS reads:

In the creation of the DB to detect the exogenous sequences of our reads, we take as an example possibly from any taxon of Cnidaria.

    cd project-tutorial
    mkdir -p 05-DECON
    cd 05-DECON
    mkdir -p 01-KRAKEN2
    cd 01-KRAKEN2_DB
    mkdir -p contaminant_kraken2

    #step #1, taxonomy (download and setup): PATIENCE, THIS WOULD TAKE A FEW HOURS.
    REMEMBER: (i) --use-ftp command would fail! (ii) you will need a x amount of gigabases to build your own DB to identify (ID) your original NGS reads:
    
    DIR_kraken2="/project-tutorial/05-DECON/01-KRAKEN2_DB/contaminant_kraken2"
    cpus=48
    
    kraken2-build --download-taxonomy --no-masking --use-ftp  --threads ${cpus} --db ${DIR_kraken2}/
    
    #step #2, library
    
    kraken2-build --download-library human --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/
    kraken2-build --download-library bacteria --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/
    kraken2-build --download-library viral --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/
    kraken2-build --download-library UniVec --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/
    kraken2-build --download-library archaea --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/
    kraken2-build --download-library plasmid --no-masking --use-ftp --threads ${cpus} --db ${DIR_kraken2}/


    #step #3, food
    mkdir -p genome_DB
    cd genome_DB
Download the genome of possible species associated with the species under study (such as food and symbiotic animals) in fasta/fna format at the NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/)

    #e.g.
    #1-Artemia franciscana (GCA_032884065.1)
    #2-Breviolum minutum (CA_000507305.1)
    #3-Symbiodinium microadriaticum (GCA_001939145.1)
    #4-Symbiodinium sp. clade A Y106 (GCA_003297005.1)
    #5-Symbiodinium sp. clade C Y103 (GCA_003297045.1)
    #6-Symbiodinium kawagutii (GCA_009767595.1)
    #7-Effrenium voratum (GCA_963377175.1)
    #8-Symbiodinium natans (GCA_905221605.1)
    #9-Symbiodinium sp. (GCA_905221615.1)
    #10-Cladocopium goreaui (GCA_947184155.1)
    #11-Symbiodinium pilosum (GCA_905231905.1)
    #12-Symbiodinium necroappetens (GCA_905231915.1)
    #13-Symbiodinium sp (GCA_905221635.1)
    #14-Symbiodinium microadriaticum (GCA_018327485.1)
    #15-Effrenium voratum (GCA_963377275.1)
    #16-Effrenium voratum (GCA_963377065.1)

add the genome in the "contaminant_kraken2" librery

    DIR_kraken2="/project-tutorial/05-DECON/01-KRAKEN2_DB/contaminant_kraken2"
    GENOME_DB="/project-tutorial/05-DECON/01-KRAKEN2_DB/genome_DB"
    cpus=48
    
    kraken2-build --add-to-library ${GENOME_DB}/*.fna --no-masking --use-ftp --threads ${cpus} --db ${DIR_krakes2}

    # step #4, Build your DB
    #Build database from LIBRARIES. PATIENCE, THIS WOULD TAKE SEVERAL HOURS. REMEMBER: you will need an x amount of gigabases to build the DB. In my example, the DB size is ~85 Gb.
    # Kraken2 inspect
    kraken2-build --build  --threads ${cpus} --db ${DIR_krakes2}

to check that the DB library has been created
    
    project-tutorial
    ├──00-raw-fastq
    ├──01-Quality_Control
    ├──02-trimming
    ├──03-ECR
    ├──04-DNA-mit
    └──05-DECON
        └──01-KRAKEN2_DB
           ├──contaminant_krake2
           |  ├──archaea
           |  ├──bacteria
           |  ├──human
           |  ├──library
           |  ├──plasmid
           |  ├──taxonomy
           |  ├──UniVec
           |  └──viral
           └──genome_DB
    #nota= verify that all folders are ok????

If everything is ok, you can delete files no longer required (those downloaded from NCBI and personal libraries). As commented in kraken2 manual:
"After building a database, if you want to reduce the disk usage of the database, you can use the --clean option for kraken2-build to remove intermediate files from the database

    kraken2-build --clean --threads ${cpus} --db ${DIR_krakes2}

    #nota= verify that all folders are ok????

## 02-decon_paired_unpaired
How to visualize result in a much more complete way (using Krona): https://github.com/marbl/Krona
First, prepare the input using KrakenTools: https://github.com/jenniferlu717/KrakenTools
and you have to install biopython: https://github.com/biopython/biopython

    cd project-tutorial/05-DECON
    mkdir -p 02-DB_PAIRED_UNPAIRED
    cd 02-DB_PAIRED_UNPAIRED
    mkdir -p 01-READS_DECONTAMINATED
    
    READ1_PE_CORRECTED="/project-tutorial/03-ECR/species.paired.A.fastq"
    READ2_PE_CORRECTED="/project-tutorial/03-ECR/species.paired.B.fastq"
    READ_SE_CORRECTED="/project-tutorial/03-ECR/species.unpaired.fastq"
    DIR_kraken2="/project-tutorial/05-DECON/01-KRAKEN2_DB/contaminant_kraken2"
    OUTPUT_DB="/project-tutorial/05-DECON/02-DB_PAIRED_UNPAIRED"
    OUTPUT_READ_DECON="/project-tutorial/05-DECON/02-DB_PAIRED_UNPAIRED/01-READS_DECONTAMINATED"
    cpus=48
    ID=species
    
    #paired reads
    #identify the exogenous sequences in the paired reads
    kraken2 \
      --db ${DIR_kraken2} \
      --threads ${cpus} \
      --use-names \
      --paired \
      --minimum-base-quality 33 \
      --report ${OUTPUT_DB}/DB_report_DECON_PE.txt \
      ${READ1_PE_CORRECTED} ${READ2_PE_CORRECTED} \
      --output ${OUTPUT_DB}/DB_DECON_PE.out
    
    #check the proportion of exogenous sequences in the paired reads
    kreport2krona.py \
      -r ${OUTPUT_DB}/DB_report_DECON_PE.txt \
      -o ${OUTPUT_DB}/MYSAMPLE_PE.krona
    ktImportText \
      ${OUTPUT_DB}/MYSAMPLE_PE.krona \
      -o ${OUTPUT_DB}/MYSAMPLE_PE.krona.html
    
    ##decontaminate the paired reads
    extract_kraken_reads.py \
      -t ${cpus} \
      --include-children \
      -k ${OUTPUT_DB}/DB_1_decon_PE.out \
      --report ${OUTPUT_DB}/DB_1_report_decon_PE.txt 
      --fastq-output -s1 ${READ1_PE_CORRECTED} -s2 ${READ2_PE_CORRECTED} \
      -o ${OUTPUT_READ_DECON}/${ID}_R1_DECON_PE.fastq -o2 ${OUTPUT_READ_DECON}/${ID}_R2_DECON_PE.fastq \
       1> ${OUTPUT_READ_DECON}/${ID}_extract_DECON_PE.log 2>&1
    
    #unpaired reads
    #identify the exogenous sequences in the unpaired read
    kraken2 \
      --db ${DIR_kraken2} \
      --threads ${cpus} \
      --use-names \
      --minimum-base-quality 33 \
      --report ${OUTPUT_DB}/DB_report_DECON_SE.txt \
      ${READ_SE_CORRECTED} \
      --output ${OUTPUT_DB}/DB_DECON_SE.out
    
    #check the proportion of exogenous sequences in the unpaired read
    kreport2krona.py \
      -r ${OUTPUT_DB}/DB_report_DECON_SE.txt \
      -o ${OUTPUT_DB}/MYSAMPLE_SE.krona
    ktImportText \
      ${OUTPUT_DB}/MYSAMPLE_SE.krona \
      -o  ${OUTPUT_DB}/MYSAMPLE_SE.krona.html
    
    #decontaminate the unpaired read
    extract_kraken_reads.py \
      -t ${cpus} \
      --include-children \
      -k ${OUTPUT_DB}/DB_1_decon_SE.out \
      --report ${OUTPUT_DB}/DB_report_DECON_SE.txt \
      --fastq-output -s1 ${READ_SE_CORRECTED} \
      -o ${OUTPUT_READ_DECON}/${ID}_DECON_SE.fastq \
      > ${OUTPUT_READ_DECON}/${ID}_extract_DECON_SE.log 2>&1
    
    #check the stadistics of the decontaminated reads
    seqkit \
      stats ${OUTPUT_READ_DECON}/*.fastq 1> ${OUTPUT_READ_DECON}/report_DECON_seqkit.txt 2>&1

look "MYSAMPLE_SE.krona.html", "MYSAMPLE_PE.krona.html", and "report_DECON_seqkit.txt" files to check the stadistics

    project-tutorial
    ├──00-raw-fastq
    ├──01-Quality_Control
    ├──02-trimming
    ├──03-ECR
    ├──04-DNA-mit
    └──05-DECON
        ├──01-KRAKEN2_DB
        └──02-DB_PAIRED_UNPAIRED
            ├──DB_report_DECON_PE.txt
            ├──DB_report_DECON_SE.txt
            ├──DB_DECON_PE.out
            ├──DB_DECON_SE.out
            ├──MYSAMPLE_PE.krona
            ├──MYSAMPLE_SE.krona
            ├──MYSAMPLE_PE.krona.html
            ├──MYSAMPLE_SE.krona.html
            └──01-READS_DECONTAMINATED  
               ├──species_R1_DECON_PE.fastq
               ├──species_R2_DECON_PE.fastq
               ├──species_DECON_SE.fastq
               ├──species_extract_DECON_PE.log
               ├──species_extract_DECON_SE.log
               └──report_DECON_seqkit.txt

## 03-Remove mitogenome reads
requires: https://github.com/josephryan/FastqSifter

    cd project-tutorial/05-DECON
    mkdir -p 03-READS_CLEANED
    READ1_PE_CLEAN_I="project-tutorial/05-DECON/02-DB_PAIRED_UNPAIRED/01-READS_DECONTAMINATED/${ID}_R1_DECON_PE.fastq"
    READ2_PE_CLEAN_I="project-tutorial/05-DECON/02-DB_PAIRED_UNPAIRED/01-READS_DECONTAMINATED/${ID}_R2_DECON_PE.fastq"
    READ_SE_CLEAN_I="project-tutorial/05-DECON/02-DB_PAIRED_UNPAIRED/01-READS_DECONTAMINATED/${ID}_DECON_SE.fastq"
    OUTPUT_READS_CLEANED="project-tutorial/05-DECON/03-READS_CLEANED"
    MITOGENOME="project-tutorial/04-DNA-mit/species_mtDNA/extended_spades/scaffolds.fasta"
    cpus=48
    ID=species
    
    FastqSifter \
      --fasta=${MITOGENOME} \
      --left=${READ1_PE_CLEAN_I} \
      --right=${READ2_PE_CLEAN_I} \
      --unp=${READ_SE_CLEAN_I} \
      --threads ${cpus}
      --savereads \
      --out=${OUTPUT_READS_CLEANED}/${ID} 

the reads decontaminated: "species.filtered.A.fq", "species.filtered.B.fq", and "species.filtered.unp.fq". 

    project-tutorial
        ├──00-raw-fastq
        ├──01-Quality_Control
        ├──02-trimming
        ├──03-ECR
        ├──04-DNA-mit
        └──05-DECON
           ├──01-KRAKEN2_DB
           ├──02-DB_PAIRED_UNPAIRED
           └──03-READS_CLEANED
               ├──species.bwa.A.sai
               ├──species.bwa.B.sai
               ├──species.bwa.sampe.sam
               ├──species.bwa.samse.sam
               ├──species.bwa.unp.sai
               ├──species.filtered.A.fq
               ├──species.filtered.B.fq
               ├──species.filtered.unp.fq
               ├──species.matched.A.fq
               ├──species.matched.B.fq
               └──species.matched.unp.fq
