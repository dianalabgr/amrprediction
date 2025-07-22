options(warn = -1)


args <- commandArgs(TRUE)
id <- args[1]
output <- args[2]
log <- args[3]
threads <- args[4]


# If the ID starts with "PNUS", then it is an assembly and we need to download it through the NCBI API (Entrez Direct System)
# If the ID start with "SAM", then it is a BioSample and we need to download the raw DNA-seq reads,
# assemble them with SPAdes, clean up the intermediate files and rename the final assembly file
if (substring(id, 1, 4) == "PNUS") {
    # system(paste0(
    #     "url=$(esearch -db assembly -query ", id, " | efetch -format docsum | ",
    #     "xtract -pattern DocumentSummary -element FtpPath_GenBank | ",
    #     "awk 'END{ split($1, a, \"/\"); print $1\"/\"a[10]\"_genomic.fna.gz\" }') && ",
    #     "curl --output /dev/null --head --silent --fail \"$url\" && ",
    #     "wget -nc -c \"$url\" -O ", output, ".gz 2>> ", log, " && ",
    #     "gunzip -c ", output, ".gz > ", output, " && rm ", output, ".gz || exit 0"
    # ))

    system(paste0(
        "url=$(esearch -db assembly -query ", id, " | efetch -format docsum | ",
        "xtract -pattern DocumentSummary -element FtpPath_GenBank | ",
        "awk 'END{ split($1, a, \"/\"); print $1\"/\"a[10]\"_genomic.fna.gz\" }'); ",
        "if [ $(echo $url | grep -c '_genomic.fna.gz') -eq 0 ]; then ",
        "curl --output /dev/null --head --silent --fail \"$url\" && ",
        "wget -nc -c \"$url\" -O ", output, ".gz 2>> ", log, " && ",
        "gunzip -c ", output, ".gz > ", output, " && rm ", output, ".gz; ",
        "else ",
        "ACC=$(esearch -db sra -query ", id, " | efetch -format docsum | ",
        "xtract -pattern Runs -element Run@acc) && ",
        "fastq-dump --split-files --gzip -O $(dirname ", output, ")/", id, "/ $ACC", " > ", log, " 2>> ", log, " && ",
        "spades.py --only-assembler -1 $(dirname ", output, ")/", id, "/${ACC}_1.fastq.gz ",
        "-2 $(dirname ", output, ")/", id, "/${ACC}_2.fastq.gz -m 50 -t ", threads, " -o $(dirname ", output, ")/", id, "/ >> ", log, " && ",
        "mv $(dirname ", output, ")/", id, "/contigs.fasta ", output, " && ",
        "rm -r $(dirname ", output, ")/", id, "/; ",
        "fi"
    ))
                
} else if (substring(id, 1, 3) == "SAM") {
    system(paste0(
        "ACC=$(efetch -db sra -id $(esearch -db biosample -query ", id, " | efetch -format docsum | ",
        "xtract -pattern Ids -element Id | awk '{ for (i=1; i<=NF; i++) if($i ~ /^SR/) print $i }') -format docsum | xtract -pattern Runs -element Run@acc) && ",
        "fastq-dump --split-files --gzip -O $(dirname ", output, ")/", id, "/ $ACC", " > ", log, " 2>> ", log, " && ",
        "spades.py --only-assembler -1 $(dirname ", output, ")/", id, "/${ACC}_1.fastq.gz ",
        "-2 $(dirname ", output, ")/", id, "/${ACC}_2.fastq.gz -m 50 -t ", threads, " -o $(dirname ", output, ")/", id, "/ >> ", log, " && ",
        "mv $(dirname ", output, ")/", id, "/scaffolds.fasta ", output, " && ",
        "rm -r $(dirname ", output, ")/", id, "/"
    ))
}
