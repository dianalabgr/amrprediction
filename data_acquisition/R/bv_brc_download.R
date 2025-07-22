options(warn = -1)


args <- commandArgs(TRUE)
id <- args[1]
output <- args[2]
log <- args[3]
threads <- args[4]


# If the ID starts with "GC", then it is an assembly and we need to download it through the NCBI API (Entrez Direct System)
# If the ID does not start with "GC", then it is a BioSample and we need to download the raw DNA-seq reads, 
# assemble them with SPAdes, clean up the intermediate files and rename the final assembly file
if (substring(id, 1, 2) == "GC") {
  system(paste("efetch -db assembly -id", id, "-format docsum | \
                xtract -pattern DocumentSummary -element FtpPath_GenBank | \
                awk 'END{{ split($1, a, /\\//); print $1\"/\"a[10]\"_genomic.fna.gz\"}}' | \
                xargs wget -nc -c -O", paste0(output, ".gz"), "2>", log, "&& \
                gunzip -c", paste0(output, ".gz"), ">", output, "&& rm", paste0(output, ".gz")))
} else if (substring(id, 1, 2) != "GC") {
  system(paste("fastq-dump --split-files --gzip -O $(dirname", paste0(output, ")/", id, "/"), id, ">", log, "2>>", log, 
               " && spades.py --only-assembler -1", paste0("$(dirname ", output, ")/", id, "/", id, "_1.fastq.gz"),
               "-2", paste0("$(dirname ", output, ")/", id, "/", id, "_2.fastq.gz"),
               "-m 50 -t ", threads, " -o", paste0("$(dirname ", output, ")/", id, "/"), ">>", log, 
               "&& mv", paste0("$(dirname ", output, ")/", id, "/scaffolds.fasta"), output, 
               "&& rm -r", paste0("$(dirname ", output, ")/", id, "/")
         ))
}
