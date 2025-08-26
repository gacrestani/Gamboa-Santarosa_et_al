renv::status()

renv::install("bioc::biomaRt")
renv::install("bioc::clusterProfiler")
renv::install("bioc::org.Dm.eg.db")
renv::install("bioc::TxDb.Dmelanogaster.UCSC.dm6.ensGene")
renv::install("MartaPelizzola/ACER")
renv::install("ThomasTaus/poolSeq")


# Also, run on debian:
sudo apt install libcurl4-openssl-dev \
libxml2-dev \
libfontconfig1-dev \
libglpk-dev 

