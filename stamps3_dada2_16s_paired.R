######################################################################
#####  GOAL: Choose and validate appropriate dada2 parameters    #####
#####         for processing  paired-end Iluumina 16S data       #####
######################################################################
#####   DATA is taken from the following paper                   #####
#
# Pyrethroid exposure alters internal and cuticle surface bacterial communities in Anopheles albimanus, ISMEJ, 2019.
# https://doi.org/10.1038/s41396-019-0445-5
# Sequencing data: First 18 samples from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA512122
#
######################################################################

rm(list = ls())
#installing dada2

#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.14")


library(dada2); packageVersion("dada2") # 1.16 or later
library(ShortRead)
library(Biostrings)
getwd()
# First download the data being used for this lab
#download.file("https://figshare.com/ndownloader/files/36394545", "data/data_IlluminaPE.zip")
#unzip("data/data_IlluminaPE.zip")

# Where is the freshly downloaded data?
        list.files()
        path <- "data_IluminaPE" # REPLACE PATH with the path to the unzipped directory containing the gzipped fastqs
# Challenge: Define the path both as a relative path and as an absolute path.
        abolute_path <- paste(getwd(), "data_IluminaPE", sep = "/")

# Read in the forward and reverse fastq filenames
        fnFs <- list.files(path, pattern="_1.fastq.gz", full.names=TRUE)
        fnRs <- list.files(path, pattern="_2.fastq.gz", full.names=TRUE)
        head(fnFs)
# Do those filenames look like what we expect?

        plotQualityProfile(fnFs[1:2]) #seqeucnign set of forward is similar
        plotQualityProfile(fnRs[1:2]) #seqeucnign set of forward is similar
# Good quality/bad quality? What are the read lengths? 
        #Quality is good. Length is longer than we need
# How might these quality profiles inform our choice of truncation lengths?
        #We need 300 length reads. 
        
# Define the paths to filtered files we are going to create
        filtFs <- file.path(path, "filtered", basename(fnFs)) 
        filtRs <- file.path(path, "filtered", basename(fnRs))
# The filtered files will be in the `filtered/` subdirectory within `path`
        
###################################################################
######  Are primers on these reads that need to be removed?  ######
######  How long is the sequenced amplicon? SEE PAPER        ######
###################################################################
        #Primer : 
        #FWD <- "CCTACGGGNGGCWGCAG"  
        #FWD <- "CCTACGGGNGGCWGCAG"  
        #REV <- "GACTACHVGGGTATCTAATCC"
        
        #allOrients <- function(primer) {
                # Create all orientations of the input sequence
        #        require(Biostrings)
        #        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        #        orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        #                     RevComp = reverseComplement(dna))
        #        return(sapply(orients, toString))  # Convert back to character vector
        #}
        #FWD.orients <- allOrients(FWD)
        #REV.orients <- allOrients(REV)
        #FWD.orients
        
        #fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
        #fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
        #filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
        
        #primerHits <- function(primer, fn) {
                # Counts number of reads in which the primer is found
        #        nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        #        return(sum(nhits > 0))
        #}
        #rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
        #      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
        #      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
        #      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
        
        #The primer is sequenced
                #Howver, the pirmer target sequence is the part of the actual seqeuence.
                #Mihai Pop - we don't need to remove the region.

        #Sequence length
                #16s v3-v4 regions were amplified - 443bp.
        
# Perform filtering and trimming
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=2, 
                             trimLeft=c(0, 0), # REPLACE XXX/YYY with proper parameter choices
                             truncLen=c(249, 249)) # REPLACE XXX/YYY with proper parameter choices
        head(out)
# Were most reads retained during filtering? If not, would any parameter tweaks help?
# How might the depth of sequencing in this data affect the questions that can be addressed?

# Learn the error model from the filtered data.
        errF <- learnErrors(filtFs, multi=TRUE) # `multi=TRUE` activated multithreading to reduce computation time
        errR <- learnErrors(filtRs, multi=TRUE)

# Visualize the error model. Points are observations, black line is fitted error model.`
        plotErrors(errF)
        plotErrors(errR)
# Do the fitted error models look reasonable?
        #Fitting looks reasonable
# Run the DADA2 method using the fitted error model.
        
        #ddF <- dada(filtFs, errF, pool=F, multi=TRUE)
        ddF <- dada(filtFs, errF, pool="pseudo", multi=TRUE)
        #ddR <- dada(filtRs, errF, pool=F, multi=TRUE)
        ddR <- dada(filtRs, errR, pool="pseudo", multi=TRUE)
        
# What pooling option makes sense? FALSE (default), "pseudo", or TRUE? See ?dada for more info
        #Prooling resuted in slightly more unique sequences relative to the no pooling.
        
# For more pseudo-pooling detail: https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling
# Challenge: Try different pooling options and compare them.
        ddF[[1]]
        ddR[[2]]
        
# Merge the denoised forward and reverse reads together.
        mm <- mergePairs(ddF, filtFs, ddR, filtRs, verbose=TRUE)
        head(mm[[1]])
        
# Were most reads retained during merging. If not, why not?

# Construct a sequence table: rows are samples, columns are ASVs, values are abundances.
        sta <- makeSequenceTable(mm)
        dim(sta)
# How many samples and ASVs are in this table?

# Remove chimeric ASVs and construct a new chimera-free sequence table.
        st <- removeBimeraDenovo(sta, multi=TRUE, verbose=TRUE)
        sum(st)/sum(sta)
# Were most reads retained during chimera removal? How about ASVs? 
# How does the fraction of chimeric ASVs compare to the fraction of chimeric reads?
# Why is that?

####################################################################
######  Inspect the number of reads passing through each step ######
######  THIS IS THE SINGLE MOST IMPORTANT SANITY CHECK!!      ######
####################################################################
# Code derived from the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
        getN <- function(x) sum(getUniques(x))
        track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mm, getN), rowSums(st))
        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track) <- basename(fnFs)
        head(track)
# In this case, most reads should make it through the entire pipeline!
# Most importantly, a large majority (>80% of reads) should merge successfully,
#  and almost all (>95%)  reads should pass chimera filtering.
# IF THAT ISN'T THE CASE, you have a problem, and need to revisit your truncation lengths
#   (merging problem) or primer removal (trimLeft, chimera problem).

# The dada2 taxonomic reference page https://benjjneb.github.io/dada2/training.html has links to a 
#   number of reference databases formatted to work with `assignTaxonomy`.I recommend Silva for 16S.
#   However, here we will use RDP in this lab for file size reasons.
        download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz?download=1",
                      "rdp_train_set_18.fa.gz")
# Assign taxonomy down to the genus level to these 16S ASVs.
        tax <- assignTaxonomy(st, "rdp_train_set_18.fa.gz", multi=TRUE)
        unname(head(tax))
        str(tax)
        #692 taxa
# Are reasonable taxonomic assignments being made to the abundant taxa?
# Challenge: Do this again using the Silva database, and compare results.
        tax2 <- assignTaxonomy(st, "/Users/MK1059/Dropbox (Personal)/Macbook Pro 16 2021/Downloads/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
        unname(head(tax2))
        str(tax2)
        
# Challenge 2: Use `addSpecies` or `assignSpecies` to do species-level assignemtn, where appropriate.
        tax2 <- addSpecies(tax2, "/Users/MK1059/Dropbox (Personal)/Macbook Pro 16 2021/Downloads/silva_species_assignment_v138.1.fa")
        #species matchign - exact match
        ##sometimes 99% is used in publications
        view(tax2)
        #No match for exact sequences...no species assigned
        # should change the variables
        
# Phyloseq is a package for the manipulation and analysis of microbiome data.
# Here we use it briefly to produce an ordination of our sequenced communities.
        library(phyloseq); library(ggplot2)

# We define a very simple data.frame that records the 3 experimental groups these samples 
#   came from (see paper for more info)
        samdf <- data.frame(row.names=rownames(st), 
                            sampleID=rownames(st), 
                            treatment=rep(c("Unexposed", "Susceptible", "Resistant"), each=6))

# Use phyloseq to plot a bray-curtis NMDS odrination of our samples, colored by treatment.
        ps <- phyloseq(sample_data(samdf), otu_table(st, taxa_are_rows=FALSE), tax_table(tax))
        plot_ordination(ps, ordinate(ps, method="NMDS", distance="bray"), color="treatment") + 
                aes(size=4) + theme_bw()

# What do you see?
# Note that this is a small subset of the total data in this study...