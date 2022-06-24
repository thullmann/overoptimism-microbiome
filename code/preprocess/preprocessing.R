# Preprocessing the OTU table #######################################################################

library(phyloseq)

# The datasets were downloaded from ftp://ftp.microbio.me/AmericanGut/ag-2017-12-04/
# from the folder 11-packaged\split\notrim\raw
# As we couldn't open the file otu_table__BODY_HABITAT_UBERON_feces__.biom with the phyloseq::import_biom command, 
# we first converted the file to txt and then to json as follows:

# # Load the reticulate package and install biom-format 
# library(reticulate)
# system("echo $PATH") # path might be without python interpreter
# reticulate::py_run_file("hello.py") # throws error but necessary to add interpreter path to PATH variable
# system("echo $PATH") # should now have additional python path
# py_install("biom-format", pip = TRUE)
# # Now convert the biom file to txt
# system("biom convert -i otu_table__BODY_HABITAT_UBERON_feces__.biom -o otu_table.txt --to-tsv --header-key taxonomy")
# # Open the txt file (e.g. with Notepad++), delete the header "Constructed from biom file" and save the txt as otu_table2.txt
# # Now convert the txt file to json
# system("biom convert -i otu_table2.txt -o otu_table__BODY_HABITAT_UBERON_feces_json.biom --to-json --table-type='OTU table' --process-obs-metadata taxonomy")
# This corresponds to the biom file that is available at Zenodo.  

ag = import_biom("data/otu_table__BODY_HABITAT_UBERON_feces_json.biom")
 
sampledata = read.delim("data/metadata__BODY_HABITAT_UBERON_feces__.txt", sep = "\t", header = TRUE, row.names = 1)

sample_data(ag) = sampledata 

# We used the same three filter steps as in Badri et al. (2020). Our code for filtering is thus similar to Badri et al.'s code,
# which can be accessed at https://github.com/MichelleBadri/NormCorr-manuscript

# remove samples with a sequencing depth of less than 10000 counts

counts = colSums(ag@otu_table@.Data)
ag.filt = prune_samples(counts >= 10000, ag)

# remove OTUs which are present in less than 30% of the remaining samples

prev = rowSums(sign(ag.filt@otu_table@.Data))
ag.filt = prune_taxa(prev > 0.3*nsamples(ag.filt), ag.filt)

counts = colSums(ag.filt@otu_table@.Data)

# remove samples with a sequencing depth under the 10\%-percentile

ag.filt = prune_samples(counts > quantile(counts, probs=seq(0, 1, .1))[2], ag.filt)


# agglomerate to genus level

# assign numbers to missing genera
taxtab = ag.filt@tax_table@.Data
missing_genus = which(taxtab[, "Rank6"] == "g__")
taxtab[missing_genus, "Rank6"] = paste0("g__", 1:length(missing_genus))
unique(taxtab[taxtab[,"Rank5"] == "f__","Rank4"]) # all species with unknown family have order Clostridiales
taxtab[taxtab[,"Rank5"] == "f__", "Rank5"] = "f__(o__Clostridiales)"
ag.filt@tax_table@.Data = taxtab
ag.genus = tax_glom(ag.filt, "Rank6")
# keep in mind that the taxa in the OTU table are still named by the OTU IDs
# attributes(ag.genus@otu_table)$taxa_are_rows

saveRDS(ag.genus, "data/ag.genus.rds")
