
# make a junction list for a thaliana tandems
# david's old tandems, not sure if they are appropriate... I may have to do this again.
tandems = read.table("/Users/grantdejong/Downloads/Arabidopsis_thaliana_gene.dis.cluster.txt", sep = "\t", header = F)

tandem_dubs = tandems[tandems$V2==2,] %>%
  mutate(V3=gsub("\\.[0-9]", "", tandem_dubs$V3)) %>%
  separate(V3, c("gene1", "gene2"), sep = " ")

write.table(tandem_dubs, "~/Desktop/scp/tandem_dubs.txt", quote = F, row.names = F)


tandem_exons = read.table("/Users/grantdejong/Desktop/scp/homologous_exons.formatted.txt", sep = "\t", header = F)

# create strandedness df = gene +/-
# incoporate old code - add gsub instead of substring