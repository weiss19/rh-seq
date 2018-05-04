args<-commandArgs(TRUE)
myfile_filtered_gene_inserts=args[1]
myfile_mwu_test_results=args[2]
insert_data <- read.table(myfile_filtered_gene_inserts, sep="\t", head=T, as.is=T)

# make a separate plot for each gene; should rewrite to stick them all on
# one plot
Wilcoxon_data <- read.table(myfile_mwu_test_results, sep="\t", head=T, as.is=T)
cutoff <- 0.01
which_signif <- which(as.numeric(Wilcoxon_data[,"pval"]) < cutoff)
for(i in 1:length(which_signif)) {
  my_gene <- as.character(Wilcoxon_data[i,"gene"])

  # get thermotolerance data for the two alleles into separate vectors
  w <- which(insert_data[,"gene"] == my_gene)
  insert_data_my_gene <- insert_data[w,]
  w <- which(insert_data_my_gene[,"allele"] == "sp")
  sp_thermotolerances <- insert_data_my_gene[w,"X39_28_log2"]
  w <- which(insert_data_my_gene[,"allele"] == "sc")
  sc_thermotolerances <- insert_data_my_gene[w,"X39_28_log2"]

  # calculate density (smoothed histogram) but do not plot
  sp_density <- density(sp_thermotolerances)
  sc_density <- density(sc_thermotolerances)

  # when you plot two traces on the same scale, you have to set the x and y limits
  # to the same for each
  xma <- max(c(sp_density$x, sc_density$x))
  xmi <- min(c(sp_density$x, sc_density$x))
  yma <- max(c(sp_density$y, sc_density$y))
  ymi <- min(c(sp_density$y, sc_density$y))

  testxmin <- 0
  testxmax <- -8.0
  testymin <- 0
  testymax <- 0.6

  # make the plot
  filename <- paste(my_gene, ".pdf", sep="")
  pdf(filename)
  par(lwd=5)
  plot(sp_density, col="#00A5CD", xlim=c(xmi, xma), ylim=c(ymi, yma), xlab="log2(39/28) per insert", ylab="density of inserts", main="")
  par(new=T)
  plot(sc_density, col="#FAB900", xlim=c(xmi, xma), ylim=c(ymi, yma), xlab="", ylab="", main=my_gene)
 legend(x=0.3*xma, y=0.9*yma, legend=c("insert in S. par","insert in S. cer"), fill=c("orange","blue"))

  dev.off()
}

