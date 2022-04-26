abundance_d <- read.table(file = 'abundance.tsv', sep = '\t', header = TRUE)
myvars <- c("target_id","est_counts")
abundance_f <- abundance_d[myvars]

gene_table <- read.table(file = 'mart_export.txt', sep = '\t', header = TRUE)

my_result <- merge(abundance_f,
                   gene_table,
                   by.x = "target_id",
                   by.y = "Transcript.stable.ID.version",
                   all = FALSE);
dim(my_result)

my_vars  <- c("Gene.name","est_counts")
table_f <- my_result[my_vars,]
a <- length(unique(table_f$Gene.name))

library(dplyr)

new_table <- table_f %>%
group_by('Gene.name') %>%
summarize(est_counts = sum(est_counts))

counts_read <- read.table(file = 'TCGA_BRCA_Gene_ReadCounts.txt', sep = '\t', header = TRUE)
merging_counts <- merge(counts_read,
                        new_table,
                        by.x = "Gene",
                        by.y = "Gene.name",
                        all = FALSE);

column_data_rd <- c("TCGA.C8.A138.01","est_counts")
plot_scatter <- merging_counts[column_data_rd]

row_sub = apply(plot_scatter, 1, function(row) all(row !=0 ))

plot_data <-plot_scatter[row_sub,]

library(ggplot2)

dev.new() 

plotI <- ggplot(plot_data, aes(x = TCGA.C8.A138.01, y = est_counts)) + 
  
  geom_point(alpha = 0.1) + 
ggtitle('Comparison of Counting Methods')

print(plotI)

dev.new()

testing_data = apply(plot_data[, 1:2] < 10^5, 1, all)
testing_data2 <-plot_data[testing_data,]
testing_data = apply(testing_data2[, 1:2] > 10^3, 1, all)
testing_data2 <-testing_data2[testing_data,]

plotII <- ggplot(testing_data2, aes(x = TCGA.C8.A138.01, y = est_counts,color=est_counts)) + 
    coord_trans(x = "log2", y = "log2") +
    geom_point(alpha = 0.1)+ 
    geom_smooth(method=lm) +
ggtitle('Comparison of Counting Methods (detail)')

print(plotII)

 
