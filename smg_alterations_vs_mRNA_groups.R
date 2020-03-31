library(magrittr);
library(ggplot2);
library(cowplot);

# clinical data (supplementary table 1)
clin <- read.delim('./data/dat_clin.tsv', stringsAsFactors = F);
rownames(clin) <- clin$MAF_sample_id;

clin <- clin[clin$study == 'TCGA' & !is.na(clin$rna_subtype), ]

# Binary matrix of missense and silent mutations
alterMat <- read.delim('./data/altMatBin_mut_cna_SMGs.tsv', stringsAsFactors = F, check.names = F);

# Gene annotation
anno <- read.delim('./data/ensemblAnnoBuild38Version86.tsv', stringsAsFactors = F);

# intersect clinical and mutation samples
samp.int <- Reduce('intersect', list(rownames(clin), colnames(alterMat)));

# Perform fisher's exact test on 2x3 table (alteration status x mRNA subgroup)
fisher.out <- apply(
    alterMat[ , samp.int], 1,
    function(x) {
        cnt <- table(factor(x, 0:1), clin[samp.int, ]$rna_subtype)
        fisher.out <- fisher.test(cnt)

        c(100 * cnt['1', ] / colSums(cnt), p = fisher.out$p.value)
    }
) %>% t %>% as.data.frame

fisher.out <- fisher.out[order(fisher.out$p), ]
fisher.out$fdr <- p.adjust(fisher.out$p, method = 'fdr')

fisher.out$gene <- anno[match(rownames(fisher.out), anno$EnsemblID), 'GeneSymbol']
fisher.out$gene <- factor(fisher.out$gene, levels = fisher.out$gene)

# Plot results
p1 <- ggplot(mapping = aes(x = gene, y = -log10(fdr)), data = fisher.out) + geom_bar(stat = 'identity', colour = 'black') +
    ylab('-Log10 FDR') + geom_hline(yintercept = -log10(0.2), color = 'red', linetype = 2) +
    scale_x_discrete(
        expand = c(0, 0)
    ) +
    theme_bw(base_size = 15) +
    theme(
        panel.border = element_rect(colour = 'black', fill = NA, size = 1), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.text = element_text(colour = 'black'), axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), 'cm')
    );

to.plot <- reshape2::melt(
    data = fisher.out,
    id.vars = c('gene', 'p', 'fdr')
);

to.plot$variable <- factor(
    to.plot$variable,
    c('OxPhos', 'Common', 'MITF-Low')
)

plot.rect <- data.frame(
    ymin = -Inf, ymax = Inf,
    xmin = seq(1, unique(to.plot$gene) %>% length) - 0.5,
    xmax = seq(1, unique(to.plot$gene) %>% length) + 0.5
)[seq(1, unique(to.plot$gene) %>% length) %% 2 == 1, ];

clust.col <- setNames(
    c('#74C476', 'lightskyblue', 'tan3'),
    c('OxPhos', 'Common', 'MITF-Low')
)

p2 <- ggplot(mapping = aes(x = as.numeric(gene), y = value, fill = variable), data = to.plot) +
    geom_rect(
        mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = plot.rect,
        fill = 'grey', alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_bar(stat = 'identity', colour = 'black', position = 'dodge', size = 0.3) +
    scale_x_continuous(
        breaks = 1:(unique(to.plot$gene) %>% length), labels = levels(to.plot$gene), expand = c(0, 0)
    ) +
    scale_fill_manual(name = NULL, values = clust.col) + ylab('Mutation frequency\nper group (%)') +
    theme_bw(base_size = 15) +
    theme(
        panel.border = element_rect(colour = 'black', fill = NA, size = 1), axis.title.x = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(colour = 'black', angle = 90, vjust = 0.5, hjust = 1), plot.margin = unit(c(1, 1, 1, 1), 'cm')
    );

the.legend <- get_legend(p2);

final.plot <- plot_grid(
    p1, p2 + theme(legend.position = 'none'),
    nrow = 2, ncol = 1, rel_widths = c(0.4, 0.6),
    align = 'hv'
);

final.plot <- plot_grid(
    final.plot, the.legend, NULL,
    nrow = 1, ncol = 3, rel_widths = c(0.9, 0.1, 0.01)
);

print(final.plot)