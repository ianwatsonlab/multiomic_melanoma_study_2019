setwd('/mnt/d/Documents/mac/Documents/c_mounts/g7_wl/COMMON/analysis/tcga_skcm/');

library(ggplot2);

# clinical data (supplementary table 1)
clin <- read.delim('./supp_tables/dat_clin.tsv', stringsAsFactors = F);
rownames(clin) <- clin$MAF_sample_id;

# Binary matrix of missense and silent mutations
syn.mis <- read.delim('./data/mutMatBin_syn_mis_allGenes.tsv', stringsAsFactors = F, check.names = F);

# Mutation consequences for SMGs
cons <- read.delim('./data/conseqMat_SMGs.tsv', check.names = F);
cons <- ifelse(cons != '0', 1, 0); # Binarize missense, LoF, and in-frame indel strings

# Gene annotation
anno <- read.delim('./data/ensemblAnnoBuild38Version86.tsv', stringsAsFactors = F);

anno$type <- NA;
anno[anno$chr %in% 1:22, 'type'] <- 'autosomal';
anno[anno$chr == 'X', 'type'] <- 'X';

# intersect clinical and mutation samples
samp.int <- Reduce('intersect', list(rownames(clin), colnames(syn.mis), colnames(cons)));

clin <- clin[samp.int, ];
syn.mis <- syn.mis[ , samp.int];
cons <- cons[ , samp.int];

# Compute expected ORs for autosomal and X-linked genes
gene.or <- apply(
    syn.mis, 1,
    function(x) {
        ct <- table(factor(x, levels = 0:1), clin$sex);
        out <- (ct['1', 'Female'] / ct['0', 'Female']) / (ct['1', 'Male'] / ct['0', 'Male']);
        return(out);
    }
);

gene.or <- gene.or[rowSums(syn.mis) >= 10]; # use gene with 10 or more mutations for reliable OR estimates
auto.x.median.or <- sapply(split(gene.or, anno[match(names(gene.or), anno$EnsemblID), 'type']), median); # NAs are lost due to split

to.plot <- lapply(
    rownames(cons),
    function(i) {

        ct <- table(
            factor(cons[i , ], 0:1),
            factor(
                clin$sex,
                levels = c('Male', 'Female')
            )
        );

        out <- data.frame(
            gene = i,
            or.manual  = (ct['1', 'Female'] / ct['0', 'Female']) / (ct['1', 'Male'] / ct['0', 'Male']),
            or = fisher.test(ct, alternative = 'two.sided', or = 1)$estimate,
            p = fisher.test(ct, alternative = 'two.sided', or = auto.x.median.or[anno[anno$EnsemblID == i, 'type']])$p.value,
            row.names = NULL
        );

        return(out);
    }
);

to.plot <- do.call(rbind, to.plot);
to.plot <- to.plot[order(to.plot$p), ];

to.plot$FDR <- p.adjust(to.plot$p, method = 'fdr');
to.plot$gene <- anno[match(to.plot$gene, anno$EnsemblID), 'GeneSymbol'];

print(to.plot);

p <- ggplot(mapping = aes(x = log2(or), y = -log10(FDR), label = gene), data = to.plot) + geom_text() +
    theme_bw(base_size = 15);

print(p);
