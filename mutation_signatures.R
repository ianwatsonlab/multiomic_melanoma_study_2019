set.seed(777);

lib <- lapply(c('magrittr', 'caret', 'NNLM', 'pheatmap', 'ggplot2', 'reshape2', 'extrafont', 'cowplot'), library, character.only = T);

perform.cv <- FALSE;
if (!perform.cv) k <- 3;

n.threads <- 1L;

to.nnmf <- read.delim('./data/trinucleotideMutCounts.tsv', stringsAsFactors = F, check.names = F);
to.nnmf <- as.matrix(to.nnmf);

if (perform.cv) {

    # For cross validation, we will perform 5 repeats of 10-fold CV. Here we are getting the indices of observations to mask per fold/repeat
    # NMF will be used to impute the values of masked observations, and the difference between the imputed and
    # true observations will be used to evaluate the NMF decomposition
    n.cells <- length(to.nnmf);
    cv.folds <- lapply(
        createMultiFolds(y = 1:n.cells, k = 10, times = 5),
        function(x) (1:n.cells)[!(1:n.cells) %in% x]
    );

    ### For each k
    k.res <- lapply(
        setNames(nm = 10:1),
        function(k) {

            ### For each fold
            fold.res <- lapply(
                cv.folds,
                function(i) {
                    to.train <- to.nnmf;
                    to.train[i] <- NA;

                    dcomp <- nnmf(
                        A = to.train, k = k, alpha = rep(0, 3), beta = rep(0, 3), method = c('scd', 'lee')[1],
                        loss = c('mse', 'mkl')[2], init = NULL, mask = NULL, W.norm = -1L, check.k = TRUE,
                        max.iter = 5000L, rel.tol = 1e-04, n.threads = n.threads
                    );

                    train.fit <- (dcomp$W %*% dcomp$H);
                    loss <- mse.mkl(obs = to.nnmf[i], pred = train.fit[i], na.rm = TRUE, show.warning = TRUE);

                    out <- data.frame(
                        # training data
                        spearman.training = cor(as.numeric(to.nnmf), as.numeric(train.fit), use = 'pairwise.complete.obs', method = 'spearman'),
                        MKL.training = min(dcomp$mkl), MSE.training = min(dcomp$mse),
                        # test data
                        spearman.test = cor(to.nnmf[i], train.fit[i], method = 'spearman'),
                        MKL.test = loss['MKL'], MSE.test = loss['MSE']
                    );

                    return(out);
                }
            ) %>% do.call(rbind, .);

            return(fold.res);
        }
    );

    mean.mse <- sapply(
        k.res,
        function(x) mean(x$MSE.test)
    );

}

### Final model
dcomp <- nnmf(
    A = to.nnmf, k = if (!perform.cv) k else which.min(mean.mse), alpha = rep(0, 3), beta = rep(0, 3), method = c('scd', 'lee')[1],
    loss = c('mse', 'mkl')[2], init = NULL, mask = NULL, W.norm = -1L, check.k = TRUE,
    max.iter = 5000L, rel.tol = 1e-04, n.threads = n.threads
);

H <- dcomp$H;
W <- dcomp$W;

colnames(W) <- rownames(H) <- paste0('NNLMsig.', 1:ncol(W));

### Plot signatures
to.plot <- t(t(W) / colSums(W));

tn.order <- rownames(to.plot)[
    do.call(
        rbind,
        strsplit(rownames(to.plot), '\\(|\\)')
    )[ , c(2, 1, 3) ] %>% as.data.frame %>% order(.)
];

to.plot <- reshape2::melt(to.plot);
to.plot$Var1 <- factor(to.plot$Var1, tn.order);

p <- ggplot(mapping = aes(x = Var1, y = to.plot$value), data = to.plot) + geom_bar(stat = 'identity') +
    theme_bw(base_size = 20) + facet_wrap( ~ Var2, nrow = 3, scales = 'free_y') + ylab('percent contribution') +
    theme(
        axis.text.x = element_text(angle = 90, family = 'mono', vjust = 0.5),
        axis.title.x = element_blank()
        );

### plot percent contribution of each tn mutation to signature
print(p);

### Compare with cosmic signatures
load('./data/signatures.exome.cosmic.v3.may2019.rda');

known.sig <- signatures.exome.cosmic.v3.may2019;

colnames(known.sig) <- gsub('\\[', '(', colnames(known.sig)) %>% gsub('\\]', ')', .);
known.sig <- t(known.sig);

tn.int <- intersect(
    rownames(known.sig),
    rownames(W)
);

to.plot <- cor(known.sig[tn.int, ], W[tn.int, ], method = 'pearson');
to.plot <- reshape2::melt(to.plot);

plt.list <- lapply(
    split(to.plot, to.plot$Var2),
    function(x) {
        x <- x[order(-x$value), ];
        x$Var1 <- factor(x$Var1, as.character(x$Var1) %>% unique);
        out <- ggplot(mapping = aes(x = Var1, y = value, fill = value), data = x) + geom_bar(stat = 'identity') +
            theme_bw(base_size = 20) + facet_wrap( ~ Var2, nrow = 3, scales = 'fixed') + ylab('Spearman cor.') +
            scale_fill_gradient2(
                low = 'blue', mid = 'white', high = 'red',
                name = NULL, guide = F
            ) +
            theme(
                axis.text.x = element_text(angle = 90, family = 'mono', vjust = 0.5, hjust = 1),
                axis.title.x = element_blank()
            );

        return(out);
    }
);

### plot correlation with cosmic signatures
p <- plot_grid(plotlist = plt.list, nrow = 3);
print(p);
