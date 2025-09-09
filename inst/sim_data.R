### simulate example data
set.seed(123);
n <- 50;
p <- 100;

data.a <- sapply(1:p, rnorm, n = n);
data.b <- sapply(1:p, rnorm, n = n);
# data.c is categorical, taking on values -1, 0, 1 (CNA loss, neutral, gain)
# note all categorical features still need to be coded numerically
data.c <- sapply(
    X = 1:p,
    FUN = function(x) {
        cna <- round(runif(n = n, min = -1, max = 1));
        return(cna);
        }
    );

rownames(data.a) <- paste0('patient.', 1:n);
colnames(data.a) <- paste0('gene.', 1:p);
rownames(data.b) <- paste0('patient.', 1:n);
colnames(data.b) <- paste0('gene.', 1:p);
rownames(data.c) <- paste0('patient.', 1:n);
colnames(data.c) <- paste0('gene.', 1:p);

# features in rows, samples in cols
data.a <- t(data.a);
data.b <- t(data.b);
data.c <- t(data.c);

simple.data <- list(data.a = data.a, data.b = data.b, data.c = data.c);
lapply(simple.data, dim);

usethis::use_data(simple.data);
