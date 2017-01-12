input <- "/Users/mcdonald/Dropbox/michor/projects/SIApop/SIApop/examples/constant-rate/passenger/inputdata.txt"
anc <- "/Users/mcdonald/Dropbox/michor/projects/SIApop/SIApop/examples/constant-rate/passenger/ancestors.txt"
out <- "/Users/mcdonald/Dropbox/michor/projects/SIApop/SIApop/examples/constant-rate/passenger/"

siapopr::siapopConstant(input = input, output_dir = out, ancestor_file = anc)
siapopr::siapopSimple(c(input, anc, out))

.Call('siapopr_siapopConstant', PACKAGE = 'siapopr', input, out, anc)

library(siapopr)
siapopConstant(observation_frequency = 1)


test <- siapopConstant(max_pop = 10000, birth_rate = 1.2, death_rate = 1, mutation_prob = 0.01,
               allow_extinction = F, num_samples = 1, sample_size = 100,
               observation_frequency = 1, num_sims = 3)

dat <- import_siapop('./')
head(dat)
tail(dat$clone_data[[1]])

samp <- dat$data[[1]]$sample_data

adj_matrix <- create_sample_adj_matrix(samp)

library(phangorn)
samptree <- as.phyDat(adj_matrix, type = "USER", levels = c(T, F))
plot(NJ(dist.hamming(samptree)))


library(ape)
samptree <- nj(dist(adj_matrix)^2)
plot(samptree)




