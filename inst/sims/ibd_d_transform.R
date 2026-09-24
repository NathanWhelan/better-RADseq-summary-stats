###############################################################################
#
#  inst/sims/ibd_d_transform.R -- should Jost's D be used as it is, or as
#  D/(1 - D), in a test for isolation by distance?
#
#  Rousset (1997) showed that FST/(1 - FST) rises in a straight line with
#  distance (1-D habitat) or with ln(distance) (2-D habitat). For a pair of
#  populations, write Q0 for the chance that two genes from one population
#  are identical and Qd for two genes from populations d steps apart. Then
#
#    FST/(1 - FST) = (Q0 - Qd)/(1 - Q0)
#    D             = (Q0 - Qd)/Q0
#    D/(1 - D)     = (Q0 - Qd)/Qd
#
#  so D = FST/(1 - FST) x Hs/(1 - Hs), with Hs = 1 - Q0 the same for every
#  pair. This script checks that with EXACT identity probabilities -- no
#  sampling noise -- in stepping-stone models, and measures how straight each
#  statistic is: `bend` is the slope over the far half of the distances
#  divided by the slope over the near half (1 = a straight line).
#
#  The model: a circular lattice of L demes (L x L in 2-D), N diploids per
#  deme, a fraction m of genes moving to a neighbouring deme each generation,
#  and infinite-alleles mutation at rate u. The lattices are large, so the
#  circle hardly bends the curves at distances 2-20. The stationary identity
#  probabilities are solved exactly in Fourier space (Malecot's recursion):
#  with phi the Fourier transform of one gene's step and k = (1 - u)^2,
#    Q(d) = c g(d),  g = inverse transform of k phi^2 / (1 - k phi^2),
#    Q0 = g(0) / (2N + g(0)),  c = (1 - Q0) / (2N).
#
#  Results are quoted in vignette("rationale"), "Isolation by distance".
#  Run from the package source directory:  Rscript inst/sims/ibd_d_transform.R
#  (a few seconds; nothing random).
#
###############################################################################

identity_by_distance <- function(dims, L, N, m, u) {
  theta <- 2 * pi * (0:(L - 1)) / L
  phi <- if (dims == 1) (1 - m) + m * cos(theta)
         else (1 - m) + (m / 2) * outer(cos(theta), cos(theta), "+")
  k <- (1 - u)^2
  g <- Re(stats::fft(k * phi^2 / (1 - k * phi^2), inverse = TRUE)) / length(phi)
  Q0 <- g[1] / (2 * N + g[1])
  g * (1 - Q0) / (2 * N)
}

bend <- function(x, y) {
  half <- floor(length(x) / 2)
  near <- seq_len(half)
  far <- (half + 1):length(x)
  unname(stats::coef(stats::lm(y[far] ~ x[far]))[2] / stats::coef(stats::lm(y[near] ~ x[near]))[2])
}

d <- 2:20
settings <- expand.grid(N_m = c("25/0.2", "5/0.05", "2/0.02"), u = c(1e-7, 1e-5),
                        habitat = c("1-D", "2-D"), stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(settings)), function(i) {
  s <- settings[i, ]
  N <- as.numeric(sub("/.*", "", s$N_m))
  m <- as.numeric(sub(".*/", "", s$N_m))
  if (s$habitat == "1-D") {
    Q <- identity_by_distance(1, 20000, N, m, s$u)
    Q0 <- Q[1]; Qd <- Q[d + 1]; x <- d
  } else {
    Q <- identity_by_distance(2, 1024, N, m, s$u)
    Q0 <- Q[1, 1]; Qd <- Q[d + 1, 1]; x <- log(d)
  }
  lin <- (Q0 - Qd) / (1 - Q0)
  D <- (Q0 - Qd) / Q0
  data.frame(habitat = s$habitat, N = N, m = m, u = s$u, Hs = round(1 - Q0, 3),
             D_from = round(min(D), 3), D_to = round(max(D), 3),
             bend_FST_linearized = round(bend(x, lin), 3), bend_D = round(bend(x, D), 3),
             bend_D_over_1_minus_D = round(bend(x, D / (1 - D)), 3),
             D_over_FST_linearized = signif(stats::median(D / lin), 4),
             Hs_over_1_minus_Hs = signif((1 - Q0) / Q0, 4))
})
print(do.call(rbind, rows), row.names = FALSE)
