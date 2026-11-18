require(here)
library(ggplot2)
source(here("src/Packages.R"))
source(here("src/Utils.R"))
load(here("Simulazioni/Hearts/no warp/sym/no_warp_sym_01.RData"))

mu     <- sam$mu
k      <- nrow(mu)
S_mat  <- S_mat_constructor(mu)
D_p    <- D_p_constructor(mu, k); k_free <- ncol(D_p)
K      <- crossprod(D_p, S_mat %*% D_p)
V      <- median(diag(S_mat %*% D_p %*% solve(K) %*% t(D_p) %*% S_mat))

# warp levels: sigma2 = 0 is the original shape
sigma2s <- c(0, 7e-4, 7e-3, 7e-2)
# or choose the levels on the prior scale sigma2/V:
# ratio_levels <- c(0, 1e-5, 1e-4, 1e-3); sigma2s <- ratio_levels * V
data.frame(sigma2 = sigma2s, sigma2_over_V = sigma2s / V)

# one fixed draw, so the panels are the same warp at increasing amplitude
set.seed(12345)
psi_unit <- backsolve(chol(K), matrix(rnorm(k_free * 2), k_free, 2))   # psi ~ N(0, K^-1) per column

# TPS kernel r^2 log r between two point sets (same as S_mat_constructor)
kappa_mat <- function(P, Q) {
  d2 <- outer(P[, 1], Q[, 1], "-")^2 + outer(P[, 2], Q[, 2], "-")^2
  0.5 * d2 * log(d2 + (d2 == 0))
}

# folding check: does the closed polygon cross itself?
seg_cross <- function(p1, p2, p3, p4) {
  o <- function(a, b, c) (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])
  o(p1,p2,p3) * o(p1,p2,p4) < 0 && o(p3,p4,p1) * o(p3,p4,p2) < 0
}
self_intersects <- function(P) {
  n <- nrow(P)
  for (i in 1:n) for (j in 1:n) {
    if (abs(i - j) <= 1 || (i == 1 && j == n) || (j == 1 && i == n)) next
    if (seg_cross(P[i, ], P[i %% n + 1, ], P[j, ], P[j %% n + 1, ])) return(TRUE)
  }
  FALSE
}

# regular grid over the bounding box
pad <- 0.25 * max(diff(range(mu[, 1])), diff(range(mu[, 2])))
nx <- ny <- 25
gx <- seq(min(mu[, 1]) - pad, max(mu[, 1]) + pad, length.out = nx)
gy <- seq(min(mu[, 2]) - pad, max(mu[, 2]) + pad, length.out = ny)
grid_pts <- as.matrix(expand.grid(x = gx, y = gy))               # x varies fastest
ij <- expand.grid(ix = 1:nx, iy = 1:ny)


make_level <- function(sigma2, lab) {
  W  <- D_p %*% (sqrt(sigma2) * psi_unit)
  SW <- S_mat %*% W
  gw <- grid_pts + kappa_mat(grid_pts, mu) %*% W
  m  <- mu + SW
  folded <- self_intersects(m)
  lev <- sprintf("%s\nsigma2 = %.3g | sigma2/V = %.3g%s", lab, sigma2, sigma2 / V,
                 if (folded) " | FOLDED" else "")
  list(
    grid = rbind(data.frame(x = gw[, 1], y = gw[, 2], line = paste0("h", ij$iy), level = lev),
                 data.frame(x = gw[, 1], y = gw[, 2], line = paste0("v", ij$ix), level = lev)),
    shape = data.frame(x = c(m[, 1], m[1, 1]), y = c(m[, 2], m[1, 2]), level = lev),
    pts   = data.frame(x = m[, 1], y = m[, 2], xo = mu[, 1], yo = mu[, 2], lm = 1:k, level = lev))
}


# inside make_level, in place of the old lev line:
lev <- sprintf("atop(sigma^2 == %s, sigma^2 / omega(eta) == %s)", sigma2s, sigma2s / V)


lv <- Map(make_level, sigma2s, ifelse(sigma2s == 0, "Original", "Warped"))
lev_order <- vapply(lv, function(z) z$shape$level[1], "")
bind <- function(nm) { d <- do.call(rbind, lapply(lv, `[[`, nm)); d$level <- factor(d$level, lev_order); d }
df_grid <- bind("grid")
levels(df_grid$level) <- lev
df_shape <- bind("shape")
levels(df_shape$level) <- lev
df_pts <- bind("pts")
levels(df_pts$level) <- lev

shape0 <- data.frame(x = c(mu[, 1], mu[1, 1]), y = c(mu[, 2], mu[1, 2]))

p <- ggplot() +
  geom_path(data = df_grid, aes(x, y, group = interaction(level, line)),
            colour = "grey70", linewidth = 0.25) +
  geom_path(data = shape0, aes(x, y), colour = "grey40", linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = df_pts, aes(x = xo, y = yo, xend = x, yend = y), colour = "red", linewidth = 0.4) +
  geom_path(data = df_shape, aes(x, y), colour = "forestgreen", linewidth = 0.9) +
  geom_point(data = df_pts, aes(x, y), shape = 21, fill = "darkgreen", size = 1.6) +
  facet_wrap(~ level, nrow = 1, labeller = label_parsed) +
  coord_equal(xlim = range(gx), ylim = range(gy)) +
  labs(title = "") +
  theme_void(base_size = 14) +
  theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') )



p

ggsave(here("Simulazioni/warp_levels.pdf"), p,
       device = cairo_pdf, bg = "transparent",
       width = 23, height = 10.5, units = "cm")


###### Plot slides ###### 
# ---- figure for the slide "The warping": plane before / after ----
s2 <- 3e-2     # warping level: pick one from sigma2s (3e-2 is visible but not folded; try 3e-1 to see folding)

# same warp as in make_level
W   <- D_p %*% (sqrt(s2) * psi_unit)
cat("sigma2 =", s2, " sigma2/V =", s2 / V, " folded:", self_intersects(mu + S_mat %*% W), "\n")

# dense closed outline through the k landmarks (periodic spline), to show the continuous curve
t_id  <- c(1:k, k + 1)
dense <- cbind(spline(t_id, c(mu[, 1], mu[1, 1]), n = 400, method = "periodic")$y,
               spline(t_id, c(mu[, 2], mu[1, 2]), n = 400, method = "periodic")$y)

# the warp acts on ANY point of the plane: x -> x + kappa(x, mu) W
warp <- function(P) P + kappa_mat(P, mu) %*% W
pad <- 0.25 * max(diff(range(mu[, 1])), diff(range(mu[, 2])))
nx <- ny <- 25
gx <- seq(min(mu[, 1]) - pad, max(mu[, 1]) + pad, length.out = nx)
gy <- seq(min(mu[, 2]) - pad, max(mu[, 2]) + pad, length.out = ny)
grid_pts <- as.matrix(expand.grid(x = gx, y = gy))               # x varies fastest
ij <- expand.grid(ix = 1:nx, iy = 1:ny)
gw  <- warp(grid_pts)       # warped grid
dw  <- warp(dense)          # warped outline
lmw <- mu + S_mat %*% W     # warped landmarks (= warp(mu))


lev <- c("Original plane", "Warped plane")
df_grid <- rbind(
  data.frame(x = grid_pts[, 1], y = grid_pts[, 2], line = paste0("h", ij$iy), level = lev[1]),
  data.frame(x = grid_pts[, 1], y = grid_pts[, 2], line = paste0("v", ij$ix), level = lev[1]),
  data.frame(x = gw[, 1],       y = gw[, 2],       line = paste0("h", ij$iy), level = lev[2]),
  data.frame(x = gw[, 1],       y = gw[, 2],       line = paste0("v", ij$ix), level = lev[2]))
df_curve <- rbind(data.frame(x = dense[, 1], y = dense[, 2], level = lev[1]),
                  data.frame(x = dw[, 1],    y = dw[, 2],    level = lev[2]))
df_lm    <- rbind(data.frame(x = mu[, 1],  y = mu[, 2],  level = lev[1]),
                  data.frame(x = lmw[, 1], y = lmw[, 2], level = lev[2]))
for (nm in c("df_grid", "df_curve", "df_lm")) {
  d <- get(nm); d$level <- factor(d$level, lev); assign(nm, d)
}

p2 <- ggplot() +
  geom_path(data = df_grid, aes(x, y, group = interaction(level, line)),
            colour = "grey70", linewidth = 0.25) +
  geom_path(data = df_curve, aes(x, y), colour = "#B8472F", linewidth = 1.1) +
  geom_point(data = df_lm, aes(x, y), shape = 21, fill = "#1B2230", colour = "white", size = 1.8) +
  facet_wrap(~ level, nrow = 1) +
  coord_equal(xlim = range(gx), ylim = range(gy)) +
  theme_void(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        panel.background  = element_rect(fill = "transparent", colour = NA),
        plot.background   = element_rect(fill = "transparent", colour = NA))

p2
ggsave(here("figs/warping_plane.pdf"), p2,
       device = cairo_pdf, bg = "transparent", width = 20, height = 8.5, units = "cm")
