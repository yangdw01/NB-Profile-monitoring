###################################################################################################
# Basic setting
# -------------------------------------------------------------------------------------------------

# directory ---------------------------------------------------------------------------------------
setwd("C:\\Users\\yangd\\OneDrive\\0. Research\\2. JQT\\Revision\\JQT\\")
# -------------------------------------------------------------------------------------------------

# library -----------------------------------------------------------------------------------------
library(dplyr)
library(reshape2)
library(tidyr)
library(dqrng)
library(LaplacesDemon)
library(fields)
library(foreach)
library(doParallel)
library(ggplot2)
library(patchwork)
library(gridExtra)  
library(viridis)    
library(grid)
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# 
###################################################################################################



##############################################################################################################
# figure 2
# ------------------------------------------------------------------------------------------------------------
# data
set.seed(5)
cl1_dat <- rmvnorm(20, mean = c(0, 0), sigma = 0.2 * diag(2))
cl2_dat <- rmvnorm(10, mean = c(2, 1), sigma = 0.2 * matrix(c(1, 0.1, 0.1, 1), nrow = 2))
all_dat <- rbind(cl1_dat, cl2_dat)
all_dat[14, 1] <- 0.9  # 14번째 행의 x값 수정

# data frame
df <- data.frame(x = all_dat[, 1],
                 y = all_dat[, 2],
                 Cluster = factor(c(rep("black", 20), rep("red", 10))))

# basic ggplot object
p_base <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = Cluster), size = 4, alpha = 0.5) +
  scale_color_manual(values = c("black" = adjustcolor("black", alpha.f = 0.5),
                                "red"   = adjustcolor("red", alpha.f = 0.5))) +
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 30),      
    axis.title.x = element_text(size = 16),                 
    axis.title.y = element_text(size = 16),                 
    legend.title = element_text(size = 20),                 
    legend.text = element_text(size = 12)                   
  )

# plot 1
p1 <- p_base +
  labs(title = "(a)")
print(p1)

# plot 2
annot_df1 <- data.frame(index = c(20, 24, 14),
                        label = c("1", "2", "3"),
                        color = c("black", "red", "red"))
annot_df1$x <- df$x[annot_df1$index]
annot_df1$y <- df$y[annot_df1$index]
annot_black1 <- subset(annot_df1, color == "black")
annot_red1   <- subset(annot_df1, color == "red")

p2 <- p_base +
  geom_text(data = annot_black1, aes(x = x, y = y, label = label),
            color = "black", vjust = -1, size = 4) +
  geom_text(data = annot_red1, aes(x = x, y = y, label = label),
            color = "red", vjust = -1, size = 4) +
  labs(title = "(b)")
print(p2)

# plot 3
annot_df <- data.frame(index = c(20, 24, 14,
                                 2, 5, 13, 11, 12, 18, 17,
                                 25, 23, 29, 21, 26),
                       label = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                              11, 12, 13, 14, 15)),
                       color = c("black", "red", "red",
                                 rep("black", 7),
                                 rep("red", 5)))
annot_df$x <- df$x[annot_df$index]
annot_df$y <- df$y[annot_df$index]
annot_black <- subset(annot_df, color == "black")
annot_red   <- subset(annot_df, color == "red")

p3 <- p_base +
  geom_text(data = annot_black, aes(x = x, y = y, label = label),
            color = "black", vjust = -1, size = 4) +
  geom_text(data = annot_red, aes(x = x, y = y, label = label),
            color = "red", vjust = -1, size = 4) +
  labs(title = "(c)")
print(p3)


# plot 4
p4 <- p_base +
  geom_text(data = annot_black, aes(x = x, y = y, label = label),
            color = "black", vjust = -1, size = 4) +
  geom_text(data = annot_red, aes(x = x, y = y, label = label),
            color = "red", vjust = -1, size = 4) +
  labs(title = "(c)")

annot_df_not3 <- annot_df[annot_df$label != "3", ]
annot_df3     <- annot_df[annot_df$label == "3", ]

p4 <- p4 + annotate("point", x = annot_df3$x-0.03, y = 0.78, 
                    shape = 21, size = 10, colour = "blue", fill = NA, stroke = 1.5)

p4 <- p4 + annotate("curve", 
                    x = annot_df3$x, y = annot_df3$y, 
                    xend = annot_df3$x - 0.4, yend = annot_df3$y + 0.4, 
                    arrow = arrow(length = unit(0.2, "cm")), 
                    colour = "blue", curvature = -0.2)
p4 <- p4 + annotate("text", 
                    x = annot_df3$x, y = annot_df3$y + 0.75, 
                    label = "Backward elimination\n is required!", 
                    colour = "blue", size = 4.5, hjust = 1)
print(p4)

# combine plots
combined <- (p1 + p2 + p4) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
print(combined)

# ggsave("figure2.png", plot = combined, width = 16, height = 4, dpi = 1000)
# ------------------------------------------------------------------------------------------------------------
# 
##############################################################################################################



 
##############################################################################################################
# figure 4
# ------------------------------------------------------------------------------------------------------------

# Zernike polynomial 1
zernike_radial <- function(n, m, rho) {
  R <- rep(0, length(rho))
  for (s in 0:((n - abs(m)) %/% 2)) {
    coeff <- (-1)^s * factorial(n - s) /
      ( factorial(s) *
          factorial((n + abs(m)) / 2 - s) *
          factorial((n - abs(m)) / 2 - s) )
    R <- R + coeff * rho^(n - 2 * s)
  }
  return(R)
}

# Zernike polynomial 2
zernike <- function(n, m, rho, theta) {
  R <- zernike_radial(n, m, rho)
  if (m >= 0) {
    return(R * cos(m * theta))
  } else {
    return(R * sin(abs(m) * theta))
  }
}

# grid
N <- 200  
x <- seq(-1, 1, length.out = N)
y <- seq(-1, 1, length.out = N)
grid_df <- expand.grid(x = x, y = y)
grid_df$rho <- sqrt(grid_df$x^2 + grid_df$y^2)
grid_df$theta <- atan2(grid_df$y, grid_df$x)
grid_df <- grid_df[grid_df$rho <= 1, ]

modes <- list()
for (n in 0:4) {
  for (m in seq(-n, n)) {
    if ((n - abs(m)) %% 2 == 0) {
      modes[[length(modes) + 1]] <- c(n, m)
    }
  }
}
print(length(modes))  

# ggplot
plot_list <- list()
for (i in seq_along(modes)) {
  mode <- modes[[i]]
  n <- mode[1]
  m <- mode[2]
  
  # Zernike
  z_val <- zernike(n, m, grid_df$rho, grid_df$theta)
  df <- data.frame(x = grid_df$x, y = grid_df$y, z = z_val)
  
  # subtitles
  label <- sprintf("(%s) n=%d, m=%d", letters[i], n, m)
  
  if(n == 0 && m == 0) {
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_raster(fill = "#F05650") +
      coord_equal() +
      theme_minimal() +
      labs(subtitle = label) +
      theme(
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
      )
  } else {
    
    p <- ggplot(df, aes(x = x, y = y, fill = z)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      coord_equal() +
      theme_minimal() +
      labs(subtitle = label) +
      theme(
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
      )
  }
  
  plot_list[[i]] <- p
}

# 4 x 4 
combined_plot <- arrangeGrob(grobs = plot_list, nrow = 4, ncol = 4)

# save
output_file <- "figure2.png"
png(filename = output_file, width = 3000, height = 3000, res = 300)
grid.newpage()
grid.draw(combined_plot)
dev.off()
# ------------------------------------------------------------------------------------------------------------
# 
##############################################################################################################
