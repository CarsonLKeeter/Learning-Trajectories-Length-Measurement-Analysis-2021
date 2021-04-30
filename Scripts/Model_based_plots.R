#######################################
## Run "Model_construction.R" First! ##
#######################################

### From "Model_construction.R"
library(rstan)
library(brms)
library(tidyverse)
library(tidybayes)

vectors <- vectors

int_1 = vectors$Intercept[, 1]
int_2 = vectors$Intercept[, 2]
int_3 = vectors$Intercept[, 3]

b_LT = vectors$b[, 2]
b_REV = vectors$b[, 3]

intercepts <- cbind(
  int_1,
  int_2,
  int_3
)

int_slope <- data.frame(
  thres1_LT = int_1 - b_LT,
  thres2_LT = int_2 - b_LT,
  thres3_LT = int_3 - b_LT,
  thres1_REV = int_1 - b_REV,
  thres2_REV = int_2 - b_REV,
  thres3_REV = int_3 - b_REV,
  thres1_BAU = int_1,
  thres2_BAU = int_2,
  thres3_BAU = int_3
)

### Sophistication Thresholds x Condition
int_slope_long <- int_slope %>% 
  gather(
    key = "Info",
    value = "Value",
    thres1_LT:thres3_BAU,
    factor_key = T
  )

int_slope_long <- int_slope_long %>% 
  mutate(
    Thres_num = as.numeric(str_extract(Info, "(\\d)+")),
    Condition = (str_sub(Info, start = -3))
  )

int_slope_long$Condition <- factor(int_slope_long$Condition, levels = c("_LT", "REV", "BAU"))

sum_label <- int_slope_long %>% 
  group_by(Condition, Thres_num) %>% 
  summarise(
    med = median(Value),
    mean = mean(Value),
    sd = sd(Value),
    L = quantile(Value, probs = 0.025),
    U = quantile(Value, probs = 0.975)
  )

sum_label$Thres_lab <- rep(c("T1", "T2", "T3"), 3)

soph_thres_p <- ggplot(
  data = int_slope_long,
  aes(
    x = Value,
    y = Condition,
    group = Thres_num
  )
) + 
  stat_interval(
    .width = c(0.5, 0.8, 0.95),
    size = 15
  ) + 
  stat_pointinterval(  
    .width = 0,
    size = 10
  ) + 
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  geom_text(
    data = sum_label,
    aes(
      x = med, 
      y = Condition,
      label = Thres_lab
    ),
    size = 3.5*2, 
    vjust = 1.5
  ) +
  scale_color_brewer(
    palette = "Blues"
  ) + 
  theme_bw(
    base_size = 14
  )  + 
  scale_y_discrete(
    labels = c("LT", "REV", "BAU")
  ) + 
  scale_x_continuous(
    n.breaks = 10
  ) +
  labs(
    x = "Latent Sophistication",
    col = "Credible Interval"
  ) + 
  theme(
    legend.position = "bottom"
  )


## Item Random Intercepts/Effects 
item_re = rstan::extract(Model)$r_2_1

item_re_col_mean <- data.frame(
  Item = c(
    "B", 
    "C", 
    "Da",    
    "Db",    
    "E",     
    "G",    
    "H",    
    "I",   
    "J",   
    "K",    
    "L", 
    "15",
    "33", 
    "M",
    "37",
    "31.1",
    "N",  
    "P",
    "Q", 
    "50",
    "32", 
    "U", 
    "21", 
    "V",    
    "X",  
    "Z" 
  ), 
  Means = colMeans(item_re)
)

item_re_col_mean <- item_re_col_mean %>% 
  mutate(
    over = case_when(
      Means < 0 ~ 0, 
      Means > 0 ~ 1
    )
  )

item_re_p <- ggplot(
  data = item_re_col_mean,
  aes(
    x = reorder(Item, Means),
    y = Means,
    fill = as.factor(over)
  )
) + 
  geom_bar(
    stat = "identity"
  ) + 
  theme_bw(
    
  ) + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + 
  labs(
    x = "Item",
    y = "Random Intercept"
  ) + 
  scale_fill_brewer(
    palette = "Paired"
  )


## Class Random Intercept/Effects
class_re <- rstan::extract(Model)$r_1_1

class_re_col_mean = data.frame(
  Class = as.factor(1:16),
  Means = colMeans(class_re)
)

class_re_col_mean <- class_re_col_mean %>% 
  mutate(
    over = case_when(
      Means < 0 ~ 0, 
      Means > 0 ~ 1
    )
  )

class_re_p <- ggplot(
  data = class_re_col_mean,
  aes(
    x = reorder(Class, Means),
    y = Means,
    fill = as.factor(over)
  )
) + 
  geom_bar(
    stat = "identity"
  ) + 
  theme_bw(
    
  ) + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + 
  labs(
    x = "Class",
    y = "Random Intercept"
  ) + 
  scale_fill_brewer(
    palette = "Paired"
  )


## Student Random Intercepts/Effects
student_re <- rstan::extract(Model)$r_3_1

student_re_means <- data.frame(
  Student.no = 1:177,
  Means = colMeans(student_re)
)

student_re_means <- student_re_means %>% 
  mutate(
    over = case_when(
      Means < 0 ~ 0, 
      Means > 0 ~ 1
    )
  )

filtered <- student_re_means %>% 
  filter(
    Means > -0.25 & Means < 0.025
  )


student_re_p <- ggplot(
  data = student_re_means,
  aes(
    x = reorder(Student.no, Means),
    y = Means,
    fill = as.factor(over)
  )
) + 
  geom_bar(
    stat = "identity"
  ) + 
  theme_bw(
    
  ) + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + 
  labs(
    x = "Student #",
    y = "Random Intercept"
  ) + 
  scale_fill_brewer(
    palette = "Paired"
  )



## Hurdle Process Plot (Stylized)

dens <- seq(0, 1, length = 10000)
db <- dbeta(dens, 5, 5)

set.seed(0)
filtered.y <- rbeta(n = 60, 5, 5)
filtered.x <- seq(0, 1, length = 60)

dots.x = rnorm(n = 400, mean = 0, 1.15)
dots.y = rnorm(n = 400, mean = 4, .25)
dots.lab = sample(x = c(0, 1 , 2), size = 200, replace = T)

dots = as.data.frame(cbind(dots.x, dots.y, dots.lab))

soph_labs = data.frame(
  lab = c("Y = 0", "Y = 1", "Y = 2", "Y = 3"),
  x = c(-2, -1, 0, 1),
  y = c(0, 0, 0, 0)
)

library(RColorBrewer)

ggplot(
  data = NULL
) + 
  geom_line(
    aes(
      x = (dens -.5)*5,
      y = db
    ),
    size = 1.5
  ) + 
  geom_point(
    data = dots,
    aes(
      x = dots.x,
      y = dots.y,
      col = as.factor(dots.lab)
    ),
    size = 4,
    alpha = .8
  ) + 
  geom_point(
    aes(
      x = (filtered.x-.5)*5,
      y =  filtered.y + 2.3
    ),
    col = "#1B9E77",
    size = 4,
    alpha = .8
  ) + 
  geom_hline(
    yintercept = 3.05,
    linetype = "dotted",
    col = "#1B9E77"
  ) + 
  geom_hline(
    yintercept = 3.06,
    linetype = "dotdash",
    col = "#1B9E77"
  ) + 
  geom_hline(
    yintercept = 3.07,
    linetype = "dotted",
    col = "#1B9E77"
  ) + 
  geom_hline(
    yintercept = 3.08,
    linetype = "dotdash",
    col = "#1B9E77"
  ) + 
  geom_hline(
    yintercept = 3.09,
    linetype = "dotted",
    col = "#1B9E77"
  ) + 
  geom_hline(
    yintercept = 3.1,
    linetype = "dotdash",
    col = "#1B9E77"
  ) + 
  geom_segment(
    aes(
      x = 0, xend = 0,
      y = -Inf, yend = max((dens -.5)*5)
    ),
    linetype = "dotted"
  ) +
  geom_segment(
    aes(
      x = -1, xend = -1,
      y = -Inf, yend = 1.1
    ),
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = 1, xend = 1,
      y = -Inf, yend = 1.1
    ),
    linetype = "dashed"
  ) +
  scale_x_continuous(
    n.breaks = 6
  ) +
  geom_segment(
    aes(x = -2, y = 3.5, xend = -2, yend = 3.15),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_segment(
    aes(x = 0, y = 3.5, xend = 0, yend = 3.15),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_segment(
    aes(x = 2, y = 3.5, xend = 2, yend = 3.15),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +  
  geom_segment(
    aes(x = -1.5, y = 2.5, xend = -1, yend = 2),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +  
  geom_segment(
    aes(x = 1.5, y = 2.5, xend = 1, yend = 2),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text(
    data = soph_labs,
    aes(
      x = x + .5,
      y = y,
      label = lab
    ),
    size = 5
  ) + 
  annotate(
    geom = "label",
    x = -2.25,
    y = 4.5,
    label = "STEP 1"
  ) + 
  annotate(
    geom = "label",
    x = -2.25,
    y = 2.5,
    label = "STEP 2"
  ) +
  theme_classic(
    base_size = 14
  ) + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  ) +
  scale_shape_discrete(
    guide = F
  ) + 
  labs(
    x = "Latent Sophistication",
    y = "",
    col = ""
  ) + 
  ylim(
    c(0, 4.5)
  ) +
  xlim(
    c(-2.5, 2.5)
  ) + 
  scale_color_brewer(
    palette = "Dark2",
    labels = c(
      "Detectable 
Strategies", 
      "Non-Codable 
Strategies", 
      "Non-Detectable 
Strategies with 
Correct Response"
    )
  )




####################################################################################









