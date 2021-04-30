### Exploratory Plots 

fig.dat <- read.csv("LT_length_data.csv")
fig.dat <- fig.dat[, -c(1, 3, 9, 10)]
fig.dat$Soph_post <- as.factor(fig.dat$Soph_post)
fig.dat$SID <- as.factor(fig.dat$SID)
fig.dat$Private <- as.factor(fig.dat$Private)
fig.dat$Class <- as.factor(fig.dat$Class)
fig.dat$Sex <- as.factor(fig.dat$Sex)
fig.dat$Sex <- revalue(fig.dat$Sex, replace = c("M" = "Boys", "F" = "Girls"))
fig.dat$Condition <- factor(x = fig.dat$Condition, levels = c("LT", "REV", "BAU"))
fig.dat$Item <- as.factor(fig.dat$Item)
fig.dat$Soph_post <- as.factor(as.ordered(fig.dat$Soph_post))
fig.dat$Correct_post <- as.factor(fig.dat$Correct_post)
fig.dat$Item <- factor(
  x = fig.dat$Item, 
  levels = c(
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
    "X15",
    "X33", 
    "M",
    "X37",
    "X31.1",
    "N",  
    "P",
    "Q", 
    "X50",
    "X32", 
    "U", 
    "X21", 
    "V",    
    "X",  
    "Z" 
  )
)

item_labs = c(
  "B" = "B", 
  "C" = "C", 
  "Da"= "Da",    
  "Db" = "Db",    
  "E" = "E",     
  "G" = "G",    
  "H" = "H",    
  "I" = "I",   
  "J" = "J",   
  "K" = "K",    
  "L" = "L", 
  "X15" = "15",
  "X33" = "33", 
  "M" = "M",
  "X37" = "37",
  "X31.1" = "31.1",
  "N" = "N",  
  "P" = "P",
  "Q" = "Q", 
  "X50" = "50",
  "X32" = "32", 
  "U" = "U", 
  "X21" = "21", 
  "V" = "V",    
  "X" = "X",  
  "Z" = "Z" 
)


## Post-assessment Sophistication x Condition 

p0 <- ggplot(
  data = na.omit(fig.dat),
  aes(
    x = Condition,
    fill = Soph_post
  )
) + 
  geom_bar(
    position="fill"
  ) + 
  labs(
    x = "Intervention condition",
    y = "Relative frequency",
    fill='Post-assessment sophistication'
  )  + 
  theme_bw(
    
  )+
  theme(
    legend.position="top"
  )  + 
  scale_fill_brewer(
    palette = "BuPu"
  ) 

## Post-assessment Sophistication x Item x Condition 

p1 <- ggplot( 
  data = na.omit(fig.dat),
  aes(
    x = Condition,
    fill = Soph_post
  )
) + 
  geom_bar(
    position="fill"
  ) + 
  facet_wrap(
    ~Item,
    as.table=TRUE,
    nrow=2,
    labeller = as_labeller(item_labs)
  ) +
  labs(
    x = "Intervention condition by Item",
    y = "Relative frequency",
    fill='Post-assessment sophistication'
  )  + 
  theme_bw(
    
  )+
  theme(
    legend.position="top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )  + 
  scale_fill_brewer(
    palette = "BuPu"
  ) 

##  Post-assessment Sophistication x Boy/Girl x Condition

p2 <- ggplot(
  data = na.omit(fig.dat),
  aes(
    x = Condition,
    fill = Soph_post
  )
) + 
  geom_bar(
    position="fill"
  ) + 
  facet_wrap(
    ~Sex,
    as.table=TRUE,
    nrow=2
  ) +
  labs(
    x = "Intervention condition",
    y = "Relative frequency",
    fill='Post-assessment sophistication'
  )  + 
  theme_bw(
    
  )+
  theme(
    legend.position="top"
  )  + 
  scale_fill_brewer(
    palette = "BuPu"
  )


## Post-assessment Sophistication x Public/Private x Condition

p3 <- ggplot(
  data = na.omit(dat),
  aes(
    x = Condition,
    fill = Soph_post
  )
) + 
  geom_bar(
    position="fill"
  ) + 
  facet_wrap(
    ~Private,
    as.table=TRUE,
    nrow=2
  ) +
  labs(
    x = "Intervention condition",
    y = "Relative frequency",
    fill='Post-assessment sophistication'
  )  + 
  theme_bw(
    
  )+
  theme(
    legend.position="top"
  )  + 
  scale_fill_brewer(
    palette = "BuPu"
  )

## Post-assessment Sophistication x Classroom x Condition

p4 <- ggplot(
  data = na.omit(fig.dat),
  aes(
    x = Condition,
    fill=Soph_post
  )
) + 
  geom_bar(
    position="fill"
  ) + 
  facet_wrap(
    vars(Class),
    nrow=2,
    as.table=TRUE
  ) +
  labs(
    x = "Intervention condition by Classroom",
    y = "Relative frequency", 
    fill = "Post-assessment sophistication"
  ) +
  theme_bw(
    base_size = 14
  ) + 
  theme(
    legend.position="top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )  + 
  scale_fill_brewer(
    palette = "BuPu"
  ) 


####################################################################################

