# Load data, show transitions
combined_data = readRDS(rds_file("msm_data")) %>% as_tibble
ids = combined_data$SWANID[1:50] %>% unique

transitions <-
  combined_data %>%
  dplyr::select(SWANID,VISIT,AGE_mod,STATUS) %>%
  group_by(SWANID) %>%
  mutate(from = STATUS,
         to = lead(STATUS)) %>%
  filter(!is.na(to)) %>%
  group_by(from,to) %>%
  tally 

nodes <- data.frame(
  state = c("1 Pre-menopausal", "2 Early Perimenopause", "3 Late Perimenopause", "4 Menopause"),
  x = c(0.5, 0, 1, 0.5),
  y = c(0, 0.5, 0.5, 1)
)

df <- transitions %>%
  left_join(nodes,  by = c("from" = "state")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(nodes,  by = c("to" = "state")) %>%
  rename(x_to = x, y_to = y)

df = df %>%
  mutate(y_from = case_when(from == to ~ y_from + 0.2,
                            TRUE ~ y_from))

df <- df %>%
  mutate(
    xm = (x_from + x_to) / 2,
    ym = (y_from + y_to) / 2
  )


ggplot() +
  # arrows
  geom_curve(
    data = df,
    aes(x = x_from, y = y_from,
        xend = x_to, yend = y_to,
        size = I(0.002*n)),
    colour = "blue",
    curvature = 0.05,
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  # nodes
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    size = 10
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = state),
    size = 5,
    vjust = c(2,-2,-2,-2),
    hjust = c(0.5,0.5,0.5,0.5)
  ) +
  scale_size_continuous(range = c(2, 3), name = "Number of transitions") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(-0.2,1.2),
                  ylim = c(-0.1,1.1)) +
  geom_text(
    data = df,
    aes(x = xm, y = ym, label = n),
    vjust = -0.3,
    size = 4
  )

ggsave(plot_file("Transition_Intensities.png"),width = cm(3),height = cm(3))

### Longitudinal data plots
ggplot(combined_data) +
  geom_line(aes(x = AGE_mod,y = FSH,group = SWANID),alpha = 0.2) +
  geom_smooth(aes(x = AGE_mod,y = FSH)) +
  theme_minimal() +
  labs(x = "Age at follow up",
       y = "FSH level")

ggsave(plot_file("FSH_timeplot.png"),height = cm(2),width = cm(3))

ggplot(combined_data) +
  geom_line(aes(x = AGE_mod,y = E2AVE,group = SWANID),alpha = 0.2) +
  geom_smooth(aes(x = AGE_mod,y = E2AVE)) +
  theme_minimal() +
  labs(x = "Age at follow up",
       y = "E2 level") +
  coord_cartesian(ylim = c(0,300))

ggsave(plot_file("E2_timeplot.png"),height = cm(2),width = cm(3))
