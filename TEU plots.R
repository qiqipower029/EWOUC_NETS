##------------------------------- Generating plots for T/E/U -------------------------------##
##------------------------------- Jieqi Tu -------------------------------##


# Read in datasets
library(tidyverse)
scenario1 = readxl::read_excel(path = "./Results/EWOUC-NETS-s1-newW2.xlsx")
scenario2 = readxl::read_excel(path = "./Results/EWOUC-NETS-s2-newW2.xlsx")
scenario3 = readxl::read_excel(path = "./Results/EWOUC-NETS-s3-newW2.xlsx")
scenario4 = readxl::read_excel(path = "./Results/EWOUC-NETS-s4-newW2.xlsx")
scenario5 = readxl::read_excel(path = "./Results/EWOUC-NETS-s5-newW2.xlsx")


# -------------------- Scenario 1 -------------------- #
s1 = scenario1 %>% 
  ggplot2::ggplot(aes(x = std.dose, y = s1.pe)) + geom_line(alpha = 0.6, color = "red") + geom_point(alpha = 0.3, color = "red") + geom_line(aes(x = std.dose, y = s1.pt), alpha = 0.6, color = "blue") + 
  geom_point(aes(x = std.dose, y = s1.pt), alpha = 0.3, color = "blue") + geom_line(aes(x = std.dose, y = d.u), alpha = 0.6, color = "purple") + 
  geom_point(aes(x = std.dose, y = d.u), alpha = 0.3, color = "purple") + geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + geom_text(x = 0.4, y = -0.02, label = "MED") + geom_text(x = 1, y = -0.01, label = "MTD") + 
  geom_text(x = 0.6, y = 0.8017, label = "Target dose") +
  theme_bw() + labs(
    x = "Dose",
    y = "T/E/U",
    title = "Scenario 1"
  )

# -------------------- Scenario 2 -------------------- #
s2 = scenario2 %>% 
  ggplot2::ggplot(aes(x = std.dose, y = s1.pe)) + geom_line(alpha = 0.6, color = "red") + geom_point(alpha = 0.3, color = "red") + geom_line(aes(x = std.dose, y = s1.pt), alpha = 0.6, color = "blue") + 
  geom_point(aes(x = std.dose, y = s1.pt), alpha = 0.3, color = "blue") + geom_line(aes(x = std.dose, y = d.u), alpha = 0.6, color = "purple") + 
  geom_point(aes(x = std.dose, y = d.u), alpha = 0.3, color = "purple") + geom_vline(xintercept = 0.8, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + geom_text(x = 0.4, y = -0.63, label = "MED") + geom_text(x = 0.8, y = -0.63, label = "MTD") + 
  geom_text(x = 0.6, y = 0.7012, label = "Target dose") +
  theme_bw() + labs(
    x = "Dose",
    y = "T/E/U",
    title = "Scenario 2"
  )

# -------------------- Scenario 3 -------------------- #
s3 = scenario3 %>% 
  ggplot2::ggplot(aes(x = std.dose, y = s1.pe)) + geom_line(alpha = 0.6, color = "red") + geom_point(alpha = 0.3, color = "red") + geom_line(aes(x = std.dose, y = s1.pt), alpha = 0.6, color = "blue") + 
  geom_point(aes(x = std.dose, y = s1.pt), alpha = 0.3, color = "blue") + geom_line(aes(x = std.dose, y = d.u), alpha = 0.6, color = "purple") + 
  geom_point(aes(x = std.dose, y = d.u), alpha = 0.3, color = "purple") + geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0.6, linetype = "dashed") + geom_text(x = 0.6, y = -0.03, label = "MED") + geom_text(x = 1, y = -0.03, label = "MTD") + 
  geom_text(x = 0.8, y = 0.3, label = "Target dose") +
  theme_bw() + labs(
    x = "Dose",
    y = "T/E/U",
    title = "Scenario 3"
  )

# -------------------- Scenario 4 -------------------- #
s4 = scenario4 %>% 
  ggplot2::ggplot(aes(x = std.dose, y = s1.pe)) + geom_line(alpha = 0.6, color = "red") + geom_point(alpha = 0.3, color = "red") + geom_line(aes(x = std.dose, y = s1.pt), alpha = 0.6, color = "blue") + 
  geom_point(aes(x = std.dose, y = s1.pt), alpha = 0.3, color = "blue") + geom_line(aes(x = std.dose, y = d.u), alpha = 0.6, color = "purple") + 
  geom_point(aes(x = std.dose, y = d.u), alpha = 0.3, color = "purple") + geom_vline(xintercept = 0.8, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + geom_text(x = 0.8, y = -1.88, label = "MED") + geom_text(x = 0.4, y = -1.88, label = "MTD") + 
  geom_text(x = 0.6, y = 0.3, label = "No target dose") +
  theme_bw() + labs(
    x = "Dose",
    y = "T/E/U",
    title = "Scenario 4"
  )


# -------------------- Scenario 5 -------------------- #
s5 = scenario5 %>% 
  ggplot2::ggplot(aes(x = std.dose, y = s1.pe)) + geom_line(alpha = 0.6, color = "red") + geom_point(alpha = 0.3, color = "red") + geom_line(aes(x = std.dose, y = s1.pt), alpha = 0.6, color = "blue") + 
  geom_point(aes(x = std.dose, y = s1.pt), alpha = 0.3, color = "blue") + geom_line(aes(x = std.dose, y = d.u), alpha = 0.6, color = "purple") + 
  geom_point(aes(x = std.dose, y = d.u), alpha = 0.3, color = "purple") + geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + geom_text(x = 1, y = -1.93, label = "MED") + geom_text(x = 0.4, y = -1.93, label = "MTD") + 
  geom_text(x = 0.6, y = 0.3, label = "No target dose") +
  theme_bw() + labs(
    x = "Dose",
    y = "T/E/U",
    title = "Scenario 5"
  )

# -------------------- Joint plot -------------------- #
library(ggpubr)

ggarrange(s1, s2, s3, s4, s5,
          ncol = 2, nrow = 3)
