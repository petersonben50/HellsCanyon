#### code/HCC_plotting_needs.R ####
# Benjamin D. Peterson


#### Colorblind friendly colors ####
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Generate needed vectors ####
color.vector <- cb.translator[c("bluishgreen", "reddishpurple", "orange", "black", "blue")]
names(color.vector) <- c("oxic", "suboxic", "no_nitrate_no_sulfide", "no_nitrate_possible_sulfide", "sulfidic")

shape.vector <- c(18, 16, 17, 15)
names(shape.vector) <- c(2016, 2017, 2018, 318)

renaming.vector <- c("DO > 0.5 mg/L",
                     "DO < 0.5 mg/L, nitrate > 0.05 mgN/L",
                     "Nitrate < 0.05 mgN/L, sulfide < 0.01 mg/L",
                     "Nitrate < 0.05 mgN/L, sulfide not measured",
                     "Sulfide > 0.01 mg/L")
names(renaming.vector) <- names(color.vector)
