
PT_Pallette <- list(
  Colors = c("#AA4643", "#4572A7", "#89A54E", "#71588F", "#4198AF", "#DB843D",
             "#93A9CF", "#D19392", "#B9CD96", "#A99BBD", "black", "grey", "white"),
  Tumor_size = setNames(c("#C0504D", "#f2dcdb", "#868686"), c(">5cm", "<5cm", "Unknown")),
  Age = setNames(c("#6360c3", "#cdccf0", "#868686"), c("<40", ">40", "Unknown")),
  Batch = setNames(c("#AA4643", "#DB843D", "#89A54E", "#4198AF", "#4572A7", "#71588F"),
                   c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6")),
  Death = setNames(c("#b3a2c7", "#e6e0ec"), c("Yes", "No")),
  Metastasis = setNames(c("#31859c", "#b7dee8"), c("Yes", "No")),
  Recurrence = setNames(c("#2685c5", "#AED6F1"), c("Yes", "No")),
  Pathological_classification = setNames(c("#e0b506", "#4198AF", "#A99BBD", "#89A54E"),
                                         c("Malignant", "Borderline", "Benign", "Normal")),
  RNA_supersubtype = setNames(c("#AA4643", "#4572A7", "#89A54E"), c("MM","MB", "Normal")),
  H3K27_subtype = setNames(c("#AA4643", "#D19392"), c("H3K27_H", "H3K27_L")),
  CNA_subtype = setNames(c("#8c564bFF", "#e377c2FF"), c("CNA_H", "CNA_L")),
  SNV_subtype=setNames(c("#006A8E","#7A3A9A","#B1283A","#3F86BC","#28ADA8"),c( "C1","C2","C3","C4","C5")),
  RNA_subtype = setNames(c("#AA4643", "#DB843D", "#71588F", "#4572A7", "#89A54E"),
                         c("HI","MI","LR","NL","Normal"))
)



PT_Pallette$Phase=setNames(c("#46B8DA","#6B58EF","#9632B8"),c("G1","S","G2M"))

PT_Pallette$Mut=setNames(c("lightgrey","steelblue"),c("Wildtype","Mut"))
PT_Pallette$Survival=setNames(c("gray40","#225EA8", "#EC7014", "#6A51A3" ,"#CB181D"),c("No","R","M-R","D-M","D-M-R"))
PT_Pallette$Progress=setNames(c("lightgrey","#225EA8", "#EC7014", "#6A51A3" ,"#CB181D"),c("No","R","M-R","D-M","D-M-R"))

library(paletteer)