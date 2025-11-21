library(seqinr)
library(readr)
library(readxl)
setwd('c:/Users/Alvaro/Documents/Oropouche/6_Mices_filt/')

#Crear un df con los epítopos anidados de los HLA
HLA <- data.frame(read_delim('Input/Human/Epitopes_Global_Nested.csv', 
                             delim= ',', show_col_types = FALSE))
HLA$Peptide_H2_1 <- rep('-', nrow(HLA))
HLA$Mices_H2_1 <- rep('-', nrow(HLA))
HLA$Mices_H2_2 <- rep('-', nrow(HLA))

H2_1 <- data.frame(read_delim('Input/Mices/Merged_mice_MHC1.txt', delim = '\t',
                              show_col_types = FALSE))

H2_2 <- data.frame(read_delim('Input/Mices/Merged_mice_MHC2.txt', delim = '\t',
                              show_col_types = FALSE))

#Indicar los péptidos que tienen afinidad por H2 de tipo 1
for (i in 1:nrow(HLA)) {
  Epitopes <- c()
  Alleles <- c()
  for (j in 1:nrow(H2_1)) {
    if (grepl(H2_1$Peptide[j], HLA$Peptide[i]) > 0) {
      for (k in 5:ncol(H2_1)) {
        Epitopes <- c(Epitopes, H2_1$Peptide[j])
        Epitopes <- unique(Epitopes)
        HLA$Peptide_H2_1[i] <- toString(Epitopes)
        if (!is.na(H2_1[j,k]) == TRUE) {
          Alleles <- c(Alleles, colnames(H2_1[k]))
          Alleles <- unique(Alleles)
          HLA$Mices_H2_1[i] <- toString(Alleles)
        }
      }
    }
  }
}

#Indicar los péptidos que tienen afinidad por H2 de tipo 2
for (i in 1:nrow(HLA)) {
  for (j in 1:nrow(H2_2)) {
    if (HLA$Peptide[i] == H2_2$Peptide[j]) {
      HLA$Mices_H2_2[i] <- 'H2-IAd'
    }
  }
}

#Filtrar péptidos que no son reconocidos por MHC de ratón
remove <- c()

for (i in 1:nrow(HLA)) {
  if (HLA$Mices_H2_1[i] == '-' & HLA$Mices_H2_2[i] == '-') {
    remove <- c(remove, i)
  }
}

HLA <- HLA[-remove, ]

write.csv(HLA, 'Output/Epitopes_Global_micefilt.csv', row.names = FALSE)
