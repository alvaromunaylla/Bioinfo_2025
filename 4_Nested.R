library(seqinr)
library(readr)
library(readxl)
setwd('c:/Users/Alvaro/Documents/OROV_project/4_Nested/')

#Generar Dataframes con las tablas del HLA1 y  2(colocar lista de alelos)
HLA1 <- data.frame(read_delim('Input/1_HLA1/Epitopes_HLA1_Global_conserved.csv', 
                              delim = ',', show_col_types = FALSE))

HLA2 <- data.frame(read_delim('Input/2_HLA2/Epitopes_HLA2_Global_conserved.csv',
                              delim = ',', show_col_types = FALSE))

#Añadir a HLA2 columna con número y string de alelos HLA2 de cada epítopo
Alelos_HLAII <- rep(0, nrow(HLA2))
HLA2$Alelos_HLAII <- Alelos_HLAII
HLA2$Num_Alelos_HLAII <- HLA2$Promiscuity

for (i in (1:nrow(HLA2))) {
  AlelosII <- c()
  for (j in (8:(ncol(HLA2)-3))) {
    if (!is.na(HLA2[i,j])) {
      AlelosII <- c(AlelosII, colnames(HLA2)[j])
      HLA2$Alelos_HLAII[i] <- toString(AlelosII)
    }
  }
}

#Añadir la información sobre los epítopos de HLA1 anidados
Epitopes_HLAI <- rep(0, nrow(HLA2))
HLA2$Epitopes_HLAI <- Epitopes_HLAI
HLA2$Num_Epitopes_HLA1 <- Epitopes_HLAI
HLA2$Alelles_HLAI <- Epitopes_HLAI
HLA2$Num_Alelos_HLAI <- Epitopes_HLAI
HLA1$Peptide <- as.character(HLA1$Peptide)
HLA2$Peptide <- as.character(HLA2$Peptide)

for (i in (1:nrow(HLA2))) {
  EpitopesI <- c()
  Num_EpitopesI <- c()
  AlelosI <- c()
  Num_AlelosI <- c()
  for (j in (1:nrow(HLA1))) {
    if (HLA2$ID[i] == HLA1$ID[j]) {
      if (length(grep(HLA1$Peptide[j], HLA2$Peptide[i]) > 0)) {
        #Completar vector con epítopos
        EpitopesI <- c(EpitopesI, HLA1$Peptide[j])
        #Completar vector con alelos 
        for (k in (8:(ncol(HLA1)-1))) {
          if (!is.na(HLA1[j,k]) == TRUE) {
            AlelosI <- c(AlelosI, colnames(HLA1[k]))
          }
        }
      }
    }
    #Añadir y contar Epítopos de HLA1
    HLA2$Epitopes_HLAI[i] <- toString(EpitopesI)
    HLA2$Num_Epitopes_HLA1[i] <- length(EpitopesI)
    #Añadir y contar alelos que tienen afinidad con los epítopos 
    AlelosI <- unique(AlelosI)
    HLA2$Alelles_HLAI[i] <- toString(AlelosI)
    HLA2$Num_Alelos_HLAI[i] <- length(AlelosI)
  }
}

#Filtrar epítopos HLA2 que no contienen epítopos anidados
No_nested <- c()

for (i in (1:nrow(HLA2))) {
  if (HLA2$Num_Alelos_HLAI[i] == 0) {
    No_nested <- c(No_nested, i)
  }
}

HLA2 <- HLA2[-No_nested, ]

#Generar el archivo output
Output <- HLA2[, c(7, 1, 4:5, 3, ncol(HLA2)-3, ncol(HLA2)-2, ncol(HLA2)-1, 
                   ncol(HLA2)-5, ncol(HLA2), ncol(HLA2)-4)]

#Modificar los nombres de HLA1 a la nomenclatura correcta
replace_name_HLAI <- function(tochange) {
  Actual_names <- unlist(strsplit(tochange, ", "))
  
  replace_one <- function(nom) {
    parts <- strsplit(nom, '')[[1]]
    New_name <- paste(c(parts[1:3], '-', parts[5], '*', parts[6:7], ':', 
                        parts[9:10]), collapse = '')
    return(New_name)
  }
  
  Changed_names <- sapply(Actual_names, replace_one)
  HLAI_correct <- paste(Changed_names, collapse = ",")
  return(HLAI_correct)
}

Output$Alelles_HLAI <- sapply(Output$Alelles_HLAI, replace_name_HLAI)

#Modificar los nombres de HLA2 a la nomenclatura correcta
replace_name_HLAII <- function(tochange) {
  Actual_names <- unlist(strsplit(tochange, ", "))
  
  replace_one <- function(nom) {
    parts <- strsplit(nom, '')[[1]]
    New_name <- paste(c('HLA-', parts[1:4], '*', parts[6:7], ':', 
                        parts[8:9]), collapse = '')
    return(New_name)
  }
  
  Changed_names <- sapply(Actual_names, replace_one)
  HLAII_correct <- paste(Changed_names, collapse = ",")
  return(HLAII_correct)
}

Output$Alelos_HLAII <- sapply(Output$Alelos_HLAII, replace_name_HLAII)

#Generar el archivo con el output
write.csv(Output, 'Output/Epitopes_Global_Nested.csv', row.names = FALSE)

      