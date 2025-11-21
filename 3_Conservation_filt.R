library(seqinr)
library(readr)
library(readxl)
setwd('c:/Users/Alvaro/Documents/OROV_project/3_Conservation/')

#Hacer un dataframe con las posiciones variables (ID - Position)
data_excel <- read_excel('Proteínas-variaciones.xlsx')
data_excel$Posición <- gsub(",", ", ", as.character(data_excel$Posición))
write.table(data_excel, 'Proteínas-variaciones.txt', sep = "\t", quote = FALSE, 
            row.names = FALSE)

Mutations <- data.frame(read_delim('Proteínas-variaciones.txt', delim = '\t', 
                                   show_col_types = FALSE))
colnames(Mutations) <- c('ID', 'Position', 'Length', 'Num', 'Mutation')

Mutattions_df <- data.frame(ID = character(), Position = character(), 
                            stringsAsFactors = FALSE)

for (i in (1:nrow(Mutations))) {
  Positions <- unlist(strsplit(Mutations$Position[i], ','))
  Positions <- trimws(Positions)
  Name_vector <- Mutations$ID[i]
  Mutattions_df <- rbind(Mutattions_df, data.frame(ID = rep(Name_vector, 
                                                            length(Positions)),
                                                   Position = Positions))
}

#Añadir posición de inicio y fin a los peptidos
Merged_results <- data.frame(read_delim('Input/2_HLA2/Merged_Global_HLA2.txt', 
                                        delim = '\t', show_col_types = FALSE))
Pept_length <- rep(0, nrow(Merged_results))
Merged_results$Pept_length <- Pept_length
Merged_results$Start <- Pept_length
Merged_results$End <- Pept_length

for (i in (1:nrow(Merged_results))) {
  Merged_results$Pept_length[i] <- nchar(Merged_results$Peptide[i])
  Merged_results$Start[i] <- Merged_results$Position[i]
  Merged_results$End[i] <- (Merged_results$Position[i] + 
                              Merged_results$Pept_length[i] - 1)
}

#Para verificar si el Pept incluye una posición mutada
#Crear un vector con las filas que incluyan pept con posiciones no cons.
Merged_results$ID <- as.character(Merged_results$ID)
Mutattions_df$ID <- as.character(Mutattions_df$ID)
Mutattions_df$Position <- as.numeric(Mutattions_df$Position)

Remove_rows <- c()

for (i in (1:nrow(Merged_results))) {
  for (j in (1:nrow(Mutattions_df))) {
    if (Merged_results$ID[i] == Mutattions_df$ID[j]) {
      if (Merged_results$Start[i] <= Mutattions_df$Position[j] & 
          Merged_results$End[i] >= Mutattions_df$Position[j]) {
        Remove_rows <- c(Remove_rows, i)
      }
    }
  }
}

Remove_rows <- unique(Remove_rows)

#Remover las filas según lo que indique el vector creado
Merged_results_conserved <- Merged_results[-Remove_rows, ]

#Ordenar las columnas
Merged_results_conserved <- Merged_results_conserved[, c(1:2, 
                                                         (ncol(Merged_results_conserved)-2):ncol(Merged_results_conserved), 
                                                         3:(ncol(Merged_results_conserved)-3))]

#Añadir columna con promiscuidad de los epítopos
Promiscuity <- rep(0, nrow(Merged_results_conserved))
Merged_results_conserved$Promiscuity <- Promiscuity
Merged_results_conserved$Promiscuity <- rowSums(!is.na(Merged_results_conserved[, 8:(ncol(Merged_results_conserved)-1)]))


write.csv(Merged_results_conserved, 
          'Output/2_HLA2/Epitopes_HLA2_Global_conserved.csv', row.names = FALSE)
