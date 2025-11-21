library(seqinr)
library(readr)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
setwd('c:/Users/Alvaro/Documents/OROV_project/3_Epitopes_T/6_Promiscuity/')

Data <- data.frame(read_xlsx('Output_Global/Epitopes_T_Global.xlsx'))

#Establecer el número total de alelos HLA1 y 2
Alelos1 <- unique(strsplit(paste(Data$Alelles_HLAI, collapse = ','), 
                          split = ',')[[1]])
Num_Alelos1 <- length(Alelos1)
Alelos2 <- unique(strsplit(paste(Data$Alelos_HLAII, collapse = ','),
                           split = ',')[[1]])
Num_Alelos2 <- length(Alelos2)

#Agregar la promiscuidad "combinada" de cada epítopo
for (i in 1:nrow(Data)) {
  Data$Promiscuity[i] <- sum(Data$Num_Alelos_HLAI[i], 
                             Data$Num_Alelos_HLAII[i])
}

#Gráfica de distribución de promiscuidad
plot <- ggplot(Data, aes(x = Promiscuity)) +
  geom_bar(color = 'blue4', alpha = 0.7) +
  xlim(0, sum(Num_Alelos1, Num_Alelos2)) +
  ggtitle('Distribución de la promiscuidad de epítopos') +
  xlab('Promiscuidad HLA1/2') +
  ylab('Número de epítopos') +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 15))
pdf('Global_promiscuity_distribution.pdf')
print(plot)
dev.off()

quantile(Data$Promiscuity)
IQR(Data$Promiscuity)

#Filtrar epítopos de baja promiscuidad
Data <- subset(Data, Promiscuity > 25)

#Tabla para cobertura
Data <- Data %>% arrange(desc(Promiscuity))

Coverage <- setNames(data.frame(matrix(ncol = 1 + Num_Alelos1 + Num_Alelos2, 
                                       nrow = 0)), c('Peptide', Alelos1, Alelos2))

for (i in (1:nrow(Data))) {
  Alelos1_2 <- c(strsplit(Data$Alelles_HLAI[i], split = ',')[[1]], 
                 strsplit(Data$Alelos_HLAII[i], split = ',')[[1]])
  New_row <- rep('-', Num_Alelos1 + Num_Alelos2)
  for (j in 1:length(Alelos1_2)) {
    for (k in 1:length(New_row)) {
      if (Alelos1_2[j] == names(Coverage)[[k+1]]) {
        New_row[k] <- 1
      }
    }
  }
  New_row <- c(Data$Peptide[i], New_row)
  Coverage[nrow(Coverage) +1,] <- New_row
}

write.csv(Data, 'Output_Global/Epitopes_Global_promiscuity_list.csv', row.names = FALSE)

#Hacer el gráfico hotspot
Data_test <- pivot_longer(Coverage, cols = starts_with("HLA"), 
                          names_to = "Alelo", values_to = "Afinidad")
Data_test$Peptide <- factor(Data_test$Peptide, levels = unique(Data_test$Peptide))

plot <- ggplot(Data_test, aes(x = Alelo, y = Peptide, fill = factor(Afinidad))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("blue4", "skyblue"), name = "Afinidad") +
  theme_minimal() +
  labs(title = "Afinidad entre péptidos y Alelos HLA",
       x = "Alelo HLA", y = "Péptido") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 5))

pdf("Global_Test_coverage.pdf", width = 12, height = 45)
print(plot)
dev.off()

#Procedimiento para facilitar la selección de conjunto de péptidos

#Para crear un vector a partir de cada fila y ponerlos en una lista
List_Affi <- list()
Names_list <- c()
for (i in 1:nrow(Coverage)) {
  Row2vector <- unname(as.vector(as.character(Coverage[i,])))
  List_Affi <- c(List_Affi, list(Row2vector))
  Names_list <- c(Names_list, Row2vector[1])
}
names(List_Affi) <- Names_list

#Función para identificar los péptidos para que combinados con un péptido dado permiten mayor covertura HLA
Identify_pept <- function(Peptido_prob) {
  #Función para emparejar un péptido con otro
  Select_pept <- function(Peptide_cov) {
    #Crear la tabla donde se colocarán los otros pépt y diferencia de covertura
    Compare <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c('Peptide', 'Contador'))
    
    for (vector in List_Affi) {
      contador <- c(0)
      for (i in 2:length(Peptide_cov)) {
        if (Peptide_cov[i] == '-' && Peptide_cov[i] != vector[i]) {
          contador <- c(contador + 1)
        }
      }
      new_row <- c(vector[1], contador)
      Compare[nrow(Compare) + 1,] <- new_row
    }
    
    #Seleccionar el pept con el que se lograría una mejor covertura
    Compare <- Compare %>% arrange(desc(Contador))
    Best_combinator <- List_Affi[[Compare$Peptide[1]]]
    
    #Hacer la "combinación" entre el vector anallizado y el vector seleccionado
    Combination <- Peptide_cov
    for (i in 1:length(Peptide_cov)) {
      if (Peptide_cov[i] == '-' && Best_combinator[i] == '1') {
        Combination[i] <- Best_combinator[i]
      }
    }
    Combination_test <- Combination
    return(list(Best_combinator[1], Combination_test))
  }
  
  VectorToSearch <- List_Affi[Peptido_prob][[1]]
  
  First_merging <- Select_pept(VectorToSearch)
  Second_merging <- Select_pept(First_merging[[2]])
  Third_merging <- Select_pept(Second_merging[[2]])
  Fourth_merging <- Select_pept(Third_merging[[2]])
  Fifth_merging <- Select_pept(Fourth_merging[[2]])
  Sixth_merging <- Select_pept(Fifth_merging[[2]])
  
  Pept_combinators <- c(First_merging[[1]], Second_merging[[1]], 
                        Third_merging[[1]], Fourth_merging[[1]],
                        Fifth_merging[[1]], Sixth_merging[[1]])
                        
  return(list(Pept_combinators, Sixth_merging[[2]]))
}

#Crear una tabla donde se pueda ver la combinación de péptidos "idónea" (?) para cada péptido
Pept_selection <- data.frame(Peptide1 = Coverage$Peptide)
Pept_selection[, c('Peptide2', 'Peptide3', 'Peptide4', 'Peptide5', 'Peptide6', 
                   'Peptide7', 'Final_coverage')] <- NA

for (i in 1:nrow(Pept_selection)) {
  #Identificar péptidos
  Complem_pepts <- Identify_pept(Pept_selection$Peptide1[i])[[1]]
  Pept_selection[i, c(2:7)] <- Complem_pepts
  
  #Remover la duplicación de primer péptido (si se forma)
  if (Pept_selection[i,7] == Pept_selection[i,1]) {
    Pept_selection[i,7] <- '-'
  }
  
  #Identificar si se logra una covertura total con 5 péptidos
  Final_Coverage <- Identify_pept(Pept_selection$Peptide1[i])[[2]]
  if (all(Final_Coverage[-1] == '1')) {
    Pept_selection$Final_coverage[i] <- 'Total'
  }
}

write.csv(Pept_selection, 'Epitopes_Global_posible_combination.csv', 
          row.names = FALSE)

write.csv(Data, 'Epitopes_Global_promiscuity_list.csv', 
          row.names = FALSE)

#PROBAR HACER UN GRÁFICO CON LOS EPÍTOPOS QUE ELEGISTE
SudAm_op1 <- Data[c(3,2,1,24,14,28,55,13,36,5,6,33,34,12,11,76,44,30,125),]

Coverage <- setNames(data.frame(matrix(ncol = 1 + Num_Alelos1 + Num_Alelos2, 
                                       nrow = 0)), c('Peptide', Alelos1, Alelos2))

for (i in (1:nrow(SudAm_op1))) {
  Alelos1_2 <- c(strsplit(SudAm_op1$Alelles_HLAI[i], split = ',')[[1]], 
                 strsplit(SudAm_op1$Alelos_HLAII[i], split = ',')[[1]])
  New_row <- rep('-', Num_Alelos1 + Num_Alelos2)
  for (j in 1:length(Alelos1_2)) {
    for (k in 1:length(New_row)) {
      if (Alelos1_2[j] == names(Coverage)[[k+1]]) {
        New_row[k] <- 1
      }
    }
  }
  New_row <- c(SudAm_op1$Peptide[i], New_row)
  Coverage[nrow(Coverage) +1,] <- New_row
}

Data_test_op1 <- pivot_longer(Coverage, cols = starts_with("HLA"), 
                          names_to = "Alelo", values_to = "Afinidad")
Data_test_op1$Peptide <- factor(Data_test_op1$Peptide, levels = unique(Data_test_op1$Peptide))

plot <- ggplot(Data_test_op1, aes(x = Alelo, y = Peptide, fill = factor(Afinidad))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("skyblue", "darkblue"), name = "Afinidad") +
  theme_minimal() +
  labs(title = "Afinidad entre péptidos y Alelos HLA",
       x = "Alelo HLA", y = "Péptido") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5))  

pdf("Output_SudAm/Coverage_Op1.pdf", width = 12, height = 4)
print(plot)
dev.off()




