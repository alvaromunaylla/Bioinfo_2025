library(dplyr)
library(readr)
library(reshape)
library(seqinr)
library(tidyr)
setwd('c:/Users/Alvaro/Documents/OROV_project/7_Multiepitopic_Protein_link_adj_sign/3_Multiepitopic_Protein_PADRE/')
Parse_NetMHCpan <- function(prediction_file, fasta_file, score_type) {
  
  #Crear dataframe
  data <- data.frame(read_delim(prediction_file, delim = '\t', 
                                escape_double = FALSE, show_col_types = FALSE))
  
  #Indicar longitud de secuencia
  fasta <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(fasta, as.string = TRUE))-8
  data$length <- rep(lengths, nrow(data))
  
  #Obtener los alelos HLA
  headers <- colnames(data)
  alleles_hla <- headers[grep("HLA", headers)]

  #Indicar de qué columnas se va a extraer el score de cada epítopo para cada alelo
  if (score_type == "Rank") {
    col_extracted <- seq(from = 7, by = 4, length.out = length(alleles_hla))
  } else {
    col_extracted <- seq(from = 6, by = 4, length.out = length(alleles_hla))
  }
  
  #Generar el dataframe output
  NetMHC_parsed <- data[-1, c(2, 1, ncol(data), 3, col_extracted)]
  colnames(NetMHC_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_hla)
  NetMHC_parsed$Position <- as.numeric(NetMHC_parsed$Position)+1
  NetMHC_parsed[c(5:ncol(NetMHC_parsed))] <- 
    sapply(NetMHC_parsed[c(5:ncol(NetMHC_parsed))], as.numeric)
  
  return(NetMHC_parsed)
}

Neoepi_scan <- Parse_NetMHCpan('1_Neoepitopes_test/NetMHCpan_results.xls', 
                               'Multiepitope_protein_PADRE.fasta', 'Rank')

#Filtrar los péptidos que no tienen afinidad para ningún alelo
for (i in 1:nrow(Neoepi_scan)) {
  for (j in 5:ncol(Neoepi_scan)) {
    if (Neoepi_scan[i,j] > 2) {
      Neoepi_scan[i,j] <- NA
    }
  }
}

Neoepi_scan <- Neoepi_scan[!rowSums(is.na(Neoepi_scan[, 5:ncol(Neoepi_scan)])) 
                           == ncol(Neoepi_scan[, 5:ncol(Neoepi_scan)]), ]

#Verificar si los péptidos con afinidad se formaron al unir los epítopos
Selected_epitopes <- read_table('Selected_epitopes.txt', col_names = TRUE)
Selected_epitopes <- Selected_epitopes$Selectes_epitopes

Original_epitope <- c()

for (i in 1:nrow(Neoepi_scan)) {
    if (any(grepl(Neoepi_scan$Peptide[i], Selected_epitopes))) {
      Original_epitope <- c(Original_epitope, i)
    }
}

Neoepi_scan <- Neoepi_scan[-Original_epitope, ]

#Verificar cuántas de las uniones son weak y strong
Neoepi_scan$Strong_bind <- rowSums(Neoepi_scan[, 5:(ncol(Neoepi_scan)-2)] <= 0.5, 
                                   na.rm = TRUE)
Neoepi_scan$Weak_bind <- rowSums(Neoepi_scan[, 5:(ncol(Neoepi_scan)-2)] > 0.5, 
                                 na.rm = TRUE)
Neoepi_scan <- Neoepi_scan[, c(1:4, ncol(Neoepi_scan) -1, ncol(Neoepi_scan))]
Neoepi_scan <- Neoepi_scan %>% arrange(desc(Strong_bind))

write.csv(Neoepi_scan, 'Multiepi_Prot_PADRE_Neoepi_1_scan.csv', row.names = FALSE)

#Detección de neoepítopos de células T CD4
NetMHCII_parse <- function(data_file, fasta_file, score_type) {
  
  #Generar dataframe
  data <- data.frame(read_delim(data_file, delim = '\t', show_col_types = FALSE, 
                                escape_double = FALSE))
  
  #Obtener longitud y nombre de secuencias
  seqs <- seqinr::read.fasta(fasta_file)
  lenghts <- nchar(getSequence(seqs, as.string = TRUE))-8
  tabla <- data.frame(ID = getName(seqs))
  tabla$length <- lenghts
  data$length <- tabla$length[match(data[, 3], tabla$ID)]
  
  #Obtener nombre de los alelos HLA
  headers <- colnames(data)
  HLA_alleles <- headers[grep("DRB", headers)]
  
  #Indicar qué columnas se extraen según el score_type
  if (score_type == "Rank") {
    cols_extract <- seq(from = 7, by = 3, length.out = length(HLA_alleles))
  } else {
    cols_extract <- seq(from = 6, by = 3, length.out = length(HLA_alleles))
  }
  
  #Crear el archivo de salida
  NetMHCII_parsed <- data[-1, c(2, 1, ncol(data), 3, cols_extract)]
  colnames(NetMHCII_parsed) <- c("Peptide", "Position", "Length", "ID", 
                                 HLA_alleles)
  NetMHCII_parsed$Position <- as.numeric(NetMHCII_parsed$Position)
  NetMHCII_parsed[c(5:ncol(NetMHCII_parsed))] <- 
    sapply(NetMHCII_parsed[c(5:ncol(NetMHCII_parsed))], as.numeric)
  NetMHCII_parsed$Position <- as.numeric(NetMHCII_parsed$Position)
  
  return(NetMHCII_parsed)
}

Neoepi2_scan <- NetMHCII_parse('1_Neoepitopes_test/NetMHCIIpan_results.xls', 
                                 'Multiepitope_protein_PADRE.fasta', 'Rank')

for (i in 1:nrow(Neoepi2_scan)) {
  for (j in 5:ncol(Neoepi2_scan)) {
    if (Neoepi2_scan[i,j] > 5) {
      Neoepi2_scan[i,j] <- NA
    }
  }
}

Neoepi2_scan <- Neoepi2_scan[!rowSums(is.na(Neoepi2_scan[, 5:ncol(Neoepi2_scan)])) 
                             == ncol(Neoepi2_scan[, 5:ncol(Neoepi2_scan)]), ]

Selected_epitopes <- read_table('Selected_epitopes.txt', col_names = TRUE)
Selected_epitopes <- Selected_epitopes$Selectes_epitopes

Original_epitope <- c()

for (i in 1:nrow(Neoepi2_scan)) {
  if (any(grepl(Neoepi2_scan$Peptide[i], Selected_epitopes))) {
    Original_epitope <- c(Original_epitope, i)
  }
}

Neoepi2_scan <- Neoepi2_scan[-Original_epitope, ]

#Verificar cuántas de las uniones son weak y strong
Neoepi2_scan$Strong_bind <- rowSums(Neoepi2_scan[, 5:(ncol(Neoepi2_scan)-2)] <= 1, 
                                    na.rm = TRUE)
Neoepi2_scan$Weak_bind <- rowSums(Neoepi2_scan[, 5:(ncol(Neoepi2_scan)-2)] > 1, 
                                  na.rm = TRUE)
Neoepi2_scan <- Neoepi2_scan[, c(1:4, ncol(Neoepi2_scan) -1, ncol(Neoepi2_scan))]
Neoepi2_scan <- Neoepi2_scan %>% arrange(desc(Strong_bind))

write.csv(Neoepi2_scan, 'Multiepi_Prot_PADRE_Neoepi_2_scan.csv', row.names = FALSE)
