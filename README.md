# I. Planteamiento del problema
La fiebre de Oropouche es una enfermedad causada por el virus Oropouche (Oropouche Virus u OROV) que prevalece en regiones amazónicas. Desde su primera detección, se han reportado más de 30 epidemias en ciertos países de Sudamérica, principalmente Brasil, donde, además, es recurrente la aparición de nuevos brotes del virus (Romero-Alvarez et al., 2023). En 2021, el Instituto Nacional de Salud (INS) emitió el último informe del estado de esta patología en Perú, notificando un total 105 casos diagnosticados ese año (Instituto Nacional de Salud, 2021), y destacando la persistencia de la enfermedad en la región amazónica en la última década (Alvarez-Falconi & Ríos, 2010; Castro et al., 2013; García et al., 2016), con casos emergentes incluso en regiones costeras (Martins-Luna et al., 2020). En conjunto, estos datos reflejan el potencial epidémico de este virus y resaltan la importancia de un programa de vigilancia.

Si bien la fiebre de Oropouche no es una enfermedad mortal, en casos particulares como en niños, adultos o personas inmunocomprometidas, la gravedad de los síntomas puede ser mayor, afectando el sistema nervioso e incluso poniendo en riesgo la vida del paciente. En todos los casos, los síntomas son tratados con medicamentos paliativos, pues no existe un tratamiento específico para esta enfermedad (Sakkas et al., 2018). Esto convierte a la prevención en un aspecto fundamental para evitar cualquier tipo de complicación. Aunque OROV ha ido adquiriendo mayor relevancia y han incrementado las investigaciones en torno a esta, fundamentalmente sobre sus aspectos moleculares, epidemiológicos y médicos, hasta la fecha no se ha aprobado una vacuna contra este virus lo cual dificulta la posibilidad de que la enfermedad sea completamente controlada y deja expuesto a los sectores de población vulnerables y a las regiones donde existe mayor propensión al surgimiento de un nuevo brote.

Considerando la recurrencia y la ventana de probabilidad para el surgimiento de nuevos brotes de OROV en el Amazonas, se resalta la necesidad imperante de desarrollar una vacuna para un mejor control epidemiológico. Dada la confianza actual en diseños basados en vacunología reversa y la eficacia probada de los virus adeno-asociados como vector en vacunas, surge la siguiente pregunta: ¿Un enfoque inmunoinformático permitiría la identificación de epítopos de OROV altamente inmunogénicos que sirvan como base para el diseño de una vacuna multiepitópica?



# II. Objetivos
1. Objetivo general

    Desarrollar un candidato a vacuna contra virus Oropouche mediante herramientas in silico.

2. Objetivos específicos

    -	Predecir epítopos T de clase I y de clase II de virus Oropouche restringidos a los alelos HLA predominantes en la población sudamericana y a nivel mundial.
    -	Predecir epítopos B lineales y conformacionales de virus Oropouche.
    -	Construir candidatos vacunales multiepitópicos con base en los epítopos predichos seleccionados.
    -	Seleccionar el mejor candidato multiepitópico mediante análisis de biofísica computacional.


# III. Metodología

1.	Obtención de secuencias del OROV

Las secuencias nucleotídicas de los segmentos S, M y L de OROV fueron obtenidas de la base de datos BV-BRC. Se seleccionaron únicamente las secuencias completas.

2.	Obtención de un genoma y proteoma de referencia

Tomando como referencia los genes anotados de la cepa BeAn19991 en GeneBank, se realizó un alineamiento para mapear las secuencias codificantes de las seis proteínas de todos los linajes de OROV. Las secuencias nucleotídicas fueron transformadas a secuencias de aminoácidos en el programa Jalview para finalmente generar un consenso de las seis proteínas de OROV. Además, se ubicaron los residuos sujetos a variaciones entre las variantes y su tasa de mutación. Se consideró como sitios no conservados o variables aquellos cuyo porcentaje de variación entre cepas fue mayor al 5%. Esta información fue añadida a un documento con las mutaciones y sus posiciones dentro de las proteínas del virus. 

3.	Búsqueda de alelos HLA

La predicción de epítopos T se realizó con los alelos HLA más comunes en países amazónicos afectados por OROV. Además, de manera complementaria se realizó una segunda predicción empleando los alelos más frecuentes a nivel mundial. Estos grupos de alelos HLA se denominaron “SudAm” y “Global”, respectivamente. Para el conjunto de datos “SudAm” se consideró a Brasil, Perú, Ecuador, Colombia y Bolivia como los países con mayor incidencia o riesgo de surgimiento de brotes de OROV.

4.	Predicción de epítopos T CD8 (epítopos CTL)

Se predijo la afinidad de unión de los péptidos de las proteínas consenso de OROV al HLA-I utilizando los programas NetMHCpan v4.1 y MHCFlurry v2.0. Los péptidos evaluados tuvieron una longitud de 9 a 11 aminoácidos. Se realizaron las predicciones tanto para el conjunto de alelos HLA-I “SudAm”, como para el conjunto “Global”. En los resultados generados por NetMHCpan v4.1, se consideró como péptidos con una afinidad elevada (strong binders) y media (weak binders) a un alelo HLA aquellos cuyo valor Rank0 era de 0.5% y 2%, respectivamente. Por otro lado, en MHCFlurry v2.0, se consideró que un péptido presentaba afinidad por un alelo HLA si la puntuación de afinidad era menor a 500 nM.

5.	Predicción de epítopos T CD4 (epítopos HTL)

Se predijo la afinidad de unión de los péptidos de las proteínas consenso de OROV al HLA-I utilizando los programas NetMHCIIpan v4.0 y MixMHC2Pred. Se estableció una longitud de péptidos de 15 aminoácidos. Las predicciones fueron realizadas empleando el conjunto de alelos HLA-II “SudAm” y “Global”. Se consideró como péptidos con una afinidad elevada (strong binders) y media (weak binders) a un alelo HLA aquellos cuyo valor Rank era de 2% y 5%, respectivamente. De forma similar, en los resultados de MixMHC2Pred se consideró que un péptido tenía afinidad hacia un alelo HLA si su valor %Rank era menor a 5.

6.	Selección de epítopos T

6.1.	Filtro de consenso entre predictores

Los resultados obtenidos de los programas para cada predicción fueron combinados. Para ello se elaboraron códigos en R: “HLA1_merging_script.R”, que combina los resultados de NetMHCpan y MHCFlurry, y “HLA2_merging_script.R” para combinar los resultados de NetMHCIIpan y MixMHC2Pred. De esta forma, se seleccionaron como candidatos a epítopos únicamente a los péptidos que, en consenso de los dos programas empleados, tenían afinidad por al menos un alelo HLA.

6.2.	Filtro de conservación del epítopo

A continuación, se filtraron los péptidos ubicados en regiones conservadas. La línea de códigos “Conservation_filt_script.R" elaborada en R se utilizó para identificar si un péptido contenía un residuo sujeto a mutaciones, usando como referencia el archivo con los residuos variables entre cepas de OROV.

6.3.	Generación de epítopos anidados

Posteriormente, se generaron epítopos anidados, es decir, epítopos CD4 que contenían epítopos CD8. Esto se realizó con el objetivo de maximizar el número de epítopos que se incluirían en la proteína multiepitópica diseñada. Para esto se generó la línea de códigos “Nested.R”, la cual da como resultado una tabla con los epítopos CD4, los epítopos CD8 que incluye y los alelos HLA de tipo 2 y 1 para los cuales estos epítopos presentan afinidad. 

6.4.	Filtro de epítopos con homología en humanos

Los epítopos fueron ingresados a BLAST para realizar una búsqueda de similitudes con el proteoma humano y el proteoma de las bacterias más abundantes de la microbiota intestinal. Se descartaron los epítopos en los que se encontraron similitudes. De esta forma se evitaría que la proteína multiepitópica pueda generar autoinmunidad o una respuesta inmune desviada hacia bacterias de la microbiota.

6.5.	Selección de epítopos según su promiscuidad

Finalmente, los epítopos que pasaron todos los filtros, fueron clasificados según su promiscuidad, es decir, número de alelos HLA con los que presenta afinidad. Luego se realizó la selección de los epítopos, de manera que se logre una cobertura total de los alelos HLA utilizados en el análisis, tanto los del conjunto de datos “SudAm”, como el “Global”. Además, se incluyó al menos un epítopo de todas las proteínas de OROV para que la vacuna diseñada tenga potencial de generar inmunidad hacia todas estas proteínas.

7.	Predicción y selección de epítopos B

La predicción de epítopos de células B se realizó en base a la secuencia consenso de la glicoproteína Gc de OROV. Se utilizaron los programas BepiPred v3.0, EpiDope (predictores por residuo), ABCPred, LBTope y SVMTrip (predictores por péptido). En los programas que realizan la predicción por péptido, se estableció una longitud de 16 aminoácidos. Posteriormente, se diseñó la línea de códigos “B_epitopes_align” para visualizar y mapear los epítopos B dentro de la proteína. Para considerar un péptido como candidato a epítopo, este tendría que haber sido predicho por Bepipred v3.0 y al menos dos de los otros cuatro predictores empleados. Se filtraron los candidatos con ausencia de glicosilaciones y ubicados en regiones conservadas de la proteína. Además, dado que en la disposición de la estructura tipo trípode conformada por las glicoproteínas Gn y Gc, solo el dominio “head” de Gc es accesible para anticuerpos, se eligieron únicamente los péptidos ubicados en dicho dominio. 

8.	Construcción de la proteína multiepitópica
Los epítopos seleccionados fueron concatenados a través de linkers para generar la proteína multiepitópica. Se realizaron distintas permutaciones de los epítopos y linkers, hasta obtener una en la cual no se hayan formado neoepítopos entre las uniones de dos epítopos. Para verificar esto, se predijeron epítopos en la proteína multiepitópica diseñada con los programas NetMHCpanI y NetMHCIIpan, y se evaluó la presencia de neoepítopos con la línea de códigos “Neoepitope_scanning”. Además, se añadieron dos secuencias en el extremo amino-terminal de la proteína: la primera correspondiente a un adyuvante, y la segunda, a un péptido-señal. Se generaron tres proteínas multiepitópicas diferenciadas por el adyuvante añadido, ya sea β-defensina 3, toxoide tétano-difteria (TpD) o el epítopo de unión pan-HLA DR (PADRE). Por otro lado, el péptido señal empleado fue la señal del activador tisular del plasminógeno (tPA).

9.	Modelización de la proteína multi-epítopo (Aún no terminado)

Las estructuras de las tres proteínas diseñadas fueron modeladas con el servidor Robetta, donde se realiza la predicción con el método ab initio mediante el algoritmo RoseTTAfold. Los modelos generados fueron visualizados con el UCSF Chimera v17.1.3. 
La calidad estructural de estos modelos fue validada empleando herramientas de evaluación estructural. En primer lugar, se utilizó ERRAT para analizar la confiabilidad de los modelos a partir del patrón de interacciones no enlazadas entre átomos. Posteriormente, se aplicó PROCHECK, que permite evaluar la estereoquímica de la estructura tridimensional mediante un diagrama de Ramachandran. Ambas herramientas fueron ejecutadas a través del servidor en línea SAVES v6.1. Entre las cinco estructuras modeladas por Robetta para cada candidato vacunal, se eligieron aquellas con los más óptimos valores de ERRAT y del gráfico Ramachandran.

10.	Evaluaciones mediante biofísica computacional (Aún no hecho)

A través de herramientas de modelado molecular por AlphaFold2 se obtendrá la estructura de los candidatos vacunales multiepitópicos, y posteriormente mediante simulaciones de dinámicas moleculares con el programa GROMACS se determinará su estabilidad. Empleando el programa ClusPro, se realizarán ensayos de acoplamiento molecular contra receptores inmunes TLR para determinar la afinidad por estos y su potencial para inducir respuesta inmune. De igual manera la estabilidad de la interacción se determinará por simulaciones de dinámica molecular. Mediante estos análisis se espera seleccionar el mejor candidato multiepitópico, el cual será producido para la evaluación experimental.
