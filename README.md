# Bioinfo_2025

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

1.	Obtención de un proteoma consenso de OROV.

Se realizará una búsqueda del genoma de todos los linajes de OROV identificados hasta la fecha en Brasil, Perú, Colombia y Ecuador, los países más afectados por la fiebre de Oropouche según reportes (9). Tomando como referencia los genes identificados en la cepa BeAn19991 (81), se extraerá la secuencia codificante de las seis proteínas de todos los linajes de OROV, y estas serán transformadas a secuencias aminoacídicas. Utilizando el programa Jalview, las secuencias serán alineadas para finalmente obtener un consenso de las seis proteínas.

2.	Búsqueda de los alelos HLA más comunes en los países más afectados por OROV y a nivel global.

Para la predicción de epítopos T se trabajará con dos grupos de alelos HLA: uno correspondiente a los alelos más comunes en Brasil, Perú, Colombia y Ecuador, y otro con los alelos más comunes a nivel mundial. Estos datos serán obtenidos mediante una revisión bibliográfica y en la base de datos “Allele Frequency Net Database”.

3.	Predicción de epítopos para células T.

Con las secuencias consenso de las proteínas se realizará la predicción de epítopos de clase I mediante los programas NetMHCpan v4.1  y MHCFlurry v2.0.4 (66), y la de epítopos de clase II a través de los programas NetMHCIIpan v4.0 (65) y MixMHC2Pred v2.0 (67). Como se menciona, se utilizarán dos programas por tipo de epítopo T, ya que ambos resultados se combinarán para tener una mayor fiabilidad en las predicciones. Los péptidos con afinidad fuerte a alelos HLA serán considerados candidatos a epítopos de células T.

4.	Predicción de epítopos para células B.

Se utilizarán los programas BepiPred v3.0 (70) y DiscoTope v3.0 (82) para la predicción de epítopos en las secuencias consenso de las proteínas. Se combinará y comparará los resultados de ambos programas para tener mayor fiabilidad en la predicción, y posteriormente se seleccionará como candidatos los péptidos con mayor probabilidad de ser epítopos de células B.

5.	Selección de epítopos.

Se aplicará una serie de filtros sobre los epítopos predichos para identificar los candidatos más óptimos. Los aspectos que se considerará en los filtros son: la conservación del péptido entre los distintos linajes del virus, para lo cual se utilizará como referencia los sitios de variación identificados al realizar el alineamiento de las proteínas de las variantes de OROV; la ausencia del péptido en el proteoma humano y el de bacterias de la microbiota intestinal,
Para seleccionar los péptidos ubicados en regiones conservadas, se utilizará un script en R que descarte los péptidos que incluyan algún residuo para el cual se haya identificado una tasa de variación significativa.
Algunos de los aspectos que se considerará en los filtros son: la conservación del péptido entre los distintos linajes del virus, ausencia del péptido en el proteoma humano o el de las bacterias de la microbiota intestinal (con el fin de evitar el riesgo de autoinmunidad), y bajo potencial de toxicidad y alergenicidad del péptido.

6.	Generación de candidatos multiepitópicos.

Los epítopos seleccionados como candidatos serán conectados entre sí mediante linkers. Se emplearán algoritmos y scripts de Python para descartar las construcciones que contengan neoepítopos en las uniones. De esta forma, se espera obtener varios candidatos multiepitópicos.

7.	Evaluaciones mediante biofísica computacional.

A través de herramientas de modelado molecular por AlphaFold2 se obtendrá la estructura de los candidatos vacunales multiepitópicos, y posteriormente mediante simulaciones de dinámicas moleculares con el programa GROMACS se determinará su estabilidad. Empleando el programa ClusPro, se realizarán ensayos de acoplamiento molecular contra receptores inmunes TLR para determinar la afinidad por estos y su potencial para inducir respuesta inmune. De igual manera la estabilidad de la interacción se determinará por simulaciones de dinámica molecular. Mediante estos análisis se espera seleccionar el mejor candidato multiepitópico, el cual será producido para la evaluación experimental.
