# Unid3 sesion 4

## Tutorial 1

### Identificación de variantes en una muestra tumoral

Antecedente: 

    La caracterización de variantes genéticas en el tejido tumoral es crucial para la oncología de precisión, ya que permite comprender la dinámica del tumor y aplicar enfoques terapéuticos dirigidos.

    Este sistema bioinformático procesa datos de NGS utilizando módulos de alineamiento, identificación y anotación de variantes. El presente documento describe el protocolo de análisis utilizado para el cribado de variantes somáticas en una muestra de interés. Como ejemplo ilustrativo, se emplearán exclusivamente las lecturas correspondientes a una fracción menor del cromosoma 5

Posterior a descargar CLC Genomicos Worckench 25.0.3 y configurar los plugins procedemos a : 

#### Paso 1 descargar e importar la data

(https://resources.qiagenbioinformatics.com/testdata/Example_data_tumor_25.zip.)

![Insertar imagen x] (https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/1.jpg)

Los archivos target_regions_chr5 y tumor_reads_chr5 fueron importados correctamente al CLC Genomics Workbench mediante la herramienta Standard Import. Estos corresponden a las regiones objetivo del exoma y a las lecturas Illumina del tumor, respectivamente

#### Paso 2  Identificación de Variantes

para eso vamos a seguir el siguiente workflow → Template Workflows → Biomedical Workflows → Whole Exome Sequencing → Somatic Cancer (WES) → Identify Variants (WES)

![Insertar imagen x] (https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/2.jpg)

Parametrizamos segun las recomendaciones del tutorial. 

destacando:

Minimum coverage setting to the value 10.

Minimum frequency is set to 5.0%

Identification of Variants in a Tumor Sample Reference Data Set.

En la etapa de Low Frequency Variant Detection, se configuraron los parámetros recomendados en el tutorial: frecuencia mínima de 5% y cobertura mínima de 10 lecturas para asegurar un llamado más confiable y de referencia usaremos "Variants in a Tumor Sample Reference Data Set"

![Insertar imagen x] (https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/3.jpg)

#### Paso 3 comprobación y analisis del QC

##### RESULTADOS

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/4.jpg)

Tras ejecutar el flujo de trabajo se generaron los archivos principales: mapeo de lecturas, reporte de cobertura, variantes sin filtrar, variantes filtradas y visualización integrada en Genome Browser. Dentro de los resultados del analisis encontramos el informe con el nombre  Target_region_coverage_report-tumor_reads_chr5 

#### Tabla de resumen

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/5.jpg)

El promedio de cobertura de las regiones objetivo fue 22.5x, superando el mínimo de 10x recomendado para el llamado de variantes. El 82.6% de los nucleótidos de las regiones objetivo alcanzaron una cobertura ≥10, indicando una adecuada profundidad para este conjunto reducido de datos

##### Tabla: Fracciones con cobertura

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/6.jpg)

La mayoría de las regiones objetivo muestran una cobertura adecuada, con un 70.97% de ellas cubiertas en más de un 80% de su longitud con ≥10 lecturas. Esto es consistente con un enriquecimiento satisfactorio en este subconjunto del cromosoma 5.

##### Histograma de fracciones con cobertura ≥10

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/8.jpg)

El histograma de cobertura por fracción del target confirma que los extremos se distribuyen principalmente hacia coberturas altas, aunque algunas regiones muestran baja cobertura, como es esperable en datasets reducidos

##### Visualización de variantes en Genome Browser

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/11.jpg)

Esta figura muestra en una sola vista las distintas capas de información del cromosoma 5: los genes anotados, la cobertura de las regiones objetivo y las variantes detectadas por el análisis. Al observarlas juntas, pude comprobar si cada variante estaba realmente respaldada por lecturas y si aparecía en una zona con buena cobertura. Para mí, como estudiante que recién está aprendiendo este tipo de análisis y area de estudio, esta visualización fue  útil para entender cómo interpretar una variante dentro de su contexto y para reconocer cuándo un resultado parece confiable y cuándo podría ser un error técnico

#### Paso 4 Revisar variantes

![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/9.jpg)
![Insertar imagen x](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion4/Imagenes/10.jpg)

Se identificaron 16 variantes en el subconjunto analizado del cromosoma 5. La frecuencia de la mayoría supera el 90%, aunque algunas presentan una frecuencia baja (25%), lo que sugiere potenciales artefactos que deben ser interpretados con cautela. Ademas podemos decir que el número de lecturas de soporte (‘Count’) y la calidad promedio permiten discriminar entre variantes confiables y posibles falsos positivos

### Apreciaciones finales y resumen del trabajo

En este trabajo pude aprender cómo se realiza un análisis de variantes utilizando el flujo *Identify Variants (WES)* del CLC Genomics Workbench. Como los datos que se entregaron corresponden solo a una parte del cromosoma 5, el objetivo principal no fue obtener resultados clínicos reales, sino entender cómo funciona el proceso completo, desde importar las lecturas hasta revisar las variantes en el Genome Browser.

Una de las primeras cosas que revisé fue la calidad de las regiones objetivo. La cobertura promedio fue de 22.5x, lo que está por encima del mínimo de 10 lecturas recomendado para poder llamar variantes. También pude ver que más del 80% de las regiones tenían una cobertura igual o superior a 10, lo que me dio más confianza en que el análisis era razonable. Aun así, había algunas zonas con muy baja cobertura, lo cual entendí que puede deberse a que el dataset está reducido y no corresponde a un exoma completo. Esto me ayudó a comprender que la calidad de los datos es un factor clave para interpretar cualquier análisis genómico.

En cuanto a las variantes, el software identificó 16 en total. La mayoría tenía frecuencias altas, por sobre el 90%, pero también aparecieron algunas con frecuencias mucho más bajas y con pocas lecturas de soporte. Durante la revisión en el Genome Browser pude darme cuenta de que este tipo de variantes deben tratarse con cautela, ya que pueden ser errores técnicos. Esto fue útil para entender que no todas las variantes que aparecen en una tabla son necesariamente reales, y que es importante siempre revisar detalles como la frecuencia, la cantidad de lecturas que la apoyan y la región donde aparece.

En general, este ejercicio me permitió entender mejor cómo se organiza y se interpreta un análisis bioinformático de exoma. Antes de hacerlo, no tenía claridad sobre qué significaban conceptos como cobertura, llamadas de variantes o cómo revisar visualmente una región del genoma. Ahora logro identificar mejor qué elementos revisar para evaluar si un resultado es confiable. Aunque todavía me falta mucho por aprender en este tipo de análisis, siento que este trabajo fue una muy buena introducción para familiarizarme con las herramientas y con la lógica del análisis de variantes en cáncer.
