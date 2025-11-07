# Tarea unid3sesion1

#### autor: Roberto Naranjo Partarrieu

## Objetivo: Evaluar la calidad de las lecturas NGS crudas y podadas mediante FastQC, y comparar R1 vs R2.

Usamos 

     module load FastQC

y ejecutamos FastQC en las lectura cruda de nuestra muestra S11

    fastqc ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R1.fastq.gz -o .
    fastqc ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R2.fastq.gz -o .

Lo cual nos genera: S11_R1_fastqc.html, S11_R1_fastqc.zip, S11_R2_fastqc.html , S11_R2_fastqc.zip

Insertar imagen 1

Ejecutamos FastQC en las lecturas filtradas y podadas 

    fastqc ~/181004_curso_calidad_datos_NGS/fastq_filter/S11_R1_filter.fastq.gz -o .
    fastqc ~/181004_curso_calidad_datos_NGS/fastq_filter/S11_R2_filter.fastq.gz -o .

S11_R1_filter_fastqc.html,  S11_R1_filter_fastqc.zip, S11_R2_filter_fastqc.html, S11_R2_filter_fastqc.zip

Insertar imagen 2

Usando los comandos Unix :

##### Contaremos el numero de lecutras en un archivo fastq usando el comando

    wc -l ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R1.fastq.gz
    wc -l ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R2.fastq.gz

Lo cual nos informa que hay ![Insertar imagen 1](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/1.png)

R1 21294 lecturas

R2 26725 lecturas

##### Previsualizar las primeras 40 líneas del mismo archivo fastq

usamos el comando: 

    zcat ~/181004_curso_calidad_datos_NGS/fastq_raw/S3_R1.fastq.gz | head -n 40
    zcat ~/181004_curso_calidad_datos_NGS/fastq_raw/S3_R2.fastq.gz | head -n 40

![Insertar imagen 2](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/2.png)

![Insertar imagen 3](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/3.png)

###### Ubicar la lectura 3 e identificar la información disponible. Describe en detalle la información entregada. ¿Donde se entrega la calidad del read?, ¿Cuál es el ID (identificador) del read? Etc. Utilice fechas y etiquetas para identificar cada parte.

El ID es `@M03564:2:000000000–D29D3:1:1101:14451:1389 1:N:0:ACAGTGG+TAGACCTA`

    El cual tiene la siguiente información: 

`M03564`: nombre del instrumento de secuenciación (Illumina MiSeq). naranjo

`2:000000000–D29D3`: ID del run y flowcell. azul

`1:1101:14451:1389`: lane 1, tile 1101, coordenadas (x=14451, y=1389). verde

`2`: indica read 2 (R2) rojo 

`N`: no pasa un filtro de control rojo

`0`: indica que no es una lectura control PhiX.

`ACAGTGG+TAGACCTA`: los índices (index1 e index2) usados para multiplexado. amarillo

calidad celeste 
![Insertar imagen 4](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/4.png)

##### Traducir el código de calidad para las primeras 10 bases del tercer leer a valores numéricos (Q) usando la codificación entregada en clase.

usaremos: 

    zcat ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R1.fastq.gz | sed -n '12p' \
    | perl -ne 'chomp; @c=split(//); for($i=0;$i<10 && $i<@c;$i++){ printf("%s\t%d\n", $c[$i], ord($c[$i])-33) }'

lo que nos permite traducir el codigo y obtener los valores numericos 



![Insertar imagen 5](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/5.png)

## PARTE 2

Bajamos los archivos HTML a nuestro equipo

    scp bioinfo1@genoma.med.uchile.cl:~/rnaranjo/unid3/sesion1/S11*_fastqc.html .

Ahora interpretamos nuestros datos visualizados en FastQC 

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/R1raw.png)
![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/R1filt.png)
![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/R2raw.png)
![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion1/Imagenes/R2filt.png)

respecto a eso podemos decir que la poda o trimming eliminó aproximadamente un 13–14 % de lecturas en ambos extremos del par. El %GC se mantiene casi constante (45 a 44–45 %), indicando que el filtrado no introdujo sesgos composicionales. Las longitudes variables en los archivos filtrados evidencian recorte de adaptadores y bases de baja calidad en los extremos 3 
Ahora si desglsamos respecto a la calidas de por base de la secuencia(A): 

* R1 crudo: La calidad por base se mantiene alta (Q ≥ 35) pero muestra una ligera caída en el extremo 3′.
* R1 filtrado: Los boxplots se estabilizan completamente en el rango verde (Q ≥ 33), evidenciando una mejora post-trimming.
* R2 crudo: La caída 3′ es más pronunciada, con dispersiones en el rango amarillo y rojo (bases Q < 25).
* R2 filtrado: El trimming corrige esta caída; la media de calidad se mantiene ≥ 30 hasta el final.

Lo cual nos hace ver que el trimming funciono o es mas notable para R2

Asi tambien nos lo hace ver el grafico "Per Sequence Quility Score" (B) donde el proceso de poda eliminó lecturas defectuosas y mejoró la homogeneidad general de la calidad promedio por read.

Por otra parte  el contenido de base por secuencia (C) la limpieza redujo los sesgos iniciales y estabilizó el contenido de bases por posición. Finalmente 
