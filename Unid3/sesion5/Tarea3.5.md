## Tarea 3.5 – nf-core/sarek**

## Informe

**Autor:** Roberto Naranjo Partarrieu 
**Muestra analizada:** **S11**

---

## **Introducción**

El análisis de variantes a partir de datos NGS sigue un flujo estándar compuesto por: (1) preprocesamiento de lecturas FASTQ, (2) alineamiento al genoma de referencia y (3) llamado de variantes. Ejecutar esta secuencia manualmente puede introducir variabilidad y errores, por lo que pipelines estandarizados como **nf-core/sarek** permiten realizar estos pasos de forma reproducible, parametrizable y documentada.

En esta tarea se utilizó SAREK para obtener variantes **germinales** y **somáticas** desde la muestra **S11**, luego se compararon ambas llamadas y se interpretó un subconjunto de variantes utilizando **gnomAD** (germinales) y **OncoKB** (somáticas).

---

Directorio donde estoy trabajando

    cd rnaranjo/unid3/sesion5/pipeline_sarek/code

## **Metodología**

### **1. Organización de carpetas y ambiente**

El análisis se desarrolló dentro del directorio:

`pipeline_sarek/`

Con la siguiente estructura:

- `data/` : archivos FASTQ

- `code/` : scripts `sarek_germinal.sh`, `sarek_somatic.sh`, `local_sarek_8cpus.config`
  
- `results/`: resultados del pipeline

![Figura 1. Creación de directorios y scripts de ejecución](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/1.png)


Los FASTQ originales se encontraban en:

    `~/181004_curso_calidad_datos_NGS/fastq_raw/`

Se copiaron a `data/` y renombraron para mejor manipulación:

```
cp ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R1.fastq.gz .

cp ~/181004_curso_calidad_datos_NGS/fastq_raw/S11_R2.fastq.gz .
```

![Figura 2. Copia de archivos FASTQ originales a data/](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/2.png)


```
  mv S11_R1.fastq.gz R1.fastq.gz
  mv S11_R2.fastq.gz R2.fastq.gz
```

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/3.png)

Antes de ejecutar SAREK se activó el ambiente:

    pyenv activate sarek_taller-pyenv

### **2. Ejecución del pipeline**

Desde `code/` se ejecutaron los análisis:

#### 2.1  **Análisis germinal (HaplotypeCaller)**

```
bash sarek_germinal.sh ../data/R1.fastq.gz ../data/R2.fastq.gz ../results S11
```

El paso de MultiQC falló por un error de conexión al intentar descargar la imagen Singularity correspondiente; el resto del pipeline completó correctamente, por lo que los VCF germinales y somáticos se utilizaron sin el reporte integrado de MultiQC

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/4.png)

#### **2. 2 Análisis somático (Mutect2 tumor-only)**

```
bash sarek_somatic.sh ../data/R1.fastq.gz ../data/R2.fastq.gz ../results S11
```

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/5.png)

SAREK generó los VCF filtrados en:

    results/variant_calling/haplotypecaller/S11/ 
    results/variant_calling/mutect2/S11/`

y un **MultiQC** integrado en:

`results/multiqc/multiqc_report.html`

### **3. Conteo y selección de variantes**

##### 3.1 Contar variantes

```
# Germinal
zcat results/variant_calling/haplotypecaller/S11/S11.haplotypecaller.filtered.vcf.gz \
  | grep -v '^#' | wc -l

# Somático
zcat results/variant_calling/mutect2/S11/S11.mutect2.filtered.vcf.gz \
  | grep -v '^#' | wc -l
```

**Germinales**: 132 variantes 

**Somatico**: 243 variantes 

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/6.png)

### 3.2 Limitación en la anotación con snpEff

Se intentó realizar la anotación funcional con **snpEff**

    snpEff ann GRCh38.86

Sin embargo, la anotación no pudo completarse debido a la ausencia de acceso a internet en el servidor:

    UnknownHostException: snpeff.blob.core.windows.net

Por esta razón, **no se utilizó snpEff para clasificar impacto funcional**, y el análisis posterior se basó directamente en los VCF generados por Sarek, complementado con consultas manuales a **gnomAD** y **OncoKB**

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/7.png)

### 3.3. Seleccionar variantes germinales

Se seleccionaron variantes germinales con **ID (rsID)** y **filtro PASS**:

```
 zcat ../results/variant_calling/haplotypecaller/S11/S11.haplotypecaller.filtered.vcf.gz   | awk '$0 ~ /^#/ {print; next} $3!="." && $7=="PASS"'   > S11_germinal_withID.vcf
```

Generando el archivo: `S11_germinal_candidates.vcf`

Se trabajó con el archivo `S11_germinal_selected.vcf`, que contiene variantes PASS con identificador rsID.  El número de variantes germinales seleccionadas fue calculado mediante:

```
grep -v '^#' S11_germinal_selected.vcf | wc -l
```

Posteriormente, se inspeccionaron las primeras variantes para su análisis:

     grep -v '^#' S11_germinal_selected.vcf | head -10

![](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid3/sesion5/Imagenes/9.png)

### 3.4 Seleccionar variantes somáticas

Se trabajó con el archivo `S11_somatic_selected.vcf`

Las primeras variantes somáticas fueron revisadas mediante:

y las contamos con

```
wc -l S11_somatic_PASS.vcf
```

head -n 20 S11_somatic_PASS.vcf > S11_somatic_selected_body.vcf

Obteniendo **123 variantes**, de las cuales se seleccionaron las primeras 20 para análisis manual

    head -n 20 S11_somatic_PASS.vcf > S11_somatic_selected_body.vcf
    
    zcat ../results/variant_calling/mutect2/S11/S11.mutect2.filtered.vcf.gz \
    | grep '^#' > S11_somatic_selected.vcf
    cat S11_somatic_selected_body.vcf >> S11_somatic_selected.vcf

__________________________________________________________________________________________________

# RESULTADOS

### **1. Calidad general**

El reporte MultiQC mostró lecturas de buena calidad, sin caída notable de Q-score y con alineamiento adecuado a GRCh38. No hubo advertencias críticas que comprometieran el llamado de variantes.

---

## **2. Variantes germinales**

total de variantes germinales:** **132**

Las variantes germinales seleccionadas fueron consultadas en **gnomAD**, observándose que la mayoría corresponden a **polimorfismos comunes** con altas frecuencias alélicas en múltiples poblaciones.

| rsID       | Coordenada (GRCh38) | AF gnomAD  | Interpretación                                     |
| ---------- | ------------------- | ---------- | -------------------------------------------------- |
| rs12621129 | chr2:197400626 T>C  | ~0.32–0.43 | Polimorfismo común, sin evidencia de patogenicidad |
| rs1492765  | chr4:54273811 T>C   | ~0.99      | Variante extremadamente frecuente                  |
| rs869978   | chr4:54273849 T>C   | ~0.75–0.80 | SNP común en múltiples ancestrías                  |

Interpretación germinal:
Las variantes germinales reflejan principalmente variación poblacional normal, sin evidencia de alelos raros o patogénicos según gnomAD.

# **3. Variantes somáticas**

Total de variantes somáticas: 243

Las variantes seleccionadas fueron consultadas manualmente en OncoKB, priorizando genes asociados a cáncer.

Las variantes seleccionadas fueron consultadas manualmente en OncoKB, priorizando genes asociados a cáncer.

| Gen   | Rol              | Evidencia OncoKB | Comentario                                                             |
| ----- | ---------------- | ---------------- | ---------------------------------------------------------------------- |
| TP53  | Supresor tumoral | Nivel 3A (gen)   | Gen frecuentemente mutado en cáncer; variante específica no actionable |
| BRCA1 | Supresor tumoral | Nivel 1 (gen)    | Asociado a reparación de DNA, sin implicancia terapéutica directa      |
| JAK2  | Oncogén          | Nivel 2 (gen)    | Alteraciones frecuentes en cáncer hematológico                         |

Interpretación somática:
La mayoría de las variantes somáticas no presentan anotación clínica directa, sugiriendo eventos pasajeros o mutaciones de significado clínico incierto.

# 4. Comparación germinal vs somático 

| Métrica               | Germinal          | Somático      |
| --------------------- | ----------------- | ------------- |
| Nº variantes          | 132               | 243           |
| Origen                | Constitucional    | Adquirido     |
| Predominio            | SNP poblacionales | SNV tumorales |
| Impacto clínico       | Bajo              | Incierto      |
| Variantes compartidas | 0                 | –             |


No se detectaron variantes compartidas entre los conjuntos germinal y somático. 
Este resultado es consistente con el enfoque metodológico del pipeline nf-core/sarek, ya que Mutect2 (modo tumor-only) aplica filtros poblacionales y heurísticos para excluir variantes germinales del conjunto somático.

Por lo tanto, las variantes detectadas en el análisis somático corresponden a eventos 
adquiridos, mientras que las variantes germinales reflejan polimorfismos constitucionales 
presentes en la población general.


# 5.Discusión y conclusiones

Este trabajo demuestra que nf-core/sarek permite obtener de forma confiable variantes germinales y somáticas a partir de una misma muestra. A pesar de limitaciones técnicas que impidieron la anotación automática con snpEff, la integración de gnomAD y OncoKB permitió contextualizar las variantes detectadas.

En conjunto, la muestra S11 no presenta variantes con relevancia clínica clara, y la mayoría de las alteraciones corresponden a polimorfismos germinales comunes o mutaciones somáticas de significado incierto. Esto resalta la importancia de combinar pipelines robustos con bases de datos externas para una correcta priorización e interpretación de variantes.
