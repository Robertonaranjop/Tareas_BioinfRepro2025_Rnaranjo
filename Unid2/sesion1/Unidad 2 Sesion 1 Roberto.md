# Unidad 2 Sesion 1

## **Ejercicios**

### Parte 1, VCFtools

antes de todo activamos vcftools 

`module load vcftools`

Trabajaremos con el archivo:`GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf`. El cual no descargaremos, si no llamaremos desde la nube, para responder las siguientes preguntas 

    *usaremos la función vcftools --gzvcf , la cual nos permite llamar el archivo que esta comprimido con ** gzip** y  --vcf cuando no este comprimdo. 

1. **¿Cuántos individuos y variantes (SNP) tiene el archivo?**
   
   usamos el comando: 
   
   ```
   $ vcftools  --gzvcf /datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf`
   ```
   
      R: Se reporta 18 individuos, con 4450360 Sitios. 

![](IUnid2/sesion1/Imagenes/1.png)

2. **¿Cuántos sitios del archivo no tienen datos perdidos?**
   
   usamos el comando:
   
   ```
      $ vcftools --vcf/datos/compartido/ChileGenomico/GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf --max-missing 1
   ```
   
     R: Se reporta que hay 382626 sitios observados a los que no les falta información, de un total de 4450360.

![](C:\Users\rnara\Desktop\Repositio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\2.png)

4. **Genera un archivo en tu carpeta de trabajo `sesion1/data`que contiene solo SNPs en una ventana de 2Mb en cualquier cromosoma. Nombra el archivo `CLG_Chr<X>_<Start>-<End>Mb.vcf`donde es número del cromosoma, es el inicio de la ventana genómica y es el final en megabases.**

     

    usamos el comando: 

```
  vcftools --vcf GATK_ChGdb_recalibrated.autosomes.12262013.snps.known.vcf \ --chr 12 --from-bp 5000000 --to-bp 7000000 \ --remove-indels --min-alleles 2 --max-alleles 2 \ --recode --stdout > CLG_Chr12_5-7Mb.vcf 
```

R:  Se trabajó en el cromosoma 5, se mantuvieron 18 individuos, se seleccionaron 4392 SNPs bialélicos dentro de la ventana 1–3 Mb

*El archivo resultante se guardó en: `../results/CLG_Chr5_1-3Mb.recode.vcf`

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\3.png)

5. **Reporta cuantas variantes tienen el archivo generado**

usamos el comando:  

```
grep -cv '^#' ../results/CLG_Chr5_1-3Mb.recode.vcf 
```

el cual cuenta todas las lineas que no comienzan con ` # ` 

R: 4392 SNPs bialélicos

6. **Reporta la cobertura promedio para todos los individuos del conjunto de datos**
   R: 2.41993
   usamos el comando: 
   
   ```vcftools
   --depth --out ../results/cobertura_promedio
   ```
   
   el cual nos permite almacenar la informacion en el directorio results, 

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\4.png)

verificamos el archivo usando:  head ../results/cobertura_promedio.idepth

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\5.png)

para poder calcular la cobertura promedio total usamos: 

```
awk 'NR>1 {sum+=$3; n++} END {print "Cobertura promedio total:", sum/n}' ../results/cobertura_promedio.idepth
```

 la expliación de este codigo es: 
     `NR>1` le dice a `awk`: ➤ “Ignora la primera línea (el encabezado)”.
   ` {sum+=$3; n++}``$3`representa la tercera columna del archivo `(MEAN_DEPTH)`.
   Calcula el promedio dividiendo `sum/n.`

7. **Calcula la frecuencia de cada alelo para todos los individuos dentro del archivo y guarda el resultado en un archivo.**
   
   usamos el comando: 
   
   ```
   frecuencias_bialelicas_chr5.frqvcftools --vcf ../results/CLG_Chr5_1-3Mb.recode.vcf \
   --freq --out ../results/frecuencias_alelos_chr5`
   ```
   
   para visualizar el resultadp usamos: 
   
   ```
    head ../results/frecuencias_alelos_chr5.frq`
   ```
   
   R:  los resultados se encuentran en la carpeta `frecuencias_alelos_chr5.frq`

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\6.png)

8. **Filtra el archivo de frecuencias para solo incluir variantes bialélicas (tip: awk puede ser útil para realizar esta tarea, tip2: puedes usar bcftools para filtrar variantes con más de dos alelos antes de calcular las frecuencias)**
   
   Usaremos el comando 
   
   ```
   frecuencias_bialelicas_chr5.frq 
   awk '$3==2' ../results/frecuencias_alelos_chr5.frq > ../results/frecuencias_bialelicas_chr5.frq`
   ```
   
   y para verificar:
   
   ```
    head ../results/frecuencias_bialelicas_chr5.frq
   ```
   
    R: los resultado se encuentra en results en el archivo: `frecuencias_bialelicas_chr5.frq`
   
   ![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\7.png)

9. **Llama a un guión escrito en lenguaje R que lee el archivo de frecuencias de variantes bialélicas y guarda un histograma con el espectro de MAF para las variantes bialélicas**.
   
    R:  primerto creamos nuestro scrip en code nano ../code/plot_maf_simple_chr5.R  sin embargo me equivoque y lo hice en data por lo cual tuve que usar el comando mv plot_maf_simple.R ../code/ para relocalizarlo. el contenido del script es el siguiente:  
   
   ```
      # === plot_maf_simple_chr5.R ===
   
   # Script para generar histograma del espectro de frecuencias alélicas (MAF)
   
   # a partir de un archivo de frecuencias bialélicas (.frq)
   
    # Leer archivo de frecuencias
   
   maf_file <- "../results/frecuencias_bialelicas_chr5.frq"
   
   # Verificar que el archivo existe
   
   if (!file.exists(maf_file)) {
     stop("El archivo de frecuencias no existe: ", maf_file)
   }
   
   # Leer los datos
   
   freq <- read.table(maf_file, header = TRUE, stringsAsFactors = FALSE)
   
   # Extraer las frecuencias alélicas de las columnas con formato A:0.472
   
   get_freq <- function(x) {
     as.numeric(sub(".*:", "", x))
   }
   
   # Convertir frecuencias
   
   allele_freqs <- apply(freq[, 5:ncol(freq)], 2, get_freq)
   
   # Calcular el MAF (menor frecuencia entre los dos alelos)
   
   maf <- apply(allele_freqs, 1, min, na.rm = TRUE)
   
   # Eliminar valores no válidos
   
   maf <- maf[!is.na(maf)]
   
   # Crear histograma y guardar imagen
   
   png("../results/histograma_MAF_chr5.png", width = 900, height = 650)
   hist(maf, breaks = 30, col = "lightblue", border = "black",
        main = "Espectro de Frecuencias de Alelos Menores (MAF)",
        xlab = "Frecuencia del Alelo Menor (MAF)",
        ylab = "Número de variantes")
   dev.off()
   
   # Guardar estadísticas básicas
   
   summary_file <- "../results/resumen_MAF_chr5.txt"
   sink(summary_file)
   cat("Número total de variantes:", length(maf), "\n")
   cat("Número de variantes con MAF < 0.05:", sum(maf < 0.05), "\n")
   cat("Porcentaje con MAF < 0.05:", round(sum(maf < 0.05) / length(maf) * 100, 2),    "%\n\n")
   cat("Resumen estadístico del MAF:\n")
   print(summary(maf))
   sink()
   
   cat("Histograma guardado en ../results/histograma_MAF_chr5.png\n")
   cat("Resumen guardado en ../results/resumen_MAF_chr5.txt\n")```    
   ```
   
    ahora para poder usar r usamos el comando:`  module load R` , verificamos escribiendo ` R`  una vez en el entorno ejecutamos: 
   
   ```
   source("../code/plot_maf_simple_chr5.R")
   ```
   
    revisamos si se guardo el histograma con: 
   
   ```
              ls -lh ../results/histograma_MAF_chr5.png ../results/resumen_MAF_chr5.txt
   ```

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\8.png)
![](C:/Users/rnara/Desktop/Repositorio/Tareas_BioinfRepro2025_Rnaranjo/Unid2/sesion1/histograma_MAF_chr5.png)

10.**¿Cuántos sitios tienen una frecuencia del alelo menor <0.05?**

 R: Porcentaje de variantes con MAF < 0.05: 41.8 %

10. **Calcula la heterocigosidad de cada individuo.**
    
    usamos el comando: 
    
    ```
    vcftools --vcf ../results/CLG_Chr5_1-3Mb.recode.vcf --het --out ../results/heterocigosidad_chr5
    ```
    
    R: heterocigosidad_chr5 en results
    
    que significa cada columna: 
    
    - **O(HOM)** = número de loci **homocigotos observados**
    
    - **E(HOM)** = número de loci **homocigotos esperados** (bajo equilibrio H–W)
    
    - **N_SITES** = número total de sitios genotipados por individuo
    
    - **F** = **índice de endogamia** → mide desviación de Hardy-Weinberg
    
    ![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\9.png)

11. **Calcula la diversidad nucleotídica por sitio.**
    
    usamos el  comando
    
    ```
    awk 'NR>1 {sum+=$3; n++} END {print "Diversidad nucleotídica promedio (π):", sum/n}' ../results/diversidad_nucleotidica_chr5.sites.pi
    ```
    
    R:el resultado es 0.349729
    
    ![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\10.png)

12. **Filtra los sitios que tengan una frecuencia del alelo menor <0.05**
    
    usamos el comando
    
    ```
    vcftools --vcf ../results/CLG_Chr5_1-3Mb.recode.vcf --maf 0.05 --out ../results/CLG_Chr5_maf05
    ```
    
    R: 3894 de 4392 almacenados en CLG_Chr5_maf05 

13. **Convierta el archivo `wolves_maf05.vcf`a formato plink.
    
    usamos el comando:
    
    ```
    vcftools --vcf ../results/wolves_maf05.vcf --plink --out ../results/wolves_maf05 
    ```
    
        no existia el archivo por lo cual debemos renombralos usando: 
    
    ```
    mv mv ../results/CLG_Chr5_maf05.recode.vcf ../results/wolves_maf05.vcf
    ```
    
    ahora para convertir el archivo usamos
    
    ```
    vcftools --vcf ../results/wolves_maf05.vcf --plink --out ../results/wolves_maf05
    ```
    
    ### Parte 2, PLINK

Pasamos a plink

Copia esos archivos a tu respositorio en una carpeta para la sesión `Unidad2/sesion1/data/chilegenomico` `cp /datos/compartido/ChileGenomico/chilean_all48_hg19.* ../data/`

1. Enlista los archivos plink que hay en `data`. ¿Qué tipos de archivos son cada uno?
   
   tenmos archivos bim, bed y fam de la base de datos copiada

2. Consulta el manual de [plink1.9](https://www.cog-genomics.org/plink/1.9/formats) y contesta utilizando comandos de plink lo siguiente. Deposita cualquier archivo que generes en la carpeta `Unididad2/Prac_Uni2/results`:
   
   a) Transforma de formato bed a formato ped (pista: sección Data Managment). El nombre del output debe ser igual, solo cambiando la extensión.

```
#crear directorio
mkdir -p ~/rnaranjo/unid2/sesion1/data/chilegenomico

#copiar los archivos
cp /datos/compartido/ChileGenomico/chilean_all48_hg19.* ~/rnaranjo/unid2/sesion1/data/chilegenomico/
```

a) Transforma de formato bed a formato ped (pista: sección Data Managment). El nombre del output debe ser igual, solo cambiando la extensión.

```
$ plink --bfile ../data/chilean_all48_hg19 --recode --out ../results/chilean_all48_hg19
```

b) Crea otro archivo ped (ojo PPPPed) pero esta vez filtrando los SNPs cuya frecuencia del alelo menor sea menor a 0.05 Y filtrando los individuos con más de 10% missing data. Tu output debe llamarse maicesArtegaetal2015_maf05_missing10

¿Cuántos SNPs y cuántos individuos fueron removidos por los filtros?

```
$ plink --bfile ../data/chilean_all48_hg19 --recode --maf 0.05  --mind 0.1 --out ../results/chilean_all48_hg19_maf05_missing10
```

c) Realiza un reporte de equilibrio de Hardy-Weinberg sobre el archivo `chilean_all48_hg19_maf05_missing10` creado en el ejercicio anterior. El nombre del archivo de tu output debe contener chilean_all48_hg19_maf05_missing10.

```
plink --file ../results/chilean_all48_hg19_maf05_missing10 --hardy --out ../results/chilean_all48_hg19_maf05_missing10
```

Observa el output y discute que es cada columna.

```
head ../results/chilean_all48_hg19_maf05_missing10.hwe
10.hwe
 CHR                          SNP     TEST   A1   A2                 GENO   O(HET)   E(HET)            P
   1                    rs9701055      ALL    T    C              18/0/28        0   0.4764    5.994e-14
   1                    rs9701055      AFF    T    C                0/0/0      nan      nan            1
   1                    rs9701055    UNAFF    T    C              18/0/28        0   0.4764    5.994e-14
   1                    rs9701055      ALL    T    C              0/16/28   0.3636   0.2975       0.3137
   1                    rs9701055      AFF    T    C                0/0/0      nan      nan            1
   1                    rs9701055    UNAFF    T    C              0/16/28   0.3636   0.2975       0.3137
   1                    rs2073813      ALL    A    G              0/17/28   0.3778   0.3064       0.3197
   1                    rs2073813      AFF    A    G                0/0/0      nan      nan            1
   1                    rs2073813    UNAFF    A    G              0/17/28   0.3778   0.3064       0.3197
```

En la siguiente tabla discutire que significa cada columna de mi resultado: 

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\11.png)

| Columna     | Ejemplo               | Descripción e interpretación                                                                                                                                                                                                                                                                                           |
| ----------- | --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **CHR**     | `1`                   | **Cromosoma** donde se encuentra el SNP (aquí, el cromosoma 1). PLINK lo toma del archivo `.bim`. Es útil para ubicar la variante en el genoma.                                                                                                                                                                        |
| **SNP**     | `rs9701055`           | **Identificador único del polimorfismo** (SNP ID). Los “rsID” provienen de bases de datos como *dbSNP* y permiten rastrear la variante y su posición.                                                                                                                                                                  |
| **TEST**    | `ALL`, `AFF`, `UNAFF` | Indica el **grupo de individuos** usado en el test: • `ALL`: todos los individuos del conjunto. • `AFF`: individuos “afectados” (casos). • `UNAFF`: “no afectados” (controles). En tu dataset, no hay fenotipo de caso/control, así que PLINK repite las tres pruebas pero `AFF` y `UNAFF` muestran `nan` (no aplica). |
| **A1 / A2** | `T` / `C`             | **Alelos observados** en esa posición: • `A1` = alelo menor (minor allele) o alternativo. • `A2` = alelo mayor (major allele) o referencia. Esto depende de la frecuencia en el conjunto.                                                                                                                              |
| **GENO**    | `18/0/28`             | **Conteo de genotipos** observados: → `A1A1 / A1A2 / A2A2` En este ejemplo: 18 homocigotos A1A1, 0 heterocigotos, 28 homocigotos A2A2. Esto permite calcular las frecuencias genotípicas.                                                                                                                              |
| **O(HET)**  | `0` o `0.3636`        | **Proporción observada de heterocigotos** = número de A1A2 dividido por el total de individuos. Por ejemplo, `0.3636` indica que ≈36 % de los individuos son heterocigotos.                                                                                                                                            |
| **E(HET)**  | `0.4764`              | **Proporción esperada de heterocigotos** bajo equilibrio Hardy–Weinberg, calculada como `2pq`, donde `p` y `q` son las frecuencias alélicas. Sirve para comparar lo observado vs. lo esperado.                                                                                                                         |
| **P**       | `5.994e-14`, `0.3137` | **Valor p del test exacto de Hardy–Weinberg**: mide la probabilidad de obtener una desviación tan extrema (o mayor) si la población estuviera en equilibrio. • **P > 0.05:** equilibrio aceptado. • **P < 0.05:** desviación significativa → posible error de genotipado, selección o estructura poblacional.          |

d) Observa el archivo `maicesArtegaetal2015.fam`. Consulta la documentación de plink para determinar que es cada columna. ¿Qué información hay y no hay en este archivo?

usamos el comando:

 head ~/rnaranjo/unid2/sesion1/data/chilegenomico/chilean_all48_hg19.fam

```
$ head ../data/chilean_all48_hg19.fam
CDSJ177 CDSJ177 0 0 1 1
CDSJ021 CDSJ021 0 0 1 1
ARI006 ARI006 0 0 1 1
ARI021 ARI021 0 0 1 1
ARI022 ARI022 0 0 2 1
CDSJ174 CDSJ174 0 0 1 1
CDSJ175 CDSJ175 0 0 1 1
CDSJ046 CDSJ046 0 0 1 1
CDSJ176 CDSJ176 0 0 1 1
CDSJ469 CDSJ469 0 0 2 1
```

respecto a eso: 

La **columna 1** corresponde al **identificador familiar** o *Family ID (FID)*, que permite agrupar individuos pertenecientes a una misma familia o población.

La **columna 2** es el **identificador individual** o *Individual ID (IID)*, que corresponde al nombre único asignado a cada muestra dentro del conjunto de datos.

La **columna 3** corresponde al **ID del padre** (*Paternal ID*). Si el valor es `0`, significa que no se conoce o no se dispone de información sobre el progenitor.

La **columna 4** corresponde al **ID de la madre** (*Maternal ID*). De igual forma, el valor `0` indica que se desconoce la información materna.

La **columna 5** representa el **sexo del individuo**, donde `1` indica **masculino**, `2` **femenino**, y `0` **desconocido**.

Finalmente, la **columna 6** corresponde al **fenotipo**, el cual puede tomar los siguientes valores: `1` = **control/no afectado**, `2` = **caso/afectado**, y `-9` o `0` = **dato perdido o no disponible**.    

e)Utilice la información del archivo `data/chilean_all48_hg19_popinfo.csv`y el comando `update-ids`de plink para cambiar los nombres de las muestras de `data/chilean_all48_hg19.fam`tal forma que el ID de familia corresponda a la información de la columna `Ancestry`en `chilean_all48_hg19_popinfo.csv`. Pista: este ejercicio requiere varias operaciones, puedes dividirlas en diferentes scripts de bash o de R y bash. Tu respuesta debe incluir todos los scripts (y deben estar en /code).

R:  Primero crearemos el script en code llamado prepare_update_ids.R

usando el comando en /code:

              nano prepare_update_ids.R

el cual contiene el siguiente script: 

```
# ==========================
# Script: prepare_update_ids_ancestry.R
# Crear un archivo update_ids.txt para PLINK
# usando la columna 'Ancestry' como nuevo Family ID (FID)
# ==========================

# Leer archivos de entrada -------------------
popinfo <- read.csv("../data/chilegenomico/chilean_all48_hg19_popinfo.csv", stringsAsFactors = FALSE)
fam <- read.table("../data/chilegenomico/chilean_all48_hg19.fam", stringsAsFactors = FALSE)

# Agregar nombres de columnas al .fam
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# ----------------------------
# Normalización de IDs
# ----------------------------
# Eliminar guiones o inconsistencias
popinfo$IndID_modified <- gsub("-", "", popinfo$IndID)
fam$IID_clean <- toupper(fam$IID)

# Verificaciones de columnas obligatorias
if (!"IndID" %in% names(popinfo)) stop("En popinfo no existe la columna 'IndID'.")
if (!"Ancestry" %in% names(popinfo)) stop("En popinfo no existe la columna 'Ancestry'.")

# ----------------------------
# Fusión entre popinfo y fam
# ----------------------------
merged <- merge(fam, popinfo, by.x = "IID_clean", by.y = "IndID_modified", all.x = TRUE)

# Filtrar filas válidas (que sí tienen información de Ancestry)
valid <- !is.na(merged$Ancestry) & merged$Ancestry != ""

# Crear el data.frame con IDs actualizados
update_ids <- data.frame(
  FID_old = merged$IID_clean[valid],
  IID_old = merged$IID_clean[valid],
  FID_new = merged$Ancestry[valid],
  IID_new = merged$IID_clean[valid],
  check.names = FALSE
)

# ----------------------------
# Guardar archivo para PLINK
# ----------------------------
write.table(
  update_ids,
  "../results/update_ids.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat("✅ Archivo update_ids.txt escrito en ../results con", nrow(update_ids), "filas.\n")

# ----------------------------
# Validación final
# ----------------------------
mism <- sum(!valid)
if (mism > 0) {
  cat("⚠️ Aviso:", mism, "muestras del .fam no encontraron 'Ancestry' en popinfo.\n")
}
```

Ejecutamos el script: 

Rscript ../code/prepare_update_ids_ancestry.R

donde se crea el archivo:  update_ids.txt

### Aplicar los nuevos IDs con PLINK

`plink \   --bfile ../data/chilegenomico/chilean_all48_hg19 \   --update-ids ../results/update_ids.txt \   --make-bed \   --out ../results/chilean_all48_hg19_updated_ids`

 Si sale “Duplicate sample ID … in --update-ids file”, deduplica el archivo (por si quedó alguna fila repetida):

`sort ../results/update_ids.txt | uniq > ../results/update_ids_unique.txt plink \   --bfile ../data/chilegenomico/chilean_all48_hg19 \   --update-ids ../results/update_ids_unique.txt \   --make-bed \   --out ../results/chilean_all48_hg19_updated_ids`

aparecio: Error: Duplicate sample ID 'ARI001 ARI001' in --update-ids file.

    ## Arreglarlo en 2 pasos

### 1) Vemos qué duplicados hay (por FID_old/IID_old)

`awk '{print $1,$2}' ../results/update_ids.txt | sort | uniq -c | awk '$1>1{print}'`

Si imprime filas, son los `IID_old`/`FID_old` repetidos.

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\12.png)

### 2) Crearmos un archivo “dedup” que deje **una sola línea por individuo**

Conserva la **primera** aparición por `IID_old` (columna 2):

`awk '!seen[$2]++' ../results/update_ids.txt > ../results/update_ids_dedup.txt`

(Esto elimina duplicados por **columna 2**. Si quisieras deduplicar por el par `FID_old IID_old`, usarías `!seen[$1 FS $2]++`).

Comprueba rápido que ya no hay duplicados:

`awk '{print $1,$2}' ../results/update_ids_dedup.txt | sort | uniq -c | awk '$1>1{print}' # (no debería imprimir nada)`

Volvemos a ejecturar plink 

plink \
  --bfile ../data/chilegenomico/chilean_all48_hg19 \
  --update-ids ../results/update_ids_dedup.txt \
  --make-bed \
  --out ../results/chilean_all48_hg19_updated_ids

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\13.png)

d) Realiza una cuna comparación entre el sexo y archivo `fam`y el `popinfo`y calcula la proporción de discordancias

primero crearemos un scrip que se almacene en ../code con el comando nano:

                  nano compare_fam_popinfo_sex.R

el cual contiene

```
# compare_fam_popinfo_sex.R

options(stringsAsFactors = FALSE)

# Rutas (ajústalas si usas otras)

popinfo_fp <- "../data/chilegenomico/chilean_all48_hg19_popinfo.csv"
fam_fp     <- "../results/chilean_all48_hg19_updated_ids.fam"  # usa el fam con FID actualizado
out_mis_fp <- "../results/sex_fam_vs_popinfo_mismatches.tsv"
out_sum_fp <- "../results/sex_fam_vs_popinfo_summary.csv"

# 1) Leer datos

pop <- read.csv(popinfo_fp)
fam <- read.table(fam_fp)
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")

# 2) Normalizar IDs y sexo de popinfo

# - En popinfo el ID suele estar en la columna "IndID"

# - El sexo en "Sex" como "F"/"M". Convertimos a plink: 1=masc, 2=fem

if (!("IndID" %in% names(pop))) stop("No se encontró columna 'IndID' en popinfo.")
if (!("Sex"   %in% names(pop))) stop("No se encontró columna 'Sex' en popinfo.")

pop$IndID_clean <- toupper(gsub("-", "", pop$IndID))
sex_map <- c("M"=1, "F"=2, "Male"=1, "Female"=2, "Hombre"=1, "Mujer"=2)
pop$SEX_info <- unname(sex_map[as.character(pop$Sex)])
pop$SEX_info[is.na(pop$SEX_info)] <- NA_integer_

fam$IID_clean <- toupper(fam$IID)

# 3) Merge por IID

m <- merge(fam, pop[, c("IndID_clean","SEX_info")],
           by.x="IID_clean", by.y="IndID_clean", all.x=TRUE)

# 4) Filtrar a los casos comparables (ambos sexos conocidos)

m$SEX_fam  <- as.integer(m$SEX)       # 1=masc, 2=fem, 0=desconocido
m$SEX_info <- as.integer(m$SEX_info)  # 1=masc, 2=fem, NA=desconocido en popinfo

cmp <- subset(m, !is.na(SEX_info) & SEX_fam %in% c(1,2))

total_cmp <- nrow(cmp)
discord   <- sum(cmp$SEX_fam != cmp$SEX_info)
prop_disc <- if (total_cmp > 0) discord / total_cmp else NA_real_

# 5) Guardar resultados

mis <- subset(cmp, SEX_fam != SEX_info,
              select=c(FID, IID, SEX_fam, SEX_info))
names(mis) <- c("FID","IID","SEX_fam(1=masc,2=fem)","SEX_popinfo(1=masc,2=fem)")
write.table(mis, out_mis_fp, sep="\t", row.names=FALSE, quote=FALSE)

sumdf <- data.frame(
  total_comparables = total_cmp,
  discordantes      = discord,
  proporcion        = if (is.na(prop_disc)) NA else round(prop_disc, 4)
)
write.csv(sumdf, out_sum_fp, row.names = FALSE)

# 6) Mensaje en consola

cat("Comparación sexo FAM vs POPINFO\n")
cat("Individuos comparables :", total_cmp, "\n")
cat("Discordantes           :", discord, "\n")
cat("Proporción discordante :", if (is.na(prop_disc)) "NA" else round(prop_disc,4), "\n")
cat("Mismatches guardados en:", out_mis_fp, "\n")
cat("Resumen guardado en    :", out_sum_fp, "\n")
```

e) Realiza una comparación entre el sexo y archivo `fam`y el `popinfo` y calcula la proporción de discordancias

```
plink \
  --bfile ../results/chilean_all48_hg19_updated_ids \
  --check-sex \
  --out ../results/sexcheck![](C:\Users\rnara\AppData\Roaming\marktext\images\2025-10-15-23-15-56-image.png)
```

ahora para el analisis crearemos un script para el analisis usando el comando: 

```
 nano compare_sex_discordance.R
```

el cual contiene: 

```
# --- 1. Cargar datos ---

fam <- read.table("../results/chilean_all48_hg19_updated_ids.fam", stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX_FAM", "PHENO")

sexcheck <- read.table("../results/sexcheck.sexcheck", header = TRUE, stringsAsFactors = FALSE)
popinfo <- read.csv("../data/chilegenomico/chilean_all48_hg19_popinfo.csv", stringsAsFactors = FALSE)

# --- 2. Armonizar IDs ---

fam$IID <- toupper(fam$IID)
sexcheck$IID <- toupper(sexcheck$IID)
popinfo$IndID <- gsub("-", "", toupper(popinfo$IndID))

# --- 3. Unir tablas ---

merged <- merge(fam, sexcheck[, c("FID", "IID", "SNPSEX")], by = "IID", all.x = TRUE)
merged <- merge(merged, popinfo[, c("IndID", "Sex")], by.x = "IID", by.y = "IndID", all.x = TRUE)

# --- 4. Convertir SEX reportado a numérico (1=M, 2=F) ---

merged$SEX_POPINFO <- ifelse(merged$Sex == "F", 2,
                        ifelse(merged$Sex == "M", 1, NA))

# --- 5. Comparar concordancias ---

merged$discord_fam_popinfo <- merged$SEX_FAM != merged$SEX_POPINFO
merged$discord_fam_genetic <- merged$SEX_FAM != merged$SNPSEX
merged$discord_popinfo_genetic <- merged$SEX_POPINFO != merged$SNPSEX

# --- 6. Calcular proporciones ---

prop_fam_popinfo <- mean(merged$discord_fam_popinfo, na.rm = TRUE)
prop_fam_genetic <- mean(merged$discord_fam_genetic, na.rm = TRUE)
prop_popinfo_genetic <- mean(merged$discord_popinfo_genetic, na.rm = TRUE)

cat("Discordancia FAM vs POPINFO:", round(prop_fam_popinfo, 3), "\n")
cat("Discordancia FAM vs GENÉTICO:", round(prop_fam_genetic, 3), "\n")
cat("Discordancia POPINFO vs GENÉTICO:", round(prop_popinfo_genetic, 3), "\n")

# --- 7. Guardar tabla combinada ---

write.csv(merged, "../results/sex_discordance_summary.csv", row.names = FALSE)
cat("\nArchivo guardado en ../results/sex_discordance_summary.csv\n")
```

Ejecutamos el codigo

`Rscript ../code/compare_fam_popinfo_sex.R`

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\14.png)

R:  61 individuos con información de sexo conocida en ambos archivos, se detectó 1 discordancia, correspondiente a la muestra Aymara ARI022, cuyo sexo fue registrado como femenino (2) en el archivo .fam y masculino (1) en el archivo popinfo.

La proporción de discordancia fue de 0.0164 (1.64%), lo que indica una alta concordancia entre ambas fuentes de información. Este resultado sugiere que los metadatos fenotípicos y la información genotípica presentan un adecuado nivel de consistencia, con errores mínimos atribuibles a posibles errores de registro o etiquetado de muestra.

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\15.png)

f)Realice una prueba de estimación de sexo usando plink y reporte los resultados en formato de tabla para todos los individuos con discordancia entre el sexto reportado en famy el calculado con plink.

R: para la estimación de sexo usamos: 

          plink \
      --bfile ../results/chilean_all48_hg19_updated_ids \
      --check-sex \
      --out ../results/sexcheck

Para aislar a los discordantes

creamos el script en code llamado `extract_sex_discordant.R`

     # Cargar datos del sexcheck de PLINK
    sexcheck <- read.table("../results/sexcheck.sexcheck", header = TRUE)
    
    # Filtrar solo los individuos con problemas de sexo
    discordant <- subset(sexcheck, STATUS == "PROBLEM")
    
    # Seleccionar columnas relevantes
    discordant_summary <- discordant[, c("FID", "IID", "PEDSEX", "SNPSEX", "F", "STATUS")]
    
    # Guardar resultados en tabla resumen
    write.table(discordant_summary, "../results/sexcheck_discordant_summary.csv",
                sep = ",", row.names = FALSE, quote = FALSE)
    
    # Mostrar resultados
    cat("Individuos con discordancia de sexo genético:\n")
    print(discordant_summary)
    cat("\nTabla guardada en ../results/sexcheck_discordant_summary.csv\n")

lo guardamos y ejectuamos: 

    R --vanilla -q -e "source('../code/extract_sex_discordant.R')"

![](C:\Users\rnara\Desktop\Repositorio\Tareas_BioinfRepro2025_Rnaranjo\Unid2\sesion1\Imagenes\16.png)

R: Se realizó la estimación de sexo genético utilizando PLINK v1.9 a partir de los genotipos del cromosoma X. De los 48 individuos analizados, se detectaron 2 casos con discordancia entre el sexo reportado en el archivo `.fam` y el inferido genéticamente.  
Los individuos ARI022 (Aymara) y CDSJ176 (Mapuche) fueron identificados como discordantes, con valores de F = 0.9249 y 0.7578 respectivamente, que los ubican fuera de los rangos esperados (F > 0.8 para hombres, F < 0.2 para mujeres). 

Estos resultados se resumen en el archivo `sexcheck_discordant_summary.csv` y corresponden a una proporción de discordancia del 4.1 % del total de muestras.

g) Genera una tabla de contingencia de individuos por sexo y ascendencia (pista: ver columna Ancestry en el archivo `popinfo`

R: partimos creando en el directorio /code un script llamados `sex_ancestry_contingency.R` el cual contiene: 

```
# Cargar datos ----
fam <- read.table("../results/chilean_all48_hg19_updated_ids.fam", 
                  header = FALSE, stringsAsFactors = FALSE)

colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

# Cargar información de ascendencia
popinfo <- read.csv("../data/chilegenomico/chilean_all48_hg19_popinfo.csv", 
                    stringsAsFactors = FALSE)

# Normalizar IDs
fam$IID_clean <- toupper(fam$IID)
popinfo$IndID_clean <- gsub("-", "", toupper(popinfo$IndID))

# Unir por ID
merged <- merge(fam, popinfo, by.x = "IID_clean", by.y = "IndID_clean", all.x = TRUE)

# Verificar columnas relevantes
if(!"Ancestry" %in% names(merged)) stop("No se encontró columna 'Ancestry' en popinfo.")

# Generar tabla de contingencia
tabla_cont <- table(merged$Ancestry, merged$SEX)

# Asignar etiquetas de sexo
colnames(tabla_cont) <- c("Masculino (1)", "Femenino (2)")

# Mostrar en consola
cat("\nTabla de contingencia: Individuos por sexo y ascendencia\n")
print(tabla_cont)

# Guardar tabla
write.csv(as.data.frame.matrix(tabla_cont), 
          "../results/sex_by_ancestry_contingency.csv", row.names = TRUE)
cat("\nTabla guardada en ../results/sex_by_ancestry_contingency.csv\n")
```

ejecutamos: 

    R --vanilla -q -e "source('../code/sex_ancestry_contingency.R')"

R:  la tabla que podemos ver: 

| Ancestry       | Masculino (1) | Femenino (2) |
| -------------- | ------------- |:------------:|
| Aymara         | 6             | 4            |
| FarNorth       | 3             | 3            |
| Mapuche        | 10            | 22           |
| SantiagoPublic | 1             | 0            |

Se construyó una tabla de contingencia cruzando el sexo (columna `SEX` del archivo `.fam`) con la ascendencia (`Ancestry` del archivo `popinfo`).  
La tabla muestra la distribución de individuos masculinos y femeninos dentro de cada grupo poblacional.  Los resultados se encuentran guardados en el archivo `sex_by_ancestry_contingency.csv`.
