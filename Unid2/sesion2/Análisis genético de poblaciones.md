## Análisis genético de poblaciones

Para este ejercicio trabajaremos en el directorio 

 cd ~/rnaranjo/unid2/sesion2

y usaremos los programas: 

    module load plink/1.90
    module load R/4.0.5

 export C=/datos/compartido/ChileGenomico
mkdir -p results

## Primera parte: Análisis de control de calidad

#### paso 1

Investigamos la cantida de genotipos perdidos de los datos de 1000G
usamos el comando:

     plink --bfile $C/chilean_all48_hg19 --missing

Respecto a los resultados podemos reposnder las siguientes preguntas: 

1. ¿Cómo se llaman los archivos que contienen las tasas de datos perdidos por SNP y por muestra?

R: Los archivos que se generan son los plink.imiss y plink.lmiss           
![texto alternativo](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/main/Unid2/sesion2/Imagenes/1oficial.png)

2. ¿Cuántas variantes se eliminaron por tener una tasa de datos perdidos mayor a 0.2?

Para responder esta pregunta usamos el siguiente comando: 

    plink --bfile $C/chilean_all48_hg19 --geno 0.2 --make-bed --out chilean_all48_hg19_2
    plink --bfile chilean_all48_hg19_2 --mind 0.2 --make-bed --out chilean_all48_hg19_3

.
1

Vemos que se eliminaron 4680 variantes de 813366 aplicando el filtro 

3. ¿Cuántos individuos tenían una tasa de datos perdidos mayor a 0.02?

Para responder esta pregunta usamos el siguiente comando 

         plink --bfile chilean_all48_hg19_3 --geno 0.02 --make-bed --out chilean_all48_hg19_4
    plink --bfile chilean_all48_hg19_4 --mind 0.02 --make-bed --out chilean_all48_hg19_5

en base a eso podemos responder que 
**insertar imagen 5

Se eliminaron 46808

4. Basados ​​en los histogramas y en sus cálculos, ¿qué valores umbrales de datos perdidos para muestras y SNPs sugerirían?

El criterio de usar un filtro de 0.2 es apropiado ya que permite eliminar 4680 variantes con perdida de datos, sin perder individuos

.

#### paso 2

El paso 2 tiene como proposito revisar las discrepacias entre la tabla de fenotipos archivos ped/bed y el inferido desde los tipo

Sexo de los individuos inferidos mediante el cromosoma X
Indiv. Femeninos  F<0.2
Indiv. Masculinos F>0.8

1. ¿Cuántos individuos fueron eliminados por discrepancia de sexo?
   Posterior a identificar los individuos con discrepacias  procedemos a eliminarlos. Para eso usamos:
   
       plink --bfile chilean_all48_hg19_5 --remove sex_discrepancy.txt --make-bed --out chilean_all48_hg19_6 

Se eliminaron 3 indiv. que fueron identificados con discrepancias.

imagen 3 indiv

2. ¿Qué riesgo(s) se corre(n) si no se eliminan?

.Los errores asociados es que nuestras inferencias estarian sustentadas en bases de datos mal tabuladas con errores de codificación. Cometiendo errores estadisticos tipo I y II.

#### paso 3

El paso 3 tiene como objetivo realizar un archivo SNPs autosomales solamente.

1. ¿Cuál es el nombre del primer conjunto de datos que solo contiene SNPs en autosomas?

la lista de ID´s es `snp_1_22.txt.` almacenado en` /unid2/sesion2 `y la informacion a `chilean_all48_hg19_7`

2. ¿Cuántos SNP se encontraron en cromosomas sexuales?
   Nuestro dataset presenta 574,624 SNP y quedaron 557,922 autosomicos por lo cua hay 16,702 SNP en cromosomas sexuales (574,624 − 557,922).

imagen paso 3

3. ¿Cómo calcularía el número de cromosomas que porta cada uno de los alelos para cada SNP?
   
         plink --bfile chilean_all48_hg19_7 --freq --out results/MAF_check

#### paso 4

1. ¿Cuál es el nombre del archivo con los resultados de la prueba de HWE?

Para esto usamos 

```c
plink --bfile chilean_all48_hg19_8 --hwe 1e-6 --hwe-all --make-bed --out chilean_all48_hg19_9 #guardamos lo que se filtro
plink --bfile chilean_all48_hg19_8 --hardy  #sin --out 
```

El nombre del archivos es `plink.hwe`

2. ¿Basándose en la distribución de los valores de *p* , le parece el umbral usado razonable o propondría otro valor?

Al ver la desviaciones  de HWE nos hace mantener p < 1e−6 como umbral razonable. Con este corte se eliminaron 1,281 SNP y quedaron 458,097.

#### paso 5

Elimimar padres desconocidos. Para mejorar los analisis de ancestría, las cuales asumen que los SNP no estan correlacionados, es decir que no tienen un fuerte ligamento

```c
 plink --bfile chilean_all48_hg19_9 --exclude $T/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
```

1. ¿Cuántos SNPs en aparente equilibrio de ligamiento se encontraron?

Se eliniaron 346,968 de 450182 usando `--indep-pairwise 50 5 0.2`

2. ¿Cuántos SNP se eliminaron por estar en regiones de inversiones conocidas?

Con el comando `--exclude range` 7,915 SNP en zonas de inversiones 

Ahora para expcluimos todos los individios que tengan un un valor pihat ≥ 0,2.

```
 plink --bfile chilean_all48_hg19_9 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
```

3. ¿Cuántos individuos quedaron luego del filtro de parentesco?

Se eliminaron 3 individos von el filtro de parentezco usando

    plink --bfile chilean_all48_hg19_9 \  --remove to_remove_by_relatedness.txt \ --make-bed \  --out chilean_all48_hg19_10

imagen paso 5(2)

4. ¿Cuál fue el mayor coeficiente de parentesco efectivamente aceptado?

al revisar aquello que superaron el filtro encontramos la pareja `48 ARI008 ARI019 `con un `PI_HAT` 0.2005

======================================================================

## Segunda parte

### unir datos locales con 1000G

Usaremos los archivos `chilean_all48_hg19_10.fam`, `chilean_all48_hg19_10.bim`, `chilean_all48_hg19_10.bed` y de los SNPs con bajo LD `.prune.in`

#### paso 1

Eliminamos los duplicados de chilegenomico 

     plink --bfile chilean_all48_hg19_10 --list-duplicate-vars suppress-first
     plink --bfile chilean_all48_hg19_10 --exclude plink.dupvar --make-bed --out chilean_all48_hg19_11

#### paso 2

### Homologar la versión del genoma

Hacer que los conjuntos de datos tengan la misma version que el conjunto del fenoma para usar las mismas coordenadas de SNP, para asegurarnos, cambiaremos las coordenadas en los datos de 1000 

     awk '{print $2, $4}' chilean_all48_hg19_12.bim> buildhapmap.txt
     plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
     #buildhapmap.txt contiene un ID y una posición física por SNP en cada línea.

`chilean_all48_hg19_12` y `1kG_MDS7` ahora tienen las mismas coordenadas.

#### paso 3

### Fusionar los conjuntos de datos HapMap y 1000 Genomes

    $ awk '{print $2, $5}' 1kG_MDS7.bim> 1kG_ref-list.txt
    $ plink --bfile chilean_all48_hg19_12 --reference-allele 1kG_ref-list.txt --make-bed --out chilean_all48_hg19_13
    $ grep "Impossible" chilean_all48_hg19_13.log | wc -l0
    $ awk '{print $2, $5, $6}' 1kG_MDS7.bim> 1kG_MDS7_tmp
    $ awk '{print $2, $5, $6}' chilean_all48_hg19_13.bim > chilean_all48_hg19_13_tmp
    $ sort 1kG_MDS7_tmp chilean_all48_hg19_13_tmp | uniq -u > all_differences.txt
    $ wc -l all_differences.txt
    0 all_differences.txt

   0 diferencias 

#### paso 4

### Intercambiar alelos de SNPs con potenciales problemas de hebra

Imprimimos el identificador SNP y elimina los duplicados. e intercambiar los alelos de los SNPs no coincide.

    $ awk '{print $ 1}' all_differences.txt | sort -u > flip_list.txt
    $ plink --bfile chilean_all48_hg19_13 --flip flip_list.txt --reference-allele 1kG_ref-list.txt --make-bed --out chilean_all48_hg19_14

 Ver si hay SNP problematicos

        $ awk '{print $2, $5, $6}' chilean_all48_hg19_14.bim > chilean_all48_hg19_14_tmp
     $ sort 1kG_MDS7_tmp chilean_all48_hg19_14_tmp | uniq -u> uncorresponding_SNPs.txt
     $ wc -l uncorresponding_SNPs.txt0 uncorresponding_SNPs.txt

#### paso 5

### Combinar 1000G con ChileGenomico.

    $ plink --bfile 1kG_MDS8 --bmerge chilean_all48_hg19_15.bed chilean_all48_hg19_15.bim chilean_all48_hg19_15.fam --allow-no-sex --make-bed --out MDS_merge

## ==============================================

## Tercera parte

=====================================================================

### Análisis de estructura poblacional

#### paso 1 MDS HapMap-ChileGenomico

    plink --bfile MDS_merge --extract indepSNP.prune.in --genome --out MDS_merge2
    plink --bfile MDS_merge --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

#### paso 2 Generar archivo con info poblaciones.

    #Convertir 
    
    awk '{print$1,$1,$2}' $G/20100804.ALL.panel > ethnicity_1kG.txt
    sed 's/JPT/ASN/g' ethnicity_1kG.txt>ethnicity_1kG2.txt
    sed 's/ASW/AFR/g' ethnicity_1kG2.txt>ethnicity_1kG3.txt
    sed 's/CEU/EUR/g' ethnicity_1kG3.txt>ethnicity_1kG4.txt
    sed 's/CHB/ASN/g' ethnicity_1kG4.txt>ethnicity_1kG5.txt
    sed 's/CHD/ASN/g' ethnicity_1kG5.txt>ethnicity_1kG6.txt
    sed 's/YRI/AFR/g' ethnicity_1kG6.txt>ethnicity_1kG7.txt
    sed 's/LWK/AFR/g' ethnicity_1kG7.txt>ethnicity_1kG8.txt
    sed 's/TSI/EUR/g' ethnicity_1kG8.txt>ethnicity_1kG9.txt
    sed 's/MXL/AMR/g' ethnicity_1kG9.txt>ethnicity_1kG10.txt
    sed 's/GBR/EUR/g' ethnicity_1kG10.txt>ethnicity_1kG11.txt
    sed 's/FIN/EUR/g' ethnicity_1kG11.txt>ethnicity_1kG12.txt
    sed 's/CHS/ASN/g' ethnicity_1kG12.txt>ethnicity_1kG13.txt
    sed 's/PUR/AMR/g' ethnicity_1kG13.txt>ethnicity_1kG14.txt
    
    $ awk '{if($1~/CDSJ/) pop="MAP"}{if($1~/ARI/) pop="AYM"} {print $1, $2, pop}' chilean_all48_hg19_14.fam > ethnicityfile_CLG.txt
    $ cat ethnicity_1kG14.txt ethnicityfile_CLG.txt | sed -e '1i \ FID IID ethnicity'> ethnicityfile.txt

#### paso 3 Graficar resultados de MDS

    Rscript $W/MDS_merged.R

1. En R, genere gráficos similares para las combinaciones Component 2 vs 3 y 3 vs 4. ¿Qué puede concluir de estos gráficos?
   Trabajaremos en nuestro propio scipt de R llamado MDS_C2C3_C3C4.R
   
       cat > MDS_C2C3_C3C4.R <<'EOF'
       args <- commandArgs(trailingOnly = TRUE)
       OUTPUT_DIR <- if (length(args) >= 1) args[1] else "figuras/mds"
       dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
       out <- function(name) file.path(OUTPUT_DIR, name)
       
       # Archivos de entrada (los que ya generaste con plink/R)
       mds_file <- "MDS_merge2.mds"
       eth_file <- "ethnicityfile.txt"
       
       mds <- read.table(mds_file, header = TRUE)
       eth <- read.table(eth_file, header = TRUE)
       
       # Une por IID (o ajusta si tu columna de emparejamiento difiere)
       if("IID" %in% names(mds) && "IID" %in% names(eth)){
         mds <- merge(mds, eth[,c("IID","ethnicity")], by="IID", all.x=TRUE)
       } else {
         stop("No encuentro columna IID en mds o ethnicityfile.txt")
       }
       
       # Paleta simple por grupo
       cols <- c(EUR="blue", ASN="orange", AMR="brown", AFR="green", AYM="gold", MAP="darkgreen")
       pchv <- c(EUR=1, ASN=16, AMR=16, AFR=1, AYM=3, MAP=3)
       
       plot_save <- function(ax="C2", ay="C3", file="plot.pdf"){
         stopifnot(ax %in% names(mds), ay %in% names(mds))
         pdf(file, width=7, height=6)
         on.exit(dev.off(), add=TRUE)
         plot(mds[[ax]], mds[[ay]], type="n",
              xlab=paste("MDS Component", sub("^C","", ax)),
              ylab=paste("MDS Component", sub("^C","", ay)))
         for(g in names(cols)){
           ii <- which(mds$ethnicity == g)
           if(length(ii)){
             points(mds[[ax]][ii], mds[[ay]][ii], col=cols[g], pch=pchv[g])
           }
         }
         legend("topleft", legend=names(cols), col=cols, pch=pchv, bty="n")
       }
       
       plot_save("C2","C3", out("MDS_C2_vs_C3.pdf"))
       plot_save("C3","C4", out("MDS_C3_vs_C4.pdf"))
       EOF

Descargamos desde el servidor usando: tambien se encuentran en /unid2/sesion2

    & "C:\Program Files\PuTTY\pscp.exe" `
      bioinfo1@genoma.med.uchile.cl:/home/bioinfo1/rnaranjo/unid2/sesion2/figuras/mds/MDS_C2_vs_C3.pdf `
      "C:\Users\rnara\Pictures\imagenes unid2.S2\MDS_C2_vs_C3.pdf"
    
    & "C:\Program Files\PuTTY\pscp.exe" `
      bioinfo1@genoma.med.uchile.cl:/home/bioinfo1/rnaranjo/unid2/sesion2/figuras/mds/MDS_C3_vs_C4.pdf `
      "C:\Users\rnara\Pictures\imagenes unid2.S2\MDS_C3_vs_C4.pdf"

De estos podemos decir que los MDS confirman que las referencias continentales están bien definidas, Aymara y Mapuche se ubican sobre un mismo cline amerindio con diferenciación detectable entre ambos, y muestran admixtura principalmente con europeos, sin señales relevantes de componente asiático y con muy poca africana.

#### paso 4

### Realizar un análisis de ascendencia

Usamos: 

      plink --bfile MDS_merge --extract indepSNP.prune.in --make-bed --out MDS_merge_r2_lt_0.2

1. ¿Cuántos SNP quedaron luego del filtro?

Lo cual mediante la linea de comando nos entrega que quedaron:70,534NPs (línea “70534 variants remaining”) .


2. ADMIXTURE asume que los individuos no están emparentados. Sin embargo, no realizamos ningún filtro. ¿Por qué?

.
