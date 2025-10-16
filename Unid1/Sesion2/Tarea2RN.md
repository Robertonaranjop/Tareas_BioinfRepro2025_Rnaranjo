## Tarea 2

## 1. Crear repositorio propio ["Tareas_BioinfRepro2025_Rnaranjo"](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo.git)

Se creo una cuenta en la plataforma Github, se establecio un nombre de usuario llamado "robertonaranjop", al  cual se le asocio un correo "robertonaranjopartarrieu@gmail.com"

#### Configuracion del Git bash

Para configurar el git.bash

```bash
$ git config --global user.email "r.naranjopartarrieu@gmail.com""  

$ git config --global user.name "robertonaranjop"  
```

y para visulizar que se hayan realizado los cambios se uso el siguiente codigo: 

```bash
$ git config --global user.email   
```

## 2. Clonar repositorio de la clase

Para poder copiar el repositorio solo del año 2021, utilizo el codigo Git clone con la opcion branch, single branch 

```bash
$ git clone --branch 2021-II --single-branch https://github.com/u-genoma/BioinfinvRepro.git
```

Ahora si quiero realizar la descarga de una carpeta en especifico de un single branch.

```bash
$ git clone --branch 2021-II --single-branch --filter=blob:none https://github.com/u-genoma/BioinfinvRepro.git
```

Para activar la funcion que nos permite disecsionar todo el repositorio usamos  sparse-checkout y para especificar el directorio usamos --cone

```bash
$ git sparse-checkout init --cone
```

  Ahora seteamos la carpeta a clonar (Unidad 2, en el caso de la version 2021 II)

```bash
$ git sparse-checkout set Unidad2
```

![..](https://github.com/Robertonaranjop/Tareas_BioinfRepro2025_Rnaranjo/blob/b45e19e0f71abb4243f5a07fbb59d5bdb6b2f149/script2.png)

## 3. Analizar Script pipeline

Respecto al siguiente codigo tenemos como mision responder las siguientes preguntas: 

1. ¿Cuántos pasos tiene este guión?
   
   R: Se presentan 6 pasos realizados: 
   
       1) Alineamiento con `gnsap` y conversion a `bam` con `samtools` y elimina el `.sam`
   
       2) Ejecución de `pstacks`
   
       3) Creación de la lista de muestras para `cstacks.src=$HOME/research/project
   
       4) Construccion del catálogo con `cstacks`.
   
       5) Ejeccución de `sstacks`
   
       6) Cálculo de la estadísticasdela poblaciones con `populations`

2. ¿Si quisieras correr este script y que funcionara en tu propio equipo, qué línea deberías cambiar ya qué?
   
   R: Deberia cambiar la fuente de origen  pregunta 3
   `src=$HOME/research/project`

3. ¿A qué equivale `$HOME`?
   
   R: Es la variable de entorno del directorio home del usuario actual.

4. ¿Qué paso del análisis hace el programa `gsnap`?
   
   R: `gsnap` alinea lecturas (FASTQ) contra una referencia (en este caso con gac_gen_broads1_e64). Produce un `SAM` con los alineamientos; luego se convierte a `BAM` con `samtools`. Ademas en el proceso filtra y ajusta el alineamiento con parametros de calidad detallados en el codigo `-n`, `-m`, `-i`, `--min-coverage`

5. ¿Qué hace en términos generales cada uno de los loops?
   
   R:Loop 1 Alineamiento con GNSAP y conversión a BAM 
   
    Recorre todas las muestras listadas en files.
    Para cada muestra: 
   
   1. Alinea con `GSNAP` → genera `.sam.`
   2) Convierte .sam a .bam (más compacto).
   3) Elimina el .sam para ahorrar espacio
   
   ```bash
   for file in $files
   do
    gsnap ... $src/samples/${file}.fq > $src/aligned/${file}.sam
    samtools view -b -S -o $src/aligned/${file}.bam $src/aligned/${file}.sam 
    rm $src/aligned/${file}.sam 
   done
   ```
   
   Loop 2  Ejecutar `pstacks`
   
   Itera sobre todas las muestras.
   Corre pstacks para identificar loci y SNPs por individuo.
   Asigna un ID único incremental (i) a cada muestra.
   
   ```bash
   i=1 
   for file in $files 
   do 
    pstacks ... -i $i -f $src/aligned/${file}.bam -o $src/stacks/ 
    let "i+=1"; 
   done
   ```
   
    Loop 3 Construcción de lista para `cstacks`
   
    Recorre todas las muestras.
    Va armando una cadena de texto con todos los argumentos -s que cstacks necesita.
   
   ```bash
    for file in $files 
    do 
    sstacks ... -c $src/stacks/batch_1 -s $src/stacks/${file} -o $src/stacks/ &>> $src/stacks/Log 
    done
   ```
   
      Loop 4 Ejecutar `sstack`
   
    Itera sobre todas las muestras.
    Corre sstacks para comparar cada muestra contra el catálogo generado por     `cstacks`.
    Va guardando resultados en stacks/ y registrando el log.

```bash
for file in $files 
do 
    sstacks ... -c $src/stacks/batch_1 -s $src/stacks/${file} -o $src/stacks/ &>> $src/stacks/Log 
done
```

```bash
#!/bin/bash 

src=$HOME/research/project 

files=”sample_01 
sample_02 
sample_03” 

#
# Align with GSnap and convert to BAM
#Corre gsnap para cada muestra, convierte SAM → BAM, borra el SAM.
# 
for file in $files
do
    gsnap -t 36 -n 1 -m 5 -i 2 --min-coverage=0.90 \
            -A sam -d gac_gen_broads1_e64 \
            -D ~/research/gsnap/gac_gen_broads1_e64 \
            $src/samples/${file}.fq > $src/aligned/${file}.sam
    samtools view -b -S -o $src/aligned/${file}.bam $src/aligned/${file}.sam 
    rm $src/aligned/${file}.sam 
done

#
# Run Stacks on the gsnap data; the i variable will be our ID for each sample we process.
# Corre pstacks por muestra, con IDs incrementales. 
#
i=1 
for file in $files 
do 
    pstacks -p 36 -t bam -m 3 -i $i \
              -f $src/aligned/${file}.bam \
              -o $src/stacks/ 
    let "i+=1"; 
done 

# 
# Use a loop to create a list of files to supply to cstacks.
# Construye una cadena con -s $src/stacks/$file para cada muestra. Este bloque no corre análisis aún, solo prepara argumentos.
# 
samp="" 
for file in $files 
do 
    samp+="-s $src/stacks/$file "; 
done 

# 
# Build the catalog; the "&>>" will capture all output and append it to the Log file.
# Construye el catálogo con todos los loci y guarda en $src/stacks. 
#
cstacks -g -p 36 -b 1 -n 1 -o $src/stacks $samp &>> $src/stacks/Log 
#
#Compara cada muestra contra el catálogo.
#
for file in $files 
do 
    sstacks -g -p 36 -b 1 -c $src/stacks/batch_1 \
             -s $src/stacks/${file} \ 
             -o $src/stacks/ &>> $src/stacks/Log 
done 

#
# Calculate population statistics and export several output files.
# Genera estadísticas de poblaciones y exportes 
#
populations -t 36 -b 1 -P $src/stacks/ -M $src/popmap \
              -p 9 -f p_value -k -r 0.75 -s --structure --phylip --genepop
```

