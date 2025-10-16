# Unidad 1 – Sesión 3.

---

### Ejercicio 1 Variable con `log10(50)` y suma con 5

**Objetivo**: Crear una variable con `log10(50)` y sumarla a otra variable con valor 5.

```r
# Calculamos el logaritmo en base 10 de 50
x <- log10(50)

# Creamos otra variable con valor 5
y <- 5

# Sumamos ambas variables
resultado <- x + y
```

----

### Ejercicio 2  Sumar 2 a todos los números entre 1 y 150

**Objetivo**: Sumar el número 2 a todos los números entre 1 y 150

```r
# Creamos una secuencia de 1 a 150 con el operador ":" (crea enteros consecutivos).
v <- 1:150

# Vectorizamos la suma: en R, sumar un escalar a un vector aplica elemento a elemento.
v_mas_2 <- v + 2

# Revisamos las primeras entradas para confirmar.
head(v_mas_2)
# También podríamos verificar las últimas entradas:
tail(v_mas_2)
```

----

### Ejercicio 3 Conteo condicional

**Objetivo**: ¿Cuántos números son mayores a 20 en el vector `-13432:234`?

```r
# Construimos el vector con un rango grande de enteros.
rango <- -13432:234

# Creamos un vector lógico que indica TRUE cuando el valor es > 20.
mayores_20_log <- rango > 20

# Sumamos el vector lógico: en R, TRUE vale 1 y FALSE vale 0, así contamos los TRUE.
conteo_mayores_20 <- sum(mayores_20_log)

conteo_mayores_20
```

---

### Ejercicio 4  Cargar archivo

**Objetivo**: Cargar en R `PracUni1Ses3/maices/meta/maizteocintle_SNP50k_meta_extended.txt` en un objeto llamado `meta_maiz`.

```r
# Establecemos el directorio donde se encuentra nuestro archivo "maizteocintle_SNP50k_meta_extended.txt"
setwd ("C:/Users/MSI/Desktop/Repositorio/BioinfinvRepro/Unidad1/Sesion3/PracUni1Ses3/maices/meta")
# Luego creamos el objeto meta_maiz
meta_maiz <- "maizteocintle_SNP50k_meta_extended.txt"

#En caso de querer cargar con los datos usamos: 
meta_maiz <- read.delim("maizteocintle_SNP50k_meta_extended.txt", 
                       fileEncoding = "macintosh")
```

----

### Ejercicio 5  For loops con divisiones y almacenamiento de resultados

#### 5a) Bucle incial

- Escribe un bucle for para que divida 35 entre 1:10 e imprima el resultado en la consola.Imprimir `35 / (1:10)`

```r
# Loop que divide 35 por cada númeo del 1 al 10
for (i in 1:10) {
  cat("35 /", i, "=", 35 / i, "\n")
}
```

#### 5b) Bucle condicionado

- Solo para números **impares** usando `next` (sin enumerarlos a mano)

```r
# Usamos 'next' para saltar los 'i' que no cumplan la condición.
for (i in 1:10) {
  # Si 'i' es par, saltamos esta iteración (no calculamos ni imprimimos).
  if (i %% 2 == 0) next
  cat("[impar] 35 /", i, "=", 35 / i, "\n")
}
```

#### 5c) Almacenar

- Guardar resultados en un **data frame** con dos columnas

- **Columna 1**: texto `"resultado para x"`  

- **Columna 2**: el valor numérico de la división

```r
# Inicializamos un data frame vacío fuera del loop (buena práctica).
# Usamos character(0) y numeric(0) para tipar las columnas.
res_df <- data.frame(
  etiqueta = character(0),
  valor    = numeric(0),
  stringsAsFactors = FALSE
)

# Iteramos de 1 a 10, manteniendo solo impares (como en 5b).
for (i in 1:10) {
  if (i %% 2 == 0) next
  etiqueta_i <- paste("resultado para", i) # construye el texto solicitado
  valor_i    <- 35 / i                     # cálculo
  # Agregamos una fila con rbind. Para grandes datos, preasignar es más eficiente.
  res_df <- rbind(res_df, data.frame(etiqueta = etiqueta_i, valor = valor_i))
}

# Vemos la tabla final.
res_df
```

----

### Ejercicio 6 Análisis script "IBR_testing.r"

Primero abrimos el script para poder analizarlo y responder las siguientes preguntas

```r
#Establecer repositorio de trabajo
setwd ("C:/Users/MSI/Desktop/Repositorio/BioinfinvRepro/Unidad1/Sesion3/PracUni1Ses3/mantel/bin")
```

###### ¿Qué hacen los dos `for` loops?

1. **Primer `for` (sobre `i` en `"present", "ccsm", "miroc", "flat", 1800–4000"`)**
   
   - **Lee** para cada condición el archivo de **resistencias/efective distances** producido por Circuitscape: `Balpina_<i>_resistances.out`.
   
   - **Convierte** esa matriz al orden de poblaciones correcto (`popNamesFP` → `popNames`) con `read.effdist(...)`.
   
   - **Calcula** la **distancia efectiva media por población**: `apply(eff.dist, 2, mean)`.
   
   - **Guarda** resultados en el entorno con nombres dinámicos:
     
     - `B.<i>` → matriz de distancias efectivas para la condición `i`.
     
     - `B.mean.<i>` → vector con la media de distancias efectivas por población para `i`.

2. **Segundo `for` (mismo conjunto de condiciones)**
   
   - **Linealiza** antes la Fst: `B.FstLin <- B.Fst/(1 - B.Fst)`.
   
   - Para cada `i`, **ejecuta un Mantel test** entre `as.dist(B.<i>)` y `as.dist(B.FstLin)` con **10 000 permutaciones** (`mantel.rtest(..., nrepet=10000)`).
   
   - **Grafica** la relación Distancia efectiva vs. Fst linealizado con `DistPlot(...)`.
   
   - **Extrae** y **acumula** en `IBRresults` el **p-value** y el **coeficiente r** del Mantel para cada superficie.

###### ¿Qué paquetes necesitas para correr el script?

- **ade4** → provee `mantel.rtest` (Mantel test con permutaciones).

- **ggplot2** → utilizado dentro de `DistPlot.R` para las gráficas (el script lo carga).

- **sp** → para **distancias geográficas** con `spDists` y manejo básico de coords.

###### ¿Qué archivos necesitas para correr el script?

- `read.fst_summary_fix.R` — función casera para leer/armar la matriz **Fst**.

- `read.effdist.R` — función casera para leer matrices de **distancia efectiva** (salidas de Circuitscape).

- `DistPlot.R` — función para **graficar** dista

----

### Ejercicio 7 Calculo de thetha con `calc.tetha`

**Objetivo**: Escribir una función `calc.tetha(Ne, u)` que calcule θ usando `θ = 4 * Ne * u`.

```r
# calc.tetha: calcula theta (θ) para poblaciones diploides
# θ = 4 * Ne * u

calc.tetha <- function(Ne, u) {
  stopifnot(is.numeric(Ne), is.numeric(u))
  4 * Ne * u
}

# Ejemplos de uso
calc.tetha(Ne = 10000, u = 0.001)
```

----

### Ejercicio 8 Prueba parcial de mantel usando `vegan`

**Objetivo**: Agregar al script del ejercicio de Mantel el **Partial Mantel test** entre la matriz Fst y las matrices del presente y LGM, **parcializando** por la matriz `flat`. (Requiere paquete `vegan`.)

Los resultados del ejercicio se encuentran alojados en el archivo `Robertomantel` dentro de este mismo repositorio 

ademas se detalla en el siguiente cuadro el codigo incorporado  

```r
# ==== AGREGADO DESDE AQUÍ: Mantel PARCIAL controlando 'flat' (vegan) ====
library(vegan)
# Construimos objetos 'dist' para vegan::mantel.partial
fst_dist     <- as.dist(B.FstLin)
present_dist <- as.dist(B.present)
ccsm_dist    <- as.dist(B.ccsm)
miroc_dist   <- as.dist(B.miroc)
flat_dist    <- as.dist(B.flat)

set.seed(123)
mp_present <- mantel.partial(fst_dist, present_dist, flat_dist,
                             method = "pearson", permutations = 9999)
mp_ccsm    <- mantel.partial(fst_dist, ccsm_dist,    flat_dist,
                             method = "pearson", permutations = 9999)
mp_miroc   <- mantel.partial(fst_dist, miroc_dist,   flat_dist,
                             method = "pearson", permutations = 9999)

cat("\n== Mantel Parcial: Fst ~ PRESENTE | FLAT ==\n");  print(mp_present)
cat("\n== Mantel Parcial: Fst ~ LGM (CCSM) | FLAT ==\n"); print(mp_ccsm)
cat("\n== Mantel Parcial: Fst ~ LGM (MIROC) | FLAT ==\n");print(mp_miroc)

# Guardar resultados a CSV
pm_results <- data.frame(
  modelo = c("present|flat", "ccsm|flat", "miroc|flat"),
  r      = c(unname(mp_present$statistic),
             unname(mp_ccsm$statistic),
             unname(mp_miroc$statistic)),
  p      = c(mp_present$signif, mp_ccsm$signif, mp_miroc$signif)
)
dir.create("../results", showWarnings = FALSE)
write.csv(pm_results, "../results/partial_mantel_results.csv", row.names = FALSE)

# ==== FIN AGREGADO ====
```

----

### Ejercicio 9 Crear Script `ExplorandoMaiz.R` almacenado en `PracUni1Ses3/maices/bin`

1) Que cargue en R el archivo 

`PPracUni1Ses3maices/meta/maizteocintle_SNP50k_meta_extended.txty` 

2) Que responda lo siguiente: 
- ¿Qué tipo de objeto creamos al cargar la base?

- ¿Cómo se ven las primeras 6 líneas del archivo?

- ¿Cuántas muestras hay?

- ¿De cuántos estados se tienen muestras?

- ¿Cuántas muestras fueron recolectadas antes de 1980?

- ¿Cuántas muestras hay de cada raza?

- En promedio ¿a qué altitud fueron recolectadas las muestras?

- ¿Y a qué altitud máxima y mínima fueron recogidas?

- Crea una nueva df de datos sólo con las muestras de la raza Olotillo

- Crea una nueva df de datos sólo con las muestras de la raza Reventador, Jala y Ancho

- Escribe la matriz anterior en un archivo llamado "submat.cvs" en /meta.

**LOS RESULTADOS SE ENTREGAN EN EL ARCHIVO EXPLORANDOMAIZ.R**


