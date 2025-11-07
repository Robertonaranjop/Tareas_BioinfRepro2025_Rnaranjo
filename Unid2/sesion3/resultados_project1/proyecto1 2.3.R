###===================================
#Tarea unid2/sesion3/proy1
#=====================================
# Project 1: Do all present-day populations from Europe display the same 3-way admixture?
#Roberto Naranjo-Partarrieu
##================
# Descargar Paquetes
##================

Sys.which("make")  # debe dar: "C:/rtools45/usr/bin/make.exe"
install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)


file.edit("~/.Renviron")

Sys.which("make")
pkgbuild::check_build_tools(debug = TRUE)





remotes::install_github("uqrmaie1/admixtools")


# Debe devolver una ruta válida, p.ej. "C:/rtools45/usr/bin/make.exe"
Sys.which("make")

# Diagnóstico detallado:
install.packages("pkgbuild")  # si no lo tienes
pkgbuild::check_build_tools(debug = TRUE)

cat('RTOOLS45_HOME=C:/rtools45
PATH="${RTOOLS45_HOME}/usr/bin;${PATH}"',
    file = "~/.Renviron", sep = "\n")



if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
pkgs <- c("admixtools","tidyverse","Rcpp","igraph","plotly")
for(p in pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    if(p == "admixtools"){
      devtools::install_github("uqrmaie1/admixtools")
    } else {
      install.packages(p)
    }
  }
}


#### Load packages:
library(admixtools)
library(tidyverse)

###==================== 
# 1) Definir poblaciones
###===================
#Las poblaciones que elegire es Orcadian (represtando el Norte) Estonian y Polish (centro), Greek.DG (meditarráneo), Sardunian.DG(meditiarraneo), Itaulian Sur (meditarraneo central)
#Fuentes Turkey_Marmara_Barcin_N.AG (ANF), Luxembourg_Mesolithic.DG (WHG), Russia_Samara_EBA_Yamnaya.AG (Steppe).
#Outgroup Mbuti.DG, CHB.DG, Papuan.DG, Russia_UstIshim_IUP.DG, Denisova.DG.

##====
BASE_DIR <- "C:/Users/rnara/Desktop/popgen_shared"
prefix   <- file.path(BASE_DIR, "v62.0_1240k_public")          # <-- EIGENSTRAT (SIN extensión)
outdir   <- file.path(BASE_DIR, "aadr_1000G_f2_project1")      # <-- aquí se guardarán los f2
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 2) Comprueba que existen los 3 archivos EIGENSTRAT
file.exists(paste0(prefix, c(".snp",".ind",".geno")))
# Debe devolver: TRUE TRUE TRUE

# 3) (Re)define tus grupos por si la sesión se limpió

target1 <- c("GBR.DG","French.DG","IBS.DG","TSI.DG","Italian_North.DG","Sardinian.DG","Orcadian.DG", "Russian.DG")

source1 <- c("Turkey_Marmara_Barcin_N.AG", "Luxembourg_Mesolithic.DG","Russia_Samara_EBA_Yamnaya.AG")

outgroup1 <- c("Mbuti.DG", "CHB.DG", "Papuan.DG","Russia_UstIshim_IUP.DG", "Denisova.DG")

all_pops <- unique(c(target1, source1, outgroup1))

# 4) Extrae f2 (solo una vez; si ya existen, puedes saltarte esto)
extract_f2(pref = prefix,
           outdir = outdir,
           pops = all_pops,
           overwrite = TRUE,
           blgsize = 0.05,
           verbose = TRUE)

# 5) Carga los bloques f2 recién creados
f2_blocks <- f2_from_precomp(outdir)

#### Outgroup-f3: shared drift entre target y sources
# pop1=outgroup; pop2=target; pop3=sources

f3_results <- f3(f2_blocks, pop1="Mbuti.DG", pop2=target1, pop3=source1)
print(f3_results)

## Todas las poblaciones objetivo muestran valores f3 positivos y estadisticamente
## significativos con las tres poblaciones fuente, confirmando que han contribuido 
## genéticamente a las poblaciones modernas 


#### f4 tests: Chequeos de asimetría
# ¿Están las poblaciones objetivo más cerca de alguna de las fuentes potenciales?

f4_results <- f4(f2_blocks, pop1 = target1, pop2 = "Turkey_Marmara_Barcin_N.AG",pop3 = "Luxembourg_Mesolithic.DG",pop4 = "Mbuti.DG")                 

print(f4_results)

## Según los valores Z y est positivos en todas las poblaciones objetivos podemos
## concluir que las poblaciones europeas modernas están más relacionadas con los 
## cazadores recolectores occidentales que con Anatolia EEF

f4_results_2 <- f4(f2_blocks, 
                   pop1 = target1, 
                   pop2 = "Turkey_Marmara_Barcin_N.AG", # Anatolia EEF
                   pop3 = "Russia_Samara_EBA_Yamnaya.AG", # Steppe
                   pop4 = "Mbuti.DG")
## Resultados f4 Target vs Anatolia/Steppe
print(f4_results_2)


#### qpWave

#####
# --- Getter robusto de nombres de poblaciones en f2_blocks ---
get_pops_safe <- function(f2, outdir = NULL, prefix = NULL) {
  cand <- list(
    tryCatch(f2$pops, error = function(e) NULL),
    tryCatch(dimnames(f2$afreqs)[[2]], error = function(e) NULL),
    tryCatch(attr(f2, "pops"), error = function(e) NULL)
  )
  pops <- unique(unlist(cand))
  if (!is.null(pops) && length(pops) > 0) return(pops)
  
  # Fallback 1: archivo de pops en el outdir (algunas versiones lo guardan)
  if (!is.null(outdir)) {
    f_txt <- file.path(outdir, "pops.txt")
    if (file.exists(f_txt)) return(readLines(f_txt))
  }
  
  # Fallback 2: leer del .ind del EIGENSTRAT y devolver lista completa
  if (!is.null(prefix) && file.exists(paste0(prefix, ".ind"))) {
    ind <- read.table(paste0(prefix, ".ind"), header = FALSE, stringsAsFactors = FALSE)
    return(sort(unique(ind[[3]])))
  }
  
  stop("No pude recuperar los nombres de poblaciones de f2_blocks.")
}

# --- Úsalo con tus objetos ya definidos ---
pops_in_f2 <- get_pops_safe(f2_blocks, outdir = outdir, prefix = prefix)
length(pops_in_f2)
head(pops_in_f2, 20)

# Verifica que todas las que quieres están presentes:
faltantes <- setdiff(c(target1, source1, outgroup1), pops_in_f2)
faltantes


#####
  
  
  wave_results <- lapply(target1, function(targ) {
    print(paste("Corriendo qpWave para:", targ))
    qpwave(f2_blocks,
           left = c(targ, source1),
           right = outgroup1)
  })
  names(wave_results) <- target1
  
  print(wave_results)
  

#### qpAdm: modelos de mezcla de 2 y 3 vías

##=================================================
#Correr qpAdm (2 vías y 3 vías) por cada target
##=================================================
  
  message("== qpAdm 2-vías (ANF + WHG) ==")
  admix_2way_ANF_WHG <- lapply(target1, function(tg){
    qpadm(f2_blocks,
          left   = c(tg, "Turkey_Marmara_Barcin_N.AG", "Luxembourg_Mesolithic.DG"),
          right  = outgroup1,
          target = tg)
  })
  names(admix_2way_ANF_WHG) <- target1
  
  message("== qpAdm 2-vías (ANF + Steppe) ==")
  admix_2way_ANF_Steppe <- lapply(target1, function(tg){
    qpadm(f2_blocks,
          left   = c(tg, "Turkey_Marmara_Barcin_N.AG", "Russia_Samara_EBA_Yamnaya.AG"),
          right  = outgroup1,
          target = tg)
  })
  names(admix_2way_ANF_Steppe) <- target1
  
  message("== qpAdm 3-vías (ANF + WHG + Steppe) ==")
  admix_3way <- lapply(target1, function(tg){
    qpadm(f2_blocks,
          left   = c(tg, source1),   # c(tg, ANF, WHG, Steppe)
          right  = outgroup1,
          target = tg)
  })
  names(admix_3way) <- target1
  
##==============================
## Extrae y ordena los  weigths 
##==============================
  library(dplyr); library(readr); library(purrr); library(stringr)
  
  extract_weights <- function(lst, model_label){
    bind_rows(lapply(names(lst), function(tg){
      x <- lst[[tg]]
      if (is.null(x) || is.null(x$weights)) return(NULL)
      x$weights %>%
        mutate(target = tg, model = model_label) %>%
        select(model, target, left, weight, se, z)
    }))
  }
  
  tab_weights <- bind_rows(
    extract_weights(admix_2way_ANF_WHG,   "2way_ANF+WHG"),
    extract_weights(admix_2way_ANF_Steppe,"2way_ANF+Steppe"),
    extract_weights(admix_3way,           "3way_ANF+WHG+Steppe")
  )
  
  write_csv(tab_weights, file.path(RES_DIR, "qpadm_weights_all.csv"))
  tab_weights %>% arrange(model, target, left) %>% print(n = 50)

###===========================
#Marcar modelos posibles
##============================
  
# Define (o re-define) rutas base y carpeta de resultados
  BASE_DIR <- "C:/Users/rnara/Desktop/popgen_shared"   # ajusta si tu base es otra
  RES_DIR  <- file.path(BASE_DIR, "resultados_project1")
  if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
  
  

  check_models <- tab_weights %>%
    group_by(model, target) %>%
    summarize(sum_w = sum(weight, na.rm = TRUE),
              any_neg = any(weight < -1e-6, na.rm = TRUE),
              any_over = any(weight > 1 + 1e-6, na.rm = TRUE),
              min_z = min(abs(z), na.rm = TRUE),
              .groups = "drop") %>%
    mutate(plausible = (!any_neg & !any_over & between(sum_w, 0.95, 1.05)))
  
  write_csv(check_models, file.path(RES_DIR, "qpadm_model_checks.csv"))
  check_models
  readr::write_csv(tab_weights, file.path(RES_DIR, "qpadm_weights_all.csv"))  
tab_weights

###===================
# Graficar
###==================

plot_df <- tab_weights %>%
  filter(model == "3way_ANF+WHG+Steppe") %>%
  mutate(Comp = recode(left,
                       "Turkey_Marmara_Barcin_N.AG"   = "ANF",
                       "Luxembourg_Mesolithic.DG"     = "WHG",
                       "Russia_Samara_EBA_Yamnaya.AG" = "Steppe"))

library(ggplot2)
p <- ggplot(plot_df, aes(x = target, y = weight, fill = Comp)) +
  geom_col() +
  geom_errorbar(aes(ymin = pmax(weight - 1.96*se, 0),
                    ymax = pmin(weight + 1.96*se, 1)),
                width = 0.25) +
  coord_flip() + theme_bw() +http://127.0.0.1:32797/graphics/d477c537-3ce5-4fe5-a1d9-0356410f2909.png
  labs(x = "", y = "Proporción (qpAdm)", fill = "Componente",
       title = "qpAdm 3-vías: ANF / WHG / Steppe")
ggsave(file.path(RES_DIR, "qpadm_3ways_barplot.png"), p, width = 8, height = 5, dpi = 300)
plot("qpadm_3ways_barplot.png")


##tabla

# Tabla solo del 3-vías, redondeada y ordenada
tab_3way <- tab_weights %>%
  dplyr::filter(model == "3way_ANF+WHG+Steppe") %>%
  dplyr::mutate(Componente = dplyr::recode(left,
                                           "Turkey_Marmara_Barcin_N.AG"   = "ANF",
                                           "Luxembourg_Mesolithic.DG"     = "WHG",
                                           "Russia_Samara_EBA_Yamnaya.AG" = "Steppe")) %>%
  dplyr::select(target, Componente, weight, se, z) %>%
  dplyr::arrange(target, dplyr::desc(weight))

readr::write_csv(tab_3way, file.path(RES_DIR, "qpadm_3ways_only.csv"))
print(tab_3way, n = 50)

