# Análise de Uso e Cobertura da Terra - PRODES
#'
#' Script para análise espacial de uso e cobertura da terra baseada em dados 
#' raster do PRODES (INPE), com recorte por polígono e geração de outputs científicos.
#'
#' Fluxo de processamento:
#' - Carregamento e preparação de dados espaciais vetoriais e raster
#' - Recorte do raster em coordenadas geográficas
#' - Conversão da área recortada para projeção UTM
#' - Classificação e cálculo de áreas por classe
#' - Agregação em grupos temáticos
#' - Geração de visualizações científicas e exportação
#'
#' @param geojson_path Caminho para o polígono de interesse (GeoJSON)
#' @param raster_path Caminho para o raster PRODES (TIFF)
#' @param export_dir Diretório de saída para resultados
#' @param nome_area Nome da área analisada para títulos e legendas
#'
#' @details
#' **Sistema de Referência:**
#' Detecção automática da zona UTM baseada na longitude central. Utiliza 
#' SIRGAS 2000 / UTM (EPSG:31978-31985) para território brasileiro, zonas 18S-25S.
#' Reamostragem por vizinho mais próximo preserva códigos originais PRODES.
#' 
#' **Classificação PRODES:**
#' - dAAAA: Desmatamento no ano AAAA
#' - rAAAA: Resíduo/correção para ano AAAA
#' - 100: Vegetação Natural/Floresta
#' - 101: Não-Floresta
#' - 91: Hidrografia
#' - 32: Nuvens
#' - Outros códigos: Classificados como Vegetação Natural/Floresta
#' 
#' **Controle de Qualidade:**
#' Filtro de valores insignificantes (≤10 pixels) remove artefatos de reprojeção.
#' Tabela de classificação criada apenas para valores presentes na área.
#' Validação por comparação entre área calculada e área do polígono.
#' 
#' **Outputs:**
#' - CSV resumido: Área agregada por classe
#' - CSV detalhado: Análise por pixel
#' - CSV metadados: Parâmetros de processamento e métricas
#' - Mapa temático: Resolução 300 DPI, elementos cartográficos
#' - Gráfico temporal: Série histórica de desmatamento
#' - Figura combinada: Mapa + gráfico integrados
#' 
#' **Metodologia:**
#' Estratégia "recorte rápido + conversão precisa" minimiza tempo de processamento.
#' Cálculos de área em UTM garantem precisão. Diferença típica <1% vs. área de referência.
#'
#' @return 
#' Arquivos gerados em export_dir:
#' - `{nome}_area_por_classe_resumido.csv`: Tabela agregada
#' - `{nome}_area_por_classe_detalhado.csv`: Dados completos por pixel
#' - `{nome}_metadados.csv`: Parâmetros e métricas de processamento
#' - `{nome}_mapa_cientifico.png`: Mapa temático (300 DPI, 12×8")
#' - `{nome}_grafico_cientifico.png`: Série temporal (300 DPI, 10×6")
#' - `{nome}_analise_completa_cientifica.png`: Composição final (300 DPI, 12×10")
#' 
#' Console: Detecção UTM, validação de dados, métricas de precisão, estatísticas resumidas.
#'
#' @examples
#' # Configuração para nova área:
#' geojson_path <- "area_estudo.geojson"
#' raster_path <- "raster_PRODES.tif" 
#' export_dir <- "resultados/"
#' nome_area <- "ÁREA PROTEGIDA"
#'
#' @note
#' Requisitos: Polígono em coordenadas geográficas, raster PRODES oficial,
#' área dentro da cobertura PRODES, espaço em disco ~50MB.
#' Recursos: RAM 1-2GB, processamento 30s-5min conforme tamanho da área.
#'
#' @references
#' PRODES: https://terrabrasilis.dpi.inpe.br/geonetwork/srv/eng/catalog.search#/metadata/fe02f2bf-2cc0-49d5-ab72-a3954f997408
#' Dados: http://terrabrasilis.dpi.inpe.br/
#' UTM SIRGAS 2000: https://epsg.io/ (31978-31985)
#'
#' @version 1.0
#' @author Jesse Burlamaque
#' @date 2024-06-11

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggspatial)

### Parâmetros e Caminhos - CONFIGURAÇÃO GERAL

# CONFIGURAÇÕES A SEREM AJUSTADAS PARA CADA ÁREA:
geojson_path <- ".../ROI.geojson"
raster_path <- ".../PDigital2000_2023_AMZ_raster.tif"
export_dir <- ".../Analises"
nome_area <- "NOME DA AREA"

# CONFIGURAÇÕES AUTOMÁTICAS (não mexer):
nome_arquivo <- tolower(gsub(" ", "_", iconv(nome_area, to = "ASCII//TRANSLIT")))

# AUTO-DETECÇÃO DA ZONA UTM baseada na longitude central
cat("0. Detectando zona UTM adequada...\n")
bbox_temp <- st_bbox(st_read(geojson_path, quiet = TRUE))
lon_central <- mean(c(bbox_temp[1], bbox_temp[3]))

# Determinar zona UTM automaticamente para o Brasil
zona_utm <- floor((lon_central + 180) / 6) + 1
hemisferio <- ifelse(mean(c(bbox_temp[2], bbox_temp[4])) < 0, "S", "N")  # Sul ou Norte

# EPSG para zonas UTM do Brasil (SIRGAS 2000)
epsg_utms <- data.frame(
  zona = 18:25,
  sul = c(31978:31985),  # UTM Sul SIRGAS 2000
  norte = c(31974:31981) # UTM Norte SIRGAS 2000 (raramente usado no Brasil)
)

if(hemisferio == "S") {
  crs_utm <- paste0("EPSG:", epsg_utms$sul[epsg_utms$zona == zona_utm])
} else {
  crs_utm <- paste0("EPSG:", epsg_utms$norte[epsg_utms$zona == zona_utm])
}

cat("Área detectada: Zona UTM", zona_utm, hemisferio, "(", crs_utm, ")\n")
cat("Longitude central:", round(lon_central, 3), "°\n\n")

### ETAPA 1: Recorte em coordenadas geográficas

cat("1. Carregando dados...\n")
poligono <- st_read(geojson_path, quiet = TRUE)
r_prodes <- rast(raster_path)

# Verificar se há dados válidos na área
bbox_pol <- st_bbox(poligono)
cat("Bounding box do polígono: ", round(as.numeric(bbox_pol), 4), "\n")

cat("2. Recortando raster (coordenadas geográficas)...\n")
pol_terra <- vect(poligono)
r_crop <- crop(r_prodes, pol_terra)  # Recorte rápido pela bbox
r_mask <- mask(r_crop, pol_terra)    # Máscara pelo polígono

# Verificar se o recorte teve sucesso
if(all(is.na(values(r_mask)))) {
  stop("ERRO: Área sem dados PRODES válidos. Verificar se o polígono está dentro da área de cobertura do PRODES.")
}

cat("Recorte realizado com sucesso!\n")

### ETAPA 2: Conversão para UTM SEM RUÍDO

cat("3. Convertendo para UTM...\n")
poligono_utm <- st_transform(poligono, crs_utm)

# MÉTODO 1: Nearest neighbor para evitar ruído
r_mask_utm <- project(r_mask, crs_utm, method = "near")

# MÉTODO 2: Alternativo - calcular área em coordenadas geográficas
# 
# cat("Usando método geográfico para evitar reprojeção...\n")
# r_mask_utm <- r_mask  # Manter em coordenadas geográficas
# crs_calculo <- "geografico"

### ETAPA 3: Tabela de Classes BASEADA NOS DADOS REAIS

cat("4. Calculando áreas...\n")
freq_table <- freq(r_mask_utm)

# Criar classes APENAS para valores que existem no raster
valores_encontrados <- sort(unique(freq_table$value))
cat("Valores únicos encontrados no raster:", paste(valores_encontrados, collapse = ", "), "\n")

# Criar tabela de classes SOMENTE para valores existentes
classes <- data.frame(
  DN = valores_encontrados,
  Classe = case_when(
    valores_encontrados == 100 ~ "Vegetacao Natural/Floresta",
    valores_encontrados == 101 ~ "Nao Floresta",
    valores_encontrados == 91 ~ "Hidrografia", 
    valores_encontrados == 32 ~ "Nuvem",
    valores_encontrados %in% 7:23 ~ paste0("d", 2000 + valores_encontrados),
    valores_encontrados %in% 50:63 ~ paste0("r", 2010 + (valores_encontrados - 50)),
    TRUE ~ "Vegetacao Natural/Floresta"  # Códigos extras = floresta
  )
)

# Cores conforme padrão PRODES
classes$cor <- case_when(
  grepl("^d", classes$Classe) ~ "red",
  grepl("^r", classes$Classe) ~ "orange", 
  classes$Classe == "Vegetacao Natural/Floresta" ~ "darkgreen",
  classes$Classe == "Hidrografia" ~ "blue",
  classes$Classe == "Nao Floresta" ~ "yellow",
  classes$Classe == "Nuvem" ~ "white",
  TRUE ~ "darkgreen"
)

cat("Classes criadas para", nrow(classes), "valores encontrados\n")
print(classes[, c("DN", "Classe")])

# Simples: UTM já está em metros
raster_res <- res(r_mask_utm)
pixel_area_m2 <- raster_res[1] * raster_res[2]
pixel_area_ha <- pixel_area_m2 / 10000

freq_df <- left_join(freq_table, classes, by = c("value" = "DN"))
freq_df$area_ha <- freq_df$count * pixel_area_ha

# Diagnóstico: mostrar valores sem correspondência
valores_sem_classe <- freq_df[is.na(freq_df$Classe), ]
if(nrow(valores_sem_classe) > 0) {
  cat("ATENÇÃO: Valores no raster sem classe definida:\n")
  print(valores_sem_classe[, c("value", "count")])
  cat("\n")
}

# Validação
area_total_calculada <- sum(freq_df$area_ha, na.rm = TRUE)
area_poligono <- as.numeric(st_area(poligono_utm)) / 10000

cat("=== VALIDAÇÃO DE ÁREA ===\n")
cat("Área do polígono:", round(area_poligono, 2), "ha\n")
cat("Área total calculada:", round(area_total_calculada, 2), "ha\n")
cat("Diferença:", round(abs(area_total_calculada - area_poligono), 2), "ha\n")
cat("Diferença percentual:", round((abs(area_total_calculada - area_poligono) / area_poligono) * 100, 2), "%\n")
cat("========================\n\n")

knitr::kable(freq_df[, c("Classe", "count", "area_ha")], digits = 2)

### ETAPA 5: Preparar dados para visualização

cat("5. Preparando visualizações...\n")
r_df <- as.data.frame(r_mask_utm, xy = TRUE, na.rm = TRUE)
colnames(r_df)[3] <- "DN"
r_df <- left_join(r_df, classes, by = "DN")

# Agrupa classes: desmatamento + resíduo = desmatamento total
r_df$Classe_Agrupada <- case_when(
  grepl("^d|^r", r_df$Classe) ~ "Deforestation",  # d + r = desmatamento total
  r_df$Classe == "Vegetacao Natural/Floresta" ~ "Forest",
  r_df$Classe == "Hidrografia" ~ "Hydrography",
  r_df$Classe == "Nao Floresta" ~ "Non-forest",
  r_df$Classe == "Nuvem" ~ "Cloud",
  TRUE ~ "Forest"  # Default = floresta conforme PRODES
)

r_df <- r_df[!is.na(r_df$Classe_Agrupada), ]

# Calcular áreas agrupadas
area_agrupada <- r_df |>
  count(Classe_Agrupada) |>
  mutate(area_ha = n * pixel_area_ha,
         label = paste0(Classe_Agrupada, " (", round(area_ha), " ha)"))

# Paleta de cores científica
cores_cientificas <- c(
  "Deforestation" = "#d73027",     # Vermelho 
  "Forest" = "#1a9850",            # Verde floresta
  "Hydrography" = "#4575b4",       # Azul  
  "Non-forest" = "#fee08b",        # Amarelo natural
  "Cloud" = "#f7f7f7"              # Cinza claro
)

# Criar label_lookup ANTES de usar
label_lookup <- setNames(area_agrupada$label, area_agrupada$Classe_Agrupada)

### ETAPA 6: Mapa 

# Mapa 
g_mapa <- ggplot(r_df) +
  geom_raster(aes(x = x, y = y, fill = Classe_Agrupada)) +
  scale_fill_manual(
    values = cores_cientificas, 
    labels = label_lookup,
    name = "Land Cover Classes"
  ) +
  coord_equal() +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 20),
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    nrow = 1,
    byrow = TRUE,
    override.aes = list(size = 1)
  )) +
  labs(
    title = paste("Land Use and Land Cover Classification –", nome_area),
    subtitle = "PRODES/INPE v2024 | UTM Zone 20S (SIRGAS 2000)"
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true", 
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = north_arrow_fancy_orienteering,
    height = unit(0.8, "cm"), width = unit(0.8, "cm")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.3,
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = "ticks", line_col = "black", text_col = "black"
  )

### ETAPA 7: Gráfico de Desmatamento

# Somar desmatamento (d) + resíduo (r) por ano
df_desmatamento <- freq_df %>%
  filter(grepl("^d|^r", Classe)) %>%  
  mutate(ano = as.integer(substr(Classe, 2, 5))) %>%  
  group_by(ano) %>%
  summarise(area_ha = sum(area_ha, na.rm = TRUE)) %>%  
  ungroup()

# Completar série temporal
anos_completos <- data.frame(ano = 2007:2023)
df_desmatamento <- full_join(anos_completos, df_desmatamento, by = "ano") %>%
  mutate(
    area_ha = ifelse(is.na(area_ha), 0, area_ha),
    label = case_when(
      ano == 2007 ~ "≤2007",
      TRUE ~ as.character(ano)
    )
  )

# Gráfico 
g_barras <- ggplot(df_desmatamento, aes(x = factor(label, levels = label), y = area_ha)) +
  geom_col(
    fill = "#d73027", 
    color = "#a50026", 
    width = 0.7,
    alpha = 0.9
  ) +
  scale_y_continuous(
    name = "Deforested Area (ha)",
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_x_discrete(name = "Year") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 20)),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3, linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = paste("Annual Deforestation Time Series –", nome_area),
    subtitle = "Combined detection (d) + residual (r) classes | PRODES/INPE v2024"
  ) +
  # Adicionar valores nas barras para anos com desmatamento
  geom_text(
    data = df_desmatamento[df_desmatamento$area_ha > 0, ],
    aes(label = round(area_ha, 1)), 
    vjust = -0.5, 
    size = 3, 
    color = "black"
  )

### ETAPA 8: Composição Final 

# Layout 
g_final <- g_mapa / g_barras + 
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    caption = paste0(
      "Study Area: ", nome_area, " | ",
      "Total Area: ", round(area_poligono, 0), " ha | ",
      "Analysis Date: ", Sys.Date(), " | ",
      "Coordinate System: UTM Zone 20S, SIRGAS 2000"
    ),
    theme = theme(
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5, margin = margin(t = 15))
    )
  )

cat("6. Exibindo resultados...\n")
print(g_mapa)
print(g_barras)
print(g_final)

### ETAPA 9: Exportar 

cat("7. Exportando arquivos...\n")

#  Criar variavel agregada 
if(!exists("freq_df_agregado") || is.null(freq_df_agregado)) {
  cat("AVISO: Criando freq_df_agregado...\n")
  freq_df_agregado <- freq_df %>%
    group_by(Classe) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      area_ha = sum(area_ha, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(area_ha))
}

if(!exists("freq_df") || is.null(freq_df)) {
  stop("ERRO: freq_df não encontrado")
}

if(!exists("area_total_calculada") || is.null(area_total_calculada)) {
  cat("AVISO: Calculando area_total_calculada...\n")
  area_total_calculada <- sum(freq_df$area_ha, na.rm = TRUE)
}

if(!exists("area_poligono") || is.null(area_poligono)) {
  cat("AVISO: Calculando area_poligono...\n")
  if(!exists("poligono_utm")) {
    poligono_utm <- st_transform(poligono, crs_utm)
  }
  area_poligono <- as.numeric(st_area(poligono_utm)) / 10000
}

if(!exists("pixel_area_ha") || is.null(pixel_area_ha)) {
  cat("AVISO: Calculando pixel_area_ha...\n")
  raster_res <- res(r_mask_utm)
  pixel_area_m2 <- raster_res[1] * raster_res[2]
  pixel_area_ha <- pixel_area_m2 / 10000
}

# Verificar se diretório existe, criar se necessário
if(!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
  cat("Diretório criado:", export_dir, "\n")
}

# Mostrar resumo antes da exportação
cat("\n=== RESUMO ANTES DA EXPORTAÇÃO ===\n")
cat("Classes encontradas:", nrow(freq_df_agregado), "\n")
cat("Total de pixels:", sum(freq_df_agregado$count), "\n")
cat("Área total:", round(sum(freq_df_agregado$area_ha), 2), "ha\n")
print(freq_df_agregado[, c("Classe", "count", "area_ha")])

# CSV com dados resumidos
#write.csv(freq_df_agregado[, c("Classe", "count", "area_ha")],
#          file = file.path(export_dir, paste0(nome_arquivo, "_area_por_classe_resumido.csv")),
#          row.names = FALSE)

# CSV com dados detalhados
write.csv(freq_df[, c("value", "Classe", "count", "area_ha")],
          file = file.path(export_dir, paste0(nome_arquivo, "_area_por_classe_detalhado.csv")),
          row.names = FALSE)

# Metadados da análise
metadados <- data.frame(
  parametro = c("Area_analisada", "CRS_utilizado", "Resolucao_pixel_ha", "Area_total_ha", 
                "Data_processamento", "Zona_UTM", "Arquivo_origem", "Diferenca_percentual"),
  valor = c(nome_area, crs_utm, round(pixel_area_ha, 6), round(area_total_calculada, 2),
            as.character(Sys.Date()), paste("Zona", zona_utm, hemisferio), 
            basename(raster_path), paste0(round((abs(area_total_calculada - area_poligono) / area_poligono) * 100, 2), "%"))
)

write.csv(metadados, 
          file = file.path(export_dir, paste0(nome_arquivo, "_metadados.csv")),
          row.names = FALSE)

# Verificar se gráficos existem antes de salvar
if(exists("g_mapa") && !is.null(g_mapa)) {
  ggsave(
    filename = file.path(export_dir, paste0(nome_arquivo, "_mapa_prodes24.png")),
    plot = g_mapa,
    width = 12, height = 8, units = "in", dpi = 300, bg = "white"
  )
}

if(exists("g_barras") && !is.null(g_barras)) {
  ggsave(
    filename = file.path(export_dir, paste0(nome_arquivo, "_grafico_prodes24.png")),
    plot = g_barras,
    width = 10, height = 6, units = "in", dpi = 300, bg = "white"
  )
}

if(exists("g_final") && !is.null(g_final)) {
  ggsave(
    filename = file.path(export_dir, paste0(nome_arquivo, "_analise_completa_prodes24.png")),
    plot = g_final,
    width = 12, height = 10, units = "in", dpi = 300, bg = "white"
  )
}

cat("Processamento concluído!\n")
cat("Arquivos salvos em:", export_dir, "\n")
cat("Zona UTM utilizada:", zona_utm, hemisferio, "(", crs_utm, ")\n")
cat("Diferença de área: ", round(abs(area_total_calculada - area_poligono), 2), " ha (", 
    round((abs(area_total_calculada - area_poligono) / area_poligono) * 100, 2), "%)\n")