---
title: "Tutorial - NDVI Temporal da Amazônia"
author: "Cássio Moquedace, Khaterine Martinez, Lara Lima, Rugana Imbana"
date: "13/09/2020"
---
## Objetivo
Criar uma animação espaço temporal do Índice de Vegetação por Diferença Normalizada (NDVI) da Amazônia Legal no ano de 2019

<p>&nbsp;</p>

## Índice de Vegetação por Diferença Normalizada (NDVI)
O NDVI trata-se de um índice que mensura a "saúde" da cobertura vegetal. O índice é derivado do cálculo da relação entre a reflexão e absorção da luz pelas bandas. Por exemplo, a vegetação saudável reflete a luz com mais intensidade na faixa do infravermelho próximo e com menor intensidade na porção do visível do espectro. Dessa forma podemos utilizar essa relação de reflectância dessas bandas para representar o potencial de saúde da cobertura do solo.

Geralmente o índice varia de -1 a 1, sendo -1 menos saudável ou água, e 1 muito saudável, provavelmente floresta densa. No entanto utilizaremos NDVI oriundo da coleção Terra MODIS da NASA, o produto MOD13A2. Essas imagens com NDVI já calculado a faixa é de -2.000 a 10.000, seguindo a mesma lógica de saúde da vegetação supracitada.

<p>&nbsp;</p>

## NDVI MODIS
O produto MOD13A2 Versão 6 que utilizaremos fornece valores de NDVI por pixel em resolução espacial de 1 quilômetro (km). Esse NDVI é derivado do *National Oceanic and Atmospheric Administration-Advanced Very High Resolution Radiometer* (NOAA-AVHRR).

### Características
- O NDVI com resolução temporal de 16 dias é gerado usando os dois grânulos de refletância de superfície composta pela média de 8 dias (MOD09A1).
- Esta entrada de refletância de superfície é baseada na abordagem de composição de azul mínima usada para gerar o produto de refletância de superfície de 8 dias.
- Um banco de dados de produtos do Índice de Vegetação Médio do *Climate Modeling Grid* (CMG) global atualizado frequentemente é usado para preencher as lacunas no conjunto de produtos CMG.

<p>&nbsp;</p>

## Execução dos códigos
### Carregando pacotes necessários e definindo local de trabalho
Limpando memória não utilizada no R
```{r message=FALSE}
gc()
```

Carregando pacotes
```{r message=FALSE}
pkg <- c("raster", "geobr", "dplyr", "dplyr", "sf", "MODIStsp", "stringr", 
         "lubridate", "ggplot2", "RColorBrewer", "scales", "ggsn", "animation",
         "grDevices", "ggpubr")

sapply(pkg, require, character.only = T)
```

Limpando todos os objetos gravados no R
```{r message=F}
rm(list = ls())
```

Definindo local de trabalho
```{r message=FALSE}
setwd("C:/R/tuto")
```

<p>&nbsp;</p>

### Baixando imagens MODIS
Criando arquivo para salvar o texto das mensagens de download, para posteriormente utilizar as datas das imagens
```{r message=FALSE, eval=FALSE, echo=TRUE}
con <- file("./saida.txt")
sink(con, append = T)
sink(con, append = T, type = "message")
```

Baixando dados
Maiores detalhes de como baixar imagens usando o pacote `MODIStsp` acesse [MODIStsp v2.0.2](https://docs.ropensci.org/MODIStsp/articles/interactive_execution.html)
```{r message=FALSE, eval=FALSE, echo=TRUE}
MODIStsp()
```

Fechando arquivo das mensagens e salvando em `.txt` na pasta de trabalho
```{r message=FALSE, eval=FALSE, echo=TRUE}
sink()
sink(type = "message")
```

<p>&nbsp;</p>

### Preparando dados
#### Criando e modificando texto das datas
Lendo arquivo
```{r message=FALSE}
data <- read.table("saida.txt", sep = "\n")
```

Remover as primeiras linhas desnecessárias
```{r message=FALSE}
data <- data.frame(data[-c(1:6), ])
```

Selecionando uma linha de cada data
```{r message=FALSE}
data <- data[seq(1, 322, 14), ]
```

Separando só as datas, retirando texto

(Utilizando `pacote::funcao`, no lugar de somente a função é para evitar conflito de funções com nomes iguais)
```{r message=FALSE}
data <- as.data.frame(stringr::str_split(data, ":"), stringsAsFactors = F)
data <- as.character(data[4, ]) %>% 
  gsub(" ", "", .)

print(data)
```

Datum e projeção que utilizaremos no mapa, mais informações das projeções de classe CRS PROJ.4, compatível com pacote `raster` do R acesse [EPSG.io](https://epsg.io/)
```{r message=FALSE}
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
```

<p>&nbsp;</p>

#### Tratando imagens
Lendo nome das imagens
```{r message=FALSE}
lista_nomes <- list.files(path = "raster/VI_16Days_1Km_v6/NDVI",
                          pattern = ".tif$", full.names = T)
```

Criando lista vazia
```{r message=FALSE}
lista_rasters <- list()
```

Criando `for` para leitura, reprojeção de cada imagem para o datum WGS84 e coordenadas geográficas e aumento dos tamanhos dos pixels das imagens para acelerar execução do código

O `fact = 5` utilizado indica que o tamanho do pixel será aumentado 5 vezes. No nosso caso que o tamanho do pixel original da imagem é de aproximadamente de 1 km, a resolução final terá 5 km cada pixel. A depender do seu objetivo com a imagem ou o tamanho da área de estudo, você pode variar esse valor ou ainda não aplicar a função `aggregate`

```{r error=TRUE, message=FALSE, warning=FALSE}
for (i in 1:length(lista_nomes)) {
  
  lista_rasters[[i]] <- raster(lista_nomes[i]) %>% 
    projectRaster(crs = wgs84, method = "bilinear") %>% 
    aggregate(fact = 3)
  
}
```

Empilhando rasters
```{r message=FALSE, error=TRUE}
r.emp <- raster::stack(lista_rasters)
print(r.emp)
```

Carregando vetor da Amazônia Legal e reprojetando para o mesmo sistema de coordenadas das imagens
```{r error=TRUE, message=FALSE, warning=FALSE}
amazonia_shp <- read_state(code_state = "all", 2018) %>% 
  st_transform(crs(r.emp)) %>% 
  dplyr::filter(abbrev_state %in% c("RO", "AC", "MT", "AM", "AP", "RR", "TO", 
                                    "MA", "PA")) 
```

Plotando imagens e vetor para visualizar a sobreposição
```{r error=TRUE, message=FALSE, fig.align='center', dpi=300}
plot(r.emp[[1]])
plot(amazonia_shp$geom, add = T)
```

Recortando imagens nos limites da Amazônia
```{r error=TRUE, message=FALSE}
mask_amz <- r.emp %>% 
  raster::crop(amazonia_shp) %>% 
  raster::mask(amazonia_shp)
```

Atribuindo `NA` em valores ausentes
```{r error=TRUE, message=FALSE}
mask_amz[mask_amz == 99999] = NA
```

Visualizando as imagens cortadas e com `NA`
```{r error=TRUE, message=FALSE, fig.align='center', dpi=300}
plot(mask_amz[[1]])
```

Atribuindo nome em cada raster dentro da pilha
```{r error=TRUE, message=FALSE}
names(mask_amz) <- data
```

<p>&nbsp;</p>

### Extraindo valores médios NDVI
Criando `data.frame` a partir das imagens
```{r error=TRUE, message=FALSE}
graphics.off()
gc()
ndvi.df  <-  rasterToPoints(mask_amz) %>% # Convertendo pixels em pontos
  data.frame() %>% # Convertendo pixels em pontos
  tidyr::gather(value = "NDVI", key = "data", contains("2019")) %>% # Organizando data.frame
  dplyr::mutate(data = str_sub(data, -10)) %>% # Removendo characterer excedente das datas
  dplyr::mutate(data = ymd(data)) # Convertendo para classe data

gc()
```

Extraindo a média do `data.frame`
```{r error=TRUE, message=FALSE}
ndvi.medio <- ndvi.df %>%
  dplyr::group_by(data) %>% # Agrupando por data
  dplyr::summarise(ndvi.m = mean(NDVI, na.rm = T)) %>% # Extraindo a média e summarizando por data
  as.data.frame() # Convertendo em data.frame
```

Criando um vetor com as datas
```{r error=TRUE, message=FALSE}
datas <- ndvi.df$data %>%
  unique()
```

Mudando a saída decimal do gráfico de `.` para `,`
```{r error=TRUE, message=FALSE}
options(OutDec = ",")
```

Criando limites para extensão do mapa, a partir do vetor da Amazônia
```{r error=TRUE, message=FALSE}
box.ndvi <- st_bbox(amazonia_shp) 
xrange <- box.ndvi$xmax - box.ndvi$xmin
yrange <- box.ndvi$ymax - box.ndvi$ymin
box.ndvi[1] <- box.ndvi[1] - (0.1 * xrange) 
box.ndvi[3] <- box.ndvi[3] + (0.1 * xrange) 
box.ndvi[2] <- box.ndvi[2] - (0.05 * yrange) 
box.ndvi[4] <- box.ndvi[4] + (0.05 * yrange)
box.ndvi <- box.ndvi %>%
  st_as_sfc()
```

Criando a legenda de texto para inserir no mapa
```{r error=TRUE, message=FALSE}
dados.sat <- "Sistemas de coordenadas geográficas\nDatum: WGS84\nColeção: Terra MODIS\nProduto: MOD13A2"
```

<p>&nbsp;</p>

### Criando animação
```{r error=TRUE, message=FALSE, eval=F, echo=T}
gc()

saveGIF(ani.height = 1000, ani.width = 1000, # Definindo resolução e tamanho da animação
        ani.res = 175, {
          
          for(i in 1:length(datas)){ # Criando sequência do plot
            
            ndvi.mes <- ndvi.df %>% # Criando data.frame com a data [i]
              dplyr::filter(data == datas[i]) %>% # Filtrando dados com data [i]
              dplyr::select(x, y, NDVI) %>% # Selecionando colunas a serem utilizadas
              rasterFromXYZ(crs = wgs84) %>% # Definindo datum
              rasterToPolygons() %>% # Convertndo data.frame para polígonos
              st_as_sf() # Convertendo para objeto classe sf
            
            gc()
            
            mes.atual <- month(datas[i], label = T, abbr = F) %>%
              str_to_title() # Filtrando mês e atribuindo a objeto para título do plot
            
            ndvi.plot <- ggplot(ndvi.mes) + # Criando plot do mapa
              geom_sf(aes(fill = NDVI), col = NA) + # Plot NDVI
              geom_sf(data = amazonia_shp, fill = NA) + # Plot contorno
              coord_sf(xlim = st_coordinates(box.ndvi)[c(1, 2), 1], # Definindo Extensão do mapa
                       ylim = st_coordinates(box.ndvi)[c(2, 3), 2],
                       expand = F) + 
              scale_fill_gradientn(colours = brewer.pal(n = 11, #Definindo cores e intervalo da paleta
                                                        name = "RdYlGn"), 
                                   limits = c(-2000, 10000), 
                                   breaks = seq(-2000, 10000, 2000)) + 
              theme_minimal() + # Tema do plot
              labs(title = mes.atual, x = "", y = "") + # Adicionando títulos
              theme(panel.border = element_rect(fill = NA), # Inserindo painel no mapa
                    panel.background = element_blank(), #Removendo cor de fundo do painel
                    plot.title  = element_text(hjust = 0.5, face = "bold")) + # Definindo estilo de escrita para título
              north(data = ndvi.mes, symbol = 1, location = "topright", # Adionando rosa dos ventos ao mapa
                    scale = 0.15) +
              scalebar(data = ndvi.mes, dist = 400, dist_unit = "km", # Adicionando escala ao mapa
                       location = "bottomright", transform = T,
                       model = "WGS84", st.size = 2, border.size = 0.1) +
              geom_text(x = -70, y = -17, label = dados.sat, size = 2.2) # Inserindo legenda de texto no mapa
            
            gc()
            
            cor.preench <- grDevices::adjustcolor(col = "grey70", # Criando uma cor para o preencdimento do gráfico de linha
                                                  alpha.f = 0.4)
            
            sazonal.plot <- ggplot(ndvi.medio, aes(x = data, y = ndvi.m)) + # Criando plot de linha
              geom_ribbon(aes(ymax = ndvi.m, ymin = min(ndvi.m) - 100), # Criando fita preenchida para o gráfico
                          fill = cor.preench) + 
              geom_line(col = "grey50") + # Inserindo linha ao gráfico
              geom_point(data = ndvi.medio[i, ], aes(x = data, y = ndvi.m,
                                                     fill = ndvi.m), pch = 21, # Inserindo o ponto
                         size = 2.5, col = "grey40") + 
              scale_fill_gradientn(colours = brewer.pal(n = 11,
                                                        name = "RdYlGn"), # Definindo escala de cor para o plot de linha
                                   limits = c(-2000, 10000), 
                                   breaks = seq(-2000, 10000, 3000)) + 
              scale_y_continuous(limits = c(min(ndvi.medio$ndvi.m) - 100, # Definindo extensão e intervalo do eixo y
                                            max(ndvi.medio$ndvi.m) + 100)) + 
              scale_x_date(date_breaks = "1 month", # Definindo extensão e intervalo do eixo x
                           labels = date_format("%b")) + 
              coord_cartesian(expand = F) + # Não expandir eixos dos gráficos
              labs(y = "NDVI", x = NULL) + # Inserindo títulos nos eixos
              theme_minimal() + #Definindo tema do plot
              theme(legend.position = "none", # Desativando legenda comum
                    panel.grid.minor = element_blank(), # Inserindo grid
                    axis.text.x = element_text(angle = 45, hjust = 0.5, # Definindo estilo de texto ao eixo x
                                               vjust = 0.5))
            
            map.linha <- ggarrange(ndvi.plot, sazonal.plot, nrow = 2, # Juntando mapa e gráfico de linha
                                   common.legend = T, legend = "right",
                                   heights = c(3, 1))
            
            gc()
            
            print(map.linha) # Imprimindo GIF
            
          }
        })
```

<p>&nbsp;</p>

### Animação pronta
<p align="center">
<img src="animation.gif" width="600" height="600">
</p>
