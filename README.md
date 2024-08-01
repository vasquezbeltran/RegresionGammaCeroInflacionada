
# Regresion gamma cero–inflacionada

## Objetivo del repositorio

Presentar el código que replica los resultados de la tesis "*Modelos de regresión gamma generalizada cero–inflacionada para la media con aplicación a gastos en educación*", PUCP, 2018, https://tesis.pucp.edu.pe/repositorio/handle/20.500.12404/12999 

## Introducción de la tesis

Una **variable semicontinua o cero–inflacionada** es aquella que puede tomar valores continuos y no negativos, incluyendo el valor cero con probabilidad no nula. En el análisis de regresión, el **modelo de dos partes (MDP)** sirve para explicar una variable respuesta semicontinua. Con el MDP se asume que la variable respuesta, $Y$, presenta una distribución mixta de probabilidades: (1) una distribución de Bernoulli para explicar si toma el valor cero o no y (2) una distribución continua positiva para cuando no es cero. Además, los parámetros son expresados de tal manera que posibilite estimar los efectos de covariables sobre (1) la media de la respuesta condicionada a que tome valores positivos, $E[Y|Y>0]$, y (2) la probabilidad de que la respuesta tome el valor cero, $P[Y=0]$. 

El objetivo de la tesis es estudiar un modelo alternativo al MDP, denominado **modelo de regresión cero–inflacionada a la media (MCIM)**, cuya parametrización permite estimar e interpretar efectos de covariables sobre la media total de la respuesta, $E[Y]$, en lugar de la media condicionada a valores positivos, $E[Y|Y>0]$. Además, optamos por la **distribución gamma generalizada (MCIM–GG)** para modelar ciertas características de los valores positivos de la respuesta, tales como la asimetría positiva y la curtosis pronunciada. 

Con el estudio de simulación, encontramos un adecuado desempeño de los **estimadores de máxima verosimilitud** del MCIM–GG bajo diferentes escenarios definidos según porcentajes de valores ceros de la respuesta y tamaños de muestra. Por último, con el estudio de aplicación, utilizamos MCIM–GG y MDP–GG para estimar los efectos de ciertas covariables sobre **gastos en educación** en adolescentes participantes del estudio Niños del Milenio en el Perú. 

## 1_Datos

- Boyden, J. (2022). *Young Lives: an International Study of Childhood Poverty: Round 3, 2009*. [data collection]. 4th Edition. UK Data Service. SN: 6853, https://beta.ukdataservice.ac.uk/datacatalogue/studies/study?id=6853 :
  - `pe_oc_householdlevel.dta` 
  - `pe_oc_childlevel.dta`
  - `pe_oc_householdmemberlevel.dta`
- `data.xlsx`: Excel con datos finales para el estudio de aplicación.

## 2_Codigo

| n   | Code file              | Type file       | Descripción   |
| --- | ---                    | ---             | ---           |
| 1   | [`masterfile.m`        ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/masterfile.m)                  | MATLAB script   | Master file que ejecuta de forma ordenada todos los códigos utilizados en la tesis, incluye tanto los códigos del capítulo 4, estudio de simulación, como los del capítulo 5, estudio de aplicación. | 
| 2   | [`gg.m`                ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/gg.m)                 | MATLAB function | Función de densidad de la distribución gamma generalizada de acuerdo a la parametrización propuesta por Manning (2005). |
| 3   | [`kfun_mcim_gg.m`      ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/kfun_mcim_gg.m)       | MATLAB function | Función de log–verosimilitud del MCIM–GG, su vector gradiente y su matriz hessiana. |
| 4   | [`kfun_mcim_g.m`       ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/kfun_mcim_g.m)        | MATLAB function | Función de log–verosimilitud del MCIM–G, su vector gradiente y su matriz hessiana. |
| 5   | [`kfun_mdp_gg.m`       ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/kfun_mdp_gg.m)        | MATLAB function | Función de log–verosimilitud del MDP–GG. |
| 6   | [`kfun_mdp_g.m`        ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/kfun_mdp_g.m)         | MATLAB function | Función de log–verosimilitud del MDP–G. |
| 7   | [`ysim_mcim_gg.m`      ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/ysim_mcim_gg.m)       | MATLAB function | Generación de valores simulados de la variable respuesta del MCIM–GG. |
| 8   | [`simulacion_generar.m`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/simulacion_generar.m) | MATLAB script   | Estudio de simulación (capítulo 4). Generación de bases de datos simulados de la variable respuesta y de las covariables del MCIM–GG. |
| 9   | [`simulacion_estimar.m`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/simulacion_estimar.m) | MATLAB script   | Estudio de simulación (capítulo 4). Estimación del modelo MCIM–GG con cada una de las bases de datos simulados. |
| 10  | [`aplicacion.m`        ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/aplicacion.m)         | MATLAB script   | Estudio de aplicación (capítulo 5). Estimación de los modelos MCIM–GG y MDP–GG. |
| 11  | `variables.do`          | STATA do–file   | Construcción de las variables del estudio de aplicación. Output file: `data.xlsx` |
| 12  | `analisisexp.R`         | R script        | Análisis exploratorio de datos. Output files: `figura51.jpg`, `cuadro51.csv` y `figura52.jpg` |
| 13  | [`aplicacion.sas`      ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/2_Codigo/aplicacion.sas)       | SAS program     | Estudio de aplicación (capítulo 5). Estimación del modelo MDP–GG. Output: `cuadro53.pdf` |

## 3_Outputs

“El reporte sube el valor agregado del algoritmo”

Capítulo 4. Estudio de simulación. 
- [`cuadro41.csv`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro41.csv) [`cuadro41.mat`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro41.mat) : &nbsp; Cuadro 4.1 Resultados de simulación MCIM–GG donde porcentaje de ceros 10%
- [`cuadro42.csv`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro42.csv) [`cuadro42.mat`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro42.mat) : &nbsp; Cuadro 4.2 Resultados de simulación MCIM–GG donde porcentaje de ceros 20%
- [`cuadro43.csv`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro43.csv) [`cuadro43.mat`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro43.mat) : &nbsp; Cuadro 4.3 Resultados de simulación MCIM–GG donde porcentaje de ceros 40%
- [`figuraa.jpg` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figuraa.jpg) [`figuraa.fig` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figuraa.fig) : &nbsp; Resultados de sesgo relativo (%) de $\omega$
- [`figurab.jpg` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurab.jpg) [`figurab.fig` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurab.fig) : &nbsp; Resultados de sesgo relativo (%) de $\beta$
- [`figurac.jpg` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurac.jpg) [`figurac.fig` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurac.fig) : &nbsp; Resultados de RECM de $\omega$
- [`figurad.jpg` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurad.jpg) [`figurad.fig` ](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/figurad.fig) : &nbsp; Resultados de RECM de $\beta$

Capítulo 5. Estudio de aplicación.
- `figura51.jpg` : &nbsp; Figura 5.1 Histograma de gasto en educación
- `cuadro51.csv` : &nbsp; Cuadro 5.1 Características de los adolescentes según decisión de gastar
- `figura52.jpg` : &nbsp; Cuadro 5.2 Gráficos de cajas y dispersión de gasto en educación
- [`cuadro52.csv`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro52.csv) [`cuadro52.mat`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro52.mat) : &nbsp; Cuadro 5.2 Estimación de coeficientes de MCIM–GG
- [`cuadro53.pdf`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro53.pdf)                                                                                                                       : &nbsp; Cuadro 5.3 Estimación de coeficientes de MDP–GG
- [`cuadro54.csv`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro54.csv) [`cuadro54.mat`](https://github.com/vasquezbeltran/RegresionGammaCeroInflacionada/blob/master/3_Outputs/cuadro54.mat) : &nbsp; Cuadro 5.4 Criterios de información de los modelos

## 4_Diseminacion

- `vasquezbeltran_tesis.pdf` : documento final
- `vasquezbeltran_tesispresentacion.pdf`

## Bibliografía

- Bayes, C. L. y Valdivieso, L. H. (2016). A beta inflated mean regression model for fractional response variables, Journal of Applied Statistics 43(10): p. 1814-1830. https://www.tandfonline.com/doi/abs/10.1080/02664763.2015.1120711
- Smith, V., Preisser, J. S., Neelon, B. y Maciejewski, M. L. (2014). A marginalized two-part model for semicontinuous data, Statistics in Medicine 33: p. 4891-4903. https://onlinelibrary.wiley.com/doi/10.1002/sim.6263
- Manning, G. W., Basu, A. y Mullahy, J. (2005). Generalized modeling approaches to risk adjustment of skewed outcomes data, Journal of Health Economics 24: p. 465-488. https://www.sciencedirect.com/science/article/abs/pii/S0167629605000056?via%3Dihub
