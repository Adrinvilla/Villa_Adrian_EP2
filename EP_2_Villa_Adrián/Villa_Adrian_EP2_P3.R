
library (BoolNet) ## Cargar libreria

Red <- loadNetwork ("Reglas.txt") ## Cargar reglas
plotNetworkWiring (Red) ## Plotear

## 1

atractores <- getAttractors (Red) ## Ver atractores de la red
atractores


## 2

plotAttractors (Red)


## 4

plotAttractors (Red)


## 5

plotStateGraph (Red)
