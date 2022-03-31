
BiocManager::install(c("IRanges", "rtracklayer", "GenomicRanges", "Rsamtools")) ## Instalar paquetes necesarios

library(Biobase) ## Cargar libreria
library(IRanges) ## Cargar libreria

x <- IRanges (start = c (11,35,40), end = c (20,50,63)) ## Crear una instancia de IRange, o sea, un objeto con 3 rangos
x ## Ver el los 3 rangos 

start (x) # Todos los inicios
end (x)   # Todos los finales
width (x) # Ancho de cada rango
range (x) # Rango total de todos


coverage(x) # La suma de la cobertura en cada posición

reduce(x)   # Une los rangos que estan encimados


exons <- reduce (x) 


reads <- IRanges (start = c (1,21,30,50,80), width = 20)
reads


countOverlaps(exons, reads)

## Honestamente, no tengo idea de que pretendian hacer en esta parte, entiendo que querian hacer rangos, y lo que nos dice "reads", es que tienen 5 rangos, divididos en intervalos de 20
## Pero no entiendo qué significa el Overlap "1 3"


## Obteniendo la anotación de un genoma

library (rtracklayer) ## Cargar libreria 

load ("human.Rdata") ## Datos resumidos del humano
human ## Bueno, ni tan resumidos, son 852534 filas

seqnames (human) ## Nombres de las secuencias del archivo
ranges (human)
strand (human) 
mcols (human)

table (mcols (human)$gene_biotype) ## De forma resumida junto con lo anterior, nos pone en una tabla de que tipo de secuencia se tienen datos 

mcols(human) <- mcols (human)[,c("source","gene_id","gene_name","gene_biotype")] ## Cambia el nombre de las columnas 
mcols(human)

## Ejercicios 

# Cómo le harian para quedarse exclusivamente con las anotaciones de "miRNA"?



# y solamente aquellas anotaciones de la cadena "-"?

negativa <- subset (human, strand (human) == "-") ## Subset para seleccionar la cadena especifica
negativa

## Anotación de secuencias mapeadas

library (Rsamtools) ## Cargar libreria

## 

what <- c ("rname", "strand", "pos", "qwidth")
param <- ScanBamParam (what = what)

bam <- scanBam ("human_mapped_small.bam", param = param)

class (bam)
lapply (bam, names)


mapGR = GRanges(
  seqnames = bam[[1]]$rname,
  ranges   = IRanges(start=bam[[1]]$pos, width=bam[[1]]$qwidth),
  strand   = bam[[1]]$strand
)

mapGR

mcols(human)$counts <- countOverlaps(human, mapGR)
mcols(human)

typeCounts <- aggregate (mcols(human)$counts, by=list("biotype"=mcols(human)$gene_biotype), sum)
typeCounts


geneCounts <- aggregate(mcols(human)$counts, by=list("id"=mcols(human)$gene_name), sum)
head(geneCounts)


minCount = 40000
typeCountsHigh = typeCounts[typeCounts$x > minCount,]
typeCountsHigh = typeCountsHigh[order(typeCountsHigh$x),]
typeCountsHigh = rbind(data.frame("biotype"="other",
                                  "x"=sum(typeCounts$x[typeCounts$x <= minCount])),
                       typeCountsHigh)

pie(typeCountsHigh$x, labels=typeCountsHigh$biotype, col=rev(rainbow(nrow(typeCountsHigh))),
    main="Number of aligned reads per biotype")

## Biotipos del mapeo

