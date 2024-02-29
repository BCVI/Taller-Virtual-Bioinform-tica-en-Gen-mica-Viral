# Ejercicio comparacion de proteomas

## Objetivo
- Familiarizarse con el programa VirClust desde el servidor web
- Comprender los archivos de entrada y salida de la herramienta
- Introducir los conceptos de grupos ortologos de proteinas, proteinas core y sintenia. 

## Para Comenzar
1. Dirijase al servidor web de VirClust ([https://rhea.icbm.uni-oldenburg.de/viridic/](https://rhea.icbm.uni-oldenburg.de/virclust/))
2. Descargar el archivo phages_arthrobacter.fasta 
3. En el tab I.DEFINE PROJECT, cargue el archivo con los genomas al servidor y digite un nombre de usuario y de proyecto
4. En el tab II.RUN CACULATIONS, corra cada uno de los pasos, uno por uno.



## Analisis 
### Step 1A = from genomes to proteins.
En este paso se usan las tablas de codones para la especie deseada para 1) Predecir donde hay genes y 2) Traducir de nucleotidos a genes.
  - Cuantos genes tiene cada genoma?
### Step 2A= Proteins to Protein Clusters (PCs)
En este paso se generan los grupos ortólogos entre las proteínas de los genomas de phages_arthrobacter.fasta . Se usa BLASTp para comparar las secuencias de proteínas y si un hit no cumple con las especificaciones se elimina. 
  - Cuantos PCs se generaron? El PC12 lo comparten cuales genomas ?
### Step 3A = Order genomes hierarchically. 
Similar a lo que hace VIRIDIC, VirClust calcula la distancia (en vez de similitud) entre todos los genes/proteínas de todos los genomas. Si un PC no esta en ningún otro genoma su distancia es 1. Luego se calcula la suma de todas las distancias (entre mayor sea la suma, menos similar es ese genoma a los otros). 
  - Cual es la distancia basada en PCs entre los genomas Yavru y Berka?
  - Genere el arbol.
      - Descargue el archivo
      - Dirijase a itol (https://itol.embl.de/)
      - Inicie sesion
      - Click en UploadTree
      - Este resultado coincide con lo que habíamos visto antes comparando los genomas a nivel de nucleotidos? En caso de que no, por que?
  - Que puede inferir del heatmap?
### Step 4A. Calculate stats and split in genome clusters (VGCs). 
Este paso crea una tabla de presencia/ausencia de los clusters proteicos en cada genoma. Luego calcula el numero de proteínas que comparten con los otros genomas, cuantas son únicas. Con estas tablas se pueden hacer graficas interesantes en R. Ex. Diagrama de venn. Este paso también genera un heatmap de PCs compartidos y su clustering
### Step 5A. Calculate core proteins for each VGC, based on PCs. 
Este paso genera las proteínas que son core para cada cluster de genomas viral.
  - Cuantos clusters se generaron?
Ahora, como vimos en la teoria no solo cuenta la similitud entre los genomas. También es importante la organización de los genomas (sintenia). 
- Ir a la pagina web (https://phamerator.org/).
  - Click en GenomeMaps
  - buscar cluster FE
  - seleccionar Berka, CabbageMan_Draft, Corgi, Yavru.
  - Clik View Maps.
- Es la sintenia (organizacion de los genes) similar entre los genomas ?  
