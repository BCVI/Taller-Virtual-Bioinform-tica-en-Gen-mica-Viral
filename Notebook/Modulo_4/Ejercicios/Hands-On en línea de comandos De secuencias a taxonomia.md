# Hands-On en línea de comandos: De secuencias a taxonomia

A partir de las secuencias de fagos que ensamblaron en el Modulo 3: Identificación de Secuencias Virales (Metagenómica viral y profagos), se va a realizar la clasificación, clustering y anotación taxonómica. 

En su directorio de `home` deben tener el directorio `Final_Viral_Contigs` donde se encuentra el archivo `viral_sequences.fna`. Este archivo contiene los contigs ensamblados ayer y dos genomas de fagos de referencia.

## 1.Realizar clasificación y clustering usando el genoma (VIRIDIC)

la herramienta **VIRIDIC** no solo se puede cirrer usando la interfaz web, en caso de tener una gran cantidad de contigs/genomas es mejor correrlo usando un servidor de HPC (High-Performance Computing).  

Ahora en su `home`cree un diretorio para la el modulo 4 clusterización y taxonomia

```bash
#Entrar al directorio de home, $HOME es una variable predefinida que contiene el path hacia su home. 
cd $HOME 
#crear directorio para modulo 4
mkdir dia4_clustering_taxonomia
#copiar los contifs en el nuevo directorio
cp /home/gama/Workshop_2024/Final_Viral_Contigs/viral_sequences.fna $HOME/dia4_clustering_taxonomia
```

Para correr VIRIDIC en el servidor, se necesita como entrada las secuencias fasta en un solo archivo (igual que en el servidor web).

```bash
#Atención: este comando no deja autocompletar así que deben estar muy seguros de como estan escritos los directorios. 
#projdir = nombre del directorio donde quieren guardar la salida de VIRIDIC
#in = archivo con las secuencias fasta.
/App/viridic_v1.0_r3.6/viridic.bash projdir=/home/came/dia4_clustering_taxonomia/final_contigs_VIRIDIC.out in=/home/came/dia4_clustering_taxonomia/viral_sequences.fna
```

 Ahora tiene un directorio llamado `final_contigs_VIRIDIC.out/`donde se encuentran todos los archivos de salida de VIRIDIC. Dentro de este directorio el resultado final se encuentra en `04_VIRIDIC_out`. 

```bash
#Para ver la matriz de similitud 
cat $HOME/dia4_clustering_taxonomia/final_contigs_VIRIDIC.out/04_VIRIDIC_out/sim_MA_genCol.csv
#Para ver los clusters generados
cat $HOME/dia4_clustering_taxonomia/final_contigs_VIRIDIC.out/04_VIRIDIC_out/clusters.csv
```

¿Cuántos clusters se generaron? ¿Coincide con la matriz de similitud?

## 2. Generación de un dendrograma/árbol filogénetico a partir de la matriz de similitud. 

Para completar este paso vamos a usar las librerias de R, ape y phagorn. Para abir R: 

```bash
cd $HOME/dia4_clustering_taxonomia
mkdir trees
cd trees
R 
```

En R debe leer el archivo de matriz de similitud, pasarlo a matriz de distancias y generar un cluster. en R la variable $HOME no funciona. 

```R
#Importar la matriz de similitud y volverla una matriz numerica
#RECUERDE REMPLAZAR MI USUARIO (came) POR EL SUYO 
mat_similitud=read.table("/home/came/dia4_clustering_taxonomia/final_contigs_VIRIDIC.out/04_VIRIDIC_out/sim_MA_genCol.csv",sep="\t",h=T)
rownames(mat_similitud)=mat_similitud$genome
mat_similitud=mat_similitud[,-1]
#Pasarla a distancias de 0 a 1
mat_distancias=1-(mat_similitud/100)
#Tansformarla a un objeto de distancias en R
distancias <- as.dist(mat_distancias)

#Cargue las librerias
library(ape)
library(phangorn)

#Crear un arbol con Neighbour-Joining (ape)
tree_nj <- nj(distancias)

#Crear un arbol con UPGMA (phangorn)
tree_upgma <- upgma(distancias)

#Guardar los arboles en un archivo de texto. 
write.nexus(tree_nj, file='NJ_tree.nex')
write.nexus(tree_upgma, file='UPGMA_tree.nex')

#Cerrar sesion
q() #presionar Y seguido de Enter. 
```

Ahora, puede descargar los archivos para visualizarlos en iTOL (para que el taller sea mas rápido descarguelos del github Modulo4/Ejercicios). Para visualizar dirijase al sitio web iTOL (interactive Tree Of Life https://itol.embl.de/ ), cree un usuario o inicie sesión. 

Ingrese a "MyTrees > Tree upload", seleccione los archivos .nex. Luego Ingrese al link con el nombre del archivo que aparece en iTOL y puede ver el árbol. 

## 3. Análisis del proteoma 

Ahora para analizar el proteoma a partir del genoma, se debe realizar primero una anotación, es decir, se deben predecir los genes en la secuencia. Para esto vamos a usar el programa pharokka (diseñado para fagos - https://github.com/gbouras13/pharokka). 

```bash
#Cree un directorio para este nuevo ejercicio
mkdir $HOME/dia4_clustering_taxonomia/proteoma

#Para correr pharokka, debe seguir las siguientes instrucciones (dependen del servidor, no siempre se necesita)
cd $HOME
cp .bashrc.conda .bashrc
# presione la tecla y seguida de enter
#salir y volver a iniciar sesión (Ctrl+D)

#luego correr pharokka (este si es el comando de siempre), -f solo sirve para que corra rapido, si desean hacer un análisis mas exahustivo deben quitarselo
#ATENCION. Este programa demora en correr al menos 3-4min
conda activate pharokka_env
cd $HOME/dia4_clustering_taxonomia/proteoma
pharokka.py -i $HOME/dia4_clustering_taxonomia/viral_sequences.fna -o annotation_out -d /home/gama/Workshop_2024/Final_Viral_Contigs/pharokka/ -t 3 -f

#Revisar la salida. 
cat $HOME/dia4_clustering_taxonomia/proteoma/annotation_out/pharokka.gbk
```

Puede ver la sintenia de las proteinas usando herramientas como https://cagecat.bioinformatics.nl/tools/clinker.

## 4. Todo en 1: PhageGCN 

Esta herramienta corre desde la anotación hasta la clasificación taxonomica. 

```bash
cp .bashrc.mamba .bashrc
# presione la tecla y seguida de enter
#salir y volver a iniciar sesión (Ctrl+D)

#Cargar phagcn
mamba activate phagcn
#Entrar al directorio PhaGCN_newICTV y ejecuta el siguiente comando ahí:
cd $HOME/PhaGCN_newICTV
python /App/PhaGCN_newICTV/run_Speed_up.py --contigs /home/came/dia4_clustering_taxonomia/viral_sequences.fna --len 8000

#En el directorio out, se encuentra la salida 
ls $HOME/PhaGCN_newICTV/out 
```

