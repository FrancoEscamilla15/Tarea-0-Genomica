###Tarea 1 - Genomica funcional 24/01/22###
### Enrique Franco Garcia ###

secuency <- readRNAStringSet("first.fasta")
secuency
secuencia <- translate(secuency)
secuencia

####Problema 1 - Counting DNA Nucleotides####
####Sin librerias especializadas#### 
sac <-readDNAStringSet("saccharomyces.txt")#Cargue mi secuencia 
sac #La visualizo 
A <- grepRaw("A",sac, all=TRUE)#Use la función grepRaw que me va a permitir ver donde
#estan posicionadas todas las adeninas dentro de mi secuencia y lo guarde en un objeto llamado A
A #Lo visualizo
adeninas <- length(A)#A partir del objeto A que cree donde estaban todas las posiciones usando la funcion 
#lenght vi que cuantos caracteres tenía A y dado que contiene todas las posiciones de la adenina me va a arojjar
#el total de adeninas de mi secuencia
adeninas #Lo visualizo para compararlo con el numero que me dio usando libreria especializada

T <- grepRaw("T",sac, all = TRUE)#Repito el mismo procedimiento que con las adeninas
T #Lo visualizo
timinas<-length(T)#Repito el mismo procedimiento que con las adeninas
timinas #Lo visualizo

G <- grepRaw("G",sac, all = TRUE)#Repito el mismo procedimiento que con las adeninas
G #Lo visualizo
guaninas <- length(G)#Repito el mismo procedimiento que con las adeninas
guaninas #Lo visualizo

C <-grepRaw("C",sac, all = TRUE)#Repito el mismo procedimiento que con las adeninas
C #Lo visualizo
citocinas<-length(C)#Repito el mismo procedimiento que con las adeninas
citocinas #Lo visualizo

####Con librerias especializadas#### 
saccharomyces <- readDNAStringSet("saccharomyces.txt") #Cargue mi secuencia 
saccharomyces #La visualizo 
oligonucleotideFrequency(saccharomyces,1) #Usando esta función de biostrings puedo ver las
#veces que aparece cada base dentro de mi secuencia, el 1 indica que solamente voy a ver las bases
#por separado, si lo cambio por otro numero puedo ver también las diferentes combinaciones posibles entre 
#bases, como los codones abajo
oligonucleotideFrequency(saccharomyces,3)#De esta manera puedo ver las veces que aparece un codon

####Problema 2 - Calculating Protein Mass####
####Sin librerias especializadas#### 
library(Biostrings)
ecoli <- readAAStringSet("e coli.txt") #Cargue una secuencia de aminoacidos 
ecoli
alphabetFrequency(ecoli) #Visualice la frecuencia de amioacidos para después usarla en una función

aminoacidos <- function(){
  A <- readline(prompt = "¿Cuantas A tiene tu secucencia: ")
  A <- as.numeric(A)
  C <- readline(prompt = "¿Cuantas C tiene tu secucencia: ")
  C <- as.numeric(C)
  D <- readline(prompt = "¿Cuantas D tiene tu secucencia: ")
  D <- as.numeric(D)
  E <- readline(prompt = "¿Cuantas E tiene tu secucencia: ")
  E <- as.numeric(E)
  F <- readline(prompt = "¿Cuantas F tiene tu secuencia: ")
  F <- as.numeric(F)
  G <- readline(prompt = "¿Cuantas G tiene tu secuencia: ")
  G <- as.numeric(G)
  H <- readline(prompt = "¿Cuantas H tiene tu secuencia: ")
  H <- as.numeric(H)
  I <- readline(prompt = "¿Cuantas I tiene tu secuencia: ")
  I <- as.numeric(I)
  K <- readline(prompt = "¿Cuantas K tiene tu secuencia: ")
  K <- as.numeric(K)
  L <- readline(prompt = "¿Cuantas L tiene tu secuencia: ")
  L <- as.numeric(L)
  M <- readline(prompt = "¿Cuantas M tiene tu secuencia: ")
  M <- as.numeric(M)
  N <- readline(prompt = "¿Cuantas N tiene tu secuencia: ")
  N <- as.numeric(N)
  P <- readline(prompt = "¿Cuantas P tiene tu secuencia: ")
  P <- as.numeric(P)
  Q <- readline(prompt = "¿Cuantas Q tiene tu secuencia: ")
  Q <- as.numeric(Q)
  R <- readline(prompt = "¿Cuantas R tiene tu secuencia: ")
  R <- as.numeric(R)
  S <- readline(prompt = "¿Cuantas S tiene tu secuencia: ")
  S <- as.numeric(S)
  T <- readline(prompt = "¿Cuantas T tiene tu secuencia: ")
  T <- as.numeric(T)
  V <- readline(prompt = "¿Cuantas V tiene tu secuencia: ")
  V <- as.numeric(V)
  W <- readline(prompt = "¿Cuantas W tiene tu secuencia: ")
  W <- as.numeric(W)
  Y <- readline(prompt = "¿Cuantas Y tiene tu secuencia: ")
  Y <- as.numeric(Y)
  aminoacido1 <- A * 71.03711
  aminoacido2 <- C * 103.00919
  aminoacido3 <- D * 115.02694
  aminoacido4 <- E * 129.04259
  aminoacido5 <- F * 147.06841
  aminoacido6 <- G * 57.02146
  aminoacido7 <- H * 137.05891
  aminoacido8 <- I * 113.08406
  aminoacido9 <- K * 128.09496
  aminoacido10 <- L * 113.08406
  aminoacido11 <- M * 131.04049
  aminoacido12 <- N * 114.04293
  aminoacido13 <- P * 97.05276
  aminoacido14 <- Q * 128.05858
  aminoacido15 <- R * 156.10111
  aminoacido16 <- S * 87.03203
  aminoacido17 <- T * 101.04768
  aminoacido18 <- V * 99.06841
  aminoacido19 <- W * 186.07931
  aminoacido20 <- Y * 163.06333 
  sum <- aminoacido1 + aminoacido2 + aminoacido3 + aminoacido4 +
    aminoacido5 + aminoacido6 + aminoacido7 + aminoacido8 + aminoacido9 +
    aminoacido10 + aminoacido11 + aminoacido12 + aminoacido13 + 
    aminoacido14 + aminoacido15 + aminoacido16 + aminoacido17 +
    aminoacido18 + aminoacido19 + aminoacido20
  print(return(print(paste("La suma total es de: ", sum)))) 
}
#Cree una funcion que lo que hace es preguntarte la cantidad de aminoacidos que 
#tiene tu secuencia, esto se puede conocer facilmente usando la funcion "alphabetFrecuency"
#la función ya tiene el peso de cada aminoácido por lo cual al momento de que el 
#usuario ingrese su cantidad de A (por ejemplo) automaticamente la función lo multiplicara por
#71.03711 que es su peso, al final la función sumara el resultado obtenido de cada calculo 
#de aminoacidos de manera que se obtendrá el tamaño total de la secuencia como lo solicita el 
#ejercicio. Se que es tedioso tener que introducir manualmente cada una de los aminoacidos
#pero no se me salión de otra forma :(

####Con librerias especializadas####
install.packages("Peptides")
library(Peptides)
mw(seq = "MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRFAWFPAPEAV
PESWLECDLPEADTVVVPSNWQMHGYDAPIYTNVTYPITVNPPFVPTENPTGCYSLTFNVDESWLQEGQT
RIIFDGVNSAFHLWCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWSDGSYLEDQDMWRMSGIF
RDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYLRVTVSLWQGETQVASGTAPFGGEII
DERGGYADRVTLRLNVENPKLWSAEIPNLYRAVVELHTADGTLIEAEACDVGFREVRIENGLLLLNGKPL
LIRGVNRHEHHPLHGQVMDEQTMVQDILLMKQNNFNAVRCSHYPNHPLWYTLCDRYGLYVVDEANIETHG
MVPMNRLTDDPRWLPAMSERVTRMVQRDRNHPSVIIWSLGNESGHGANHDALYRWIKSVDPSRPVQYEGG
GADTTATDIICPMYARVDEDQPFPAVPKWSIKKWLSLPGETRPLILCEYAHAMGNSLGGFAKYWQAFRQY
PRLQGGFVWDWVDQSLIKYDENGNPWSAYGGDFGDTPNDRQFCMNGLVFADRTPHPALTEAKHQQQFFQF
RLSGQTIEVTSEYLFRHSDNELLHWMVALDGKPLASGEVPLDVAPQGKQLIELPELPQPESAGQLWLTVR
VVQPNATAWSEAGHISAWQQWRLAENLSVTLPAASHAIPHLTTSEMDFCIELGNKRWQFNRQSGFLSQMW
IGDKKQLLTPLRDQFTRAPLDNDIGVSEATRIDPNAWVERWKAAGHYQAEAALLQCTADTLADAVLITTA
HAWQHQGKTLFISRKTYRIDGSGQMAITVDVEVASDTPHPARIGLNCQLAQVAERVNWLGLGPQENYPDR
LTAACFDRWDLPLSDMYTPYVFPSENGLRCGTRELNYGPHQWRGDFQFNISRYSQQQLMETSHRHLLHAE
EGTWLNIDGFHMGIGGDDSWSPSVSAEFQLSAGRYHYQLVWCQK",monoisotopic = FALSE)

#Con ayuda de la libreria llamada Peptides use la funcion mw (molecular weight), esta funcion se encarga de 
#calcular la suma de la masa de cada aminoacido presente en la secuencia con base en la 
#escala pI/Mw, tengo entendido que se basa en la pagina web expesy la cual también permite calcular
#el peso molecular 
#En el argumento monoisotopic indica si se deben usar pesos monoisotopicos de aminoacidos (T o F)
#Observacion: en el ejemplo de Rosalind me fije que varia un poco el resultado y buscando mas...
#descubri que Expasy colaboro en esta libreria, asi que supongo que por eso varia un poco.

