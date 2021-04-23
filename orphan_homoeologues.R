install.packages("rapportools")
library(rapportools)
library(plyr)
library(tidyr)
library(dplyr)

homo <- read.csv("~/Documents/Sobia/orphans/Homoeologues/homo.csv")
# 
# homo$all<-paste(homo$A, homo$B, homo$D)
# homo<-homo[ !duplicated(homo$all), ]

homo<-homo%>%
  mutate(Cardinality=case_when(
Cardinality == "00:00:01" ~	"0:0:1",
Cardinality == "00:01:00" ~	"0:1:0",
Cardinality == "00:01:01" ~	"0:1:1",
Cardinality == "01:00:00" ~	"1:0:0",
Cardinality == "01:00:01" ~	"1:0:1",
Cardinality == "01:01:00" ~	"1:1:0",
Cardinality == "01:01:01" ~	"1:1:1"
  ))
# 
# homo<-homo%>% 
#   mutate(Cardinality_type=case_when(
#     Cardinality_type == "Dyad"  ~ "Pair",
#     Cardinality_type == "Triad" ~ "Triplet",
#     Cardinality_type == "Singlet" ~ "Single"
#   ))
# 

{
  setwd("nr/")
filesA<-list.files(pattern = "A")
filesB<-list.files(pattern = "B")
filesD<-list.files(pattern = "D")

list<-list()
for(file in filesA){
  cat<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,1]
  chrom<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,3]
  data<-read.csv(file, header=F)
  data<-as.character(data$V1)
  data<-data[-grep("LC", data)]
  data<-as.data.frame(data)
  data$Category<-paste(cat)
  data$Chromosome<-paste(chrom)
  colnames(data)<-c("Gene", "designation", "Chromosome")
  list[[length(list)+1]]<-data
}

all_A<-as.data.frame(do.call(rbind.data.frame, list))

list<-list()
for(file in filesB){
  cat<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,1]
  chrom<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,3]
  data<-read.csv(file, header=F)
  data<-as.character(data$V1)
  data<-data[-grep("LC", data)]
  data<-as.data.frame(data)
  data$Category<-paste(cat)
  data$Chromosome<-paste(chrom)
  colnames(data)<-c("Gene", "designation", "Chromosome")
  list[[length(list)+1]]<-data
}

all_B<-as.data.frame(do.call(rbind.data.frame, list))

list<-list()
for(file in filesD){
  cat<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,1]
  chrom<-do.call('rbind', strsplit(file,'-',fixed=TRUE))[,3]
  data<-read.csv(file, header=F)
  data<-as.character(data$V1)
  data<-data[-grep("LC", data)]
  data<-as.data.frame(data)
  data$Category<-paste(cat)
  data$Chromosome<-paste(chrom)
  colnames(data)<-c("Gene", "designation", "Chromosome")
  list[[length(list)+1]]<-data
}

all_D<-as.data.frame(do.call(rbind.data.frame, list))
}

{
#  A genome 
poaceae<-all_A[(all_A$designation=="Pooeceae"),]
pooideae<-all_A[(all_A$designation=="Pooideae"),]
triticeae<-all_A[(all_A$designation=="Triticeae"),]
wheat<-all_A[(all_A$designation=="SSG"),]
poaceae$score<-4
pooideae$score<-3
triticeae$score<-2
wheat$score<-1
poaceae<-poaceae[ !duplicated(poaceae$Gene), ]
pooideae<-pooideae[ !duplicated(pooideae$Gene), ]
triticeae<-triticeae[ !duplicated(triticeae$Gene), ]
wheat<-wheat[ !duplicated(wheat$Gene), ]
A<-rbind(poaceae, pooideae, triticeae, wheat)
A$Genome<-"A"

#  B genome
poaceae<-all_B[(all_B$designation=="Pooeceae"),]
pooideae<-all_B[(all_B$designation=="Pooideae"),]
triticeae<-all_B[(all_B$designation=="Triticeae"),]
wheat<-all_B[(all_B$designation=="SSG"),]
poaceae$score<-4
pooideae$score<-3
triticeae$score<-2
wheat$score<-1
poaceae<-poaceae[ !duplicated(poaceae$Gene), ]
pooideae<-pooideae[ !duplicated(pooideae$Gene), ]
triticeae<-triticeae[ !duplicated(triticeae$Gene), ]
wheat<-wheat[ !duplicated(wheat$Gene), ]
B<-rbind(poaceae, pooideae, triticeae, wheat)
B$Genome<-"B"

# D genome
poaceae<-all_D[(all_D$designation=="Pooeceae"),]
pooideae<-all_D[(all_D$designation=="Pooideae"),]
triticeae<-all_D[(all_D$designation=="Triticeae"),]
wheat<-all_D[(all_D$designation=="SSG"),]
poaceae$score<-4
pooideae$score<-3
triticeae$score<-2
wheat$score<-1
poaceae<-poaceae[ !duplicated(poaceae$Gene), ]
pooideae<-pooideae[ !duplicated(pooideae$Gene), ]
triticeae<-triticeae[ !duplicated(triticeae$Gene), ]
wheat<-wheat[ !duplicated(wheat$Gene), ]
D<-rbind(poaceae, pooideae, triticeae, wheat)
D$Genome<-"D"

}

sorted<-A[order(A$Gene, abs(A$score) ), ]
A<-sorted[ !duplicated(sorted$Gene), ]
colnames(A)<-c("GeneA", "DesignationA", "Chromosome","score","Genome")

sorted<-B[order(B$Gene, abs(B$score) ), ]
B<-sorted[ !duplicated(sorted$Gene), ]
colnames(B)<-c("GeneB", "DesignationB", "Chromosome","score","Genome")

sorted<-D[order(D$Gene, abs(D$score) ), ]
D<-sorted[ !duplicated(sorted$Gene), ]
colnames(D)<-c("GeneD", "DesignationD", "Chromosome","score","Genome")


mergeA<-merge(homo, A, by.x="A", by.y="GeneA", all=T)
mergeB<-merge(mergeA, B, by.x="B", by.y="GeneB", all=T)
mergeD<-merge(mergeB, D, by.x="D", by.y="GeneD", all=T)
All_merged<-mergeD[!with(mergeD,is.na(DesignationA)& is.na(DesignationB) & is.na(DesignationD)),]
All_merged<-All_merged[,c(3, 2, 1, 4,8, 12, 16, 5, 6)]
colnames(All_merged)<-c("A", "B", "D", "Group", "Taxonomy_A", "Taxonomy_B", "Taxonomy_D", "Cardinality_original", "Cardinality_type")
All_merged$Cardinality_Orphans<-paste(All_merged$Taxonomy_A, ":", All_merged$Taxonomy_B, ":", All_merged$Taxonomy_D, sep="")
table(All_merged$Cardinality_Orphans)
poaceae_triads<-All_merged[(All_merged$Cardinality_Orphans == "Poaceae:Poaceae:Poaceae"),]
pooideae_triads<-All_merged[(All_merged$Cardinality_Orphans == "Pooideae:Pooideae:Pooideae"),]
triticeae_triads<-All_merged[(All_merged$Cardinality_Orphans == "Triticeae:Triticeae:Triticeae"),]
wheat_triads<-All_merged[(All_merged$Cardinality_Orphans == "Wheat:Wheat:Wheat"),]

catagorys<- as.data.frame(data.frame(do.call('rbind', strsplit(as.character(unique(All_merged$Cardinality_Orphans)),':',fixed=TRUE))))
colnames(catagorys)<-c("A", "B", "D")
# next thing to do is change cardianlity column to 1:0:1 to compare use gsub or grep maybe?

catagorys<-catagorys%>% 
  mutate(groupA=case_when(
    A == "Pooeceae" | A == "Pooideae" | A == "Triticeae" | A == "SSG"   ~ 1,
    A == "NA"   ~ 0
  ))
catagorys<-catagorys%>% 
  mutate(groupB=case_when(
    B == "Pooeceae" | B == "Pooideae" | B == "Triticeae" | B == "SSG"   ~ 1,
    B == "NA"   ~ 0
  ))
catagorys<-catagorys%>% 
  mutate(groupD=case_when(
    D == "Pooeceae" | D == "Pooideae" | D == "Triticeae" | D == "SSG"   ~ 1,
    D == "NA"   ~ 0
  ))

catagorys$Cardinality_Orphans<-paste(catagorys$A, ":", catagorys$B, ":", catagorys$D, sep="")
catagorys$Group2<-paste(catagorys$groupA, ":", catagorys$groupB, ":", catagorys$groupD, sep="")
catagorys<-catagorys[,c(7,8)]
All_merged<-left_join(All_merged, catagorys)
All_merged$Match<-All_merged$Cardinality_original == All_merged$Group2
table(All_merged$Match)
write.csv(All_merged, "~/Documents/Sobia/orphans/Homoeologues/merged.csv")

plot(table(homo$Cardinality_type))
plot(table(All_merged$Cardinality_type))


# prepare for hotspot analysis
colnames(A)<-c("Gene", "Designation", "Chromosome","score","Genome")
colnames(B)<-c("Gene", "Designation", "Chromosome","score","Genome")
colnames(D)<-c("Gene", "Designation", "Chromosome","score","Genome")

all<-rbind(A, B, D)
all$TRG<-1
all<-all[,c(1,6)]
write.csv(all, "~/Documents/Sobia/orphans/hotspots/all_orphans.csv")
