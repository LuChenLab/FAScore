# Prediction

library(FAScore)
library(dplyr)

## HD
SPECIES <- c("HomSap", "MusMus")
LINEAGE <- c("Erythroid", "Lymphoid", "Myeloid")

meta <- readRDS("metainfo.Rds")
my_gtf <- list(HomSap = "Homo_sapiens.GRCh38.101.gtf",
               MusMus = "Mus_musculus.GRCm38.101.gtf",
               MacMul = "Macaca_mulatta.Mmul_10.101.gtf",
               RatNor = "Rattus_norvegicus.Rnor_6.0.101.gtf",
               DanRer = "Danio_rerio.GRCz11.93.gtf")
model = readRDS("finalModel.Rds")


### calculate dynamic score
FAScore_HD <- list()
for( s in SPECIES[2]){
  for(l in LINEAGE){
    df_AS <- read.table(paste0(s,"_",l,"_AS.txt"), header = T)
    colnames(df_AS) <- gsub("[.]","-",colnames(df_AS))
    df_Gene <- read.table(paste0(s,"_",l,"_Gene.txt"), header = T)
    colnames(df_Gene) <- gsub("[.]","-",colnames(df_Gene))
    df_Iso <- read.table(paste0(s,"_",l,"_Iso.txt"), header = T)
    colnames(df_Gene) <- gsub("[.]","-",colnames(df_Gene))
    
    df_AS <- df_AS[apply(df_AS, 1, function(x) {sum(x >= 0.05 & x <= 0.95, na.rm = T) >= 2}), ]
    
    MyObj <-FAScoreDataSet(colData = meta[[s]][[l]], 
                           AS = df_AS, 
                           GENE = df_Gene,
                           ISOFORM = df_Iso)
    
    MyObj <- CalcuFeature(MyObj, group.by = "CellType", cores = 10)
    MyObj <- CalcuDyScore(MyObj, maxRange = 10, maxSlope = 1, type = "Gene")
    MyObj <- CalcuDyScore(MyObj, maxRange = 1, maxSlope = 1, type = "AS")
    MyObj <- ASmapIso(MyObj, gtf = my_gtf[[s]], AStype = "exonic", cores = 10)
    MyObj <- ChooseIso(MyObj)
    MyObj <- matchAppris(MyObj, gtf = my_gtf[[s]], species = s)
    MyObj <- CalcuFAScore(MyObj, model = model)
    MyObj <- calcuGMM(MyObj, LOG = 2)
    FAScore_HD[[s]][[l]] <- MyObj
  }
}



## FHO

Species <- c("HomSap", "MacMul", "MusMus", "RatNor", "DanRer")
FHOmeta <- readRDS("FL.meta.Rds")



FAScore_FHO <- list()
for( s in Species){
  df_AS <- read.table(paste0(,s,"_FHO_AS.txt"), header = T)
  colnames(df_AS) <- gsub("[.]","-",colnames(df_AS))
  df_Gene <- read.table(paste0(,s,"_FHO_Gene.txt"), header = T)
  colnames(df_Gene) <- gsub("[.]","-",colnames(df_Gene))
  df_Iso <- read.table(paste0(,s,"_FHO_Iso.txt"), header = T)
  colnames(df_Gene) <- gsub("[.]","-",colnames(df_Gene))
  
  df_AS <- df_AS[apply(df_AS, 1, function(x) {sum(x >= 0.05 & x <= 0.95, na.rm = T) >= 2}), ]
  
  MyObj <-FAScoreDataSet(colData = FHOmeta[[s]], 
                         AS = df_AS, 
                         GENE = df_Gene,
                         ISOFORM = df_Iso)
  
  MyObj <- CalcuFeature(MyObj, group.by = "StageAbbreviation", cores = 10)
  MyObj <- CalcuDyScore(MyObj, maxRange = 10, maxSlope = 1, type = "Gene")
  MyObj <- CalcuDyScore(MyObj, maxRange = 1, maxSlope = 1, type = "AS")
  MyObj <- ASmapIso(MyObj, gtf = my_gtf[[s]], AStype = "exonic", cores = 10)
  MyObj <- ChooseIso(MyObj)
  MyObj <- matchAppris(MyObj, gtf = my_gtf[[s]], species = s)
  MyObj <- CalcuFAScore(MyObj, model = model)
  MyObj <- calcuGMM(MyObj, LOG = 2)
  FAScore_FHO[[s]] <- MyObj
}





saveRDS(FAScore_HD, file = "FAScore_HD.Rds")
saveRDS(FAScore_FHO, file = "FAScore_FHO.Rds")








