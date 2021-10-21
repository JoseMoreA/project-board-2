
# RCBD model --------------------------------------------------------------
pheno_train <- read.csv("../data/2021-09-17T174543phenotype_download.csv",header = T)
df.tr <- split(pheno_train,pheno_train$studyYear)
geno.ID.tr <- colnames(pheno_train)[62]
traits.tr <- colnames(pheno_train)[c(42,51,52,61)]
blues <- list()
for (h in names(df.tr)){
  for (i in traits.tr){
    frm <- paste(i,"~ 0 +", geno.ID.tr)
    fix_frm <- formula(paste(frm, collapse = " "))
    mod.name <- paste("model",h,i,sep="_")
    df.tr[[h]][,geno.ID.tr] <- as.factor(df.tr[[h]][,geno.ID.tr])
    mm1 <- mmer(fixed = fix_frm,
                random = ~  replicate  ,
                data = df.tr[[h]])
    blues[[mod.name]] <- mm1$Beta
    blues[[mod.name]]$Effect <-
      substr(blues[[mod.name]]$Effect,
             start = 14, stop=1000000L) 
  }
}


# Cross validation --------------------------------------------------------


geno.mat <- readRDS("../../Downloads/geno_mat.rds")
# blues <- blues2 # BLUEs from RCBD model
traits <- colnames(blues[[1]])[-1]
geno.ID <- "Effect"

for(df in names(blues)){
  blues.cv <- blues
  blues.cv[[df]] <- droplevels(blues.cv[[df]][ which(blues.cv[[df]][,geno.ID] %in% rownames(geno.mat)), ])
  geno.mat.cv <- geno.mat[rownames(geno.mat) %in% 
                         as.character(blues.cv[[df]][,geno.ID]),  ]
  G_mat <- rrBLUP::A.mat(X = as.matrix(geno.mat.cv), 
                         min.MAF = .1,shrink = T)
  k=7;cyc=50
  for (c in 1:cyc){
    flds <- createFolds(seq(1:nrow(blues.cv[[df]])),
                        k = k, list = TRUE, returnTrain = FALSE)
    for (e in traits[4]){
      acc.0 <- NULL
      for (f in 1:k){
        blues.cv <- blues
        blues.cv[[df]] <- droplevels(blues.cv[[df]][ which(blues.cv[[df]][,geno.ID] %in% rownames(geno_df3)), ])
        blues.cv[[df]][flds[[f]],e] <- NA
        blues.cv[[df]][,geno.ID] <- as.factor(blues.cv[[df]][,geno.ID])
        frm.cv <- paste(e,"~ 1 ",collapse = " ")
        fix_frm.cv <- formula(paste(frm.cv, collapse = " "))
        mm <- mmer(fixed = fix_frm.cv,
                   random = ~vs(Effect,Gu=G_mat),
                   rcov = ~vs(units),
                   data = blues.cv[[df]])
        PRED0 <- data.frame(ID=names(mm$U$`u:Effect`[[e ]])
                            ,gebv=as.numeric(mm$U$`u:Effect`[[e ]]))
        PRED.vp0 <- merge(blues2[[df]][flds[[f]],c(1,5)],PRED0,by.x =  geno.ID,by.y = "ID") 
        acc.0  <- c(acc.0,cor(PRED.vp0[,e],
                              PRED.vp0[,"gebv"],use = "complete.obs"))
      }
      tr.tr <- paste(e,df,sep=".")
      acc.1.all[[tr.tr]] <- append(acc.1.all[[tr.tr]],values = acc.0)
    }
  }
}