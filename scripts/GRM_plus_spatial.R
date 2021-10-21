
# GRM + spatial fit -------------------------------------------------------
pheno_train <- read.csv("../data/2021-09-17T174543phenotype_download.csv",header = T)
geno.mat <- readRDS("../data/geno_mat.rds")
list.tr <- split(pheno_train,pheno_train$studyYear)
traits <- colnames(pheno_train)[c(42,51,52,61)]

acc.all<- list()
for(df in names(list.tr)){
  list.tr.cv <- list.tr
  list.tr.cv[[df]] <- droplevels(list.tr.cv[[df]][list.tr.cv[[df]]$replicate ==2,])
  geno.mat.cv <- geno.mat[rownames(geno.mat) %in% 
                       as.character(list.tr.cv[[df]]$germplasmName),  ]
  G_mat <- rrBLUP::A.mat(X = as.matrix(geno.mat.cv), 
                         min.MAF = .1,shrink = T)
  k=7;cyc=5
  for (c in 1:cyc){
    flds <- caret::createFolds(seq(1:nrow(list.tr.cv[[df]])),
                        k = k, list = TRUE, returnTrain = FALSE)
    for (e in traits[4]){
      acc.k <- NULL
      for (f in 1:k){
        df.tr <- list.tr.cv[[df]]
        df.tr <- droplevels(df.tr[ which(df.tr$germplasmName %in% rownames(geno.mat.cv)), ])
        df.tr[flds[[f]],e] <- NA
        df.tr$germplasmName <- as.factor(df.tr$germplasmName)
        frm.cv <- paste(e,"~ 1 ",collapse = " ")
        fix_frm.cv <- formula(paste(frm.cv, collapse = " "))
        mm <- mmer(fixed = fix_frm.cv,
                   random = ~  vs(germplasmName,Gu=G_mat) + 
                     vs(rowNumber) + vs(colNumber)+
                     vs(spl2D(rowNumber,colNumber)),
                   rcov = ~vs(units),
                   data = df.tr)
        PRED0 <- data.frame(ID=names(mm$U$`u:germplasmName`[[e ]])
                            ,gebv=as.numeric(mm$U$`u:germplasmName`[[e ]]))
        PRED.vp0 <- merge(list.tr.cv[[df]][flds[[f]],c(19,61)],PRED0,by.x =  "germplasmName",by.y = "ID")
        acc.k  <- c(acc.0,cor(PRED.vp0[,e],
                              PRED.vp0[,"gebv"],use = "complete.obs"))
      }
      tr.tr <- paste(e,df,sep=".")
      acc.all[[tr.tr]] <- append(acc.all[[tr.tr]],values = acc.k)
    }
  }
}