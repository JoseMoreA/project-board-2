

# GRM + spatial fit -------------------------------------------------------

df18 <- droplevels(tr$`2018`[tr$`2018`$replicate ==2,])

geno18 <- geno.mat[rownames(geno.mat) %in% 
                     as.character(tr$`2018`$germplasmName),  ]
G_mat <- rrBLUP::A.mat(X = as.matrix(geno18), 
                       min.MAF = .1,shrink = T)
  k=7;cyc=5
  acc.1.all<- list()
  for (c in 1:cyc){
    flds <- caret::createFolds(seq(1:nrow(df18)),
                        k = k, list = TRUE, returnTrain = FALSE)
    # for (e in trait.cv){
      acc.0 <- NULL
      for (f in 1:k){
        blues.cv <- df18
        blues.cv <- droplevels(blues.cv[ which(blues.cv$germplasmName %in% rownames(geno18)), ])
        blues.cv[flds[[f]],trait[4]] <- NA
        blues.cv$germplasmName <- as.factor(blues.cv$germplasmName)
        mm <- mmer(fixed = Yield.LSU_01.0000138 ~ 1,
                   random = ~  vs(germplasmName,Gu=G_mat) + 
                     vs(rowNumber) + vs(colNumber)+
                     vs(spl2D(rowNumber,colNumber)),
                   rcov = ~vs(units),
                   data = blues.cv)
        PRED0 <- data.frame(ID=names(mm$U$`u:germplasmName`[[e ]])
                            ,gblup.trait= 
                              as.numeric(mm$U$`u:germplasmName`[[e ]]))
        # PRED.vp0 <- droplevels(data.frame(ID=blues.cv[flds[[f]],"germplasmName" ],
                                # gblup=PRED0$gblup.trait[PRED0$ID %in% blues.cv[flds[[f]],"germplasmName"] ] ))
        PRED.vp0 <- merge(df18[flds[[f]],c(19,61)],PRED0,by.x =  "germplasmName",by.y = "ID")
        acc.0  <- c(acc.0,cor(PRED.vp0$Yield.LSU_01.0000138,PRED.vp0$gblup.trait,use = "complete.obs"))
      }
      tr.tr <- paste(e,"df18",sep=".")
      acc.1.all[[tr.tr]] <- append(acc.1.all[[tr.tr]],values = acc.0)
    # }
  }
# }