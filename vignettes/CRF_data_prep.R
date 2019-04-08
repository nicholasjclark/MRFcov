## ----eval=FALSE----------------------------------------------------------
#  temp <- tempfile()
#  download.file('https://ndownloader.figshare.com/files/2075362',
#                temp)
#  dataset <- read.csv(temp, as.is = T)
#  unlink(temp)

## ----eval=FALSE----------------------------------------------------------
#  dataset$dipping_round <- as.factor(dataset$dipping_round)
#  dataset$field_site <- as.factor(dataset$field_site)
#  dataset[,c(1,2,5,6)] <- NULL

## ----eval=FALSE----------------------------------------------------------
#  levels(dataset$dipping_round)[1]
#  levels(dataset$field_site)[1]

## ----eval=FALSE----------------------------------------------------------
#  library(dplyr)
#  analysis.data = dataset %>%
#    cbind(.,data.frame(model.matrix(~.[,'field_site'],
#                                    .)[,-1])) %>%
#    cbind(.,data.frame(model.matrix(~.[,'dipping_round'],
#                                    .)[,-1])) %>%
#    dplyr::select(-field_site,-dipping_round) %>%
#    dplyr::rename_all(funs(gsub("\\.|model.matrix", "", .)))
#  

## ----eval=FALSE----------------------------------------------------------
#  analysis.data[, 1:16] <- ifelse(analysis.data[, 1:16] > 0, 1, 0)
#  analysis.data[, 17:20] <- scale(analysis.data[, 17:20], center = T, scale = T)

