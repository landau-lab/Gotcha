
MitoClassifier <- function(mitomutations){

   wholedataset <- mitomutations %>%
     dplyr::filter(perc_filter=="pass") %>%
     dplyr::filter(wilcox=="informative") %>%
     dplyr::select(MutationCall,OriginalWhiteListMatch,Site,Heteroplasmy) %>%
     tidyr::pivot_wider(names_from=Site,values_from=Heteroplasmy, values_fill=NA)
     #!!!!! Select the number of variants I want to use

   mydata_genotyped <- wholedataset[complete.cases(wholedataset),] %>% # remove rows with missing data
     dplyr::select(-OriginalWhiteListMatch) %>%
     dplyr::filter(MutationCall!="Non-genotyped" & MutationCall!="HET") %>%
     dplyr::mutate(MutationCall=factor(MutationCall, levels=c("MUT","WT")))

   # !!!!!!!!!!!!!!!!!!!!!!!
   # take the same number of MUT than WT cells
   number_wt <- table(mydata_genotyped$MutationCall)["WT"]
   mydata_genotyped_mut <- mydata_genotyped %>%
     dplyr::filter(MutationCall=="MUT")
   mydata_genotyped_wt <- mydata_genotyped %>%
     dplyr::filter(MutationCall=="WT")
   mydata_genotyped <- rbind.data.frame(mydata_genotyped_mut,mydata_genotyped_wt)
   category_withless <- names(which.min(table(mydata_genotyped$MutationCall)))
   min_cells <- table(mydata_genotyped$MutationCall)[category_withless]
   mydata_genotyped <- rbind.data.frame(mydata_genotyped_mut[1:min_cells,],mydata_genotyped_wt[1:min_cells,])

   #### Create the dataset for which I want to predict the genotypes
   mydata_nongenotyped <- wholedataset[complete.cases(wholedataset),] %>%
     dplyr::filter(MutationCall=="Non-genotyped")

   #### Sampling
   #!!!!!! prop 75% is the default but it can be other sizes
   mydata_split <- initial_split(mydata_genotyped, prop = 0.75)
   recipe <- mydata_split %>%
     training() %>%
     recipe(MutationCall ~.) %>%
     prep()
   # apply same recipe to the testing data
   testing_set <- recipe %>%
     bake(testing(mydata_split))
   # Performing the same operation over the training data is redundant,
   # because that data has already been prepped. To load the prepared training data into a variable, we use juice()
   training_set <- juice(recipe)
   glimpse(training_set)


   # Build the random forest model
   mydata_randomforest <-  rand_forest(trees = 1000L, mode = "classification") %>%
     set_engine("randomForest") %>% #     set_engine("ranger") %>%
     fit(MutationCall ~ ., data = training_set) # MutationCall vs everything else (.)

   # Use the model to make genotyping predictions and assess accuracy  #

   myPROB = 0.5

   PredictAndComputeAccuracy <- function(myPROB){

     metrics <- mydata_randomforest %>%
       predict(testing_set) %>%
       cbind.data.frame(predict(mydata_randomforest, testing_set, type="prob")) %>%
       bind_cols(testing_set) %>%
       # !!!!! remove those cells with low assignment probability
       dplyr::filter(`.pred_MUT` > myPROB | `.pred_WT` > myPROB) %>%
       dplyr::select(-`.pred_MUT`,-`.pred_WT`) %>%
       metrics(truth = MutationCall, estimate = .pred_class)

     accuracy <- as.character(metrics[1,3])
     kappa <- as.character(metrics[2,3])

     # Predict in the non-genotyped data
     mydata_nongenotyped_predicted <- mydata_randomforest %>%
       predict(mydata_nongenotyped) %>%
       cbind.data.frame(predict(mydata_randomforest, mydata_nongenotyped, type="prob")) %>%
       cbind.data.frame(mydata_nongenotyped) %>%
       dplyr::select(-contains(">"),-MutationCall) %>%
       dplyr::rename(Predicted=.pred_class) %>%
       dplyr::rename(mut_prob=.pred_MUT) %>%
       dplyr::rename(wt_prob=.pred_WT) %>%
       dplyr::mutate(Predicted=ifelse(as.numeric(mut_prob)<myPROB & as.numeric(wt_prob)<myPROB,paste("Lowconf",Predicted,sep=" "),as.character(Predicted)))

     # Get the whole dataset
     # Add the gotcha genotyped cells (training+testing) to the predictions
     wholedataset_withpredictions <- wholedataset %>%
       dplyr::left_join(mydata_nongenotyped_predicted) %>%
       #dplyr::select(MutationCall,OriginalWhiteListMatch,Predicted) %>%
       dplyr::select(-contains(">")) %>%
       dplyr::mutate(Classifier=MutationCall) %>%
       dplyr::mutate(Predicted=as.character(Predicted)) %>%
       dplyr::mutate(Classifier=if_else(is.na(Predicted),Classifier,Predicted))


     # Calculate genotyping efficiency

     # Efficiency GOTCHA genotyping
     efficiency_gotcha <- table(wholedataset_withpredictions$MutationCall) %>% as.list()
     efficiency_gotcha_result <- (as.numeric(efficiency_gotcha["MUT"])+as.numeric(efficiency_gotcha["WT"]))/sum(as.numeric(efficiency_gotcha))
     # 18.94%

     # Efficiency classifier
     efficiency_gotchaplusmito <- table(wholedataset_withpredictions$Classifier) %>% as.list()
     efficiency_gotchaplusmito_result <- (as.numeric(efficiency_gotchaplusmito["MUT"])+as.numeric(efficiency_gotchaplusmito["WT"]))/sum(as.numeric(efficiency_gotchaplusmito))
     # 86.55%

     out <- wholedataset_withpredictions %>%
       dplyr::mutate(probability_threshold=myPROB) %>%
       dplyr::mutate(genotyping_efficiency_gotcha=efficiency_gotcha_result) %>%
       dplyr::mutate(genotyping_efficiency_classifier=efficiency_gotchaplusmito_result) %>%
       dplyr::mutate(genotyping_accuracy_classifier=accuracy) %>%
       dplyr::mutate(genotyping_kappa_classifier=kappa)

     return(out)

   }

 #  tables_classifier <- lapply(seq(0.5,1,0.01),PredictAndComputeAccuracy)
 #  alltables_classifier <- do.call(rbind,tables_classifier)
 #
 #
 #  alltables_classifier[,7:ncol(alltables_classifier)] %>%
 #    unique() %>%
 #    ggplot(aes(as.numeric(probability_threshold),as.numeric(genotyping_accuracy_classifier), col=genotyping_efficiency_classifier)) +
 #    geom_point() +
 #    theme_classic() +
 #    labs(x="Genotype probability\nclassifier confidence",y="Genotype accuracy") +
 #    scale_color_distiller(direction = "-1") +
 #    ylim(0.97,1)
 #
 #
 #  wholedataset_withpredictions <- PredictAndComputeAccuracy(myPROB = 0.5) %>%
 #    dplyr::mutate(Predicted=ifelse(is.na(Predicted),"Non-genotyped",Predicted))
 #
 #  # barplot
 #
 #  gotcha_efficiency <- round(wholedataset_withpredictions$genotyping_efficiency_gotcha[1]*100,2)
 #  mito_efficiency <- round(wholedataset_withpredictions$genotyping_efficiency_classifier[1]*100,2)
 #
 #
 #  gotcha <- table(wholedataset_withpredictions$MutationCall)
 #  mitoclassifier <- table(wholedataset_withpredictions$Classifier)
 #
 #  barplotinput <- rbind(gotcha,mitoclassifier) %>%
 #    as.data.frame() %>%
 #    tibble::rownames_to_column(var = "technique") %>%
 #    dplyr::mutate(efficiency=ifelse(technique=="gotcha",gotcha_efficiency,mito_efficiency )) %>%
 #    tidyr::pivot_longer(cols = c("HET","MUT","Non-genotyped","WT")) %>%
 #    dplyr::rename(number=value) %>%
 #    dplyr::rename(genotype=name)  %>%
 #    #dplyr::mutate(genotype=as.factor(genotype)) %>%
 #    dplyr::mutate(genotype= gsub("HET","ambiguous",gsub("MUT","mutant",gsub("WT","wild type",gsub("Non-genotyped",NA,genotype))))) %>%
 #    dplyr::mutate(technique=gsub("gotcha","gotcha",gsub("mitoclassifier","gotcha + mito", technique)))
 #
 #
 #  barplot_ggplot <- ggplot(barplotinput, aes(technique, number, fill=genotype)) +
 #    geom_col(width = 0.5) +
 #    scale_fill_manual(values = c("darkorchid4","firebrick2","cornflowerblue"),na.value="grey") +
 #    labs(x="",y="Number of cells", fill="Genotype") +
 #    geom_text(aes(label = efficiency), size = 3, hjust = 0.5, vjust = 2, position = "stack") +
 #    theme_classic()


   #ggsave(barplot_ggplot,
   #       file = paste(plotdirshare,"Barplot.gotcha.mito.genoefficiency.pdf",sep=""), dpi = 300)

}
