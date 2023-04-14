# 1) Single protein Cox models for dementia risk (visits 3-6)
# Full follow-up
# Long-term follow-up (>=15 yr)
# Near-term follow-up (<15 yr)
ds<-c("soma_v3_SMP_log2",
      "soma_v3_SMP_log2_fu15less","soma_v3_SMP_log2_fu15more")
dss<-c("",
       "_fu15less","_fu15more")


for(j in 1:length(dss)) {
  
  d=get(ds[j])
  
  if(ds[j] %in% c("soma_v3_SMP_log2")) {
    Model<-c("",
             "+age+female+racecen+edu+apoe",
             "+age+female+racecen+edu+apoe+egfrcr",
             "+age+female+racecen+edu+apoe+egfrcr+cursmk+bmi+diabetes+hypertens")
    model<-c(0,1,2,3)
    event<-c("event_dementia","event_dementia_end_v5")
    fu<-c("fu_dementia","fu_dementia_end_v5")
    outcome<-c("Incident Dementia","Incident Dementia censored at v5")
  } else if(ds[j] %in% c("soma_v3_SMP_log2_fu15less","soma_v3_SMP_log2_fu15more")) {
    Model<-c("+age+female+racecen+edu+apoe+egfrcr+cursmk+bmi+diabetes+hypertens")
    model<-c(3)
    event<-c("event_dementia")
    fu<-c("fu_dementia")
    outcome<-c("Incident Dementia")
  }
  
  
  for(k in 1:length(outcome)) {  
    for(m in 1:length(Model)) {
      h<-data.frame(seq_id=character(),beta=double(),se=double(),p=double(),HR=double(),
                    lci=double(),uci=double(),stringsAsFactors=FALSE)
      
      for(i in 1:length(seq_ids)) {
        h[i,"seq_id"]<-seq_ids[i]
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",seq_ids[i],Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        h[i,"beta"]<-summary(f)$coefficient[1,1]
        h[i,"se"]<-summary(f)$coefficient[1,3]
        h[i,"p"]<-summary(f)$coefficient[1,5]
        h[i,"HR"]<-summary(f)$coefficient[1,2]
        h[i,"lci"]<-exp(confint(f)[1,1])
        h[i,"uci"]<-exp(confint(f)[1,2])
      }
      
      h$p_FDR<-p.adjust(h$p, method = "fdr")
      h<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","targetfullname")],h,by.x="seqid_in_sample",by.y="seq_id")
      h<-h[order(h$p),]
      write.table(h,file=paste(out.dir,"cox-incident dementia//Cox_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}





# 2) Netboost module eigengene Cox models for dementia risk (visits 3-6)
# Full follow-up
# Long-term follow-up (>=15 yr)
# Near-term follow-up (<15 yr)
ds<-c("MEs",
      "MEs_fu15less","MEs_fu15more")
dss<-c("",
       "_fu15less","_fu15more")

for(j in 1:length(dss)) {
  
  d=get(ds[j])
  
  if(ds[j] %in% c("MEs")) {
    Model<-c("",
             "+age+female+racecen+edu+apoe",
             "+age+female+racecen+edu+apoe+egfrcr",
             "+age+female+racecen+edu+apoe+egfrcr+cursmk+bmi+diabetes+hypertens")
    model<-c(0,1,2,3)
    event<-c("event_dementia","event_dementia_end_v5")
    fu<-c("fu_dementia","fu_dementia_end_v5")
    outcome<-c("Incident Dementia","Incident Dementia censored at v5")
  } else if(ds[j] %in% c("MEs_fu15less","MEs_fu15more")) {
    Model<-c("+age+female+racecen+edu+apoe+egfrcr+cursmk+bmi+diabetes+hypertens")
    model<-c(3)
    event<-c("event_dementia")
    fu<-c("fu_dementia")
    outcome<-c("Incident Dementia")
  }
  
  
  for(k in 1:length(outcome)) {  
    for(m in 1:length(Model)) {
      h<-data.frame(ME=character(),beta=double(),se=double(),p=double(),HR=double(),
                    lci=double(),uci=double(),stringsAsFactors=FALSE)
      
      for(i in 1:length(MEvalues)) {
        h[i,"ME"]<-MEvalues[i]
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~",MEvalues[i],Model[m],sep=""))
        f<-coxph(fmla,data=d, na.action=na.exclude)
        h[i,"beta"]<-summary(f)$coefficient[1,1]
        h[i,"se"]<-summary(f)$coefficient[1,3]
        h[i,"p"]<-summary(f)$coefficient[1,5]
        h[i,"HR"]<-summary(f)$coefficient[1,2]
        h[i,"lci"]<-exp(confint(f)[1,1])
        h[i,"uci"]<-exp(confint(f)[1,2])
      }
      
      h$p_FDR<-p.adjust(h$p, method = "fdr")
      h<-h[order(h$p),]
      write.table(h,file=paste(out.dir,"cox-incident dementia//Cox_",outcome[k],"_Model",model[m],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    }
  }
}




# 3) Elastic net dementia prediction (visits 3-6)
# Full follow-up
# Long-term follow-up (>=15 yr)
# Near-term follow-up (<15 yr)
### Elastic net with cox regression
elasticNetCox <- function(candiateVars, forcedVars = forced_in_Vars, data = d.temp,
                          endpoint, time){
  
  if(length(candiateVars)>0 & length(c(candiateVars,forcedVars))>=2) {
    ## Cross validation for hyperparamter tuning
    set.seed(2524)
    cv.fit <- cv.glmnet(data.matrix(data[, c(forcedVars, candiateVars)]),
                        Surv(time, endpoint), family="cox", alpha = 0.5, maxit = 2000, nfolds=5,
                        penalty.factor=c(rep(0, length(forcedVars)), rep(1, length(candiateVars))))
    ## Fit the model
    fit <- glmnet(as.matrix(data[, c(forcedVars, candiateVars)]),
                  Surv(time, endpoint), family="cox", alpha = 0.5, maxit = 2000,
                  penalty.factor=c(rep(0, length(forcedVars)), rep(1, length(candiateVars))))
    #browser()
    Coefficients <- coef(fit, s = cv.fit$lambda.1se)
    nonzeorcoefs <- Coefficients[which(Coefficients != 0)]
    nonzeor_protein = names(nonzeorcoefs) = rownames(Coefficients)[ which(Coefficients != 0)]
    #nonzeor_flags = sapply(nonzeor_protein, function(x) flag_status[which(names(flag_status) == x)] )
    #nonzeor_names = sapply(nonzeor_protein, function(x) names(proteList_all)[which(proteList_all == x)] )
    SurvObj<-paste("Surv(",fu[k],",",event[k],")~",sep="")
    if(length(nonzeorcoefs)>1) {
      cox_formula <- as.formula(paste(SurvObj, paste(c(nonzeor_protein, forcedVars), collapse= "+") ))
    } else {
      cox_formula <- as.formula(paste(SurvObj, "~1"))
    }
    final_cox <- coxph(cox_formula, data = data, method = "efron")
    # Base Model
    hr_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 1]
    se_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 3]
    p_fM = summary(final_cox)$coefficients[1:length(nonzeor_protein), 5]
    hr_fMs = paste0(round(hr_fM, 2), " (", round(se_fM, 2), ")" )
    p_fMs = sapply(formatC(p_fM, format = "e", digits = 1), toString)
    
    elsNet_results = cbind(unlist(nonzeor_protein), 
                           round(exp(nonzeorcoefs), 2), round(nonzeorcoefs, 2),
                           hr_fMs, p_fMs) %>% data.table()
    names(elsNet_results) <-c("seq_id", "HR(elasticNet)", "logHR(elasticNet)", "logHR(Cox)", "Pvalue(Cox)")
    return(list(elsNet_results = elsNet_results, optimalNet = Coefficients, finalCox = final_cox,fit=fit,cv.fit=cv.fit,flagx=0))
  } else if(length(forcedVars)!=0) {
    SurvObj<-paste("Surv(",fu[k],",",event[k],")~",sep="")
    cox_formula <- as.formula(paste(SurvObj, paste(c(forcedVars), collapse= "+") ))
    final_cox <- coxph(cox_formula, data = data, method = "efron")
    # Base Model
    hr_fM = summary(final_cox)$coefficients[1:length(forcedVars), 1]
    se_fM = summary(final_cox)$coefficients[1:length(forcedVars), 3]
    p_fM = summary(final_cox)$coefficients[1:length(forcedVars), 5]
    hr_fMs = paste0(round(hr_fM, 2), " (", round(se_fM, 2), ")" )
    p_fMs = sapply(formatC(p_fM, format = "e", digits = 1), toString)
    
    a<-NA
    b<-NA
    elsNet_results = cbind(unlist(forcedVars),a,b,hr_fMs, p_fMs) %>% data.table()
    names(elsNet_results) <-c("seq_id", "HR(elasticNet)", "logHR(elasticNet)", "logHR(Cox)", "Pvalue(Cox)")
    return(list(elsNet_results = elsNet_results, finalCox = final_cox, flagx=0))
  } else {
    return(list(elsNet_results = NA,  finalCox = NA,flagx=1))
  }
}


get.risk.coxph.ex <-
  function(mdl, t0, lp) {
    bash    = basehaz(mdl)
    lambda0 = approx(bash$time, bash$hazard, t0)$y
    risk    = 1 - exp( - lambda0 * exp( lp ) )
    return(risk)
  }




###################################################################################
## naive models
ds<-c("soma_v3_SMP_log2","soma_v3_SMP_log2_fu15less","soma_v3_SMP_log2_fu15more")
dss<-c("_Full-Term Dementia","_Near-Term Dementia","_Long-Term Dementia")
event<-c("event_dementia")
fu<-c("fu_dementia")
outcome<-c("Incident Dementia")
E_covar3=c("age","female","racecenF_B","racecenF_W","racecenJ_B","racecenW_W","edu","apoe_gt1","apoe_missing","egfrcr","cursmk","bmi","diabetes","hypertens")
EN<-c("EN0","EN3CF","M3F")
sheet=c("Full-Term Dementia","Near-Term Dementia 0-15","Long-Term Dementia >15")


for(j in 1:length(ds)) {
  d.temp<-get(ds[j])
  
  for(k in 1:length(outcome)) {
    
    en<-NULL
    h<-data.frame(model=character(),c_stat=double(),se_c=double(),stringsAsFactors=FALSE)
    
    for(n in 1:length(EN)) {
      
      if (grepl("EN",EN[n])) {
        num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
        t<-read_excel(path = paste(out.dir,"Midlife Prot. of Dementia - Netboost Prediction - Prot. Lists.xlsx",sep=""), sheet=sheet[j])
        tops=t$seqid_in_sample
        print(length(tops))
        vars=tops
      } else {
        vars=NULL
      }
      
      if (grepl("C",EN[n]) | grepl("M",EN[n])) {
        num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
        C=get(paste("E_covar",num,sep=""))
      } else {
        C=NULL
      }
      
      if (grepl("F",EN[n])) {
        forced_in_Vars=C
      } else {
        forced_in_Vars=NULL
      }
      
      if (!grepl("F",EN[n]) & (grepl("C",EN[n]) | grepl("M",EN[n]))) {
        vars=c(vars,C)
      }
      
      print("spot0")
      elsNet_results_all = elasticNetCox(candiateVars=vars, forcedVars = forced_in_Vars, data = d.temp,
                                         endpoint = d.temp[,event[k]], time = d.temp[,fu[k]])
      
      if(elsNet_results_all$flagx==0) {
        print("spot1")
        e<-elsNet_results_all$elsNet_results
        e$model<-EN[n]
        en<-rbind(en,e)
        print("spot2")
        # c-stat
        h[n,"model"]<-EN[n]
        h[n,"c_stat"]<-summary(elsNet_results_all$finalCox)$concordance[1]
        h[n,"se_c"]<-summary(elsNet_results_all$finalCox)$concordance[2]
        print("spot3")
        # print model
        print(EN[n])
        print(summary(elsNet_results_all$finalCox))
        print("spot4")
      } else {
        h[n,"model"]<-EN[n]
      }
    }
    
    en$n<-1:nrow(en)
    en<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","target","targetfullname")],en,by.x="seqid_in_sample",by.y="seq_id",all.y=TRUE)
    en<-en[order(en$n),]
    en<-en[ , -which(names(en) %in% c("n"))]
    write.table(en,file=paste(out.dir,"Elastic Network/cox/naive/Elastic Network_Cox (naive)_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    write.table(h,file=paste(out.dir,"Elastic Network/cox/naive/Elastic Network_Cox_cstat (naive)_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}



##################################################################################
## 10-fold cross validation
# cross-validation for cox
ds<-c("soma_v3_SMP_log2","soma_v3_SMP_log2_fu15less","soma_v3_SMP_log2_fu15more")
dss<-c("_Full-Term","_Near-Term","_Long-Term")
event<-c("event_dementia")
fu<-c("fu_dementia")
outcome<-c("Incident Dementia")
E_covar3=c("age","female","racecenF_B","racecenF_W","racecenJ_B","racecenW_W","edu","apoe_gt1","apoe_missing","egfrcr","cursmk","bmi","diabetes","hypertens")
EN<-c("EN0","EN3CF","M3F")
sheet=c("Full-Term Dementia","Near-Term Dementia 0-15","Long-Term Dementia >15")


for(j in 1:length(ds)) {
  d.temp<-get(ds[j])
  
  
  for(k in 1:length(outcome)) {
    
    c_validation<-data.frame(model=character(),c_predall=double(),se_predall=double(),sum_c1=double(),sum_c2=double(),
                             sum_c3=double(),sum_c4=double(),sum_c5=double(),sum_c6=double(),
                             sum_c7=double(),sum_c8=double(),sum_c9=double(),sum_c10=double(),
                             stringsAsFactors=FALSE)
    
    
    index<-d.temp[,event[k]]==1
    d.temp.e=d.temp[index,]
    d.temp.ne=d.temp[!index,]
    
    #set.seed(2524)
    set.seed(12345)
    folds.e <- split(sample(nrow(d.temp.e), nrow(d.temp.e),replace=FALSE), as.factor(1:10))
    #set.seed(9527)
    set.seed(23333)
    folds.ne <- split(sample(nrow(d.temp.ne), nrow(d.temp.ne),replace=FALSE), as.factor(1:10))
    
    d.temp2<-d.temp
    
    
    for(n in 1:length(EN)) {
      
      new.pred.all<-NULL
      cv_var_coef<-NULL
      
      
      for(q in 1:10) {
        
        if (grepl("EN",EN[n])) {
          num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
          t<-read_excel(path = paste(out.dir,"Midlife Prot. of Dementia - Netboost Prediction - Prot. Lists.xlsx",sep=""), sheet=sheet[j])
          tops=t$seqid_in_sample
          print(length(tops))
          vars=tops
        } else {
          vars=NULL
        }
        
        if (grepl("C",EN[n]) | grepl("M",EN[n])) {
          num=str_extract(EN[n], "\\-*\\d+\\.*\\d*")
          C=get(paste("E_covar",num,sep=""))
        } else {
          C=NULL
        }
        
        if (grepl("F",EN[n])) {
          forced_in_Vars=C
        } else {
          forced_in_Vars=NULL
        }
        
        if (!grepl("F",EN[n]) & (grepl("C",EN[n]) | grepl("M",EN[n]))) {
          vars=c(vars,C)
        }
        
        
        d1<-rbind.data.frame(d.temp.e[folds.e[[q]],],d.temp.ne[folds.ne[[q]],])
        d2<-rbind.data.frame(d.temp.e[-folds.e[[q]],],d.temp.ne[-folds.ne[[q]],])
        
        print(paste(EN[n],"__",q,sep=""))
        
        elsNet_results_part =  elasticNetCox(candiateVars=vars, forcedVars = forced_in_Vars, data = d2,
                                             endpoint = d2[,event[k]], time = d2[,fu[k]])
        
        if(elsNet_results_part$flagx==0) {
          
          print("spot1")
          
          f<-elsNet_results_part$finalCox
          
          print("spot2")
          
          new.pred = predict(f, d1)
          new.pred1=data.frame(new.pred)
          rownames(new.pred1)=d1$SampleId
          
          print("spot3")
          
          new.pred.all<-rbind.data.frame(new.pred.all,new.pred1)
          print(dim(d1))
          print(length(new.pred))
          print(dim(new.pred1))
          print(dim(new.pred.all))
          
          print("spot4")
          
          fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~ new.pred",sep=""))
          
          print("spot5")
          
          res.1 = survConcordance(fmla, d1)
          
          print("spot6")
          
          c_validation[n,"model"]=EN[n]
          c_validation[n,paste("sum_c",q,sep="")]=res.1$concordance
          
          
          # save coefficients
          coef<-data.frame(coef(elsNet_results_part$finalCox))
          colnames(coef)<-paste("coef_loop",q,sep="")
          coef$Var<-row.names(coef)
          if (q==1 | is.null(cv_var_coef)) {
            cv_var_coef=coef
          } else if(q>1) {
            cv_var_coef<-merge(cv_var_coef,coef,by="Var",all=TRUE)
          }
          
        } else {
          c_validation[n,"model"]=EN[n]
          if(q==1) {
            cv_var_coef=NULL
          }
        }
        
      }
      
      
      if(!is.null(new.pred.all))  {
        print(dim(d.temp2))
        names(new.pred.all)=paste0("new.pred_",EN[n])
        new.pred.all$SampleId=rownames(new.pred.all)
        d.temp2<-merge(d.temp2,new.pred.all,by="SampleId",all=T)
        print(dim(d.temp2))
        
        fmla<-as.formula(paste("Surv(",fu[k],",",event[k],")~ new.pred_",EN[n],sep=""))
        res.2 = survConcordance(fmla, d.temp2)
        c_validation[n,"c_predall"]=res.2$concordance
        c_validation[n,"se_predall"]=res.2$std.err
      }
      
      if(!is.null(cv_var_coef)) {
        coefs<-grep("coef_loop", names(cv_var_coef), value=TRUE)
        cv_var_coef$N_select=rowSums(!is.na(cv_var_coef[,coefs]))
        cv_var_coef<-merge(annot[,c("seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","target","targetfullname")],cv_var_coef,by.x="seqid_in_sample",by.y="Var",all.y=TRUE)
        write.table(cv_var_coef,file=paste(out.dir,"Elastic Network/cox/CV (10 fold)/Elastic net_cox_Validation_coefficients_",EN[n],"_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
      }
    }
    
    
    # C-stat compare
    cstat_compare<-data.frame(ModelA=character(),ModelB=character(),p_cstat=character(),
                              stringsAsFactors=FALSE)
    
    lc1=c("M3F")
    lc2=c("EN3CF")
    for(c in 1:length(lc1)) {
      c1=lc1[c]
      c2=lc2[c]
      compare<-compareC(d.temp2[,fu[k]], d.temp2[,event[k]], d.temp2[,paste0("new.pred_",c1)],  d.temp2[,paste0("new.pred_",c2)])
      cstat_compare[c,"ModelA"]=lc1[c]
      cstat_compare[c,"ModelB"]=lc2[c]
      cstat_compare[c,"p_cstat"]=compare$pval
    }
    write.table(cstat_compare,file=paste(out.dir,"Elastic Network/cox/CV (10 fold)/Elastic net_cox_cstat_compare_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
    
    
    sumCs<-grep("sum_c", names(c_validation), value=TRUE) 
    c_validation$mean_c=rowMeans(c_validation[,sumCs])
    #c_validation$se_c=rowSds(as.matrix(c_validation[,sumCs]))/sqrt(10)
    write.table(c_validation,file=paste(out.dir,"Elastic Network/cox/CV (10 fold)/Elastic net_cox_Validation_Cstat_",outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}



### counting how many times selected
for(j in 1:length(ds)) {
  for(k in 1:length(outcome)) {
    
    cv_var_coef_all<-NULL
    
    for(n in 1:length(EN)) {
      cv_var_coef=read.delim(file = paste(out.dir,"Elastic Network/cox/CV (10 fold)/Elastic net_cox_Validation_coefficients_", EN[n],"_",outcome[k],dss[j],".txt",sep=""), header = TRUE, sep = "\t", stringsAsFactors =FALSE)
      cv_var_coef<-cv_var_coef[order(-cv_var_coef$N_select),]
      cv_var_coef$Model<-EN[n]
      cv_var_coef<-cv_var_coef[cv_var_coef$N_select>=7,c("Model","seqid_in_sample","uniprot_id","entrezgeneid","entrezgenesymbol","target","targetfullname","N_select")]
      
      if(n==1) {
        cv_var_coef_all<-cv_var_coef
      } else {
        cv_var_coef_all<-rbind.data.frame(cv_var_coef_all,cv_var_coef)
      }
    }
    
    write.table(cv_var_coef_all,file=paste(out.dir,"Elastic Network/cox/CV (10 fold)/Elastic net_variables selected 7 times or more_", outcome[k],dss[j],".txt",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
  }
}
