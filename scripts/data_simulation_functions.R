
run_mech_level_simulation <- function(input_data = '', simulation_runs = '', variable_effect = '') {

  for( i in 1:simulation_runs) {
    settings = zigControl(maxit=1000,verbose=TRUE)
    vector_sig_diff <-c()
    base_data <- input_data
    for( f in 1:length(variable_effect )){
      var_effect <- variable_effect[f]
      for( v in sample_number_variable ) { ## This function is what makes the counts for each simulated dataset
        temp_datasets <- list()
        ## Make treatment groups and model matrix
        trt_group <- data.frame(matrix(ncol = 0, nrow = v))
        trt_group$Time <-rep(c("Arrival","Exit"),each= v/2) ## make treatment group metadata
        trt_group$Treatment <- rep(c("Ctrl","Trt"),each= v/4)
        trt_group$TT <- paste(trt_group$Time,trt_group$Treatment,sep="")
        mod = model.matrix(~ 0 + TT, trt_group)
        ### different data simulation methods, based on distribution
        if(distribution_type == 'negbinom_mu_data_simulation') {
          if( effect_type == 'fixed_sd' ) {
            temp_datasets <- negbinom_mu_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else if( effect_type == 'variable_sd' ) {
            temp_datasets <- negbinom_mu_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else if( effect_type == 'trt_effect' ) {
            temp_datasets <- negbinom_mu_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else if( effect_type == 'drug_effect' ) {
            temp_datasets <- negbinom_mu_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else {
            temp_datasets <- negbinom_mu_data_simulation(data = base_data )
          }
        }
        else if(distribution_type == 'negbinom_prob_data_simulation') {
          if( effect_type == 'fixed_sd' ) {
            temp_datasets <- negbinom_prob_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else if( effect_type == 'variable_sd' ) {
            temp_datasets <- negbinom_prob_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else if( effect_type == 'trt_effect' ) {
            temp_datasets <- negbinom_prob_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else if( effect_type == 'drug_effect' ) {
            temp_datasets <- negbinom_prob_data_simulation(data = base_data , effect_type = effect_type, effect_size =  var_effect)
          }
          else {
            temp_datasets <- negbinom_prob_data_simulation(data = base_data )
          }
        }
        else { ## This uses the normal distribution by default
          if( effect_type == 'fixed_sd' ) {
            temp_datasets <- normal_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else if( effect_type == 'variable_sd' ) {
            temp_datasets <- normal_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else if( effect_type == 'trt_effect' ) {
            temp_datasets <- normal_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else if( effect_type == 'drug_effect' ) {
            temp_datasets <- normal_data_simulation(data = base_data , effect_type = effect_type, fixed_sd =  var_effect)
          }
          else {
            temp_datasets <- normal_data_simulation(base_data)
          }
        }
        
        ### Make the data objects at gene, mechanism, and class levels
        df <- do.call("rbind",temp_datasets)
        df <- as.data.frame(df)
        df <- replace(df, is.na(df), 0) ## confirm this
        trt_group <- as.data.frame(trt_group,row.names = colnames(df))
        dt <- as.data.table(df)
        ## aggregate to class and mechanism level
        dt[, header :=( gene_norm$header ), ]
        setkey(dt, header)
        dt <- annotations[dt]  # left outer join
        dt_class <- dt[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
        dt_class_analytic <- newMRexperiment(counts=dt_class[, .SD, .SDcols=!'class'],phenoData = AnnotatedDataFrame(trt_group))
        rownames(dt_class_analytic) <- dt_class$class
        dt_mech <- dt[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
        dt_mech_analytic <- newMRexperiment(counts=dt_mech[, .SD, .SDcols=!'mechanism'],phenoData = AnnotatedDataFrame(trt_group))
        rownames(dt_mech_analytic) <- dt_mech$mech
        ## make other objects for fitZig levels at class and mechanism levels
        d.ms.class <- dt_class_analytic
        d.ms.mech <- dt_mech_analytic
        #
        ##
        ###
        ##### ZeroInflatedGaussianModel ####
        ###
        ##
        # Gene level
        d.ms = newMRexperiment(df,phenoData = AnnotatedDataFrame(trt_group),featureData = NULL)
        p = cumNormStatFast(d.ms)
        d.ms = cumNorm(d.ms,p) ### We adjust our container to contain the normalization pecentile using cumNorm()
        nf = normFactors(d.ms)   ### We find the normalization factors per sample using normFactors()
        d.ms.fltr = filterData(d.ms,present = 10,depth = 2)   #### removing taxa that might be scarcely available
        d.fit = fitZig(d.ms.fltr, mod, control= settings) #### fitting the zero inflated gaussian model (ZIG) using the above specified model
        coef.df = MRcoefs(d.fit,number=94) #### recovering the model parameters and significance tests
        coef.df = coef.df[order(coef.df$adjPvalues),]#### ordering the taxa based on their adjusted significance
        zig_mod = d.fit$fit
        contrast_mod = d.fit$fit$design
        contrast_groups =  makeContrasts(((TTArrivalCtrl+TTExitCtrl)/2)-((TTArrivalTrt+TTExitTrt)/2), levels=contrast_mod)
        res_contrasts = contrasts.fit(zig_mod , contrast_groups)
        res_contrasts = eBayes(res_contrasts)
        log_FC_fitZig = topTable(res_contrasts, coef=1, adjust.method="fdr",number = 1000)
        sig_diff = length(which(log_FC_fitZig$adj.P.Val<0.05))
        gene_results <- as.data.table(log_FC_fitZig, keep.rownames = TRUE)
        gene_results <- gene_results[adj.P.Val<0.05]
        setkey(gene_results,rn)
        # class level, d.ms.class
        p = cumNormStatFast(d.ms.class)
        d.ms.class = cumNorm(d.ms.class,p) ### We adjust our container to contain the normalization pecentile using cumNorm()
        nf = normFactors(d.ms)   ### We find the normalization factors per sample using normFactors()
        d.ms.class.fltr = filterData(d.ms.class,present = 10,depth = 2)   #### removing taxa that might be scarcely available
        d.fit = fitZig(d.ms.class.fltr, mod, control= settings) #### fitting the zero inflated gaussian model (ZIG) using the above specified model
        coef.df = MRcoefs(d.fit,number=94) #### recovering the model parameters and significance tests
        coef.df = coef.df[order(coef.df$adjPvalues),]#### ordering the taxa based on their adjusted significance
        zig_mod = d.fit$fit
        contrast_mod = d.fit$fit$design
        contrast_groups =  makeContrasts(((TTArrivalCtrl+TTExitCtrl)/2)-((TTArrivalTrt+TTExitTrt)/2), levels=contrast_mod)
        res_contrasts = contrasts.fit(zig_mod , contrast_groups)
        res_contrasts = eBayes(res_contrasts)
        log_FC_fitZig = topTable(res_contrasts, coef=1, adjust.method="fdr",number = 1000)
        class_results <- as.data.table(log_FC_fitZig, keep.rownames = TRUE)
        class_results <- class_results[adj.P.Val<0.05]
        setkey(class_results,rn)
        # mech level, d.ms.mech
        p = cumNormStatFast(d.ms.mech)
        d.ms.mech = cumNorm(d.ms.mech,p) ### We adjust our container to contain the normalization pecentile using cumNorm()
        nf = normFactors(d.ms)   ### We find the normalization factors per sample using normFactors()
        d.ms.mech.fltr = filterData(d.ms.mech,present = 10,depth = 2)   #### removing taxa that might be scarcely available
        d.fit = fitZig(d.ms.mech.fltr, mod, control= settings) #### fitting the zero inflated gaussian model (ZIG) using the above specified model
        coef.df = MRcoefs(d.fit,number=94) #### recovering the model parameters and significance tests
        coef.df = coef.df[order(coef.df$adjPvalues),]#### ordering the taxa based on their adjusted significance
        zig_mod = d.fit$fit
        contrast_mod = d.fit$fit$design
        contrast_groups =  makeContrasts(((TTArrivalCtrl+TTExitCtrl)/2)-((TTArrivalTrt+TTExitTrt)/2), levels=contrast_mod)
        res_contrasts = contrasts.fit(zig_mod , contrast_groups)
        res_contrasts = eBayes(res_contrasts)
        log_FC_fitZig = topTable(res_contrasts, coef=1, adjust.method="fdr",number = 1000)
        mech_results <- as.data.table(log_FC_fitZig, keep.rownames = TRUE)
        mech_results <- mech_results[adj.P.Val<0.05]
        setkey(mech_results,rn)
        ## to name the output
        name <- v
        ## First output: the number of significantly different genes
        simulated_datasets[[name]] <- df
        vector_sig_diff <- c(vector_sig_diff, sig_diff) # save the number of significantly different genes
        simulated_fitzig[[name]] <- sig_diff
        ## Gene level summary:
        temp_simulation_gene_summary <- simulation_gene_summary
        setkey(temp_simulation_gene_summary, header)
        temp_simulation_gene_summary <- temp_simulation_gene_summary[row_num %in% gene_results$rn, `:=` (Significant = 1, Samples = eval(name))]
        temp_simulation_gene_summary <- temp_simulation_gene_summary[!row_num %in% gene_results$rn, `:=` (Significant = 0, Samples = eval(name))]
        output_temp_gene_summary <- rbind(output_temp_gene_summary,temp_simulation_gene_summary, fill=TRUE)
        ## Class level summary:
        temp_simulation_class_summary <- simulation_class_summary
        setkey(temp_simulation_class_summary, class)
        temp_simulation_class_summary <- temp_simulation_class_summary[class %in% class_results$rn, `:=` (Significant = 1, Samples = eval(name))]
        temp_simulation_class_summary <- temp_simulation_class_summary[!class %in% class_results$rn, `:=` (Significant = 0, Samples = eval(name))]
        output_temp_class_summary <- rbind(output_temp_class_summary,temp_simulation_class_summary, fill=TRUE)
        ## Mechanism level summary:
        temp_simulation_mech_summary <- simulation_mech_summary
        setkey(temp_simulation_mech_summary, mech)
        temp_simulation_mech_summary <- temp_simulation_mech_summary[mech %in% mech_results$rn, `:=` (Significant = 1, Samples = eval(name))]
        temp_simulation_mech_summary <- temp_simulation_mech_summary[!mech %in% mech_results$rn, `:=` (Significant = 0, Samples = eval(name))]
        output_temp_mech_summary <- rbind(output_temp_mech_summary,temp_simulation_mech_summary, fill=TRUE)
      }
      ## Gene level summary:
      out1_gene_summary <- output_temp_gene_summary[, var_size := eval(var_effect)]
      output_temp1_gene_summary <- rbind(output_temp1_gene_summary ,out1_gene_summary, fill=TRUE)
      ## Class level summary:
      out1_class_summary <- output_temp_class_summary[,var_size := eval(var_effect)]
      output_temp1_class_summary <- rbind(output_temp1_class_summary,out1_class_summary, fill=TRUE)
      ## Mechanism level summary:
      out1_mech_summary <- output_temp_mech_summary[,var_size := eval(var_effect)]
      output_temp1_mech_summary <- rbind(output_temp1_mech_summary ,out1_mech_summary, fill=TRUE)
    }
    sim_name <- paste("Simulation_",i,sep='')
    # full_results[[sim_name]] <- simulated_fitzig
    # simulation_results <- cbind(simulation_results, temp_name = vector_sig_diff)
    # names(simulation_results)[names(simulation_results) == "temp_name"] <- sim_name
    # full_simulated_results[[sim_name]] <- simulated_datasets
    ## Gene level summary:
    out2_gene_summary <- output_temp1_gene_summary[,Simulation := eval(sim_name)]
    output_gene_summary <- rbind(output_gene_summary ,out2_gene_summary, fill=TRUE)
    ## Class level summary:
    out2_class_summary <- output_temp1_class_summary[,Simulation := eval(sim_name)]
    output_class_summary <- rbind(output_class_summary ,out2_class_summary, fill=TRUE)
    ## Mechanism level summary:
    out2_mech_summary <- output_temp1_mech_summary[,Simulation := eval(sim_name)]
    output_mech_summary <- rbind(output_mech_summary ,out2_mech_summary, fill=TRUE)
  }
}





normal_data_simulation <- function(data = '' , effect_type='',fixed_sd='', effect_size='') {
  if( effect_type == 'fixed_sd' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]),fixed_sd)
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), fixed_sd)
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), fixed_sd)
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]),fixed_sd)
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'variable_sd' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]), as.numeric(Ctrl_arrival_summary[gene,"sd"]) * effect_size)
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"])* effect_size)
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"])* effect_size)
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"])* effect_size)
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'trt_effect' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]) + effect_size , as.numeric(Trt_arrival_summary[gene,"sd"]))
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]) + effect_size , as.numeric(Trt_exit_summary[gene,"sd"]))
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'drug_effect' ) {
    for(gene in 1:nrow(data)) {
      if ( amr_gene_summary[gene,"group"] == enriched_group ){
        temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
        temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
        temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
        temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]) + as.numeric(effect_size), as.numeric(Trt_exit_summary[gene,"sd"]))
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_datasets[[gene]] <- temp_row
      }
      else{
        temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
        temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
        temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
        temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"]))
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_datasets[[gene]] <- temp_row
      }
    }
  }
  else {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]),as.numeric(Ctrl_arrival_summary[gene,"sd"]))
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"]))
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  return(temp_datasets)
}




## Neg binom, with size and mu
negbinom_mu_data_simulation <- function(data = '' , effect_type='',fixed_sd='', effect_size='') {
  if( effect_type == 'fixed_sd' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= fixed_sd ,mu= as.numeric(Ctrl_arrival_summary[gene,"mu"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= fixed_sd,mu= as.numeric(Trt_arrival_summary[gene,"mu"]) )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= fixed_sd,mu= as.numeric(Ctrl_exit_summary[gene,"mu"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= fixed_sd ,mu= as.numeric(Trt_exit_summary[gene,"mu"]) )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if(effect_type == 'variable_sd' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),mu= as.numeric(Ctrl_arrival_summary[gene,"mean"])* effect_size )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),mu= as.numeric(Trt_arrival_summary[gene,"mean"]) * effect_size )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),mu= as.numeric(Ctrl_exit_summary[gene,"mean"])* effect_size )
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),mu= as.numeric(Trt_exit_summary[gene,"mean"]) * effect_size )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if(effect_type == 'trt_effect' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),mu= as.numeric(Ctrl_arrival_summary[gene,"mu"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),mu= as.numeric(Trt_arrival_summary[gene,"mu"])+ effect_size )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),mu= as.numeric(Ctrl_exit_summary[gene,"mu"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),mu= as.numeric(Trt_exit_summary[gene,"mu"]) + effect_size )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'drug_effect' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      if ( amr_gene_summary[gene,"group"] == enriched_group ) {
        temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),mu= as.numeric(Ctrl_arrival_summary[gene,"mu"]) )
        temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),mu= as.numeric(Trt_arrival_summary[gene,"mu"]) )
        temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),mu= as.numeric(Ctrl_exit_summary[gene,"mu"]) )
        temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]) + as.numeric(effect_size),mu= as.numeric(Trt_exit_summary[gene,"mu"]) )
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_row[is.na(temp_row)] = 0
        temp_datasets[[gene]] <- temp_row
      } else {
        temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),mu= as.numeric(Ctrl_arrival_summary[gene,"mu"]) )
        temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),mu= as.numeric(Trt_arrival_summary[gene,"mu"]) )
        temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),mu= as.numeric(Ctrl_exit_summary[gene,"mu"]) )
        temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),mu= as.numeric(Trt_exit_summary[gene,"mu"]) )
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_row[is.na(temp_row)] = 0
        temp_datasets[[gene]] <- temp_row
      }
    }
  }
  else {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),mu= as.numeric(Ctrl_arrival_summary[gene,"mu"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),mu= as.numeric(Trt_arrival_summary[gene,"mu"]) )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),mu= as.numeric(Ctrl_exit_summary[gene,"mu"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),mu= as.numeric(Trt_exit_summary[gene,"mu"]) )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  return(temp_datasets)
}


### Neg binom, with size and probability
negbinom_prob_data_simulation <- function(data = '' , effect_type='',fixed_sd='', effect_size='') {
  if( effect_type == 'fixed_sd' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= fixed_sd,prob= as.numeric(Ctrl_arrival_summary[gene,"prob"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= fixed_sd,prob= as.numeric(Trt_arrival_summary[gene,"prob"]) )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= fixed_sd,prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= fixed_sd,prob= as.numeric(Trt_exit_summary[gene,"prob"]) )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if(effect_type == 'variable_sd' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),prob= as.numeric(Ctrl_arrival_summary[gene,"prob"])* effect_size )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),prob= as.numeric(Trt_arrival_summary[gene,"prob"]) * effect_size)
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) * effect_size)
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),prob= as.numeric(Trt_exit_summary[gene,"prob"]) * effect_size )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if(effect_type == 'trt_effect' ) {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),prob= as.numeric(Ctrl_arrival_summary[gene,"prob"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),prob= as.numeric(Trt_arrival_summary[gene,"prob"]) + effect_size)
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),prob= as.numeric(Trt_exit_summary[gene,"prob"]) + effect_size )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'drug_effect' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      if ( amr_gene_summary[gene,"group"] == enriched_group ){
        temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),prob= as.numeric(Ctrl_arrival_summary[gene,"prob"]) )
        temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),prob= as.numeric(Trt_arrival_summary[gene,"prob"]) )
        temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) )
        temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]) + as.numeric(effect_size),prob= as.numeric(Trt_exit_summary[gene,"prob"]) )
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_row[is.na(temp_row)] = 0
        temp_datasets[[gene]] <- temp_row
        } else {
        temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),prob= as.numeric(Ctrl_arrival_summary[gene,"prob"]) )
        temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),prob= as.numeric(Trt_arrival_summary[gene,"prob"]) )
        temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) )
        temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),prob= as.numeric(Trt_exit_summary[gene,"prob"]) )
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_row[is.na(temp_row)] = 0
        temp_datasets[[gene]] <- temp_row
        }
    }
  }
  else {
    for(gene in 1:nrow(data)) {
      temp_row_Arrival_Ctrl =rnbinom(v/4,size= as.numeric(Ctrl_arrival_summary[gene,"succ"]),prob= as.numeric(Ctrl_arrival_summary[gene,"prob"]) )
      temp_row_Arrival_Trt = rnbinom(v/4,size= as.numeric(Trt_arrival_summary[gene,"succ"]),prob= as.numeric(Trt_arrival_summary[gene,"prob"]) )
      temp_row_Exit_Ctrl = rnbinom(v/4,size= as.numeric(Ctrl_exit_summary[gene,"succ"]),prob= as.numeric(Ctrl_exit_summary[gene,"prob"]) )
      temp_row_Exit_Trt = rnbinom(v/4,size= as.numeric(Trt_exit_summary[gene,"succ"]),prob= as.numeric(Trt_exit_summary[gene,"prob"]) )
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  return(temp_datasets)
}


mix_binom_normal_data_simulation <- function(data = '' , effect_type='',fixed_sd='', effect_size='') {
  if( effect_type == 'fixed_sd' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]),fixed_sd)
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), fixed_sd)
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), fixed_sd)
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]),fixed_sd)
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'variable_sd' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      num_Arrival_Ctrl <- length(which(rbinom(v/4, 1, as.numeric(Ctrl_arrival_summary[gene,"prob"]))==1))
      num_0_Arrival_Ctrl <- (v/4 -num_Arrival_Ctrl )
      temp_1_row_Arrival_Ctrl = rnorm(num_Arrival_Ctrl,as.numeric(Ctrl_arrival_summary[gene,"mean"]), as.numeric(Ctrl_arrival_summary[gene,"sd"]) * effect_size)
      temp_0_row_Arrival_Ctrl = rep(0, each = num_0_Arrival_Ctrl)
      
      num_Arrival_Trt <- length(which(rbinom(v/4, 1, as.numeric(Trt_arrival_summary[gene,"prob"]))==1))
      num_0_Arrival_Trt <- (v/4 - num_Arrival_Trt )
      temp_1_row_Arrival_Trt = rnorm(num_Arrival_Trt,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]) * effect_size)
      temp_0_row_Arrival_Trt = rep(0, each = num_0_Arrival_Trt)
      
      num_Exit_Ctrl <- length(which(rbinom(v/4, 1, as.numeric(Ctrl_exit_summary[gene,"prob"]))==1))
      num_0_Exit_Ctrl <- (v/4 -num_Exit_Ctrl )
      temp_1_row_Exit_Ctrl = rnorm(num_Exit_Ctrl,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]) * effect_size)
      temp_0_row_Exit_Ctrl = rep(0, each = num_0_Exit_Ctrl)
      num_Exit_Trt <- length(which(rbinom(v/4, 1, as.numeric(Trt_exit_summary[gene,"prob"]))==1))
      num_0_Exit_Trt <- (v/4 -num_Exit_Trt )
      temp_1_row_Exit_Trt = rnorm(num_Exit_Trt,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"]) * effect_size)
      temp_0_row_Exit_Trt = rep(0, each = num_0_Exit_Trt)
      temp_row = c(temp_1_row_Arrival_Ctrl,temp_0_row_Arrival_Ctrl,temp_1_row_Arrival_Trt,temp_0_row_Arrival_Trt,temp_1_row_Exit_Ctrl,temp_0_row_Exit_Ctrl,temp_1_row_Exit_Trt,temp_0_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_row[is.na(temp_row)] = 0
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'trt_effect' ) {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]) + effect_size , as.numeric(Trt_arrival_summary[gene,"sd"]))
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]) + effect_size , as.numeric(Trt_exit_summary[gene,"sd"]))
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  else if( effect_type == 'drug_effect' ) {
    for(gene in 1:nrow(data)) {
      if ( amr_gene_summary[gene,"group"] == enriched_group ){
        temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
        temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
        temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
        temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]) + as.numeric(effect_size), as.numeric(Trt_exit_summary[gene,"sd"]))
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_datasets[[gene]] <- temp_row
      }
      else{
        temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]) ,as.numeric(Ctrl_arrival_summary[gene,"sd"]))
        temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
        temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
        temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"]))
        temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
        temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
        temp_datasets[[gene]] <- temp_row
      }
    }
  }
  else {
    for(gene in 1:nrow(data)) {
      ## If "effect_type" == ref_data
      ## now AMR counts are calculated individually for each of two treatment groups and time points
      temp_row_Arrival_Ctrl = rnorm(v/4,as.numeric(Ctrl_arrival_summary[gene,"mean"]),as.numeric(Ctrl_arrival_summary[gene,"sd"]))
      temp_row_Arrival_Trt = rnorm(v/4,as.numeric(Trt_arrival_summary[gene,"mean"]), as.numeric(Trt_arrival_summary[gene,"sd"]))
      temp_row_Exit_Ctrl = rnorm(v/4,as.numeric(Ctrl_exit_summary[gene,"mean"]), as.numeric(Ctrl_exit_summary[gene,"sd"]))
      temp_row_Exit_Trt = rnorm(v/4,as.numeric(Trt_exit_summary[gene,"mean"]), as.numeric(Trt_exit_summary[gene,"sd"]))
      temp_row = c(temp_row_Arrival_Ctrl,temp_row_Arrival_Trt,temp_row_Exit_Ctrl,temp_row_Exit_Trt)
      temp_row[temp_row<0] = 0 # change negative values to 0, ## probably have to change this
      temp_datasets[[gene]] <- temp_row
    }
  }
  return(temp_datasets)
}
