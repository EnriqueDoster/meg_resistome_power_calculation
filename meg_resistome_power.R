### Only using NCBA2 dataset
library(metagenomeSeq)
library(data.table)
require(ggplot2)
require(vegan)
set.seed(154)  

### First load the data, normalize and make data.table object ####
## Borrows code from Steven Lakin's "metagenomeseq_analytic_template.R"

setwd('/home/enrique/Dropbox/Projects/meg_resistome_power_calculation')
source('scripts/meg_utility_functions.R')
source('scripts/data_simulation_functions.R')
metadata_filepath = 'NCBA2_metadata.csv'
megares_annotation_filename = 'megaresbio_annotations_v1.01.csv'
sample_column_id = 'ID'

# Load data
amr <- newMRexperiment(read.table('NCBA2_AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys
metadata <- read.csv(metadata_filepath, header=T)
metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])

#Normalize
cumNorm(amr)

# Aggregate the normalized counts for AMR using the annotations data table, SQL
# outer join, and aggregation with vectorized lapply
amr_norm <- data.table(MRcounts(amr, norm=T))
amr_raw <- data.table(MRcounts(amr, norm=F))

amr_norm[, header :=( rownames(amr) ), ]
setkey(amr_norm, header)
amr_norm <- annotations[amr_norm]  # left outer join

amr_raw[, header :=( rownames(amr) ), ]
setkey(amr_raw, header)
amr_raw <- annotations[amr_raw]  # left outer join

# Group the AMR data by level for analysis
amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
rownames(amr_class_analytic) <- amr_class$class

amr_class_raw <- amr_raw[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_raw_analytic <- newMRexperiment(counts=amr_class_raw[, .SD, .SDcols=!'class'])
rownames(amr_class_raw_analytic) <- amr_class_raw$class

amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_analytic) <- amr_mech$mechanism

amr_mech_raw <- amr_raw[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_raw_analytic <- newMRexperiment(counts=amr_mech_raw[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_raw_analytic) <- amr_mech_raw$mechanism

amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
rownames(amr_group_analytic) <- amr_group$group

amr_group_raw <- amr_raw[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_raw_analytic <- newMRexperiment(counts=amr_group_raw[, .SD, .SDcols=!'group'])
rownames(amr_group_raw_analytic) <- amr_group_raw$group

amr_gene_analytic <- newMRexperiment(
  counts=amr_norm[, .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
amr_gene_raw_analytic <- newMRexperiment(
  counts=amr_raw[, .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])

rownames(amr_gene_analytic) <- amr_norm$header
rownames(amr_gene_raw_analytic) <- amr_raw$header

amr_melted_analytic <- rbind(melt_dt(MRcounts(amr_class_analytic), 'Class'),
                             melt_dt(MRcounts(amr_mech_analytic), 'Mechanism'),
                             melt_dt(MRcounts(amr_group_analytic), 'Group'),
                             melt_dt(MRcounts(amr_gene_analytic), 'Gene'))
amr_melted_raw_analytic <- rbind(melt_dt(MRcounts(amr_class_raw_analytic), 'Class'),
                                 melt_dt(MRcounts(amr_mech_raw_analytic), 'Mechanism'),
                                 melt_dt(MRcounts(amr_group_raw_analytic), 'Group'),
                                 melt_dt(MRcounts(amr_gene_raw_analytic), 'Gene'))


# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[, sample_column_id]), ])
setkey(metadata, ID)


############## Power calculation ###################
## This takes the amr_norm matrix and just removes the class, mechanism, and group columns
amr_melted_raw_analytic<- as.data.table(amr_melted_raw_analytic) 
norm_gene <- amr_melted_raw_analytic[Level_ID=="Mechanism"][,Level_ID :=NULL ]
setkey(norm_gene,ID)
norm_gene <- metadata[norm_gene]
norm_gene <- norm_gene[,Group := paste0(Treatment,Time) ]
norm_gene[, Round_count := round(Normalized_Count)]

mech_annotations <- as.data.table(annotations)
mech_annotations[,header := NULL]
mech_annotations[,group:=NULL]
setkey(mech_annotations, mechanism)
mech_annotations <- unique(mech_annotations, by = key(mech_annotations))

setkey(norm_gene, Name)

norm_gene <- mech_annotations[norm_gene]

### Need to add the probability of "success" or getting a count of that gene
amr_gene_summary <- norm_gene[, .(median = median(Round_count),
              mean = mean(Round_count),
              sd = sd(Round_count),
              succ = length(which(Round_count >0)), 
              prob = (length(which(Round_count >0))/.N )),
          by=.(mechanism,Group,class)]
amr_gene_summary[, mu := (succ*(1-prob))/prob, by=.(mechanism,Group)]
amr_gene_summary[,lapply(.SD,function(x){ifelse(is.na(x),0,x)})]

Ctrl_arrival_summary <- amr_gene_summary[Group=="CtrlArrival"]
Ctrl_exit_summary <- amr_gene_summary[Group=="CtrlExit"]
Trt_arrival_summary <- amr_gene_summary[Group=="TrtArrival"]
Trt_exit_summary <- amr_gene_summary[Group=="TrtExit"]


#
##
####
##### Simulation Code ##########
###
##
#

## This is used to loop the simulation code , these variables need to be edited:
# distribution_type
# sample_number_variable
# simulation_runs
# effect type 

distribution_type = 'negbinom_mu_data_simulation' # or 'negbinom_prob_data_simulation' or if empty variable it runs normal distribution
sample_number_variable <- seq(12,96,by = 4) ## number of samples per simulated dataset, has to be divisible by 4
simulation_runs <- 10
effect_type <- 'variable_sd' ## can be 'drug_effect' and enriched_group, or 'fixed_sd' , 'trt_effect', or 'variable_sd'
#effect_size <- as.numeric(2)

## For enriched genes "drug_effect"
enriched_group <- 'CTX' ## Gene group to be enriched
#fixed_sd <- seq(5,35,by = 5)

## renames the effect sizes
#variable_effect <- fixed_sd
variable_effect <- c(0.5,1,2)

# # Make empty datasets ####
gene_norm <- amr_raw[,c("class","mechanism","group") :=NULL]
simulation_gene_summary <- as.data.table(gene_norm$header)
names(simulation_gene_summary)[1] <- "header"
simulation_gene_summary$row_num <- seq(1,nrow(simulation_gene_summary))

simulation_class_summary <- as.data.table(amr_class$class)
names(simulation_class_summary)[1] <- "class"
simulation_class_summary$row_num <- seq(1,nrow(simulation_class_summary))

simulation_mech_summary <- as.data.table(amr_mech$mechanism)
names(simulation_mech_summary)[1] <- "mech"
simulation_mech_summary$row_num <- seq(1,nrow(simulation_mech_summary))

## gene level summary:
output_gene_summary <- data.table(header = '',
                                  row_num='',
                                  Significant = '',
                                  Samples = '',
                                  Simulation='',
                                  var_size='')
output_temp1_gene_summary <- data.table(header = '',
                                       row_num='',
                                       Significant = '',
                                       Samples = '',
                                       var_size='')
output_temp_gene_summary <- data.table(header = '',
                                  row_num='',
                                  Significant = '',
                                  Samples = '')
## Class level summary:
output_class_summary <- data.table(header = '',
                                  row_num='',
                                  Significant = '',
                                  Samples = '',
                                  Simulation='',
                                  var_size='')
output_temp1_class_summary <- data.table(header = '',
                                        row_num='',
                                        Significant = '',
                                        Samples = '',
                                        var_size='')
output_temp_class_summary <- data.table(header = '',
                                       row_num='',
                                       Significant = '',
                                       Samples = '')
## Mechanism level summary:
output_mech_summary <- data.table(header = '',
                                  row_num='',
                                  Significant = '',
                                  Samples = '',
                                  Simulation='',
                                  var_size='')
output_temp1_mech_summary <- data.table(header = '',
                                        row_num='',
                                        Significant = '',
                                        Samples = '',
                                        var_size='')
output_temp_mech_summary <- data.table(header = '',
                                       row_num='',
                                       Significant = '',
                                       Samples = '')
simulated_datasets <- list()
simulated_fitzig <- list()
full_results <- list()
full_simulated_results <- list()
simulation_results <- data.frame(matrix(ncol= 1,nrow = length(sample_number_variable)))
names(simulation_results)[1] <- "sample_numbers"
for( i in 1:length(sample_number_variable)) {
  v <- sample_number_variable[i]
  simulation_results$sample_numbers[i] <- paste(v) #,"_Samples",sep=''
}

#
##
###
#### Start the simulation ######
###
##
#

## add input_data

input_data <- amr_mech

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
      ### different distributions for data simulation
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
      dt[, mechanism :=( amr_mech$mechanism ), ]
      setkey(dt, mechanism)
      dt <- mech_annotations[dt]  # left outer join
      dt_class <- dt[, lapply(.SD, sum), by='class', .SDcols=!c('mechanism')]
      dt_class_analytic <- newMRexperiment(counts=dt_class[, .SD, .SDcols=!'class'],phenoData = AnnotatedDataFrame(trt_group))
      rownames(dt_class_analytic) <- dt_class$class
      dt_mech <- dt[, lapply(.SD, sum), by='mechanism', .SDcols=!c('class')]
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
    ## Class level summary:
    out1_class_summary <- output_temp_class_summary[,var_size := eval(var_effect)]
    output_temp1_class_summary <- rbind(output_temp1_class_summary,out1_class_summary, fill=TRUE)
    ## Mechanism level summary:
    out1_mech_summary <- output_temp_mech_summary[,var_size := eval(var_effect)]
    output_temp1_mech_summary <- rbind(output_temp1_mech_summary ,out1_mech_summary, fill=TRUE)
  }
  sim_name <- paste("Simulation_",i,sep='')
  ## Class level summary:
  out2_class_summary <- output_temp1_class_summary[,Simulation := eval(sim_name)]
  output_class_summary <- rbind(output_class_summary ,out2_class_summary, fill=TRUE)
  ## Mechanism level summary:
  out2_mech_summary <- output_temp1_mech_summary[,Simulation := eval(sim_name)]
  output_mech_summary <- rbind(output_mech_summary ,out2_mech_summary, fill=TRUE)
}


####### Plot for results ##########

output_class_summary <- output_class_summary[class !=''] 
output_class_summary[, Significant:=as.numeric(Significant)]
summarized_class_output <- output_class_summary[, .(Proportion = sum(Significant)/.N), by= .(class,Samples,var_size)]
str(summarized_class_output)
summarized_class_output[, Samples:=as.numeric(Samples)]

output_mech_summary <- output_mech_summary[mech !=''] 
output_mech_summary[, Significant:=as.numeric(Significant)]
summarized_mech_output <- output_mech_summary[, .(Proportion = sum(Significant)/.N), by= .(mech,Samples,var_size)]
str(summarized_mech_output)
summarized_mech_output[, Samples:=as.numeric(Samples)]



## plot results at gene level
negbinom_mu_simulation_mech <- ggplot(data=summarized_mech_output , 
  aes(x=Samples,y=Proportion,group=mech, color=mech , label = mech)) + #group=Simulation, 
  #stat_smooth(aes(y = Proportion, group=mech, colour=mech), method=lm, formula = y ~ poly(x,2), level=0.95 , fill=NA) + #tw
  geom_text(data=subset(summarized_mech_output, Samples ==52  & var_size == 2),aes(Samples,Proportion,label=mech), position =  position_dodge(15)) +
  #geom_label(aes(fill = factor(mech)), colour = "white", fontface = "bold")+
  geom_line(aes(group=interaction(mech,var_size),linetype=var_size))+ #aes(linetype=count)
  #geom_point() +   
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="none",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  ) +
  #xlim(0,96) + ylim(0,1) +
  #scale_x_continuous(limits = c(14,96)) + scale_y_continuous(limits= c(0,1))+
  xlab('\n Number of samples in simulation') +
  ylab('Proportion significantly different per 10 simulations\n') +
  ggtitle('NegBinom distribution (mu,size) and variable drug effect on mu - AMR mech level \n')
negbinom_mu_simulation_mech 



#### plot the results at the class level##

negbinom_mu_simulation_class <- ggplot(data=summarized_class_output , 
       aes(x=Samples,y=Proportion,group=class, color=class , label = class)) + #group=Simulation, 
      #stat_smooth(aes(y = Proportion, group=class, colour=class), method=lm, formula = y ~ poly(x,2), level=0.95 , fill=NA) + #tw
      geom_text(data=subset(summarized_class_output, Samples ==24  & var_size == 1),
                aes(Samples,Proportion,label=class), position =  position_dodge(15)) +
      #geom_label(aes(fill = factor(class)), colour = "white", fontface = "bold")+
      geom_line(aes(group=interaction(class,var_size),linetype=var_size))+ #aes(linetype=count)
      #geom_point() +   
       theme(
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         strip.text.x=element_text(size=24),
         strip.text.y=element_text(size=24, angle=0),
         axis.text.x=element_text(size=16, angle=20, hjust=1),
         axis.text.y=element_text(size=22),
         axis.title=element_text(size=26),
         legend.position="none",
         panel.spacing=unit(0.1, "lines"),
         plot.title=element_text(size=32, hjust=0.5),
         legend.text=element_text(size=18),
         legend.title=element_blank()
       ) +
  #xlim(0,96) + ylim(0,1) +
  #scale_x_continuous(limits = c(14,96)) + scale_y_continuous(limits= c(0,1))+
  xlab('\n Number of samples in simulation') +
  ylab('Proportion significantly different per 10 simulations\n') +
  ggtitle('NegBinom distribution (mu,size) and variable drug effect on mu - AMR Class level \n')
negbinom_mu_simulation_class
  


normal_dist_simulation <- ggplot(data=sim_long,aes(x=sample_numbers,y=value,group=Simulation, color=Simulation)) +
  geom_line()+ #aes(linetype=count)
  #geom_point() +   
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="none",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  ) +
  xlab('\n Number of samples') +
  ylab('Count of significantly different AMR gene abundances\n') +
  ggtitle('Power calculation using NCBA2 dataset \n')
normal_dist_simulation 


