#Author: Arvis Sulovari, PhD
#The scripts and functions below were used to cally out analyses presented in Sulovari et al., 2017 (See README)

##########################################################
##~!~!~!~##Converting the above into a function##~!~!~!~##
##########################################################


#Read-in the locations of probes that fall into common CNV regions
sage_common_probes <- read.table("SAGE//assoc//AF_CNVs_SAGE_commonProbeLocations.txt",header=F)

cidr_common_probes_EA <- read.table("CIDR/assoc/PLINK/Gene_level_collapse/CIDR_EA_rare_locations.txt")
cidr_common_probes_EA <- na.omit(cidr_common_probes_EA$V5)

cidr_common_probes_AA <- read.table("CIDR/assoc/PLINK/Gene_level_collapse/CIDR_AA_rare_locations.txt")
cidr_common_probes_AA <- na.omit(cidr_common_probes_AA$V5)

ozalc_common_probes <- read.table("OZ-ALC//association/PLINK/Gene_level_collapse/OZALC_rare_locations.txt")
ozalc_common_probes <- na.omit(ozalc_common_probes$V5)




##Step 1: Read-in file and create CNV_status column
sage_test <- read.delim("SAGE/sage_4090sample_phenotypes.txt",header=T)
sage_test$Probe_CNVstatus <- rep(2,nrow(sage_test))

cidr_test <- read.delim("CIDR/CIDR_1913sample_phenotypes.txt",header=T)
cidr_test$Probe_CNVstatus <- rep(2,nrow(cidr_test))
#Merge with CIDR PCA data
pca_cidr <- read.table("CIDR/assoc/CIDR_ALL_autosomes_PCA_20.eigenvec",header=T)
cidr_test <- merge(cidr_test,pca_cidr,by="SampleID")

ozalc_test <- read.delim("OZ-ALC/OZALC_1989sample_phenotypes.txt",header=T)
ozalc_test$Probe_CNVstatus <- rep(2,nrow(ozalc_test))
#Merge with OZALC PCA data
pca_ozalc <- read.table("OZ-ALC/association/OZALC_ALL_autosomes_PCA_20.eigenvec",header=T)
ozalc_test <- merge(ozalc_test,pca_ozalc,by="SampleID")


##Step 2: Remove missing phenotypes and reformat final_type column
sage_test_clean <- sage_test[sage_test$final_type=="CASE" | sage_test$final_type=="CNTL",]
sage_test_clean$final_type <- factor(sage_test_clean$final_type)
sage_test_clean$final_type <- ifelse(sage_test_clean$final_type=="CASE",2,1)

cidr_test_clean <- cidr_test[cidr_test$type=="CA" | cidr_test$type=="CO",]
cidr_test_clean$type <- factor(cidr_test_clean$type)
cidr_test_clean$final_type <- ifelse(cidr_test_clean$type=="CA",2,1)
cidr_test_clean$chip.id <- cidr_test_clean$SampleID

ozalc_test_clean <- ozalc_test[ozalc_test$casestatus== 0 | ozalc_test$casestatus==1,]
ozalc_test_clean$casestatus <- factor(ozalc_test_clean$casestatus)
ozalc_test_clean$final_type <- ifelse(ozalc_test_clean$casestatus==1,2,1)
ozalc_test_clean$chip.id <- ozalc_test_clean$SampleID


##Step 3: Separate two major ethinic groups
sage_test_clean_EA <- unique(sage_test_clean[sage_test_clean$race=="WHITE",])
sage_test_clean_EA <- sage_test_clean_EA[sage_test_clean_EA$chip.id != "NA",]

sage_test_clean_AA <- unique(sage_test_clean[sage_test_clean$race=="BLACK",])
sage_test_clean_AA <- sage_test_clean_AA[sage_test_clean_AA$chip.id != "NA",]

cidr_test_clean_EA <- unique(cidr_test_clean[cidr_test_clean$ethnicity==6,])
cidr_test_clean_EA <- cidr_test_clean_EA[cidr_test_clean_EA$SampleID != "NA",]

cidr_test_clean_AA <- unique(cidr_test_clean[cidr_test_clean$ethnicity==4,])
cidr_test_clean_AA <- cidr_test_clean_AA[cidr_test_clean_AA$SampleID != "NA",]

ozalc_test_clean_EA <- ozalc_test_clean[ozalc_test_clean$SampleID != "NA",]


#Make sure CIDR and OZALC have the corrrect column names for the functions below to execute properly
#chip.id, age_int, tissue , sex.x

cidr_test_clean_AA$chip.id <- cidr_test_clean_AA$SampleID
cidr_test_clean_EA$chip.id <- cidr_test_clean_EA$SampleID

cidr_test_clean_AA$age_int <- cidr_test_clean_AA$age
cidr_test_clean_EA$age_int <- cidr_test_clean_EA$age

#CIDR pheno files needs to be merged with _test_clean to get DNA source information
cidr_pheno2 <- read.delim("CIDR/phenotypes/CIDR_PhenoFile2_withDNAsource.txt",header=T)
cidr_pheno2$SampleID <- paste0(cidr_pheno2$Human1Mv1_C.Sentrix.ID,"_",cidr_pheno2$Human1Mv1_C.Sentrix.Position)
cidr_test_clean_AA <- merge(cidr_test_clean_AA,cidr_pheno2,by="SampleID")
cidr_test_clean_EA <- merge(cidr_test_clean_EA,cidr_pheno2,by="SampleID")

cidr_test_clean_AA$tissue <- cidr_test_clean_AA$DNA.Source
cidr_test_clean_EA$tissue <- cidr_test_clean_EA$DNA.Source

cidr_test_clean_AA$sex.x <- cidr_test_clean_AA$sex
cidr_test_clean_EA$sex.x <- cidr_test_clean_EA$sex


ozalc_test_clean_EA$age_int <- ozalc_test_clean_EA$age
ozalc_test_clean_EA$tissue <- ozalc_test_clean_EA$DNA.Source.x 
ozalc_test_clean_EA$sex.x <- ozalc_test_clean_EA$sex

#Create logistic regression output matrix where the ORs and P-values will be stored
#Probe no., probe name, OR_2by2, CI_2by2_upper, CI_2by2_lower, P_2by2, OR_2by3_del, CI_2by3_upper_del, CI_2by3_lower_del, P_2by3_del,
#...OR_2by3_dup, CI_2by3_upper_dup, CI_2by3_lower_dup, P_2by3_dup, OR_2by2_dup, CI_2by2_upper_dup, CI_2by2_lower_dup, P_2by2_dup,
#...OR_2by2_del, CI_2by2_upper_del, CI_2by2_lower_del, P_2by2_del

sage_out_mx_EA <- array(NA, dim=c(24000,24))
sage_out_mx_AA <- array(NA, dim=c(24000,24))

cidr_out_mx_EA <- array(NA, dim=c(16000,24))
cidr_out_mx_AA <- array(NA, dim=c(10000,24))

ozalc_out_mx_EA <- array(NA, dim=c(9000,24))




#Test 1

cnvAssoc_2by3 <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){

  
cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))

loc <- as.numeric(probe_file[probe_pos,]$Loci.start)

#ADD CHROMOSOME INFORMATION!!!!!!
loc_chr <- probe_file[probe_pos, ]$Chr

sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)

#*
sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
#* 

for (i in 1:length(sampleID_atProbe_loc)) {
      t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
      cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }


#Run the 2 by 3

if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0) {
  
  cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
  
}

else {
  stop("Single level CNV status in 2_by_3 DELETIONS")
}


if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus) != 0){
  
  cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus <- 3
}


else {
  stop("Single level CNV status in 2_by_3 DUPLICATIONS")
}




#COnsider CNV status to be a factor
cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
cnv_file_unique <- within(cnv_file_unique, Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
#keep only non-missing phenotype samples

#Run logistic regression model
#model_probe1 <- glm(factor(final_type)~Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
model_probe1 <- speedglm(final_type ~ Probe_CNVstatus + age_int + factor(tissue) + factor(sex.x) + PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma(log))
  
s <- summary(model_probe1)
or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
temp_pos = 1

assoc_out_mx[temp_pos,1] <- probe_pos
assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)

assoc_out_mx[temp_pos,9] <- or[2,1]
assoc_out_mx[temp_pos,10] <- or[2,2]
assoc_out_mx[temp_pos,11] <- or[2,3]
assoc_out_mx[temp_pos,12] <- as.numeric(as.character(coef(s)[2,4]))
assoc_out_mx[temp_pos,13] <- or[3,1] 
assoc_out_mx[temp_pos,14] <- or[3,2]
assoc_out_mx[temp_pos,15] <- or[3,3]
assoc_out_mx[temp_pos,16] <- as.numeric(as.character(coef(s)[3,4]))
  
return(assoc_out_mx)
}




#Test 2

cnvAssoc_2by2 <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)

  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
  
  #Run the 2 by 2
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2 | cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0){ 
            
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2 | cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
    
  }
  
  else {
    stop("Single level CNV status in 2_by_2")
  }
  
  
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int + factor(tissue) + factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  

  assoc_out_mx[temp_pos,5] <- or[2,1]
  assoc_out_mx[temp_pos,6] <- or[2,2]
  assoc_out_mx[temp_pos,7] <- or[2,3]
  assoc_out_mx[temp_pos,8] <- as.numeric(as.character(coef(s)[2,4]))

  return(assoc_out_mx)
  
}





#Test 3


cnvAssoc_2by2_DELETIONS <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
    
  #Run the 2 by 2
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0){
  cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
  }
  
  else {
    stop("Single level CNV status in 2_by_2_DELETIONS")
  }
  
  
  
  cnv_file_unique <- cnv_file_unique[cnv_file_unique$Probe_CNVstatus <= 2,]
  
  
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int + factor(tissue) + factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  
  assoc_out_mx[temp_pos,17] <- or[2,1]
  assoc_out_mx[temp_pos,18] <- or[2,2]
  assoc_out_mx[temp_pos,19] <- or[2,3]
  assoc_out_mx[temp_pos,20] <- as.numeric(as.character(coef(s)[2,4]))
  
  return(assoc_out_mx)
  
}



#Test 4

cnvAssoc_2by2_DUPLICATIONS <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
    
  }

  #Run the 2 by 2
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus) != 0){
  cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus <- 3
  }
  
  else {
    stop("Single level CNV status in 2_by_2_DUPLICATIONS")
  }
  
  
  cnv_file_unique <- cnv_file_unique[cnv_file_unique$Probe_CNVstatus >= 2,]
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  #cnv_file_unique$final_type <- factor(cnv_file_unique$final_type)
  
  #COMMENT
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int + factor(tissue) + factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  #model_probe1 <- speedglm(factor(final_type) ~ factor(Probe_CNVstatus) + age_int + factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  
  assoc_out_mx[temp_pos,21] <- or[2,1]
  assoc_out_mx[temp_pos,22] <- or[2,2]
  assoc_out_mx[temp_pos,23] <- or[2,3]
  assoc_out_mx[temp_pos,24] <- as.numeric(as.character(coef(s)[2,4]))
  
  return(assoc_out_mx)
  
}




############################
########FUNCTIONS########
############################


CNV_pop_frequency <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,threshold){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- probe_file[probe_pos,]$Loci.start
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  l <- length(sampleID_atProbe_loc)/nrow(cnv_file_unique)
  
  #DEBUG
  #ifelse(l >= threshold,return(l),return())
  #ifelse(l >= threshold,return(i),return())
  #
  
  ifelse(l <= threshold,return(i),return())
  
}




# GENERATE LOCATION OF PROBES WITHIN CNV's WITH FREQUENCY >= 0.25%

sink("AF_CNVs_SAGE_commonProbeLocations.txt")

for(i in 1:1100000){t <- CNV_pop_frequency(i,sage_probeID,cnv_calls_file = sage,cnv_file_unique = sage_test_clean,threshold = 0.0025); 
                  ifelse(t > 0, print(t),print("NOTHING"))
 }

sink()




###Second version of CNV frequency calculator (for 4 categories of rare variants)###
AA_rare_locations <- array(NA, dim=c(1100000,4))
EA_rare_locations <- array(NA, dim=c(1100000,4))

ozalc_EA_rare_locations <- array(NA, dim=c(370500,5))

cidr_AA_rare_locations <- array(NA, dim=c(1100000,5))
cidr_EA_rare_locations <- array(NA, dim=c(1100000,5))

CNV_pop_frequency_rare <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique){
  
  a <- array(NA, dim=c(1,4))
  
  loc <- probe_file[probe_pos,]$Loci.start
  
  #ADD CHROMOSOME INFORMATION
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  l <- length(sampleID_atProbe_loc)
  
  #DEBUG
  #ifelse(l >= threshold,return(l),return())
  #ifelse(l >= threshold,return(i),return())
  #
  
  
  ifelse(l > 0 & l <= 5, a[1,1] <- probe_pos,print("Nothing 1"))
  ifelse(l > 0 & l <= 20, a[1,2] <- probe_pos,print("Nothing 2"))
  ifelse(l > 0 & l <= 40, a[1,3] <- probe_pos,print("Nothing 3"))
  ifelse(l > 0 & l <= 100, a[1,4] <- probe_pos,print("Nothing 4"))
  
  return(a)
}







CNV_pop_frequency_4rare_common <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique){
  
  arr <- array(NA, dim=c(1,5))
  
  loc <- probe_file[probe_pos,]$Loci.start
  
  #ADD CHROMOSOME INFORMATION
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  l <- length(sampleID_atProbe_loc)
  
  #DEBUG
  #ifelse(l >= threshold,return(l),return())
  #ifelse(l >= threshold,return(i),return())
  #
  
  
  ifelse(l > 0 & l <= 5, arr[1,1] <- probe_pos,print("Nothing 1"))
  ifelse(l > 0 & l <= 20, arr[1,2] <- probe_pos,print("Nothing 2"))
  ifelse(l > 0 & l <= 40, arr[1,3] <- probe_pos,print("Nothing 3"))
  ifelse(l > 0 & l <= 100, arr[1,4] <- probe_pos,print("Nothing 4"))
  ifelse(l > 5, arr[1,5] <- probe_pos,print("Nothing 5"))
  
  return(arr)
}




#SAGE: EA + AA
for(i in 1:1100000){
                    
                    t <- CNV_pop_frequency_rare(i,sage_probeID,sage_EA,cnv_file_unique = sage_test_clean_EA) 
                    t -> EA_rare_locations[i,]
                    
                    t <- CNV_pop_frequency_rare(i,sage_probeID,sage_AA,cnv_file_unique = sage_test_clean_AA) 
                    t -> AA_rare_locations[i,]
                    
}


write.table(EA_rare_locations,"SAGE_EA_rare_locations.txt")

# 
# 
# for(i in 1:1100000){t <- CNV_pop_frequency_rare(i,sage_probeID,sage_AA,cnv_file_unique = sage_test_clean_AA) 
#                     t -> AA_rare_locations[i,]
# }

write.table(AA_rare_locations,"SAGE_AA_rare_locations.txt")


#CIDR: EA + AA 

for(i in 567977:1100000){
                    
                    t <- CNV_pop_frequency_4rare_common(i,cidr_probeID,cidr_EA,cnv_file_unique = cidr_test_clean_EA) 
                    cidr_EA_rare_locations[i,] <- t
                    
                    t <- CNV_pop_frequency_4rare_common(i,cidr_probeID,cidr_AA,cnv_file_unique = cidr_test_clean_AA) 
                    t -> cidr_AA_rare_locations[i,]
                    
}


write.table(cidr_EA_rare_locations,"CIDR_EA_rare_locations.txt")



# for(i in 1:1100000){t <- CNV_pop_frequency_4rare_common(i,cidr_probeID,cidr_AA,cnv_file_unique = cidr_test_clean_AA) 
#                     t -> cidr_AA_rare_locations[i,]
# }

write.table(cidr_AA_rare_locations,"CIDR_AA_rare_locations.txt")


#OZALC: EA

for(i in 1:370500){t <- CNV_pop_frequency_4rare_common(i,ozalc_probeID,ozalc_EA,cnv_file_unique = ozalc_test_clean_EA) 
                    ozalc_EA_rare_locations[i,] <- t
}

write.table(ozalc_EA_rare_locations,"OZALC_rare_locations.txt")



############################################
########Association Testing AREA############
############################################

#EA SAMPLES
# 
# for (i in 1:200000){
#   
#   tryCatch( {
#   
#     n <- sage_common_probes[i,]
#       
#     a <- suppressWarnings(cnvAssoc_2by2(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
#     b <- suppressWarnings(cnvAssoc_2by3(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
#     c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
#     d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
#     
#     sage_out_mx_EA[i,c(5:8)] <- a[n,c(5:8)]
#     sage_out_mx_EA[i,c(1:4,9:16)] <- b[n,c(1:4,9:16)]
#     sage_out_mx_EA[i,c(17:20)] <- c[n,c(17:20)]
#     sage_out_mx_EA[i,c(21:24)] <- d[n,c(21:24)]
#     
#     print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
#   )
# }
# 



#SAGE EA working version

system.time(

for (i in 1:24000){
  
  n <- sage_common_probes[i,]
  
  tryCatch( {
    
    a <- suppressWarnings(cnvAssoc_2by2(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
    sage_out_mx_EA[i,c(1:4,5:8)] <- a[1,c(1:4,5:8)]
   
    print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
  tryCatch( {  
    b <- suppressWarnings(cnvAssoc_2by3(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
    sage_out_mx_EA[i,c(1:4,9:16)] <- b[1,c(1:4,9:16)]
  
    print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
  tryCatch( {  
    c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
    sage_out_mx_EA[i,c(1:4,17:20)] <- c[1,c(1:4,17:20)]
  
    print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
  tryCatch( {
    d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,sage_probeID,sage_EA,sage_test_clean_EA,sage_out_mx_EA))
    sage_out_mx_EA[i,c(1:4,21:24)] <- d[1,c(1:4,21:24)]
    print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
}

)




#SAGE AA working version


system.time(
  
  for (i in 1:24000){
    
    n <- sage_common_probes[i,]
    
    tryCatch( {
      
      a <- suppressWarnings(cnvAssoc_2by2(n,sage_probeID,sage_AA,sage_test_clean_AA,sage_out_mx_AA))
      sage_out_mx_AA[i,c(1:4,5:8)] <- a[1,c(1:4,5:8)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      b <- suppressWarnings(cnvAssoc_2by3(n,sage_probeID,sage_AA,sage_test_clean_AA,sage_out_mx_AA))
      sage_out_mx_AA[i,c(1:4,9:16)] <- b[1,c(1:4,9:16)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,sage_probeID,sage_AA,sage_test_clean_AA,sage_out_mx_AA))
      sage_out_mx_AA[i,c(1:4,17:20)] <- c[1,c(1:4,17:20)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {
      d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,sage_probeID,sage_AA,sage_test_clean_AA,sage_out_mx_AA))
      sage_out_mx_AA[i,c(1:4,21:24)] <- d[1,c(1:4,21:24)]
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
  }
      
)



#CIDR EA

#Start at 400
system.time(
  
  for (i in 1:16000){
    
    n <- cidr_common_probes_EA[i]
    
    tryCatch( {
      
      a <- suppressWarnings(cnvAssoc_2by2(n,cidr_probeID,cidr_EA,cidr_test_clean_EA,cidr_out_mx_EA))
      cidr_out_mx_EA[i,c(1:4,5:8)] <- a[1,c(1:4,5:8)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      b <- suppressWarnings(cnvAssoc_2by3(n,cidr_probeID,cidr_EA,cidr_test_clean_EA,cidr_out_mx_EA))
      cidr_out_mx_EA[i,c(1:4,9:16)] <- b[1,c(1:4,9:16)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,cidr_probeID,cidr_EA,cidr_test_clean_EA,cidr_out_mx_EA))
      cidr_out_mx_EA[i,c(1:4,17:20)] <- c[1,c(1:4,17:20)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {
      d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,cidr_probeID,cidr_EA,cidr_test_clean_EA,cidr_out_mx_EA))
      cidr_out_mx_EA[i,c(1:4,21:24)] <- d[1,c(1:4,21:24)]
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
  }
  
)


#CIDR AA

system.time(
  
  for (i in 1:10000){
    
    n <- cidr_common_probes_AA[i]
    
    tryCatch( {
      
      a <- suppressWarnings(cnvAssoc_2by2(n,cidr_probeID,cidr_AA,cidr_test_clean_AA,cidr_out_mx_AA))
      cidr_out_mx_AA[i,c(1:4,5:8)] <- a[1,c(1:4,5:8)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      b <- suppressWarnings(cnvAssoc_2by3(n,cidr_probeID,cidr_AA,cidr_test_clean_AA,cidr_out_mx_AA))
      cidr_out_mx_AA[i,c(1:4,9:16)] <- b[1,c(1:4,9:16)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,cidr_probeID,cidr_AA,cidr_test_clean_AA,cidr_out_mx_AA))
      cidr_out_mx_AA[i,c(1:4,17:20)] <- c[1,c(1:4,17:20)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {
      d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,cidr_probeID,cidr_AA,cidr_test_clean_AA,cidr_out_mx_AA))
      cidr_out_mx_AA[i,c(1:4,21:24)] <- d[1,c(1:4,21:24)]
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
  }
  
)



#OZALC EA
system.time(
  
  for (i in 1:8500){
    
    n <- ozalc_common_probes[i]
    
    tryCatch( {
      
      a <- suppressWarnings(cnvAssoc_2by2(n,ozalc_probeID,ozalc_EA,ozalc_test_clean_EA,ozalc_out_mx_EA))
      ozalc_out_mx_EA[i,c(1:4,5:8)] <- a[1,c(1:4,5:8)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      b <- suppressWarnings(cnvAssoc_2by3(n,ozalc_probeID,ozalc_EA,ozalc_test_clean_EA,ozalc_out_mx_EA))
      ozalc_out_mx_EA[i,c(1:4,9:16)] <- b[1,c(1:4,9:16)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {  
      c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(n,ozalc_probeID,ozalc_EA,ozalc_test_clean_EA,ozalc_out_mx_EA))
      ozalc_out_mx_EA[i,c(1:4,17:20)] <- c[1,c(1:4,17:20)]
      
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
    tryCatch( {
      d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(n,ozalc_probeID,ozalc_EA,ozalc_test_clean_EA,ozalc_out_mx_EA))
      ozalc_out_mx_EA[i,c(1:4,21:24)] <- d[1,c(1:4,21:24)]
      print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
    )
    
  }
  
)



###########################################
##########Meta-Analysis Area###############
###########################################

#Idea: for each shared probe between studies create table of 5 rows (total population-specific studies) and 6 columns.
#The 6 columns will be the 4 cells of 2x2 table (ai,bi,ci,di) and 2 covariates: average age, percent male.


#Prepare data for the meta-analysis

#Step 1: remove overlap between samples (CIDR and SAGE). Keep duplicated samples in CIDR
#overlap_sage_cidr <- read.table("CIDR_SAGE_overlapped.txt",header=F)
overlap_sage_cidr <- read.table("CIDR_SAGE_overlapped_v2.txt",header=F)

#AA
arr_AA<- array(NA, dim=c(1360,1))

for (i in 1:1360){
  tryCatch( {  
  arr_AA[i,1] <- which(sage_test_clean_AA[,"chip.id"]==toString(overlap_sage_cidr[i,"V1"]))

  print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
}

arr_AA <- na.omit(arr_AA)

#EA 
arr_EA<- array(NA, dim=c(1360,1))

for (i in 1:1360){
  tryCatch( {  
    arr_EA[i,1] <- which(sage_test_clean_EA[,"chip.id"]==toString(overlap_sage_cidr[i,"V1"]))
    
    print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
}

arr_EA <- na.omit(arr_EA)

#Remove rows with duplicated samples from CIDR
#cidr_test_clean_AA_META <- cidr_test_clean_AA[-c(arr_AA),]
#cidr_test_clean_EA_META <- cidr_test_clean_EA[-c(arr_EA),]

#Remove rows with duplicated samples from SAGE
sage_test_clean_AA_META <- sage_test_clean_AA[-c(arr_AA),]
sage_test_clean_EA_META <- sage_test_clean_EA[-c(arr_EA),]




#Step 2: find overlap between probe locations across studies
#i) overlap between SAGE_AA and CIDR_AA

AA_shared_probes <- intersect(sage_probeID[sage_common_probes$V1,"ID"],cidr_probeID[cidr_common_probes_AA,"ID"])
#write.table(AA_shared_probes,"AA_shared_probes.txt") 

#ii) overlap between SAGE_EA, CIDR_EA and OZALC
EA_shared_probes <- intersect(sage_probeID[sage_common_probes$V1,"ID"],cidr_probeID[cidr_common_probes_EA,"ID"])

#**********************
EA_shared_probe_file <- array(NA,dim=c(1,16000))


for (i in 1:16000){
  tryCatch( {
EA_shared_probe_file[i,] <- sage_probeID[which(sage_probeID$ID==EA_shared_probes[i]),]

print(i)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )

}


#**********************

#write.table(EA_shared_probes,"EA_shared_probes.txt") 


#iii) overlap between SAGE, CIDR and OZALC
ThreeStudies_shared_probes_tmp <- intersect(EA_shared_probes,AA_shared_probes)
ThreeStudies_shared_probes <- intersect(ozalc_probeID[ozalc_common_probes,"ID"],ThreeStudies_shared_probes_tmp)
#write.table(ThreeStudies_shared_probes,"ThreeStudies_shared_probes.txt")


meta_out_array <- array(NA,dim=c(5,4))

#Test 1: AA samples -> sage_AA and cidr_AA 
n <- length(AA_shared_probes)

for (i in 1:n){
  #find out how many cases/controls in SAGE_AA have/don't have CNVs at probe i
    
  
  #find out how many cases/controls in CIDR_AA have/don't have CNVs at probe i
  
}


#Test 2: EA samples -> sage_EA, cidr_EA and ozalc


#Test 3: AA and EA samples



cnvAssoc_meta1 <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique_meta,meta_assoc_out_mx){
  
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <- probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  
  #*
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  #* 
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
  
  #Run the 2 by 3
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0) {
    
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
    
  }
  
  else {
    stop("Single level CNV status in 2_by_3 DELETIONS")
  }
  
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus) != 0){
    
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus <- 3
  }
  
  
  else {
    stop("Single level CNV status in 2_by_3 DUPLICATIONS")
  }
  
  
  
  
  #COnsider CNV status to be a factor
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique, Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  #keep only non-missing phenotype samples
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ Probe_CNVstatus + age_int + factor(tissue) + factor(sex.x) + PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma(log))
  
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  assoc_out_mx[temp_pos,9] <- or[2,1]
  assoc_out_mx[temp_pos,10] <- or[2,2]
  assoc_out_mx[temp_pos,11] <- or[2,3]
  assoc_out_mx[temp_pos,12] <- as.numeric(as.character(coef(s)[2,4]))
  assoc_out_mx[temp_pos,13] <- or[3,1] 
  assoc_out_mx[temp_pos,14] <- or[3,2]
  assoc_out_mx[temp_pos,15] <- or[3,3]
  assoc_out_mx[temp_pos,16] <- as.numeric(as.character(coef(s)[3,4]))
  
  return(assoc_out_mx)
}









######FOR CIDR ONLY (no factor(tissue) in the logistic regression model)######



#Test 1

cnvAssoc_2by3 <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <- probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  
  #*
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  #* 
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
  
  #Run the 2 by 3
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0) {
    
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
    
  }
  
  else {
    stop("Single level CNV status in 2_by_3 DELETIONS")
  }
  
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus) != 0){
    
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus <- 3
  }
  
  
  else {
    stop("Single level CNV status in 2_by_3 DUPLICATIONS")
  }
  
  
  
  
  #COnsider CNV status to be a factor
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique, Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  #keep only non-missing phenotype samples
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ Probe_CNVstatus + age_int +  factor(sex.x) + PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma(log))
  
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  assoc_out_mx[temp_pos,9] <- or[2,1]
  assoc_out_mx[temp_pos,10] <- or[2,2]
  assoc_out_mx[temp_pos,11] <- or[2,3]
  assoc_out_mx[temp_pos,12] <- as.numeric(as.character(coef(s)[2,4]))
  assoc_out_mx[temp_pos,13] <- or[3,1] 
  assoc_out_mx[temp_pos,14] <- or[3,2]
  assoc_out_mx[temp_pos,15] <- or[3,3]
  assoc_out_mx[temp_pos,16] <- as.numeric(as.character(coef(s)[3,4]))
  
  return(assoc_out_mx)
}




#Test 2

cnvAssoc_2by2 <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
  
  #Run the 2 by 2
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2 | cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0){ 
    
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2 | cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
    
  }
  
  else {
    stop("Single level CNV status in 2_by_2")
  }
  
  
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int +  factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  
  assoc_out_mx[temp_pos,5] <- or[2,1]
  assoc_out_mx[temp_pos,6] <- or[2,2]
  assoc_out_mx[temp_pos,7] <- or[2,3]
  assoc_out_mx[temp_pos,8] <- as.numeric(as.character(coef(s)[2,4]))
  
  return(assoc_out_mx)
  
}





#Test 3


cnvAssoc_2by2_DELETIONS <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
  }
  
  
  #Run the 2 by 2
  
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus) != 0){
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus<2,]$Probe_CNVstatus <- 1
  }
  
  else {
    stop("Single level CNV status in 2_by_2_DELETIONS")
  }
  
  
  
  cnv_file_unique <- cnv_file_unique[cnv_file_unique$Probe_CNVstatus <= 2,]
  
  
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int +  factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  
  assoc_out_mx[temp_pos,17] <- or[2,1]
  assoc_out_mx[temp_pos,18] <- or[2,2]
  assoc_out_mx[temp_pos,19] <- or[2,3]
  assoc_out_mx[temp_pos,20] <- as.numeric(as.character(coef(s)[2,4]))
  
  return(assoc_out_mx)
  
}



#Test 4

cnvAssoc_2by2_DUPLICATIONS <- function(probe_pos,probe_file,cnv_calls_file,cnv_file_unique,assoc_out_mx){
  
  cnv_file_unique$Probe_CNVstatus <- rep(2,nrow(cnv_file_unique))
  
  loc <- as.numeric(probe_file[probe_pos,]$Loci.start)
  
  #ADD CHROMOSOME INFORMATION!!!!!!
  loc_chr <-probe_file[probe_pos, ]$Chr
  
  
  sampleID_atProbe_loc <- unique(cnv_calls_file[which(cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc &  cnv_calls_file$Chr ==  paste0("chr",loc_chr)),]$chip.id)
  sampleID_atProbe_loc <- intersect(sampleID_atProbe_loc,unique(cnv_file_unique$chip.id))
  
  
  for (i in 1:length(sampleID_atProbe_loc)) {
    t <- cnv_calls_file[which(cnv_calls_file$chip.id==sampleID_atProbe_loc[i] & cnv_calls_file$Start <= loc & cnv_calls_file$End >= loc),]$Type1
    cnv_file_unique[which(cnv_file_unique$chip.id==sampleID_atProbe_loc[i]),]$Probe_CNVstatus <- t[1]
    
  }
  
  #Run the 2 by 2
  if(length(cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus) != 0){
    cnv_file_unique[cnv_file_unique$Probe_CNVstatus>2,]$Probe_CNVstatus <- 3
  }
  
  else {
    stop("Single level CNV status in 2_by_2_DUPLICATIONS")
  }
  
  
  cnv_file_unique <- cnv_file_unique[cnv_file_unique$Probe_CNVstatus >= 2,]
  cnv_file_unique$Probe_CNVstatus <- factor(cnv_file_unique$Probe_CNVstatus)
  #cnv_file_unique$final_type <- factor(cnv_file_unique$final_type)
  
  #COMMENT
  cnv_file_unique <- within(cnv_file_unique,Probe_CNVstatus <- relevel(Probe_CNVstatus,ref=2))
  
  
  #Run logistic regression model
  #model_probe1 <- glm(factor(final_type)~ Probe_CNVstatus+age_int+factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family="binomial")
  model_probe1 <- speedglm(final_type ~ factor(Probe_CNVstatus) + age_int +  factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  #model_probe1 <- speedglm(factor(final_type) ~ factor(Probe_CNVstatus) + age_int + factor(sex.x)+PC1+PC2+PC3+PC4+PC5,data=cnv_file_unique,family = Gamma())
  s <- summary(model_probe1)
  or <- exp(cbind(OR = coef(model_probe1), confint(model_probe1)))
  
  
  temp_pos = 1
  
  assoc_out_mx[temp_pos,1] <- probe_pos
  assoc_out_mx[temp_pos,2] <- probe_file[probe_pos,]$Chr
  assoc_out_mx[temp_pos,3] <- probe_file[probe_pos,]$Loci.start
  assoc_out_mx[temp_pos,4] <- as.character(probe_file[probe_pos,]$ID)
  
  
  assoc_out_mx[temp_pos,21] <- or[2,1]
  assoc_out_mx[temp_pos,22] <- or[2,2]
  assoc_out_mx[temp_pos,23] <- or[2,3]
  assoc_out_mx[temp_pos,24] <- as.numeric(as.character(coef(s)[2,4]))
  
  return(assoc_out_mx)
  
}





########################END######################################










##############TEST#2##################



for (i in 1:10){
  
    
    a <- suppressWarnings(cnvAssoc_2by2(i,sage_probeID,sage,sage_test_clean,sage_out_mx))
    b <- suppressWarnings(cnvAssoc_2by3(i,sage_probeID,sage,sage_test_clean,sage_out_mx))
    c <- suppressWarnings(cnvAssoc_2by2_DELETIONS(i,sage_probeID,sage,sage_test_clean,sage_out_mx))
    d <- suppressWarnings(cnvAssoc_2by2_DUPLICATIONS(i,sage_probeID,sage,sage_test_clean,sage_out_mx))
    
    sage_out_mx[i,c(5:8)] <- a[i,c(5:8)]
    sage_out_mx[i,c(1:4,9:16)] <- b[i,c(1:4,9:16)]
    sage_out_mx[i,c(17:20)] <- c[i,c(17:20)]
    sage_out_mx[i,c(21:24)] <- d[i,c(21:24)]
  
}
  
  


#TO DO: Read in every data table as ff object


#P_s function

P_s <- function(probe_coords, CNV_calls){
	
	#Generate sequence of CNV coordinates for each CNV region in CNV_calls object
	for (i in seq(1, nrow(CNV_calls))) {
		
		start <- which(probe_coords$V2==toString(CNV_calls[i,]$Loci.start))
		
		end <- which(probe_coords$V2==toString(CNV_calls[i,]$Loci.end))

		print(seq(start,end))
	}

	#Merge generated sequence of CNV coordinates above with those available in probe_coords


}

#Populate output matrix using type information from original CNV calling data

for (s in sage_out[,1]){
	
	for (p in sage_out[1,]) {
		
		if(sage[!is.null(which(sage$chip.id==toString(s) & sage$Loci.start==toString(p))),]$Type1 < 2) {
			sage_out[which(sage_out[,1]==toString(s)),which(sage_out[1,]==toString(p))] <- 1
			
			} else if(sage[!is.null(which(sage$chip.id==toString(s) & sage$Loci.start==toString(p))),]$Type1 > 2){
				sage_out[which(sage_out[,1]==toString(s)),which(sage_out[1,]==toString(p))] <- 3
        
			} else {sage_out[which(sage_out[,1]==toString(s)),which(sage_out[1,]==toString(p))] <- 2
		}
	}
}

