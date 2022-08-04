setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/calling/strelka2")

A549 = "A549_work/annovar/A549_PASS.hg38_multianno_k50.vcf.gz"
MCF7 = "MCF7_work/annovar/MCF7_PASS.hg38_multianno_k50.vcf.gz"
HAP1 = "HAP1_work/annovar/HAP1_PASS.hg38_multianno_k50.vcf.gz"

library(VariantAnnotation)
library(ggpubr)
library(MutationalPatterns)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(dplyr)
library(gridExtra)
source("/g/strcombio/fsupek_cancer1/indels_signatures_project/indels_sig_WGEX/bin/functions_TD.R")
source("/g/strcombio/fsupek_cancer1/indels_signatures_project/indels_sig_WGEX/bin/indel_83_classification_hg38.R") # from Ivan to re-write the classify function
indels83 = read.table("/g/strcombio/fsupek_cancer1/indels_signatures_project/83_indel_categories.csv", sep=",", h=T, stringsAsFactors = F)[,2]
snvs96 = rownames(readRDS(system.file("states/mut_mat_data.rds", package = "MutationalPatterns")  ))

nb_positive_samples = function(id, matrix_geno, matrix_filter){
  genos = matrix_geno[id, ]
  filts = matrix_filter[id, ]
  genos_kept = genos #[which(filts == "PASS")] #add this will keep a lot of germline with VAF~1
  sum ( genos_kept!="0/0" & genos_kept!="0|0" ) # better to remove if found in at least 2 different conditions?
  #sum ( genos_kept!="0/0" & genos_kept!="." & genos_kept!="0|0" ) # better to remove if found in at least 2 different conditions?
  #sum ( genos_kept!="0/0" & genos_kept!="." & genos_kept!="0|0" ) + sum( (genos_kept=="0/0" | genos_kept=="." | genos_kept=="0|0" ) & filts != "PASS")
}

type_mut = function(vect_mut){ # vector of "chr1:283492_ACT/A". returns if it is SNV, INS or DEL (one per input element)
  unlist(lapply(vect_mut, function(v){
    tmp = unlist(strsplit(v, "_"))[2]
    ref = unlist(strsplit(tmp, "/"))[1]
    alt = unlist(strsplit(tmp, "/"))[2]
    if(nchar(ref) > 1) {return("DEL")} else {
      if(nchar(alt) > 1) {return("INS")} else {return("SNV")}
    } 
  }))
}

# for the mutational profiles analysis
genome = "hg38"
ref_genome = available.genomes()[which(grepl(genome,available.genomes()) & ! grepl("masked",available.genomes()))]
human_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
               "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX", "chrY", "chrM")

remove_lowVAF = FALSE
vaf_thr = 0.25

for(cellline in c("A549", "HAP1", "MCF7")){
  print(paste(date(), " INFO: working on cell line: ", cellline, sep=""))
  vcf = open(VcfFile(get(cellline),  yieldSize=50000))
  vcf_chunk = readVcf(vcf, "hg38")
  i=1
  
  while(dim(vcf_chunk)[1] != 0) {
    seqlevels(vcf_chunk) = human_chrs
    print(paste(date(), "working on chunk ", i, sep=""))
    gg0 = geno(vcf_chunk, "GT")
    gnomd0 = info(vcf_chunk)[,"gnomAD_genome_ALL"]
    # to compute the VAF:
    ad0 = geno(vcf_chunk, "AD")
    dp0 = geno(vcf_chunk, "DP")
    ft0 = geno(vcf_chunk, "FT")
    # initialize vectors of zeros
    if(i == 1){
      allsm = colnames(gg0)
      allsm = unlist(lapply(allsm, function(s) gsub("_REP1","",s)))
      for(s in allsm) { for(type in c("SNV", "INS", "DEL")) { assign(paste(s, type, sep="_"), 0 ) } }
    }
    #ids_keep = which( apply(gg0, 1, function(r) nb_positive_samples(r)) < 2 )
    ids_keep = which( unlist(lapply(1:nrow(gg0), function(id) nb_positive_samples(id, gg0, ft0))) < 2 )
    ids_keep = intersect(ids_keep, which(gnomd0 < 0.001 | is.na(gnomd0) )) # if not found, consider MAF=0
    gg = gg0[ids_keep, ]
    ad = ad0[ids_keep, ]
    dp = dp0[ids_keep, ]
    ft = ft0[ids_keep, ]
    vcf_chunk_filt = vcf_chunk[ids_keep, ]

    # compute the vectors of indels and SNVs observed in the chunk -- for VAF and remove SNVs
    indels_id = which(nchar(unlist(lapply(rownames(gg), function(x) unlist(strsplit(x ,"_"))[2])))>3)
    indels_names = rownames(gg)[indels_id]
    snvs_id = setdiff(1:nrow(gg), indels_id)
    
    # compute VAF for SNVs in that chunk
    for(sm in allsm){
      # compute AF, DP, FT and GG vectors for the SNVs 
      adtmp = ad[snvs_id, grepl(sm, colnames(ad))] ; dptmp = dp[snvs_id, grepl(sm, colnames(ad))] 
      ggtmp = gg[snvs_id, grepl(sm, colnames(gg))] ; fttmp = ft[snvs_id, grepl(sm, colnames(ft))]
      mutsnv = names(ggtmp)[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS")] # list of snvs where the sample is mutated and PASS
      if(i==1) { assign(paste("all_SNVs",sm,sep="_"), mutsnv) } else {assign(paste("all_SNVs",sm,sep="_"), c( get(paste("all_SNVs",sm,sep="_")), mutsnv)) }
      if(length(mutsnv) > 0){ # in case the sample is not mutated in this chunk
        vaftmp = unlist(lapply(mutsnv, function(curmut) adtmp[[curmut]][[2]] / dptmp[curmut]) ) # compute the VAF of each obs snv
        names(vaftmp) = mutsnv
        if(remove_lowVAF) { 
          kept_snvs_vaf = names(vaftmp)[which(vaftmp > vaf_thr)] 
          vaftmp = vaftmp[which(vaftmp > vaf_thr)]
        } else { kept_snvs_vaf = mutsnv} # which SNVs we should kept
        if(i==1) { assign(paste("all_VAF",sm,sep="_"), vaftmp) } else {assign(paste("all_VAF",sm,sep="_"), c( get(paste("all_VAF",sm,sep="_")), vaftmp)) }
      } else { kept_snvs_vaf = c() } # if no mutations for that sample
      
      # count number of each mutation types for the chunk
      ggtmp2 = gg[c(kept_snvs_vaf, indels_names), which(grepl(sm, colnames(gg)))] # keep good SNVs + all indels
      fttmp2 = ft[c(kept_snvs_vaf, indels_names), which(grepl(sm, colnames(ft)))] # keep good SNVs + all indels
      ggtmpmut = ggtmp2[which(ggtmp2!="0/0" & ggtmp2!="." & ggtmp2!="0|0" & fttmp2 == "PASS")]
      type_mut_tmp = type_mut(names(ggtmpmut))
      if(!is.null(type_mut_tmp)){ # if no SNV or indels in that sample
        for(type in c("SNV", "INS", "DEL")) { 
          if(!is.na(table(type_mut_tmp)[type])) {assign(paste(sm, type, sep="_"), get(paste(sm, type, sep="_")) +  as.numeric(table(type_mut_tmp)[type]))}else{
            assign(paste(sm, type, sep="_"), get(paste(sm, type, sep="_")) + 0) # if zero variant for that sample and that type of mutation
          }
        }
      }
      
      # compute 96 nucleotides context for both SNVs and indels
      for(type in c("SNV", "indel")){
        if(type=="SNV") {ids_cur = snvs_id} else {ids_cur = indels_id}
        vcf_chunk_cur = vcf_chunk_filt[ids_cur,] # get SNVs or indels
        ggcur = gg[ids_cur, grepl(sm, colnames(gg))] ; ftcur = ft[ids_cur, grepl(sm, colnames(ft))] # get SNVs or indels
        mutcur_ids = which(ggcur!="0/0" & ggcur!="." & ggcur!="0|0" & ftcur == "PASS") # get mutations PASS in the sample
        if(length(mutcur_ids) > 0){
          if(type == "SNV"){
            nut3_context = type_context(rowRanges(vcf_chunk_cur[mutcur_ids,]), ref_genome) # get nut3 context for mutations PASS in the sample, SNVs or indels
            all_mut96 = paste(substr(nut3_context$context,1,1), "[", nut3_context$types, "]", substr(nut3_context$context,3,3), sep = "") # 96 nut
            all_mut96_lowVAF = all_mut96[which(vaftmp <= vaf_thr)] # 96 nut for low VAF
            all_mut96_highVAF = all_mut96[which(vaftmp > vaf_thr)] # 96 nut for high VAF
            if(i==1) { assign(paste("mut96",type,sm,sep="_"), all_mut96) } else {assign(paste("mut96",type,sm,sep="_"), c( get(paste("mut96",type,sm,sep="_")), all_mut96)) }
            if(i==1) { assign(paste("mut96",type,"lowVAF",sm,sep="_"), all_mut96_lowVAF) } else {assign(paste("mut96",type,"lowVAF",sm,sep="_"), c( get(paste("mut96",type,"lowVAF",sm,sep="_")), all_mut96_lowVAF)) }
            if(i==1) { assign(paste("mut96",type,"highVAF",sm,sep="_"), all_mut96_highVAF) } else {assign(paste("mut96",type,"highVAF",sm,sep="_"), c( get(paste("mut96",type,"highVAF",sm,sep="_")), all_mut96_highVAF)) }
          }
          # UNCOMMENT TO BUILD INDELS CLASSES
          if(type == "indel"){
            datindels = data.frame(CHROM = unlist(lapply(names(mutcur_ids), function(x) unlist(strsplit(x, ":"))[1])),
                                   POS = as.numeric(unlist(lapply(names(mutcur_ids), function(x) {tmp = unlist(strsplit(x, ":"))[2]; unlist(strsplit(tmp,"_"))[1]}))),
                                   REF = unlist(lapply(names(mutcur_ids), function(x) {tmp = unlist(strsplit(x, ":"))[2]; tmp2 = unlist(strsplit(tmp,"_"))[2]; unlist(strsplit(tmp2,"/"))[1]})),
                                   ALT = unlist(lapply(names(mutcur_ids), function(x) {tmp = unlist(strsplit(x, ":"))[2]; tmp2 = unlist(strsplit(tmp,"_"))[2]; unlist(strsplit(tmp2,"/"))[2]})),
                                   ID="rsX", QUAL=".", FILTER="PASS")
            if(nrow(datindels) > 0 ){
              #datindels$REF = as.character(datindels$REF) ; datindels$ALT = as.character(datindels$ALT)
              #indels_class = indel_83_classification(datindels)$ID83
              #if(i==1) { assign(paste("mut83",type,sm,sep="_"), indels_class) } else {assign(paste("mut83",type,sm,sep="_"), c( get(paste("mut83",type,sm,sep="_")), indels_class)) }
              if(i==1) {assign(paste("indelmatVCF_", sm, sep=""), datindels)} else {assign(paste("indelmatVCF_", sm, sep=""), rbind(get(paste("indelmatVCF_", sm, sep="")), datindels))}
            }
          }
        }
      }
    }
    
    i = i + 1
    vcf_chunk = readVcf(vcf, "hg38")
  } # here the whole VCF has been read
  print(paste(date(), " INFO: the whole VCF has been read"))
  # compute a dataframe per sample containing number of mutations per each type
  for(s in allsm){
    for(type in c("SNV", "INS", "DEL")){
      d = data.frame(clone = s, nb_mut = get(paste(s, type, sep="_")), type = type)
      if(s == allsm[1] & type == "SNV") {finald = d} else {finald = rbind(finald, d)}
    }
  }
  # assign plot of number of mutation per type to each sample
  # do not conserve the MCF7 control
  assign(paste("p_", cellline, ifelse(remove_lowVAF, "filtered", ""), sep=""),
         ggbarplot(finald[which(finald$clone != "MCF7_CONTROL_B"),], x="clone", y="nb_mut", fill="type", label=T) + ggtitle(cellline) +  rotate_x_text(45)
  )
  assign(paste("finald_", cellline, sep=""), finald)
  # build the matrix of 96 mutation types
  for(s in allsm){
    dt = table( get(paste("mut96_SNV", s, sep="_")) )
    t = rep(0, length(setdiff(snvs96, names(dt)))) ; names(t) = setdiff(snvs96, names(dt)) ; dt = c(dt, t) ; dt = dt[snvs96]
    dt_lowVAF = table( get(paste("mut96_SNV_lowVAF", s, sep="_")) ) # for low and high VAF we should add 0 when one of the 96nut is not seen
    t = rep(0, length(setdiff(names(dt), names(dt_lowVAF)))) ; names(t) = setdiff(names(dt), names(dt_lowVAF)) ; dt_lowVAF = c(dt_lowVAF, t) ; dt_lowVAF = dt_lowVAF[names(dt)]
    dt_highVAF = table( get(paste("mut96_SNV_highVAF", s, sep="_")) ) 
    t = rep(0, length(setdiff(names(dt), names(dt_highVAF)))) ; names(t) = setdiff(names(dt), names(dt_highVAF)) ; dt_highVAF = c(dt_highVAF, t) ; dt_highVAF = dt_highVAF[names(dt)]
    dd = as.numeric(dt) ; dd_lowVAF = as.numeric(dt_lowVAF) ; dd_highVAF = as.numeric(dt_highVAF)
    if (match(s, allsm) == 1) { ddf = data.frame(dd) ; rownames(ddf) = names(dt) } else { ddf = cbind( ddf, as.numeric(dd)) }
    if (match(s, allsm) == 1) { ddf_lowVAF = data.frame(dd_lowVAF) ; rownames(ddf_lowVAF) = names(dt_lowVAF) } else { ddf_lowVAF = cbind( ddf_lowVAF, as.numeric(dd_lowVAF)) }
    if (match(s, allsm) == 1) { ddf_highVAF = data.frame(dd_highVAF) ; rownames(ddf_highVAF) = names(dt_highVAF) } else { ddf_highVAF = cbind( ddf_highVAF, as.numeric(dd_highVAF)) }
  }
  colnames(ddf) = allsm ; colnames(ddf_lowVAF) = allsm ; colnames(ddf_highVAF) = allsm
  assign(paste("mat96_", cellline, sep="") , ddf)
  assign(paste("mat96_", cellline, "_lowVAF", sep="") , ddf_lowVAF)
  assign(paste("mat96_", cellline, "_highVAF",  sep="") , ddf_highVAF)
  
  # build the matrix of the 83 indel types
  # UNCOMMENT TO BUILD INDEL CLASSES
  # for(s in allsm){
  #   dt = table( get(paste("mut83_indel", s, sep="_")) )
  #   t = rep(0, length(setdiff(indels83, names(dt)))) ; names(t) = setdiff(indels83, names(dt)) ; dt = c(dt, t) ; dt = dt[indels83]
  #   dd = as.numeric(dt)
  #   if (match(s, allsm) == 1) { ddf = data.frame(dd) ; rownames(ddf) = names(dt) } else { ddf = cbind( ddf, as.numeric(dd)) }
  # }
  # colnames(ddf) = allsm
  # assign(paste("mat83indels_", cellline, sep="") , ddf)
}

save(mat83indels_A549, file="Rdata/mat83indels_A549_k50.Rdata")
save(mat83indels_HAP1, file="Rdata/mat83indels_HAP1_k50.Rdata")
save(mat83indels_MCF7, file="Rdata/mat83indels_MCF7_k50.Rdata")

save(mat96_A549, file="Rdata/mat96_A549_k50.Rdata"); save(mat96_A549_lowVAF, file="Rdata/mat96_A549_lowVAF_k50.Rdata"); save(mat96_A549_highVAF, file="Rdata/mat96_A549_highVAF_k50.Rdata")
save(mat96_MCF7, file="Rdata/mat96_MCF7_k50.Rdata"); save(mat96_MCF7_lowVAF, file="Rdata/mat96_MCF7_lowVAF_k50.Rdata"); save(mat96_MCF7_highVAF, file="Rdata/mat96_MCF7_highVAF_k50.Rdata")
save(mat96_HAP1, file="Rdata/mat96_HAP1_k50.Rdata"); save(mat96_HAP1_lowVAF, file="Rdata/mat96_HAP1_lowVAF_k50.Rdata"); save(mat96_HAP1_highVAF, file="Rdata/mat96_HAP1_highVAF_k50.Rdata")

for(cellline in c("A549", "HAP1", "MCF7")){
  alldat = list(ls(pattern = paste("all_SNVs_", cellline, sep="")))[[1]]
  allsm = unlist(lapply(alldat, function(x) unlist(strsplit(x, paste(cellline,"_",sep="")))[2]))
  for(sm in allsm){
    save(list=paste("all_SNVs", cellline, sm, sep = "_"), file=paste("analysis/Rdata/all_SNVs_", cellline, "_", sm, ".Rdata", sep = ""))
  }
}


### PLOTS ###
library(ggridges)
for(smvaf in ls(pattern="all_VAF")){
  print(smvaf)
  if(match(smvaf, ls(pattern="all_VAF")) == 1)  { vafdat = data.frame(sample = gsub("all_VAF_", "", smvaf), VAF = get(smvaf)) } else {
    vafdat = rbind(vafdat, data.frame(sample = gsub("all_VAF_", "", smvaf), VAF = get(smvaf)))
  }
}
for(cellline in c("A549", "HAP1", "MCF7")){
  assign(paste("p_vaf_", cellline, sep=""), 
         ggplot(vafdat[which(grepl(cellline, vafdat$sample) & !grepl("MCF7_CONTROL", vafdat$sample)),], aes(x = VAF, y = sample, fill = sample)) + 
           geom_density_ridges() + theme_ridges() + theme(legend.position = "none") + ggtitle(cellline) )
}

ggplot(vafdat, aes(x = VAF, y = sample, fill = sample)) + geom_density_ridges() + theme_ridges() + theme(legend.position = "none")
grid.arrange(p_vaf_HAP1, p_vaf_A549, p_vaf_MCF7, nrow=1)

grid.arrange(p_HAP1, p_A549, p_MCF7, nrow=1)
grid.arrange(p_HAP1filtered, p_A549filtered, p_MCF7filtered, nrow=1)

## build indels VCF for sigprofiler
setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/calling/strelka2/analysis/")
for(cellline in c("A549", "HAP1", "MCF7")){
  alldat = list(ls(pattern = paste("indelmatVCF_", cellline, sep="")))[[1]]
  allsm = unlist(lapply(alldat, function(x) unlist(strsplit(x, paste(cellline,"_",sep="")))[2]))
  for(sm in allsm){
    vcf = get(alldat[grepl(sm, alldat)])
    vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
    colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
    write.table(c("##fileformat=VCFv4.2"), paste0("data/indels/ours/", cellline, "_", sm,".vcf"), row.names = F,col.names = F,quote = F)
    fwrite(vcf, paste0("data/indels/ours/", cellline, "_", sm,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")
  }
}


# ratios like Behjati paper
for(cellline in c("A549", "HAP1", "MCF7")){
  finald = get(paste("finald_", cellline, sep=""))
  for(type in c("ratioSNV", "ratioINS")){
    for(sm in as.character(unique(finald$clone))){
      nbDEL = finald[which(finald$clone == sm & finald$type == "DEL"), "nb_mut"]
      nbINS = finald[which(finald$clone == sm & finald$type == "INS"), "nb_mut"]
      nbSNV = finald[which(finald$clone == sm & finald$type == "SNV"), "nb_mut"]
      if(type == "ratioSNV") ratio = (nbDEL + nbINS) / nbSNV
      if(type == "ratioINS") ratio = nbDEL / nbINS
      if(type == "ratioSNV" && sm == as.character(unique(finald$clone))[1]){dd = data.frame(sm = sm, ratio = ratio, type = type)} else {
        dd = rbind(dd, data.frame(sm = sm, ratio = ratio, type = type)) }
    }
  }
  dd$sm <- factor(dd$sm, levels = as.character(unique(dd$sm)))
  assign(paste("pr", cellline, sep="_"), ggdotchart(dd, x = "sm", y ="ratio", color = "type", palette = "jco", size = 3, group = "type",
             add = "segment", add.params = list(color = "lightgray", size = 1.5),
             position = position_dodge(0.3),
             ggtheme = theme_pubclean() )
  )
}
grid.arrange(pr_HAP1, pr_A549, pr_MCF7, nrow=1)


#### 96 MUTATIONAL PROFILE ####

pdf("profiles.pdf", 10 ,5) #8,7 when export
plot_96_profile(mat96_A549, condensed = T)
plot_96_profile(mat96_HAP1, condensed = T)
plot_96_profile(mat96_MCF7[ , which(!grepl("CONTROL", colnames(mat96_MCF7)))], condensed = T)
dev.off()

library(corrplot)
pdf("matrix_types_mut.pdf", 10,4)
par(mfrow=c(1,3))
for(cellline in c("A549", "HAP1", "MCF7")){
  mat = get(paste("mat96_", cellline, sep=""))
  if(cellline == "MCF7"){ mat = mat[ , which(!grepl("CONTROL", colnames(mat)))] }
  for(c in colnames(mat)){
    types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    res = unlist(lapply(types, function(t) sum(mat[which(grepl(t, rownames(mat))), c]) ))
    names(res) = types ; dd = data.frame(res) ; colnames(dd) = c
    if(c == colnames(mat)[1]){dat = dd} else {dat = cbind(dat, dd)}
  }
  assign(paste("mat6_", cellline, sep=""), dat)
  corrplot(as.matrix(log10(dat)), is.cor=F, xlab="sample", ylab="", tl.pos = "lt", tl.col = "black", tl.offset=1)
}
dev.off()

####  PCA  ####
for(cellline in c("A549", "HAP1", "MCF7")){
  load(paste("Rdata/mat96_",cellline,"_lowVAF_k50.Rdata", sep=""))
  matlow = get(paste("mat96_",cellline,"_lowVAF", sep=""))
  colnames(matlow) = paste(colnames(matlow), "lowVAF", sep="_")
  load(paste("Rdata/mat96_",cellline,"_highVAF_k50.Rdata", sep=""))
  mathigh = get(paste("mat96_",cellline,"_highVAF", sep=""))
  colnames(mathigh) = paste(colnames(mathigh), "highVAF", sep="_")
  mat = mathigh #cbind(matlow, mathigh)
  if(cellline == "A549") {matAll = mat} else {matAll = cbind(matAll, mat)}
  snvs.pca <- prcomp(t(mat), center = T)
  assign(paste("pca", cellline, sep=""), autoplot(snvs.pca, data.frame(sample=colnames(mat)), col="sample")) #+ scale_color_brewer(palette = "Dark2")
}
ggarrange(pcaA549, pcaHAP1, pcaMCF7, nrow=1)

matAll2 = matAll[, which(!grepl("MCF7_CONTROL", colnames(matAll)))]
colnames(matAll2) = gsub("_highVAF","",colnames(matAll2))
colnames(matAll2) = gsub("_lowVAF","",colnames(matAll2))
library(stringi)
cell=unlist(lapply(colnames(matAll2), function(x){unlist(stri_split_fixed(str = x, pattern = "_", n=2))[1]}))
type_rad0 = unlist(lapply(colnames(matAll2), function(x){unlist(stri_split_fixed(str = x, pattern = "_", n=2))[2]}))
type_rad = rep("control", length(type_rad0))
type_rad[which(grepl("PR" , type_rad0))] = "proton"
type_rad[which(grepl("HE" , type_rad0))] = "helium"

autoplot(prcomp(t(matAll2), center = T), data.frame(radiation=type_rad, cell_line=cell), col="radiation", shape="cell_line", 
         label=TRUE, label.size = 2.75, size=3) 


#### SNV SIGNATURES ####
setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/calling/strelka2/analysis")
load("Rdata/mat96_A549_k50.Rdata") ; load("Rdata/mat96_HAP1_k50.Rdata") ; load("Rdata/mat96_MCF7_k50.Rdata") #load("Rdata/MCF7_withoutcontrol/mat96_MCF7_k50.Rdata")
mat96_MCF7 = mat96_MCF7[ , which(!grepl("CONTROL", colnames(mat96_MCF7)))]
library("MutationalPatterns")
library("xlsx")

mattmp=cbind(mat96_A549, mat96_HAP1)
nmf_res <- extract_signatures(mattmp, rank = 4, nrun = 1)
p1 <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
p2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative")
grid.arrange(p1, p2, nrow=2)

### WITH MUTAGENS SAMPLES ###
# adding the Nik-Zainal samples
datNZ = read.table("data/denovo_subclone_subs_final_Kucab2019.txt", h=T, stringsAsFactors = F, sep="\t")
datNZ_suppT2 = read.xlsx("data/suppTable2_Kucab2019.xlsx", sheetIndex = 1)
#datNZ_rad = datNZ[which(datNZ$Sample.Name %in% c("MSM0.71", "MSM0.56")), ] # keep only the radiation samples
datNZ_rad = datNZ[c( which(datNZ$Sample.Name %in% as.character(unique(datNZ_suppT2[which(datNZ_suppT2$Has.Signature == "YES"), "Sample.Name"]))),
                     which(datNZ$Sample.Name == "MSM0.71") ) , ] # we keep the cosmic sample

# create the GRanges used for context
df_tmp = datNZ_rad[,c("Chrom", "Pos")] ; df_tmp$end = df_tmp$Pos ; colnames(df_tmp) = c("chr","start","end") ; df_tmp$chr=paste("chr",df_tmp$chr,sep="") 
GR_NZ_rad = makeGRangesFromDataFrame(df_tmp) ; seqlevels(GR_NZ_rad) = human_chrs[1:24] # create the Grange object with good chrs (chrX was chr23 e.g.)

#compute the mut tables
ref_genomeNZ = available.genomes()[which(grepl("hg19",available.genomes()) & ! grepl("masked",available.genomes()))]
nut3_contextNZ = type_context(GR_NZ_rad, ref_genomeNZ) # get nut3 context for mutations PASS in the sample, SNVs or indels
all_mut96NZ = paste(substr(nut3_contextNZ$context,1,1), "[", paste(datNZ_rad$Ref, datNZ_rad$Alt, sep=">"), "]", substr(nut3_contextNZ$context,3,3), sep = "") 
for(s in unique(datNZ_rad$Sample.Name)){
  dt = table( all_mut96NZ[which(datNZ_rad$Sample.Name == s)] ) # keep mutations of the current sample
  t = rep(0, length(setdiff(snvs96, names(dt)))) ; names(t) = setdiff(snvs96, names(dt)) ; dt = c(dt, t) ; dt = dt[snvs96] # merge 0 if no mut + observed
  dd = as.numeric(dt)
  if (match(s, unique(datNZ_rad$Sample.Name)) == 1) { ddf = data.frame(dd) ; rownames(ddf) = names(dt) } else { ddf = cbind( ddf, as.numeric(dd)) }
}
mat96_NZ = ddf ; colnames(mat96_NZ) = unique(datNZ_rad$Sample.Name) #colnames(mat96_NZ) = c("Solar", "Gamma")

mat96_NZgroups = mat96_NZ
colnames(mat96_NZgroups) = unlist(lapply(colnames(mat96_NZgroups), function(x) paste(datNZ_suppT2[which(as.character(datNZ_suppT2$Sample.Name) == x), "Group"], 
                                                                                     match(x,datNZ_suppT2$Sample.Name), sep="."))) # add number for non duplications
# write the ID concordance to match samples with indels data
write.table(data.frame(Sample.name = colnames(mat96_NZ), newSM = colnames(mat96_NZgroups)), file="data/Kucab_snvs_IDs_groups.txt", sep="\t", row.names=F, col.names=T, quote=F)
plot_cosine_heatmap(cos_sim_matrix(mat96_HAP1, mat96_NZgroups), col_order = colnames(mat96_NZgroups)[order(colnames(mat96_NZgroups))])
plot_cosine_heatmap(cos_sim_matrix(mat96_A549, mat96_NZgroups), col_order = colnames(mat96_NZgroups)[order(colnames(mat96_NZgroups))])

### WITH DNA-REPAIR SAMPLES ###
datNZ2 = read.table("data/denovo_subs_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t")

df_tmp = datNZ2[,c("Chrom", "Pos")] ; df_tmp$end = df_tmp$Pos ; colnames(df_tmp) = c("chr","start","end") ; df_tmp$chr=paste("chr",df_tmp$chr,sep="") 
df_tmp = df_tmp[which(df_tmp$chr %in% human_chrs), ]
GR_NZ_rad = makeGRangesFromDataFrame(df_tmp) ; seqlevels(GR_NZ_rad) = human_chrs[1:24] # create the Grange object with good chrs (chrX was chr23 e.g.)

#compute the mut tables
ref_genomeNZ = available.genomes()[which(grepl("hg19",available.genomes()) & ! grepl("masked",available.genomes()))]
nut3_contextNZ = type_context(GR_NZ_rad, ref_genomeNZ) # get nut3 context for mutations PASS in the sample, SNVs or indels
all_mut96NZ = paste(substr(nut3_contextNZ$context,1,1), "[", paste(datNZ2$Ref, datNZ2$Alt, sep=">"), "]", substr(nut3_contextNZ$context,3,3), sep = "") 
for(s in unique(datNZ2$Sample)){
  dt = table( all_mut96NZ[which(datNZ2$Sample == s)] ) # keep mutations of the current sample
  t = rep(0, length(setdiff(snvs96, names(dt)))) ; names(t) = setdiff(snvs96, names(dt)) ; dt = c(dt, t) ; dt = dt[snvs96] # merge 0 if no mut + observed
  dd = as.numeric(dt)
  if (match(s, unique(datNZ2$Sample)) == 1) { ddf = data.frame(dd) ; rownames(ddf) = names(dt) } else { ddf = cbind( ddf, as.numeric(dd)) }
}
mat96_NZ2 = ddf ; colnames(mat96_NZ2) = unique(datNZ2$Sample) #colnames(mat96_NZ) = c("Solar", "Gamma")

mat96_NZ2groups = mat96_NZ2
colnames(mat96_NZ2groups) = unlist(lapply(colnames(mat96_NZ2groups), function(x) paste(datNZ2[which(as.character(datNZ2$Sample) == x)[1], "Ko_gene"], #we have all the mut it is not the supp table so take first which
                                                                                     match(x,datNZ2$Sample), sep="."))) # add number for non duplications
# write the ID concordance to match samples with indels data
write.table(data.frame(Sample = colnames(mat96_NZ2), newSM = colnames(mat96_NZ2groups)), file="data/Zou_snvs_IDs_groups.txt", sep="\t", row.names=F, col.names=T, quote=F)
plot_cosine_heatmap(cos_sim_matrix(mat96_HAP1, mat96_NZ2groups), col_order = colnames(mat96_NZ2groups)[order(colnames(mat96_NZ2groups))])
plot_cosine_heatmap(cos_sim_matrix(mat96_A549, mat96_NZ2groups), col_order = colnames(mat96_NZ2groups)[order(colnames(mat96_NZ2groups))]) 

## everything merged 
mattmp=cbind(mat96_A549, mat96_HAP1, mat96_MCF7, mat96_NZ2groups, mat96_NZgroups)

library("NMF")
mattmp = mattmp + 0.0001
estimate <- nmf(mattmp, rank = 2:20, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")
plot(estimate)
nmf_res <- extract_signatures(mattmp, rank = 14, nrun = 5) 
p1 <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
p2 <- plot_contribution(nmf_res$contribution[,1:20], nmf_res$signature, mode = "absolute") + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
p3 <- ggplot(melt(nmf_res$contribution[,1:20]), aes(x = Var1, y = Var2, fill = value)) +  
  geom_tile(colour = "black") + scale_fill_gradient(low = "white", high = "red") +
  scale_x_continuous(breaks=1:nrow(nmf_res$contribution)) + labs(y = "Sample", x = "Signature") + guides(fill=guide_legend(title="Contribution"))

grid.arrange(p1, p2, nrow=2)

# cosine similarities between our signatures and cosmic signatures
library(SigsPack)
data(sigProfiler20190522)
csm = cos_sim_matrix(sigProfiler20190522, nmf_res$signatures) ; colnames(csm) = paste("Sig", 1:ncol(csm), sep="") # cosmic similarities
plot_cosine_heatmap(csm)
apply(csm, 2, function(c)  c[which(c>0.9)])
# cosine similarities between our signatures and NZ samples
nzsm = cos_sim_matrix(cbind(mat96_NZ2groups, mat96_NZgroups), nmf_res$signatures) ; colnames(nzsm) = paste("Sig", 1:ncol(nzsm), sep="") # NZ similarities
plot_cosine_heatmap(nzsm)
apply(nzsm, 2, function(c)  c[which(c>0.9)])
# using both cosmic and NZ signatures
csnzsm = cos_sim_matrix(cbind(sigProfiler20190522, mat96_NZ2groups, mat96_NZgroups), nmf_res$signatures) ; colnames(csnzsm) = paste("Sig", 1:ncol(csnzsm), sep="") 
plot_cosine_heatmap(csnzsm)
apply(csnzsm, 2, function(c) c[which(c>0.9)])

# write the signatures/activities/samples in order to use later the decompose from sigprofiler, that needs a tab in input
letters = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M')
nmf_sigs = nmf_res$signatures
colnames(nmf_sigs) = paste("SBS96", letters, sep="") 
nmf_sigs = cbind(data.frame(MutationType = rownames(nmf_sigs)), nmf_sigs)
write.table(nmf_sigs, "data/nmf_13_signatures.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)

nmf_conts <- t ( round(nmf_res$contribution * colSums(nmf_res$signatures)) )
colnames(nmf_conts) = paste("SBS96", letters, sep="") 
nmf_conts = cbind(data.frame(Samples = rownames(nmf_conts)), nmf_conts)
write.table(nmf_conts, "data/nmf_246_13_activities.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)

nmf_samples = cbind(mat96_A549, mat96_HAP1, mat96_MCF7, mat96_NZ2groups, mat96_NZgroups)
nmf_samples = cbind(data.frame(MutationType = rownames(nmf_samples)), nmf_samples)
write.table(nmf_samples, "data/nmf_247_samples.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)


#### INDELS NMF WITHOUT POOLED ####
indels83 = read.table("/g/strcombio/fsupek_cancer1/indels_signatures_project/83_indel_categories.csv", sep=",", h=T, stringsAsFactors = F)[,3]
load("Rdata/mat83indels_A549_k50.Rdata") ; mat83indels_A549=mat83indels_A549[indels83,]
load("Rdata/mat83indels_MCF7_k50.Rdata") ; mat83indels_MCF7=mat83indels_MCF7[indels83,]
load("Rdata/mat83indels_HAP1_k50.Rdata") ; mat83indels_HAP1=mat83indels_HAP1[indels83,]

mat83indels_MCF7 = mat83indels_MCF7[, which(!grepl("CONTROL", colnames(mat83indels_MCF7)))]
plotExchangeSpectra_indel(mat83indels_A549)
plotExchangeSpectra_indel(mat83indels_HAP1)
plotExchangeSpectra_indel(mat83indels_MCF7)

matindels = cbind(mat83indels_HAP1, mat83indels_A549, mat83indels_MCF7)
estimate <- nmf(matindels, rank = 2:10, method = "brunet", 
                nrun = 10, seed = 123456, .opt = "v-p")
plot(estimate)
nmf_res <- extract_signatures(matindels, rank = 2, nrun = 3) # also looking at cos sim matrix. more than 2 sign == high cosines == high similarities
plotExchangeSpectra_indel(as.data.frame(nmf_res$signatures))
plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(melt(nmf_res$contribution), aes(x = Var1, y = Var2, fill = value)) +  
  geom_tile(colour = "black") + scale_fill_gradient(low = "white", high = "blue") +
  scale_x_continuous(breaks=1:nrow(nmf_res$contribution)) + labs(y = "Sample", x = "Signature") + guides(fill=guide_legend(title="Contribution"))


#### INDELS WITH DNA REPAIR AND MUTAGENS SAMPLES ####
datNZ_indels_Zou = read.table("data/denovo_indels_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t")
datNZ_indels_Kucab = read.table("data/denovo_subclone_indels_final_Kucab2019.txt", h=T, stringsAsFactors = F, sep="\t")
datNZ_suppT2 = read.xlsx("data/suppTable2_Kucab2019.xlsx", sheetIndex = 1)
datNZ_indels_Kucab$newSM = unlist(lapply(datNZ_indels_Kucab$Sample.Name, function(x) paste(datNZ_suppT2[which(as.character(datNZ_suppT2$Sample.Name) == x), "Group"], 
                                                                                     match(x,datNZ_suppT2$Sample.Name), sep="."))) # add number for non duplications
datNZ_indels_Zou$newSM = unlist(lapply(datNZ_indels_Zou$Sample, function(x) paste(datNZ_indels_Zou[which(as.character(datNZ_indels_Zou$Sample) == x)[1], "Ko_gene"], #we have all the mut it is not the supp table so take first which
                                                                                       match(x,datNZ_indels_Zou$Sample), sep="."))) # add number for non duplications
for(study in c("Zou", "Kucab")){
  datNZ2_indels = get(paste("datNZ_indels_", study, sep=""))
  datNZ2_indels$Chrom = paste("chr", datNZ2_indels$Chrom, sep="")
  datNZ2_indels = datNZ2_indels[which(datNZ2_indels$Chrom %in% human_chrs), ]
  datNZ2_indels$newSM = gsub("\\.", "_", datNZ2_indels$newSM)
  datNZ2_indels$newSM = gsub("/", "_", datNZ2_indels$newSM)
  datNZ2_indels$newSM = gsub(" ", "_", datNZ2_indels$newSM)
  for(NZ2sm in unique(datNZ2_indels$newSM)){
    df_tmp = datNZ2_indels[which(datNZ2_indels$newSM == NZ2sm), ]
    vcf = data.frame(CHROM = df_tmp$Chrom, POS = df_tmp$Pos, REF = df_tmp$Ref, ALT = df_tmp$Alt, ID="rsX", QUAL=".", FILTER="PASS")
    vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
    colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
    write.table(c("##fileformat=VCFv4.2"), paste0("data/indels/",study,"/",NZ2sm,".vcf"), row.names = F,col.names = F,quote = F)
    fwrite(vcf, paste0("data/indels/",study,"/",NZ2sm,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")
  }
}


#### RATIO DEL/INS (BEHJATI et al.) ####
setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/calling/strelka2/analysis/")
indels_vcf = list.files("data/indels/ours_withoutMCF7control/", pattern = "vcf")
for(cellline in c("A549", "HAP1", "MCF7")){
  cl_vcf = indels_vcf[which(grepl(cellline, indels_vcf))]
  cl_sm = gsub(paste(cellline, "_", sep=""), "", gsub(".vcf", "", cl_vcf))
  for(sm in cl_sm){
    curdat = read.table(paste("data/indels/ours_withoutMCF7control/", cl_vcf[which(grepl(sm , cl_vcf))], sep=""), sep="\t", h=F, stringsAsFactors = T)
    curdat$type = "DEL"
    curdat[which(nchar( as.character(curdat[,5]) ) > 1), "type"] = "INS"
    curratio = table(curdat$type)["DEL"] / table(curdat$type)["INS"]
    if(grepl("HE", sm)) { cond = "HE" } else { if(grepl("PR", sm)) { cond = "PR"} else { cond = "CT"} }
    if(cellline == "A549" & sm == cl_sm[1]) { alldat = data.frame(clone=sm, cellline=cellline, ratio=curratio, condition=cond) } else {
      alldat = rbind(alldat, data.frame(clone=sm, cellline=cellline, ratio=curratio, condition=cond))
    }
  }
}
for(cellline in c("A549", "HAP1", "MCF7")){
  assign(paste("pp", cellline, sep="_"),
         ggboxplot(alldat[which(alldat$cellline==cellline),], x = "condition", y = "ratio", add = "jitter", palette = "jco") +
           theme(legend.position="none") + ggtitle(cellline) )
}
ggarrange(pp_A549, pp_HAP1, pp_MCF7, nrow=1)
