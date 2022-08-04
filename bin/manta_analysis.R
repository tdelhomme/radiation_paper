##### MANTA STRUCTURAL VARIANTS #####
library(VariantAnnotation)

setwd("/fsupek_cancer1/radiation_project/radiation_project/calling/manta/")

A549 = "manta_output/A549_diploidSV_PASS.vcf.gz"
HAP1 = "manta_output/HAP1_diploidSV_PASS.vcf.gz"
MCF7 = "manta_output/MCF7_diploidSV_PASS.vcf.gz"
A549_inv = "inversions/A549_diploidSV_PASS_convertInvertion.vcf.gz"
HAP1_inv = "inversions/HAP1_diploidSV_PASS_convertInvertion.vcf.gz"
MCF7_inv = "inversions/MCF7_diploidSV_PASS_convertInvertion.vcf.gz"

nb_positive_samples = function(id, matrix_geno, matrix_filter){
  genos = matrix_geno[id, ]
  filts = matrix_filter[id, ]
  genos_kept = genos #[which(filts == "PASS")] #add this will keep a lot of germline with VAF~1
  sum ( genos_kept!="0/0" & genos_kept!="0|0" ) # better to remove if found in at least 2 different conditions?
}

for(cellline in c("A549", "HAP1", "MCF7")){
  print(paste(date(), " INFO: working on cell line: ", cellline, sep=""))
  vcf = open(VcfFile(get(cellline),  yieldSize=100000))
  vcf_chunk0 = readVcf(vcf, "hg38")
  vcf_chunk = vcf_chunk0[which(!grepl("BND", rownames(vcf_chunk0))) , ] # remove BND because we run the convertInversions.py
  vcf_chunk = vcf_chunk[which(is.na(info(vcf_chunk)$EVENT) | !duplicated(info(vcf_chunk)$EVENT)) , ] # same event for non-inv == one SV
  if(cellline == "MCF7") { vcf_chunk = vcf_chunk[ , which(!grepl("CONTROL", colnames(vcf_chunk)))] }
  
  vcf_inv = open(VcfFile(get(paste(cellline, "_inv", sep="")),  yieldSize=100000))
  vcf_chunk_inv0 = readVcf(vcf_inv, "hg38")
  vcf_chunk_inv = vcf_chunk_inv0[which(grepl("INV", rownames(vcf_chunk_inv0))) , ] # add INV after removing BND
  events_ids = info(vcf_chunk_inv)$EVENT # compute the events IDs-- unique or NA are single inverted junctions
  vcf_chunk_inv = vcf_chunk_inv[ which( is.na(events_ids) | !events_ids %in% events_ids[duplicated(events_ids)]), ] # remove reciprocal inversions (should be artefacts)

  if(cellline == "MCF7") { vcf_chunk_inv = vcf_chunk_inv[ , which(!grepl("CONTROL", colnames(vcf_chunk_inv)))] }
  
  # merge SV without BND and the inversions
  gg0 = rbind(geno(vcf_chunk, "GT"), geno(vcf_chunk_inv, "GT"))
  ft0 = rbind(geno(vcf_chunk, "FT"), geno(vcf_chunk_inv, "FT"))
  type0 = c(info(vcf_chunk)[,"SVTYPE"], info(vcf_chunk_inv)[,"SVTYPE"])
  len0 = c(info(vcf_chunk)[,"SVLEN"], info(vcf_chunk_inv)[,"SVLEN"])
  chr0=c(as.character(seqnames(rowRanges(vcf_chunk,"seqnames"))), as.character(seqnames(rowRanges(vcf_chunk_inv,"seqnames"))))
  loc0=c(start(ranges(rowRanges(vcf_chunk,"seqnames"))), start(ranges(rowRanges(vcf_chunk_inv,"seqnames"))))
  ref0=rep("A", length(loc0)) #as.character(ref(vcf_chunk))
  alts0=ref0 #rep(NA, length(ref))
  
  ids_keep = which( unlist(lapply(1:nrow(gg0), function(id) nb_positive_samples(id, gg0, ft0))) < 2 )
  gg = gg0[ids_keep, ]
  ft = ft0[ids_keep, ]
  len = len0[ids_keep, ]
  muts_id = paste(chr0[ids_keep], ":", loc0[ids_keep], "_", ref0[ids_keep], "/", alts0[ids_keep], sep="")
  is.na(len) <- lengths(len) == 0 # replace integer(0) by NAs, otherwise the unlist() removes them
  len = abs(unlist(len))
  len[which(is.na(len))] = 1000 # large insertions more than 1000 are length NA in the VCF (CAUTION: BND are also NA)
  
  allsm = colnames(gg0)
  allsm = unlist(lapply(allsm, function(s) gsub("_REP1","",s)))
  if(cellline == "MCF7") {allsm = allsm[which(!grepl("CONTROL", allsm))]}
  for(sm in allsm){
    ggtmp = gg[, grepl(sm, colnames(gg))] ; fttmp = ft[, grepl(sm, colnames(ft))] 
    muts_id_sm = muts_id[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS")]
    mutSV = names(ggtmp)[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS")] # list of SVs where the sample is mutated and PASS

    # consider unique BND
    # bnd_ids = which(grepl("BND", mutSV))
    # bnd_bc = paste(unlist(lapply(mutSV[bnd_ids], function(s) unlist(strsplit(s, ":"))[1])),
    #                unlist(lapply(mutSV[bnd_ids], function(s) unlist(strsplit(s, ":"))[2])),
    #                unlist(lapply(mutSV[bnd_ids], function(s) unlist(strsplit(s, ":"))[3])), sep=":")
    # kept_ids_sv = c(which(!grepl("BND", mutSV)), bnd_ids[which(!duplicated(bnd_bc))] ) # keep not BND or BND unique
    kept_ids_sv = 1:length(mutSV) # because we do not look at SV we take all
    lenSV = len[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS")][kept_ids_sv] # length of SVs where the sample is mutated and PASS and remove multiple BND
    mutSV = mutSV[kept_ids_sv] # we remove multiple BND
    muts_id_sm = muts_id_sm[kept_ids_sv] # we remove multiple BND
    assign(paste("all_SVs",sm,sep="_"), muts_id_sm) ; save(list=paste("all_SVs",sm,sep = "_"), file=paste("Rdata/all_SVs_", sm, ".Rdata", sep = ""))
    
    
    mutSVtype = gsub("Manta", "", unlist(lapply(mutSV, function(x) unlist(strsplit(x, ":"))[1]))) # type of SVs
    d = data.frame(clone = sm, nb_SV = as.numeric(table(mutSVtype)), type = names(table(mutSVtype)))
    dlen = data.frame(clone = sm, SV_length = log10(lenSV), type = mutSVtype)
    if(sm == allsm[1]) {finald = d} else {finald = rbind(finald, d)}
    if(sm == allsm[1]) {finaldlen = dlen} else {finaldlen = rbind(finaldlen, dlen)}
  }
  # consider BND as INV because we took them unique
  finald$type = as.character(finald$type) # do not consider it as factor
  finald[which(finald$type == "BND"), "type"] = "INV"
  # assign data frames
  assign(paste("finald_", cellline, sep=""), finald)
  assign(paste("finaldlen_", cellline, sep=""), finaldlen)
  # assign plot of number of mutation per type to each sample
  assign(paste("p_", cellline, sep=""),
         ggbarplot(finald, x="clone", y="nb_SV", fill="type", label=T, palette = "jco") + ggtitle(cellline) +  rotate_x_text(45)
  )
  # assign plot of length of mutation per type to each sample
  assign(paste("p_length_", cellline, sep=""),
         ggviolin(finaldlen, x="type", y="SV_length", fill="clone", add = "boxplot") + ggtitle(cellline)
  )
}
grid.arrange(p_HAP1, p_A549, p_MCF7, nrow=1) #5.5 15 for pdf export
grid.arrange( ggpar(p_length_HAP1, legend = "right"),  ggpar(p_length_A549, legend = "right"), ggpar(p_length_MCF7, legend = "right"), nrow=3)

# plot ratio DEL/INS
for(cellline in c("A549", "HAP1", "MCF7")){
  dd = get(paste("finald_", cellline, sep=""))
  cl_sm = gsub(paste(cellline, "_", sep=""), "", unique(dd$clone))
  for(sm in cl_sm){
    curdat = dd[which(grepl(sm, dd$clone)), ]
    if( length(curdat[which(curdat$type == "INS"), "nb_SV"]) == 0) { nb_ins = 0 } else { nb_ins = curdat[which(curdat$type == "INS"), "nb_SV"]}
    curratio = curdat[which(curdat$type == "DEL"), "nb_SV"] / ( nb_ins + curdat[which(curdat$type == "DUP"), "nb_SV"])
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



# plot lengths with comparisons
for(clone in c("A549", "HAP1", "MCF7")){
  ddfl = get(paste("finaldlen_", clone, sep=""))
  ddfl$class = rep("control", nrow(ddfl))
  ddfl[which(grepl("_HE", ddfl$clone)), "class"] = "Helium"
  ddfl[which(grepl("_PR", ddfl$clone)), "class"] = "Proton"
  if(clone == "MCF7"){
    mycomps=list(c("Helium", "Proton"))
  } else { mycomps=list(c("control", "Helium"), c("control", "Proton"), c("Helium", "Proton")) }
  for(tt in c("DEL", "DUP", "INS", "INV")){
    assign(paste("p_", tt, sep=""), ggboxplot(ddfl[which(ddfl$type==tt) , ],
                                              x="class", y="SV_length", color = "class", palette="jco", add="point") + # jitter seems to add perturbation also in y axis 
             stat_compare_means(comparisons = mycomps) )
  }
  grid.arrange(p_DEL+ggtitle(paste(clone, " DEL", sep="")), p_INS+ggtitle("INS"), p_DUP+ggtitle("DUP"), p_INV+ggtitle("INV"),nrow=1)
}

# build table with summary statistics
celllines = c("A549", "HAP1")
conditions = c("C", "H", "P")
resdat = data.frame(stat = rep(c("sum_nb_SV", "mean_length"), each=length(celllines)*length(conditions)),
                    condition = rep(conditions, each=length(celllines)),
                    cl = rep(celllines, length(conditions)),
                    value=NA)
for(i in 1:nrow(resdat)){
  curcl = resdat[i , "cl"]
  curcond = resdat[i , "condition"]
  datnb = get(paste("finald_", curcl, sep=""))
  datl = get(paste("finaldlen_", curcl, sep=""))
  if(grepl("nb", resdat[i, "stat"])) { resdat[i, "value"] = sum(datnb[which(grepl(paste(curcl, curcond, sep="_"), datnb$clone)), 2]) } else {
    resdat[i, "value"] = median(datl[which(!datl$type == "BND" & grepl(paste(curcl, curcond, sep="_"), datl$clone)), 2])
  }
}

#### replication timing ####
setwd("/fsupek_cancer1/radiation_project/radiation_project/regional_enrichment/REA_output_merged_SV/TXT")
setwd("/fsupek_cancer1/radiation_project/radiation_project/regional_enrichment/REA_output_merged_SV_H3K36me3/TXT")
allbins = list.files(".", pattern = "allbins")
for(curbin in allbins){
  sm = gsub("_allbins.txt", "", curbin)
  curdat = read.table(curbin, sep="\t", h=T, stringsAsFactors = F)
  res = table(curdat$bin)
  names(res) = paste("bin", 1:length(res), sep="")
  if(curbin == allbins[1]){
    alldat = data.frame(sm = sm, bins = names(res), counts = as.numeric(res))
  } else { alldat = rbind(alldat, data.frame(sm = sm, bins = names(res), counts = as.numeric(res))) }
}
ggarrange( ggline(alldat[which(grepl("A549", alldat$sm)),], x="bins", y="counts", color="sm", size=1, palette = "jco"),
           ggline(alldat[which(grepl("HAP1", alldat$sm)),], x="bins", y="counts", color="sm", size=1, palette = "jco"),
           ggline(alldat[which(grepl("MCF7", alldat$sm)),], x="bins", y="counts", color="sm", size=1, palette = "jco"), nrow=1)


#### clusters ####  load the functions in analysis_clusters.R script in strelka2/analysis folder 
for(cellline in c("A549", "HAP1", "MCF7")){
  print(paste("Starting cell line: ", cellline, sep=""))
  if(exists("alldist")) rm(list="alldist") ; if(exists("allmuts")) rm(list="allmuts") 
  alldat = list(ls(pattern = paste("all_SVs_", cellline, sep="")))[[1]]
  allsm = unlist(lapply(alldat, function(x) unlist(strsplit(x, "all_SVs_"))[2]))
  if(cellline == "MCF7") allsm = allsm[which(!grepl("CONTROL", allsm))]
  for(sm in allsm){
    allmut0 = get(paste("all_SVs_", sm, sep = ""))
    if(exists(paste("dist_", sm, sep=""))) rm(list=paste("dist_", sm, sep=""))
    for(chr in paste("chr", 1:22, sep="")){
      print(chr)
      allmutchr = allmut0[which(grepl(paste(chr, ":", sep=""), allmut0))]
      if(length(allmutchr) > 0){
        res = graph_clusters(allmutchr, max_len = 50000)
        if(!is.na(res)) { if(!exists(paste("dist_", sm, sep=""))) { 
          comps = components(graph_from_data_frame(d = res$edges0, vertices = res$nodes0, directed = FALSE))$csize 
          assign(paste("comps_", sm, sep=""), comps)
          assign(paste("dist_", sm, sep=""), res$alldistances)
          assign(paste("muts_", sm, sep=""), res$muts)
        } else { 
          assign(paste("dist_", sm, sep=""), c( get(paste("dist_", sm, sep="")), res$alldistances))
          assign(paste("muts_", sm, sep=""), c( get(paste("muts_", sm, sep="")), res$muts))
          assign(paste("comps_", sm, sep=""), c( get(paste("comps_", sm, sep="")), comps))
        }
        }
      }
    }
    if(!exists("alldist") & exists(paste("dist_", sm, sep=""))) {
      alldist=data.frame(d=get(paste("dist_", sm, sep="")), sm=sm)
      allmuts=data.frame(d=get(paste("muts_", sm, sep="")), sm=sm) 
      allcomps=data.frame(d=get(paste("comps_", sm, sep="")), sm=sm) } else {
        if(exists(paste("dist_", sm, sep=""))){ # if not exists, no clustered mutations
          alldist = rbind(alldist, data.frame(d=get(paste("dist_", sm, sep="")), sm=sm) )
          allmuts = rbind(allmuts, data.frame(d=get(paste("muts_", sm, sep="")), sm=sm) )
          allcomps = rbind(allcomps, data.frame(d=get(paste("comps_", sm, sep="")), sm=sm) )
        }
      }
  }
  assign(paste("alldist_", cellline, sep=""), alldist)
  assign(paste("allmuts_", cellline, sep=""), allmuts)
  assign(paste("allcomps_", cellline, sep=""), allcomps)
}
for(cellline in c("A549", "HAP1", "MCF7")){
  dd = get(paste("allcomps_", cellline, sep=""))
  cl_sm = gsub(paste(cellline, "_", sep=""), "", unique(dd$sm))
  for(sm in cl_sm){
    curdat = dd[which(grepl(sm, dd$sm)), ]
    if(grepl("HE", sm)) { cond = "HE" } else { if(grepl("PR", sm)) { cond = "PR"} else { cond = "CT"} }
    if(cellline == "A549" & sm == cl_sm[1]) { alldat = data.frame(clone=sm, cellline=cellline, comp_size=curdat[,1], condition=cond) } else {
      alldat = rbind(alldat, data.frame(clone=sm, cellline=cellline, comp_size=curdat[,1], condition=cond))
    }
  }
}
for(cellline in c("A549", "HAP1", "MCF7")){
  curdatp = alldat[which(alldat$cellline==cellline),]
  curdatp$sign = rep("no", nrow(curdatp)) ; curdatp[which(curdatp$comp_size>3) , "sign"] = "yes" 
  if(length(which(curdatp$sign=="yes")) > 1){
    assign(paste("pp", cellline, sep="_"),
           ggboxplot(curdatp[which(curdatp$sign=="no"),], x = "condition", y = "comp_size", add = "jitter") +theme(legend.position="none") + ggtitle(cellline) +
             geom_point(d = data.frame(x=curdatp[which(curdatp$sign == "yes"), "condition"], y=curdatp[which(curdatp$sign == "yes"), "comp_size"]),
                        aes(x,y,colour="blue")) )
    } else {
      assign(paste("pp", cellline, sep="_"),
           ggboxplot(curdatp[which(curdatp$sign=="no"),], x = "condition", y = "comp_size", add = "jitter") +theme(legend.position="none") + ggtitle(cellline) )
    }
}
ggarrange(pp_A549, pp_HAP1, pp_MCF7, nrow=1)


#### clustering with SNVs ####
setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/calling/manta")
for(dist_clust in c(1000, 10000, 50000)){
  for(cellline in c("A549", "HAP1", "MCF7")){
    alldat = list(ls(pattern = paste("all_SVs_", cellline, sep="")))[[1]]
    allsm = unlist(lapply(alldat, function(x) unlist(strsplit(x, paste("all_SVs_",cellline,"_",sep="")))[2]))
    for(sm in allsm){
      all_SVs = get(paste("all_SVs_",cellline,"_",sm, sep=""))
      load(paste("../strelka2/analysis/Rdata/all_SNVs_", cellline, "_", sm, ".Rdata", sep=""))
      all_SNVs = get(paste("all_SNVs_",cellline,"_",sm, sep=""))
      for(SV in all_SVs){
        SV_chr = unlist(strsplit(SV, ":"))[1]
        SV_pos = as.numeric(unlist(strsplit(unlist(strsplit(SV, ":"))[2], "_"))[1])
        snv_chr = all_SNVs[which(grepl(paste(SV_chr, ":", sep=""), all_SNVs))]
        snv_pos = as.numeric(unlist(lapply(snv_chr, function(x) unlist(strsplit(unlist(strsplit(x, ":"))[2], "_"))[1])))
        clust_nb = length(which( abs(snv_pos - SV_pos) <= dist_clust ))
        if(SV == all_SVs[1]) { clust_nb_tmp = clust_nb } else { clust_nb_tmp = clust_nb_tmp + clust_nb }
        if(SV == all_SVs[1]) { SV_is_clust = (clust_nb>0) } else { SV_is_clust = c(SV_is_clust, clust_nb>0) }
      }
      nb_clust_rel = sum(SV_is_clust) / length(all_SVs)
      if(sm == allsm[1] & cellline == "A549") { 
        clustnbf = data.frame(nb_clust = clust_nb_tmp, sm = paste(sm, cellline, sep="_"), dist=dist_clust, nb_clust_rel = nb_clust_rel) } else {
        clustnbf = rbind(clustnbf, data.frame(nb_clust = clust_nb_tmp, sm = paste(sm, cellline, sep="_"), dist=dist_clust, nb_clust_rel = nb_clust_rel) )
      }
    }
  }
  if(dist_clust == 1000) {allclustnbf = clustnbf} else { allclustnbf = rbind(allclustnbf, clustnbf)}
}
p1 = ggplot(data=allclustnbf[which(allclustnbf$dist==1000),], aes(x=sm, y=nb_clust_rel)) + geom_bar(stat="identity", fill="steelblue") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 = ggplot(data=allclustnbf[which(allclustnbf$dist==10000),], aes(x=sm, y=nb_clust_rel)) + geom_bar(stat="identity", fill="steelblue") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 = ggplot(data=allclustnbf[which(allclustnbf$dist==50000),], aes(x=sm, y=nb_clust_rel)) + geom_bar(stat="identity", fill="steelblue") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggarrange(p1+ggtitle("dist=1Kb"), p2+ggtitle("dist=10Kb"), p3+ggtitle("dist=50Kb"), nrow=1)
