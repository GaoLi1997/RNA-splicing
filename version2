# Author: Gaoli 
# data: 2021/4
rm(list = ls())
setwd('/home/gaoli/data_sys/RNA_splicing_data/xiugai_collapse')
library(tidyverse)
library(ggplot2)
library(stringr)
library(pheatmap)
library(parallel)
cl <- makeCluster(mc <- getOption("cl.cores", 4))
#library(export)
# get file
gtf_file <- read.table("flair.collapse.isoforms.gtf", stringsAsFactors = F, 
                       sep = "\t", quote = '"')
count_file <- read.table("counts_matrix.tsv", 
                         stringsAsFactors = F, sep = "\t",quote = '"',header = T)

##step1 split file and arrange 
#1.1 split gtf file which is created by flair
gtf_file_split <- separate(data = gtf_file, col = V9, sep = ";", into = c('gene_id','transcript_id','exon_num'), remove = T)
gtf_file_split$gene_id <- substring(gtf_file_split$gene_id, 9,) #delete 'gene_id' 
gtf_file_split$transcript_id <- substring(gtf_file_split$transcript_id,16,) #delete 'transcript_id'
gtf_file_split = mutate(gtf_file_split,exon_element = paste(V1,V4,V5,sep = '_'))


#1.2 split counts_matrix file
count_file_split <- separate(data = count_file, col = ids, sep ='_', 
                             into = c('transcript_id','gene_id'),remove = T) 


gtf_file_split = gtf_file_split[which(gtf_file_split$V3 == 'exon'),]

find_exon_element = function(gtf_file_split, count_file_split){
  split_by_tf <- split(1:dim(gtf_file_split)[1], gtf_file_split$transcript_id)
  res <- parLapply(cl, split_by_tf, function(x, gtf_file_split, count_file_split){
    x <- unlist(x)
    trans_data_frame = gtf_file_split[x,]
    isoform_id = trans_data_frame[1,10]
    exon_list = paste(trans_data_frame$exon_element, collapse = '-')
    len = length(trans_data_frame$exon_element)
    sing_count_file = data.frame(transcript_id = isoform_id, exonlist = exon_list, exon_num = len)
    return(sing_count_file)
  },gtf_file_split = gtf_file_split, count_file_split = count_file_split)
  res <- do.call(rbind, res)
  return(res)
}

iso_exon = find_exon_element(gtf_file_split, count_file_split)
write.table(iso_exon, file = '/home/gaoli/data_sys/RNA_splicing_data/iso_exon.txt', sep ="\t")

count_file_split = merge(x = count_file_split, y = iso_exon, all = TRUE)

##step2 get the exon composition for every transcript isoform
 
write.table(count_file_split, file = '/home/gaoli/data_sys/RNA_splicing_data/count_file_exon.txt', sep ="\t")


#density plot to see distubution
ggplot(data = count_file_split) +
  geom_density(mapping = aes(x=log2(UCSC_conditionA_batch1+1)),alpha = 0.5,color = "black", fill = "gray",adjust = 3)+
  xlab(expression(log[2](UCSC_count))) + ylab('Density') +
  ggtitle('UCSC_data_isoform_count') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 1.7, linetype="dotted")
graph2ppt(file="UCSC.pptx", width=10, height=8)

ggplot(data = count_file_split) +
  geom_density(mapping = aes(x=log2(NA12878_conditionB_batch1+1)),alpha = 0.5,color = "black", fill = "gray",adjust = 3)+
  xlab(expression(log[2](count))) + ylab('Density') +
  ggtitle('NA12878_data_isoform_count') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 2.2, linetype="dotted")
graph2ppt(file="NA12878.pptx", width=10, height=8)

#step3  delete the genes which expressed less than 20
#PS: consider gene counts rather than isoform counts

count_293T.PLVX1 = count_file_split[, c(1:3,12:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.PLVX1_conditionA_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.PLVX1_conditionA_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 

count_293T.PLVX2 = count_file_split[, c(1:2,4,12:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.PLVX2_conditionA_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.PLVX2_conditionA_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 


count_293T.RBM18 = count_file_split[,c(1:2,5,12:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.RBM18_conditionB_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.RBM18_conditionB_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 


count_293T.RBM18_2 = count_file_split[,c(1:2,6,12:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.RBM18_2_conditionB_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.RBM18_2_conditionB_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 


count_293T_RBM23 = count_file_split[,c(1:2,7,12:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.RBM23_conditionC_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.RBM23_conditionC_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 


count_293T_RBM23_2 = count_file_split[,c(1:2,8,12:13)] %>%
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(X293T.RBM23_2_conditionC_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,X293T.RBM23_2_conditionC_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 

count_CCRF = count_file_split[,c(1:2,9,12:13)] %>%
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(CCRF_conditionD_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,CCRF_conditionD_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 

count_UCSC = count_file_split[,c(1:2,10,12:13)] %>%
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(UCSC_conditionE_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,UCSC_conditionE_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() 

count_NA12878 = count_file_split[,c(1,2,11:13)] %>% 
  group_by(., gene_id) %>%  mutate(.,gene_expression = sum(NA12878_conditionF_batch1)) %>% 
  arrange(.,gene_id) %>% ungroup() %>% filter(.,NA12878_conditionF_batch1 > 0) %>% 
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame()

#delete the isforms which have one exon_num
#consider deleting first and last exon

#step4  delete genes which have only one isform
all_file = list(count_293T.PLVX1,count_293T.PLVX2, count_293T.RBM18, count_293T.RBM18_2,
                count_293T_RBM23, count_293T_RBM23_2, count_CCRF, count_UCSC, count_NA12878)

result = lapply(all_file, function(x) { 
  data = as.data.frame(x)
  gene_transform_num  = count(data, gene_id)
  gene_id_delete <- gene_transform_num[which(gene_transform_num$n == 1),1]
  for (gene in gene_id_delete){
    data = data[-which(data$gene_id == gene),]
  }
  return(data)
})


#step5 choose 10bp as the first and last exons

modify_two_exon = lapply(result, function(x) {
  data1 = as.data.frame(x)
  gene_list = unique(data1$gene_id)
  for (b in gene_list){
    data_frame <- data1[which(data1$gene_id == b),] 
    exon_library = c()
    i = 1 
    while (i <= dim(data_frame)[1]) {
      exon_composition = data_frame$exonlist[i] %>% as.character() %>% 
        strsplit(.,'-') %>% unlist(.) %>% sort(.) 
      first_exon = exon_composition[1] %>% strsplit(., '_') %>% unlist()
      last_exon = tail(exon_composition,1) %>% strsplit(., '_') %>% unlist()
      #choose 10bp as the first and last exons
      first_exon[2] = as.numeric(first_exon[3]) - 10
      last_exon[3] = as.numeric(last_exon[2]) + 10
      
      new_first_exon = paste(first_exon, collapse = '_')
      new_last_exon = paste(last_exon, collapse = '_')
      
      exon_composition[1] = new_first_exon
      exon_composition[length(exon_composition)] = new_last_exon
      
      position = which(data1$transcript_id == data_frame$transcript_id[i])
      data1$exon_list[position] = paste(exon_composition, collapse = '-')
      
      exon_library = c(exon_library,exon_composition) %>% unique() %>% sort()
      i = i+1}
    if(length(exon_library) <=2){data1[-which(data1$gene_id == b),]}
  }
  return(data1)})



#step6 

f1 = function(new_data_frame, exon1_split, exon2_split,exon1,exon2){
  con1 = con2 = con3 = con4 = 0
  for (z in 1:nrow(new_data_frame)){
    isoform.exon = new_data_frame$exon_list[z] %>% str_split(.,'-') %>% unlist()
    isoform_scale = str_split(isoform.exon,'_') %>% unlist() %>%     
      matrix(.,length(isoform.exon), byrow = T, dimnames = list(NULL, c("chr", "start",'end'))) %>% 
      as.data.frame()
    isoform_scale$start <- as.numeric(as.character(isoform_scale$start))
    isoform_scale$end <- as.numeric(as.character(isoform_scale$end))
    
    if(min(as.numeric(exon1_split[2]), as.numeric(exon2_split[2])) >= min(isoform_scale$start) & 
       max(as.numeric(exon1_split[3]),as.numeric(exon2_split[3])) <= max(isoform_scale$end)){   
      #con1 = exon1+exon2; con2 = exon1-exon2; con3 = -exon1 + exon2; con4 = -exon1-exon2
      if (str_detect(new_data_frame$exon_list[z],exon1)==TRUE & str_detect(new_data_frame$exon_list[z],exon2)==TRUE){
        con1 = con1+new_data_frame[,3][z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==TRUE & str_detect(new_data_frame$exon_list[z],exon2)==FALSE){
        con2 = con2+new_data_frame[,3][z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==FALSE & str_detect(new_data_frame$exon_list[z],exon2)==TRUE){
        con3 = con3+new_data_frame[,3][z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==FALSE & str_detect(new_data_frame$exon_list[z],exon2)==FALSE){
        con4 = con4+new_data_frame[,3][z]}}
  }
  return(list(con1,con2,con3,con4))
}

clusterExport(cl,"f1")
clusterEvalQ(cl,{library(tidyverse)
  library(stringr)})

# step7 fisher.test & return gene
get_bed_file = parLapply(cl, modify_two_exon, function(x){
  bed_file = c() 
  sig_gene_list = c()
  data2 = as.data.frame(x)
  new_gene_list = unique(data2$gene_id)
  for (n in new_gene_list){
    new_data_frame <- data2[which(data2$gene_id == n),] 
    exon_library = new_data_frame$exon_list %>% str_split(.,'-') %>% unlist() %>% unique() %>% sort()
    exon_matrix = matrix(rep(0,times = length(exon_library)^2),length(exon_library),
                         dimnames = list(a=exon_library,b=exon_library)) %>% as.data.frame()
    two_exon <- combn(exon_library,2)
    for (number in 1:dim(two_exon)[2]){
      exon1 = two_exon[1,number]
      exon2 = two_exon[2,number]  #find partially overlap exons
      exon1_split = str_split(exon1,'_') %>% unlist()
      exon2_split = str_split(exon2,'_') %>% unlist()
      
      scale1 = c(exon1_split[2]:exon1_split[3]) 
      scale2 = c(exon2_split[2]:exon2_split[3])
      if (length(intersect(x = scale1, y = scale2)) == 0){
        constrct_matrix = f1(new_data_frame, exon1_split, exon2_split,exon1,exon2) %>% unlist() %>% matrix(.,2)
        pvalue_less <- fisher.test(constrct_matrix, alternative = 'less')[[1]]  
        pvalue_greater <- fisher.test(constrct_matrix, alternative = 'greater')[[1]]
        pvalue = ifelse(pvalue_less <= pvalue_greater, 
                        p.adjust(pvalue_less ,method = "bonferroni")%>%log(.), #exclusive 
                        p.adjust(pvalue_greater, method = "bonferroni")%>% -log(.)) #inclusive 
        if(abs(pvalue) > 5){bed_file = c(bed_file,exon1,exon2)} 
      }else{pvalue <- NA}
      exon_matrix[which(row.names(exon_matrix) == exon1), which(names(exon_matrix)==exon2)] = pvalue #将p值填写在之前构建matrix中
      exon_matrix[which(row.names(exon_matrix) == exon2), which(names(exon_matrix)==exon1)] = pvalue
    }
    if (abs(min(unlist(c(exon_matrix)),na.rm=T)) > 5){
      sig_gene_list = c(n, sig_gene_list)} 
  }
  return(list(bed_file,sig_gene_list))
})
