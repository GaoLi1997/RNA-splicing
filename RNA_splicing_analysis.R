# Author：gaoli 
# data： 2021/4
rm(list = ls())
setwd('E:/project')
library(tidyverse)
library(ggplot2)
library(stringr)
library(pheatmap)
library(export)
# get file
gtf_file <- read.table("flair.output.isoforms.gtf", stringsAsFactors = F, 
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


##step2 get the exon composition for every transcript isoform
gtf_file_split = gtf_file_split[which(gtf_file_split$V3 == 'exon'),]
for (i in unique(gtf_file_split$transcript_id)){ 
  trans_data_frame = gtf_file_split[which(gtf_file_split$transcript_id == i),] #得到单个转录本的exon信息
  exon_list = paste(trans_data_frame$exon_element, collapse = '-')
  len = length(trans_data_frame$exon_element)
  count_file_split[which(count_file_split$transcript_id == i),5] = exon_list 
  count_file_split[which(count_file_split$transcript_id == i),6] = len}  #得到单个转录本的的exon数量

colnames(count_file_split)[5] = 'exon_list'
colnames(count_file_split)[6] = 'exon_num' 
write.table(count_file_split, file = 'count_file_exon', sep ="\t")

count_file_split = read.csv('count_file_exon', sep = '\t', header = T)


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

count_file_UCSC = group_by(count_file_split, gene_id) %>%  #按照基因分组
  mutate(.,gene_expression = sum(UCSC_conditionA_batch1)) %>%  #计算基因表达量
  arrange(.,gene_id) %>% ungroup() %>% filter(.,UCSC_conditionA_batch1 > 0) %>%  #删除表达量为0的isoform
  filter(.,gene_expression > 20) %>% filter(.,exon_num > 1) %>% as.data.frame() %>% .[, -4] #删除表达量低于20的基因，删除只有一个exon的isoform，只提取UCSC的data
  #delete the isforms which have one exon_num
  #consider deleting first and last exon

#step4  delete genes which have only one isform
gene_transform_num <- count(count_file_UCSC, gene_id)
gene_id_delete <- gene_transform_num[which(gene_transform_num$n == 1),1]
for (gene in gene_id_delete){
  count_file_UCSC = count_file_UCSC[-which(count_file_UCSC$gene_id == gene),]
}

#step5 choose 10bp as the first and last exons
count_file_UCSC <-  arrange(count_file_UCSC, gene_id) ##sort dataframe by geneid
gene_list = unique(count_file_UCSC$gene_id)
for (b in gene_list){
  data_frame <- count_file_UCSC[which(count_file_UCSC$gene_id == b),] 
  exon_library = c()
  i = 1 #遍历每行
  while (i <= dim(data_frame)[1]) {
    exon_composition = data_frame$exon_list[i] %>% strsplit(.,'-') %>% 
      unlist(.) %>% sort(.) 
    first_exon = exon_composition[1] %>% strsplit(., '_') %>% unlist()
    last_exon = tail(exon_composition,1) %>% strsplit(., '_') %>% unlist()
    #choose 10bp as the first and last exons
    first_exon[2] = as.numeric(first_exon[3]) - 10
    last_exon[3] = as.numeric(last_exon[2]) + 10
    
    new_first_exon = paste(first_exon, collapse = '_')
    new_last_exon = paste(last_exon, collapse = '_')
    
    exon_composition[1] = new_first_exon
    exon_composition[length(exon_composition)] = new_last_exon
    
    position = which(count_file_UCSC$transcript_id == data_frame$transcript_id[i])
    count_file_UCSC$exon_list[position] = paste(exon_composition, collapse = '-')
    
    exon_library = c(exon_library,exon_composition) %>% unique() %>% sort()
    i = i+1}
  if(length(exon_library) <=2){count_file_UCSC[-which(count_file_UCSC$gene_id == b),]}
  }

write.table(count_file_UCSC,'count_file_USCS.txt', sep = '\t')

count_file_UCSC = read.csv('count_file_USCS.txt', sep='\t', header = T)

#step6 构建一个计算四种condition的function

f1 = function(new_data_frame){      #此处的dataframe是每个基因的所有isoform信息构成的新数据框
  for (z in 1:nrow(new_data_frame)){
    isoform.exon = new_data_frame$exon_list[z] %>% str_split(.,'-') %>% unlist()
    isoform_scale = str_split(isoform.exon,'_') %>% unlist() %>%     #遍历每行，计算每个isoform的长度范围
      matrix(.,length(isoform.exon), byrow = T, dimnames = list(NULL, c("chr", "start",'end'))) %>% 
      as.data.frame()
    if(min(exon1_split[2], exon2_split[2]) >= min(isoform_scale$start) & 
       max(exon1_split[3],exon1_split[3]) <= max(isoform_scale$end)){   #判断是否此isoform范围囊括的exon1记忆exon2的范围
      #con1 = exon1+exon2; con2 = exon1-exon2; con3 = -exon1 + exon2; con4 = -exon1-exon2
      if (str_detect(new_data_frame$exon_list[z],exon1)==TRUE & str_detect(new_data_frame$exon_list[z],exon2)==TRUE){
        con1 = con1+new_data_frame$UCSC_conditionA_batch1[z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==TRUE & str_detect(new_data_frame$exon_list[z],exon2)==FALSE){
        con2 = con2+new_data_frame$UCSC_conditionA_batch1[z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==FALSE & str_detect(new_data_frame$exon_list[z],exon2)==TRUE){
        con3 = con3+new_data_frame$UCSC_conditionA_batch1[z]
      }else if (str_detect(new_data_frame$exon_list[z],exon1)==FALSE & str_detect(new_data_frame$exon_list[z],exon2)==FALSE){
        con4 = con4+new_data_frame$UCSC_conditionA_batch1[z]}}
  }
  return(list(con1,con2,con3,con4))
}


# step7 fisher.test & heatmap
new_gene_list = unique(count_file_UCSC$gene_id)
bed_file = c() 
sig_gene_list = c()  

for (x in new_gene_list){
  new_data_frame <- count_file_UCSC[which(count_file_UCSC$gene_id == x),] 
  
  exon_library = new_data_frame$exon_list %>% str_split(.,'-') %>% unlist() %>% unique() %>% sort()
  exon_matrix = matrix(rep(0,times = length(exon_library)^2),length(exon_library),
                       dimnames = list(a=exon_library,b=exon_library)) %>% as.data.frame()  # 构建全为0的矩阵

  two_exon <- combn(exon_library,2)
  for (number in 1:dim(two_exon)[2]){
    exon1 = two_exon[1,number]
    exon2 = two_exon[2,number]  #find partially overlap exons
    exon1_split = str_split(exon1,'_') %>% unlist()
    exon2_split = str_split(exon2,'_') %>% unlist()
    scale1 = c(exon1_split[2]:exon1_split[3]) #判断两个外显子是否有重复的部分
    scale2 = c(exon2_split[2]:exon2_split[3])
    if (length(intersect(x = scale1, y = scale2)) == 0){
      con1 = con2 = con3 = con4 = 0
      constrct_matrix = f1(new_data_frame) %>% unlist() %>% matrix(.,2)
      pvalue_less <- fisher.test(constrct_matrix, alternative = 'less')[[1]]  
      pvalue_greater <- fisher.test(constrct_matrix, alternative = 'greater')[[1]]
      pvalue = ifelse(pvalue_less <= pvalue_greater, 
                      p.adjust(pvalue_less ,method = "bonferroni")%>%log(.), #exclusive 
                      p.adjust(pvalue_greater, method = "bonferroni")%>% -log(.)) #inclusive 
      #计算p值，矫正，取log，赋予正负值(exclusive&inclusive)
      if(abs(pvalue) > 5){bed_file = c(bed_file,exon1,exon2)} #  p值显著，保存两个exon信息，后续构建bed文件
    }else{pvalue <- NA}
    exon_matrix[which(row.names(exon_matrix) == exon1), which(names(exon_matrix)==exon2)] = pvalue #将p值填写在之前构建matrix中
    exon_matrix[which(row.names(exon_matrix) == exon2), which(names(exon_matrix)==exon1)] = pvalue
  }
  if (abs(min(unlist(c(exon_matrix)),na.rm=T)) > 5){
    sig_gene_list = c(x, sig_gene_list) #将p值显著的gene保存下来
    png(file = paste(x,'.png',sep = ''), width = 1000, height = 600)
    pheatmap(exon_matrix,main = x,angle_col = 45, 
             color = colorRampPalette(colors = c("blue","white","red"))(100),
             display_numbers = T, cluster_cols = F,cluster_rows = F) 
    dev.off()} 
}

bed_file1 = str_split(bed_file,'_') %>% unlist() %>% 
  matrix(.,ncol = 3,byrow = T, dimnames = list(NULL,c("chrom", "chromStart",'chromEnd'))) %>% 
  as.data.frame()


bed_file1$chromStart = as.numeric(bed_file1$chromStart)
bed_file1$chromEnd = as.numeric(bed_file1$chromEnd)
bed_file1$chrom = as.factor(bed_file1$chrom)
bed_file1 = mutate(bed_file1, vars = bed_file1$V3 - bed_file1$V2) %>% filter(., vars > 10)
write.table(bed_file1, 'UCSC1.bed', sep = '\t')

