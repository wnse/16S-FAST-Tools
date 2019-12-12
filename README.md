# 16S-FAST-Tools

16S-FAST数据分析工具包，适用于具有连接文库（L）和拼接文库（A）的微生物16S高通量测序数据分析。根据L数据中UMI信息对A数据中的序列进行组装拼接，形成16S全长序列，进行物种分析及后续分析。

##L umi  
根据测序DNA文库序列结构，利用cutadapt软件切除umi前后的引物及接头等序列，得到14bp的umi序列。
##L umiID
根据L双端数据的umi序列，进行umi间的匹配，形成umiID,即umi的分组。
##A umi
同样利用cutadapt软件，获取在A数据的R2序列中的umi，并与上述L数据的umi对应，将A序列与umiID关联。
##过滤低质量序列
使用trimmomatic软件，去除A数据的R1序列中的低质量序列。
##序列分组
根据A序列与umiID关系，将A数据中R1序列进行分组。
##组装
对每一组序列进行组装，组装软件为spades。
##去除Contigs引物
使用cutadapt软件去除16S全长序列两端扩增引物及其他序列。
##Contig筛选及聚类
根据序列长度筛选Contig，如果多条来自于同一umiID的Contig符合条件则全部保留。  
使用cd-hit在特定相似度水平对Contigs聚类形成Cluster。    
如果来自同一umiID的Contig相似度低于该水平，去掉该umiID的全部Contig，反之，则只保留序列长度最长的Contig。  
选择同一Cluter中最长的Contig作为该Cluster的代表序列。
##物种分类
使用mothur软件对聚类后的代表序列进行物种分类。
