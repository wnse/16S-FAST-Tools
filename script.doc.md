### mkdir:

- NAME  
      mkdir - Created on Fri Nov  1 10:04:10 2019

- DESCRIPTION  
      @author: yk

- FUNCTIONS  
      mkdir(dir_name)  
          建立目录

- FILE  
      script/mkdir.py



### submit_cutadapt:

- NAME
      submit_cutadapt - Created on Thu Oct 31 10:24:34 2019

- DESCRIPTION
      @author: yk

- FUNCTIONS
  - submit_cutadapt(in_file, out_file, log_file, adapter, g_a, path, threads=1, info=0)
            去除接头序列
            参数：
                in_file: 输入fastq文件
                out_file: 输出fastq文件
                log_file: 输出日志
                adapter: 接头序列
                g_a: 去除方式（g/a）
                path: cutadapt路径
                threads: 进程数
                info: 是否输出每条序列的接头信息  
            返回：
                无

- FILE
      script/submit_cutadapt.py



### submit_Trimmomatic_SE:

- NAME
      submit_Trimmomatic_SE

- FUNCTIONS
  - submit_Trimmomatic_SE(in_file, out_file, path, threads=1)
            去除低质量序列（单端）
            参数：
                in_file: 输入fastq文件
                out_file: 输出fastq文件
                path: Trimmomatic路径
                threads: 进程数
            返回：

- FILE
     script/submit_Trimmomatic_SE.py



### get_L_UMI:

- NAME
      get_L_UMI

- DESCRIPTION
      Created on Fri Nov  1 14:32:45 2019
      根据去除引物和接头的连接文库序列，统计连接文库的umi对应关系和数量
      序列长度等于14为umi，否则为*
      @author: yk

- FUNCTIONS
  - get_L_UMI(LR1, LR2)
            根据去除引物和接头的连接文库序列，统计连接文库的umi对应关系和数量
            参数：
                LR1: 连接文库去除接头和引物的R1数据（type=dict）
                LR2: 连接文库去除接头和引物的R2数据（type=dict）
            返回：
                df: 统计信息（type=pandas）

- FILE
      script/get_L_UMI.py



### define_UMI_ID:

- NAME
      define_UMI_ID

- DESCRIPTION
      Created on Thu Oct 31 14:48:08 2019
      根据连接文库序列的umi，定义umiID。
      1.排序，确定数量最多的配对关系；
      2.在剩余的配对关系中删除：
          已经确定的umiID的umi2对应的所有其他umi配对中，该umi2为最高数量的umi1，除非该配对关系数量大于第一配对关系数量的特定倍数；
          已经确定的umiID的umi1对应的所有其他umi配对中，该umi1为最高数量的umi2，除非该配对关系数量大于第一配对关系数量的特定倍数；
      3.如果配对关系数量小于指定的最小值，配对结束；
      4.如果已无配对，配对结束；
      @author: yk

- FUNCTIONS
  - define_UMI_ID(UMI_dict, UMI1_dict, UMI2_dict, cutoff=0.5, min_counts=1)
            根据连接文库序列的umi，定义umiID
            参数：
                UMI_dict: 连接文库所有umi配对关系及数量（type=dict）
                UMI1_dict: 连接文库所有umi1对应的umi2及数量
                UMI2_dict: 连接文库所有umi2对应的umi1及数量
                cutoff: 可以作为umiID的配对关系所需要大于第一配对关系数量的最小比例
                min_counts: 配对关系的最小数量
            返回：
                out_list: [umiID, umi1, umi2, counts]

- FILE
      script/define_UMI_ID.py



### reverse_complement:

- NAME
      reverse_complement - Created on Thu Oct 31 13:06:49 2019

- DESCRIPTION
      @author: yk

- FUNCTIONS
  - reverse_complement(seq_list)
            序列取反向互补
            参数：
                seq_list: 序列列表（tpye=list）
            返回：
                rc_list: 序列列表（type=list）

- DATA
      generic_dna = DNAAlphabet()

- FILE
      script/reverse_complement.py



### corresponding_seq2umi:

- NAME
      corresponding_seq2umi

- DESCRIPTION
      Created on Mon Nov  4 10:27:33 2019
      获取拼接文库中seq_id与umi对应关系
      @author: yk

- FUNCTIONS
  - corresponding_seq2umi(adapter_info_file, primer_info_file, umi_dict)
            处理拼接文库数据R2的接头，确定UMI，并序列ID与连接文库UMI进行关联
            参数：
                adapter_info_file: adapter序列文件
                primer_info_file: primer序列文件
                umi_dict: 连接文库确定的umi与umiID关系（type=dict）
            返回：
                aUMI: 序列ID与umi的对应关系（type=dict）
                sta_list: 统计信息（type=list)

- FILE
      script/corresponding_seq2umi.py



### group_umi_seq:

- NAME
      group_umi_seq - Created on Mon Nov  4 16:24:11 2019

- DESCRIPTION
      @author: yk

- FUNCTIONS
  - group_into_umis(seq, aUMI_dict)
  - group_into_umis_unpaired(seq)
  - queu_group_umi_seq(seq, umi2ID_dict, out_dir, aUMI_dict, threads=1, paired=1)
        基于umi将序列分组（多进程）
        参数：
            seq: 输入序列文件（fastq）
            umi2ID_dict: umi与umi ID的对应关系（type=dict）
            out_dir: 分组序列输出目录
            aUMI_dict: umi与序列的对应关系（type=dict)
            threads: 进程数
            paired: 默认参数
        返回：
            umi_counts_sub: 每个umi的序列数（type=dict）
            umi_id_counts_sub: 每个umi ID的序列数（type=dict）
            umi2seq_tmp: 分组过程检测的序列数据（type=list）
  - write_to_file(tmp_list)

- FILE
      script/group_umi_seq.py



### plot_seq2umi:

- NAME
      plot_seq2umi - Created on Tue Nov  5 14:08:05 2019

- DESCRIPTION
      @author: yk

- FUNCTIONS
  - plot_seq2umi(np_seq2UMI, out_file)
            对分组过程检测的序列数据作图
            参数：
                np_seq2UMI: 分组过程检测的序列数据
                out_file: 输出文件（png）
            返回：
                无

-  FILE
      script/plot_seq2umi.py





### assemble_fastq_from_list:

- NAME
      assemble_fastq_from_list

- DESCRIPTION
      Created on Tue Nov  5 14:56:20 2019
      序列组装
      @author: yk

- FUNCTIONS
  - assemble(file_list, out_dir, log_dir, spades_path, threads=1)
            多进程提交序列组装任务
            参数：
                file_list: 用于组装的序列文件列表
                out_dir: 组装结果输出目录
                log_dir: 组装日志输出目录
                spades_path: spades.py程序路径
                threads: 进程数
            返回：
                ass_success: 组装成功(结果中有contigs.fasta)数量
                ass_failed: 组装失败数量
  - spades_sub(file_list)
        提交spades组装任务
        参数：
            list: [file, outdir, logfile, spades_path]
                file: 用于组装的fastq序列文件
                outdir: 组装结果输出目录
                logfile: 组装日志输出文件
                spades_path: spades.py程序路径
        返回：
            程序执行的返回值（0,1）

- FILE
      script/assemble_fastq_from_list.py



### merge_contigs:

- NAME
      merge_contigs - Created on Wed Nov  6 00:21:48 2019

- DESCRIPTION
      @author: yk

- FUNCTIONS
  - merge_contigs(contig_list, merge_fa_name)
            合并组装后contigs
            参数：
                contig_list: contig文件路径列表（type=list）
                merge_fa_name: 输出文件
            返回：
                df: 每个umiID的contig数目（type=pandas）

- FILE
      script/merge_contigs.py



### cut_fa_by_len:

- NAME
      cut_fa_by_len

- DESCRIPTION
      Created on Mon Nov 18 14:26:48 2019
      筛选fasta文件符合长度要求的序列
      @author: yk

- FUNCTIONS
  - cut_fa_by_len(input_fa, output_fa, min_len=0, max_len=0)
            筛选fasta文件符合长度要求的序列
            参数：
                input_fa: 输入fasta文件
                output_fa: 输出fasta文件
                min_len: 序列最小长度
                max_len: 序列最大长度
            返回：
                c: 过滤后序列数量

- FILE
     script/cut_fa_by_len.py



### submit_cdhit:

- NAME
      submit_cdhit

- FUNCTIONS
  - submit_cdhit(infa, outfa, identity, path, threads=1)
            序列聚类
            参数：
                infa: 输入fasta序列
                outfa: 输出fasta代表序列
                identity: 相似度
                path: cd-hit路径
                threads: 线程数
            返回：
                无

- FILE
    script/submit_cdhit.py



### get_consensus_seq_from_cdhit:

- NAME
      get_consensus_seq_from_cdhit

- FUNCTIONS
  - get_consensus_seq_from_cdhit(clstr, fa, outfa)
            从cd-hit聚类结果中过滤一致性umiID
            参数：
                clstr: cd-hit聚类结果文件
                fa: 用于聚类的所有序列fa文件
                outfa: 输出结果文件
            返回：
                total_size: 过滤后的序列数量
                rep_seq_tab: 代表序列及对应clust的序列ID（type=pandas）

- FILE
     script/get_consensus_seq_from_cdhit.py



### remove_chimeras_uchime:

- NAME
      remove_chimeras_uchime

- FUNCTIONS
  - remove_chimeras_uchime(fa, uchimeout, chf, nonchf, usearch_path)
            去除嵌合序列
            参数：
                fa: 输入fasta文件
                uchimeout: 去除嵌合过程文件
                chf: 嵌合序列
                nonchf: 非嵌合序列
                usearch_path: usearch路径
            返回：
                无

- FILE
     script/remove_chimeras_uchime.py



### update_seqID:

- NAME
      update_seqID

- FUNCTIONS
  - update_seqID(infa, outfa, tag)
            重新命名序列
            参数：
                infa: 输入fasta文件
                outfa: 输出fasta文件
                tag: 序列标签
            返回：
                id_track: 序列ID对应关系（type=dict）

- FILE
      script/update_seqID.py



### submit_mothur:

- NAME
      submit_mothur - @author: yk

- FUNCTIONS
  - submit_mothur(infa, taxdb, tax, mothur, threads)
            物种分类
            参数：
                infa: 输入fasta文件
                taxdb: mothur数据库
                tax: mothur数据库分类文件
                mothur: mothur路径
                threads: 进程数
            返回：
                无

- FILE
      script/submit_mothur.py



### get_tax_from_mothur_with_level:

- NAME
      get_tax_from_mothur_with_level

- FUNCTIONS
  - get_tax_from_mothur_with_level(tax_file, asv_file, tax_ratio_file, tax_reads_file)
            统计物种分类结果
            参数：
                tax_file: mothur注释结果（.taxonomy文件）
                asv_file: asv或者otu表格（每个asv或者otu的数量）
                tax_ratio_file: 输出文件（物种序列比例）
                tax_reads_file: 输出文件（物种序列数量）
            返回：
                无

- FILE
      script/get_tax_from_mothur_with_level.py
