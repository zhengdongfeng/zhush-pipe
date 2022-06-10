
================================在工作目录下创建目录========================
#1 prepare  
mkdir rawdata
mkdir report
mkdir database ## to prepare genome, mRNA, pep, do blastp, and go annotation, ko annotation  and so on 
mkdir result;cd result

============================生成样本信息文件，链接原始数据==================
#s0 
mkdir 0.info
vi 0.info/sampleinfo.xls
cd ../rawdata
ls ../*.gz  | while read i; do ln -s $i ./ ;done
cd ../result

#data qc
================================测序数据质控===============================
# do qc  s1
/home/rongzhengqin/bin/fqQC/fastqQC --run-fastqc --seqlen=150 --outdir=1.QC --indir=../rawdata 0.info/sampleinfo.xls
ls *.sh | while read i; do qsub -l ncpus=4 $i ;done

# do qc s2
cd ../rawdata && ls *.zip  | while read i; do unzip $i; done && cd  -
cd result
doqsub "/home/rongzhengqin/bin/fqQC/fastqQC --run-stat --seqlen=150 --outdir=1.QC --indir=../rawdata --qual=33 0.info/sampleinfo.xls" （注意seqlen， qual的修改）


=============================数据比对=============================================
#若注释文件为gff/gff3格式的话需要先将gff/gff3转换成gtf格式
gffread -T -o gffread.gtf Ptrichocarpa_210_v3.0.gff3

#s1 genome index
/opt/software/STAR/bin/Linux_x86_64/STAR --runThreadN 12 --limitGenomeGenerateRAM 100000000000 --runMode genomeGenerate -genomeDir ./GenomeDir/ --genomeFastaFiles yourGenome.fa  --sjdbGTFfile yourGtf
#doqsub -m 大小可根据基因组实际大小选择，10X为宜

#s2 run aln
/data/BIN/zhush/bin/pipe/GenomeGuild/STAR_script --num-threads 4 --genomeDir ../database/genome-fmt-index/ --readFilesIn ../rawdata/ --sjdbGTFfile ../database/Arabidopsis_thaliana.TAIR10.27.fmt.gtf --twopassMode Basic --readFilesCommand zcat --outSAMattrIHstart 0 --outSAMunmapped 'Within KeepPairs' --outSAMtype 'BAM SortedByCoordinate' --strand_specific unstranded 0.info/sampleinfo.xls
ls STARaln_*.sh | while read i; do doqsub -t 4 -m 40G ${i};done

#If the aln job failed ,just re-run STAR with parameters "--limitBAMsortRAM 4172944425" as the error message suggested

===========================================重构转录本===============================
grep -v "#" 0.info/sampleinfo.xls | cut -f 1 | while read i ;do echo "stringtie -p 2 -G /data/database/ftp.ensembl.org/release-81/Rattus/Rattus_norvegicus.Rnor_6.0.84.fmt.gtf -o ${i}/${i}_STR.gtf ${i}/Aligned.sortedByCoord.out.bam" > ${i}.cl.sh; done
ls *.cl.sh | while read i; do doqsub -t 2 ${i} ;done

# merge gtfs
ls ./*/*_STR.gtf >> gtf_samples.txt
mkdir 3.merged_asm
stringtie --merge -G /data/database/Gencode/r24/gencode.v24.annotation.gtf -o 3.merged_asm/merged.gtf -T 1 -f 0.01 gtf_samples.txt

cd 3.merged_asm
python /home/rongzhengqin/pipe/RNAseq/merged_gtf2novel4sT.py  merged.gtf   ##  to mark novel ids，产生merged.fmt.gtf 文件
gffread ../3.merged_asm/merged.fmt.gtf  -g /data/database/Gencode/r24/hg38.fa  -w exons.fa

# to filter  rRNA, tRNAs, RuBisCO 
#doqsub -t 16 "blastn -query exons.fa -num_threads 16 -db /data/database/rRNA_database/rRNA.5S.SSU.LSU.fmt.fa -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -out blastn.rRNA.result.xls"
#python /home/rongzhengqin/pipe/LNC/filter.py merged.fmt.gtf blastn.rRNA.result.xls > merged.gtf 
#gffread ../3.merged_asm/merged.gtf  -g /data/database/Gencode/r24/hg38.fa  -w exons.fa



# run mapping QC
========================================比对结果的基本统计=====================================================
mkdir 2.mapQC && cd 2.mapQC
doqsub "python /home/rongzhengqin/bin/MappingQC/mappingQC.py ../0.info/sampleinfo.xls star ../"
cd ..


=====================================RNA 水平的质量控制========================================================
# run RNA qc
ls ./*/Aligned.sortedByCoord.out.bam | while read i; do samtools index $i ;done
ls ./*/Aligned.sortedByCoord.out.bam | while read i ; do doqsub -t 1 "python /home/rongzhengqin/bin/MappingQC/RNAseq-RNAqc.py ./3.merged_asm/merged.gtf ${i}" ;done

====wait
mkdir 5.rnaqc && cd 5.rnaqc
python /home/rongzhengqin/bin/MappingQC/RNAseq-RNAqc_plot.py ../0.info/sampleinfo.xls ../ /Aligned.sortedByCoord.out.bam.baohe_cover.dat.npz
cd ..


## to quant
====================================RNA 表达定量==========================================================
************************************************************************************* method 2 use RSEM
grep ">" 3.merged_asm/exons.fa | sed s'/>//'g | awk '{print $2"\t"$1}' | sed s'/gene=//'g > gene.transcript.table
/data/tools/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts 3.merged_asm/exons.fa --est_method RSEM --aln_method bowtie2 --gene_trans_map gene.transcript.table --prep_reference --output_dir 3.merged_asm/
python /home/rongzhengqin/pipe/LNC/s0_aln_estimate.py  0.info/sampleinfo.xls ../rawdata/ 3.merged_asm/exons.fa 1
ls *.exprs_estimate_aln.sh  | while read i; do doqsub -t 4 ${i} ;done

python /home/rongzhengqin/pipe/LNC/s1_estimate_matrix.py  0.info/sampleinfo.xls  2 && mv Exprs/DAT.matrix.anno Exprs/DAT.matrix.anno.TPM
python /home/rongzhengqin/pipe/LNC/s1_estimate_matrix.py  0.info/sampleinfo.xls  0 && mv Exprs/DAT.matrix.anno Exprs/DAT.matrix.anno.Count

mkdir norm 
==对novel 进行过滤，满足至少有3个样本，TPM > 1，2个样本时选择1，2 则表示筛选至少在1个样本中，TPM > 0.5的转录本  
python /home/rongzhengqin/pipe/RNAseq/filter_matrix.py Exprs/DAT.matrix.anno.TPM  1 0.5 > norm/screen.isoforms.tpm_table.matrix.anno.xls 
grep -v "#"  norm/screen.isoforms.tpm_table.matrix.anno.xls  | cut -f 2 | python /home/rongzhengqin/scripts/select_tabfile.py - Exprs/DAT.matrix.anno.Count 1 > norm/screen.isoforms.count_table.matrix.anno.xls 

***********************************************************************************************************


=================================================预测lncRNA================================================
mkdir 6.lncRNA_predict && cd 6.lncRNA_predict
sort -t $'\t' -k1,1 -k4,4n -k5,5n /data/database/ftp.ensembl.org/release-81/Rattus/Rattus_norvegicus.Rnor_6.0.84.fmt.gtf > sorted.ref.gtf
python /home/rongzhengqin/bin/fmt_2gtf.py  sorted.ref.gtf ../3.merged_asm/merged.fmt.gtf
cut -f 2 ../norm/screen.isoforms.tpm_table.matrix.anno.xls | python /home/rongzhengqin/scripts/select_tabfile.py - stringtie.info.xls  5 > tmp && mv tmp stringtie.info.xls
python /home/rongzhengqin/pipe/RNAseq/lnc_mRNA_2predict.py stringtie.info.xls /data/database/ftp.ensembl.org/release-81/Rattus/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fmt.fa

===use cpat
# if all, will use 4 cpu cores,  for 1 model, use 1 cpu core
python /home/rongzhengqin/bin/lnc4cpat.py unknown.fa Human ## or Fly  Mouse Zebrafish  or all
grep noncoding coding_predict.result.xls | cut -f 1 | python /home/rongzhengqin/scripts/select_fa.py - unknown.fa > predict.noncoding.fa
grep -P "\tcoding" coding_predict.result.xls | cut -f 1 | python /home/rongzhengqin/scripts/select_fa.py - unknown.fa > predict.coding.fa
rm unknown.fa *.out *.out.dat *.out.r

mkdir scan_miR && cd scan_miR
python /home/rongzhengqin/pipe/miRseq/scan_predict.py ../predict.noncoding.fa animal 1 300  #  if plant use plant
sh Scan_microRNAs.sh
cd ..
python /home/rongzhengqin/pipe/miRseq/scan_predict_result_parse.py scan_miR/*/*.out  && rm -rf scan_miR
python /home/rongzhengqin/pipe/RNAseq/lnc_stat_stringtie_predict.py stringtie.info.xls coding_predict.result.xls miRpare_predict_result.xls

=== 统计lnc 的种类，放在后面， 比如 intron ， anti 的 ，exon-intron 的 等等
=== 转录本长度分布，predict_lncRNA，predict_sRNA，predict_coding, coding , lncRNA,sRNA,imm.fa,pseudo.fa
=== exon个数分布，predict_lncRNA，predict_sRNA，predict_coding, coding , lncRNA,sRNA,imm.fa,pseudo.fa （分成1，2，3，4，5，>5）
python /data/BIN/zhush/bin/pipe/LNC/gene_exon_stat.py coding,know_nc,novel_coding,novel_nc ref_novel_predicted.info.xls coding.fa lnc.fa predict.coding.fa predict.noncoding.fa
cd ..


=================================================== 对转录组序列 进行功能注释 =============================
### coding.fa  imm.fa  lnc.fa  predict.coding.fa  predict.noncoding.fa  pseudo.fa  srna.fa
========s1. 获取序列，核酸及pep预测
mkdir 8.seq_anno && cd 8.seq_anno 
# 若存在已知的序列, 得到参考的cds 序列：
gffread /data/database/Gencode/r24/gencode.v24.annotation.gtf -x known.cds.fa  -g /data/database/Gencode/r24/hg38.fa
ls ../6.lncRNA_predict/*.fa | while read i; do ln -s ${i} ./ ;done 

# use  known to predict novel seq 's  proteins 
---------------------------------------------更新transdecoder-------------------------------------------
如果是链特异性文库，TransDecoder.LongOrfs这一步(第一步)需要加-Ｓ参数，非链特异性时如下第一步所示，后面几步不变：
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.LongOrfs -t known.cds.fa -S
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.Predict -t known.cds.fa  --single_best_orf
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.LongOrfs -t predict.coding.fa -S
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.Predict -t predict.coding.fa --train known.cds.fa.transdecoder.pep --single_best_orf
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.LongOrfs -t coding.fa -S
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.Predict -t coding.fa --train known.cds.fa.transdecoder.pep --single_best_orf
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.LongOrfs -t imm.fa -S
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.Predict -t imm.fa --train known.cds.fa.transdecoder.pep --single_best_orf
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.LongOrfs -t pseudo.fa -S
/data/tools/trinityrnaseq-2.2.0/TransDecoder-3.0.0/TransDecoder.Predict -t pseudo.fa --train known.cds.fa.transdecoder.pep --single_best_orf

cat predict.coding.fa.transdecoder.pep known.cds.fa.transdecoder.pep > known.novel.transdecoder.pep

rm known.cds.fa.transdecoder.pep && ls *.pep | while read i; do python /data/BIN/suzhencheng/bin/extract_transdecoder_revise.py ${i};done

awk '{print $1"\t"$5"\t"$6"\t"$7"\t"$8}' predict.coding.fa.transdecoder.pep.info.txt > predict.coding.pep.orf.xls

python /home/rongzhengqin/pipe/RNAseq/orf_types_plot.py predict.coding.pep.orf.xls
rm -rf *.bed *.cds *.pep *.gff3 *.mRNA *.info.txt transdecoder.tmp*
python /data/BIN/zhush/bin/pipe/LNC/gene_exon_stat.py coding,know_nc,novel_coding,novel_nc ../6.lncRNA_predict/ref_novel_predicted.info.xls coding.fa lnc.fa predict.coding.fa predict.noncoding.fa

mkdir NR eggnog ipr kegg
#可以利用 /data/BIN/zhush/bin/tools/fasta_tool 工具来对fasta文件进行分割
#blastx 预测novel coding
(根据物种信息选择,参见/data/database/nr/README.txt)
doqsub -t 16 -s novelnr "blastp -num_threads 16 -query predict.coding.fa.transdecoder.pep.sel.fa -db /data/database/nr/current/NR_PRI(根据物种信息选择,参见/data/database/nr/README.txt) -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue stitle' -max_target_seqs 1 -evalue 0.00001 -out NR/blastp.predict.coding.result.xls"
doqsub -t 16 -s novelegg "blastp -query predict.coding.fa.transdecoder.pep.sel.fa -num_threads 16 -db /data/database/eggNOG/current/knowneggnog -outfmt 6 -max_target_seqs 5 -evalue  0.00001 -out eggnog/eggnog.result.xls"
doqsub -t 16 -s novelsp "blastp -num_threads 16 -query predict.coding.fa.transdecoder.pep.sel.fa -db /data/database/ftp.ncbi.nlm.nih.gov/2015-07-14/blast/swissprot -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -out NR/blastp.sp.predict.coding.result.xls"

# GO，kegg ，interproscan 用 pep.fa 预测
cat *.pep.sel.fa > known.novel.pep.fa 
doqsub -t 16 -s go "/data/tools/interproscan-5.15-54.0/interproscan.sh -appl PfamA,TIGRFAM,SMART,SuperFamily,PRINTS -dp -f tsv,html --goterms  --iprlookup -t p -i known.novel.pep.fa"
doqsub -t 16 -s kegg "blastp -query known.novel.pep.fa -num_threads 23 -db /data/database/KEGG/Current/Plants(根据物种信息选择)  -outfmt 6  -max_target_seqs 1 -evalue 0.00001 -out kegg/blast.ko.xls"
=== wait for summary ==== 
mv known.novel.pep.fa.html.tar.gz known.novel.pep.fa.tsv ipr/
# summary 
cd eggnog && python /home/rongzhengqin/bin/AnnoDB/blast2eggNOG.py eggnog.result.xls /data/database/eggNOG/v4/eggNOG.info.db && rm eggnog.result.xls && cd ..
cd NR && python  /home/rongzhengqin/bin/AnnoDB/DBanno_stat.py  blastp.predict.coding.result.xls blastp.sp.predict.coding.result.xls /data/project/blankfile NR,Swissprot,None ../predict.coding.fa && rm -rf DBanno_stat_venn_venn* && cd ..
cd kegg && python /home/zhush/scripts/getKO.py /data/database/KEGG/Current/KO.id.list blast.ko.xls /data/database/KEGG/Current/ko_pathway_name_class.list /data/database/KEGG/Current/ko_db.info && rm blast.ko.xls && cd ..
cd ipr && python /home/rongzhengqin/bin/AnnoGoKegg/ipr_go_fmt_v2.py known.novel.pep.fa.tsv /data/database/GO/Current/obo_topology.txt /data/database/GO/Current/GO_alt_id.alt_id > known.novel.pep.fa.tsv.topo.xls
mkdir Detail_anno && mv known.novel.pep.fa.html.tar.gz Detail_anno && cd Detail_anno && tar zxvf known.novel.pep.fa.html.tar.gz && cd ..
python /home/rongzhengqin/pipe/LNC/iprscananno2html.py known.novel.pep.fa.tsv ; cd ..

python /home/lizhenzhong/bin/LNC/DB_anno_total.py NR,eggNOG,kegg,ipr NR/DBanno_detail.xls,eggnog/gene2eggNOG.annotation.xls,kegg/blast.ko.xls.KO.xls,ipr/known.novel.pep.fa.tsv.topo.xls

#py /home/rongzhengqin/pipe/LNC/iprscananno2html.py 
#python /home/rongzhengqin/bin/AnnoGoKegg/ipr_go_fmt.py know.novel.transdecoder.pep.sel.fa.tsv  /data/database/GO/current/GO_obo_flat.txt /data/database/GO/current/GO_alt_id.txt > total.iprscan.tsv.fmt.xls


=======================================================snqc ===========================================================================
cd norm
python /home/rongzhengqin/pipe/RNAseq/split_matrix_anno2non_c.py  ../6.lncRNA_predict/ref_novel_predicted.info.xls screen.isoforms.count_table.matrix.anno.xls
python /home/rongzhengqin/pipe/RNAseq/split_matrix_anno2non_c.py  ../6.lncRNA_predict/ref_novel_predicted.info.xls screen.isoforms.tpm_table.matrix.anno.xls
cd ..
mkdir 7.snqc_coding && cd 7.snqc_coding
python /home/rongzhengqin/bin/exprs_QC/exprs_qc.py ../0.info/sampleinfo.xls ../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls 1
python /home/rongzhengqin/bin/stat_exprs_insamples.py ../0.info/sampleinfo.xls ../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls  12  ## 12 为样本数目(使用实际的样本数目)
cd ..
mkdir 13.snqc_noncoding && cd 13.snqc_noncoding
python /home/rongzhengqin/bin/exprs_QC/exprs_qc.py ../0.info/sampleinfo.xls ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls 1
python /home/rongzhengqin/bin/stat_exprs_insamples.py ../0.info/sampleinfo.xls ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls 12  ## 12 为样本数目
cd ..


# do diff analysis
=====================================================For 差异表达分析====================================================
# diff analysis
mkdir 9.diff_coding && cd 9.diff_coding 
awk -F"\t" '{print $1"|"$2"\t"$0}' ../norm/screen.isoforms.count_table.matrix.anno.coding.xls |  cut -f 1,4-100 > readcounts.xls
grep -v "#" ../0.info/sampleinfo.xls | awk  '{print $5"\t"$1}'  > samples.list
## data file (readcounts.xls) :
#gene   r1      r2      r3      r4      r5      （需确认样本名的准确性）
Unigene8220_All 13.0939143585184        3.61943302983515        1.59861917485423        1.81238339505037        1.22487385249009
## contrasts.list  （生成对比的列表，需自己填写）
cond_1  cond_2
cond_1  cond_5
cond_2  cond_4
cond_2  cond_5
## samples.list
cond_1  r1  
cond_2  r2  
cond_3  r3  
###########################
/data/tools/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix readcounts.xls --method edgeR --samples_file samples.list  --contrasts contrasts.list --output edgeR.diff.dir 
#/data/tools/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix readcounts.xls --method edgeR --samples_file samples.list  --contrasts contrasts.list --output edgeR.diff.dir --dispersion 0.16
ls edgeR.diff.dir/readcounts.xls.*.edgeR.DE_results | while read i; do python /home/rongzhengqin/pipe/LNC/s2_fmt_DEresult.py ${i} 1; done ## to check logFC\A3\A8A/B\A3\A9\A3\AC\C8\F4\D0\E8ҪlogFC(B/A)\A3\AC\D4\F2-1

awk -F"\t" '{print $1"_vs_"$2}'  contrasts.list  | while read i; do mkdir $i ;done

awk -F"\t" '{print $1"_vs_"$2}'  contrasts.list  | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.Diff.sig.xls ${i}/Diff.sig.xls ;done
awk -F"\t" '{print $1"_vs_"$2}'  contrasts.list  | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.Diff.total.xls ${i}/Diff.total.xls ;done
awk -F"\t" '{print $1"_vs_"$2}'  contrasts.list  | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.MA_n_Volcano.pdf ${i}/MA_n_Volcano.pdf ;done

cat  ../9.diff_coding/*/Diff.sig.xls | grep -v "^#" |  cut -f 2 | sort |  uniq | python /home/rongzhengqin/scripts/select_fa.py - ../3.merged_asm/exons.fa > diffexprs.transcripts_coding.fa
cd ..



=================================================== 富集分析=========================================================
# to do gokegg enrich analysis
mkdir 10-12.diff_go_kegg_coding; cd 10-12.diff_go_kegg_coding 
ln -s ../8.seq_anno/ipr/known.novel.pep.fa.tsv.topo.xls ./
ln -s ../8.seq_anno/kegg/blast.ko.xls.KO.xls ./

# 2 do analysis
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do mkdir $i ;done

按组别做注释 
## 注意，对于模式生物，必须添加如下参数：
在第一步中，添加
# --taxid=4530（4530为水稻） --geneinfo=/data/database/ftp.ncbi.nlm.nih.gov/2015-07-14/gene/gene_info 
在第二步中，前面的参数继续添加，同时，如果是植物，则调用plant（拟南芥+水稻的信息，以避免出现非植物类的通路）
# --pathway-list=/data/database/KEGG/Current/plant_pathway_name_class.list
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/lihuamei/bin/go_Annotation.py --nohuman --gene2ipr=../known.novel.pep.fa.tsv.topo.xls --genesigxls=../../9.diff_coding/${i}/Diff.sig.xls --bg=../../9.diff_coding/${i}/Diff.total.xls --sampleinfo=../../0.info/sampleinfo.xls --annodata=../../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls --maxnum=30 && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/lihuamei/bin/kegg_ipr_Anotation.py --nohuman --gene2ipr=../known.novel.pep.fa.tsv.topo.xls  --ko=../blast.ko.xls.KO.xls --genesigxls=../../9.diff_coding/${i}/Diff.sig.xls --bg=../../9.diff_coding/${i}/Diff.total.xls --sampleinfo=../../0.info/sampleinfo.xls --annodata=../../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/lihuamei/bin/anno_results_fmt.py --run-stat --sig_go=./go_annotation/Diff_exprs.GO_enrich.lst --sig_kegg=./kegg_annotation/Diff_exprs.KEGG_enrich.lst --sig_ipr=./IPR_annotation/Diff_exprs.IPR_enrich.lst && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/chenghh/bin/pipe/sig_anno_nr.py ../../../8.seq_anno/NR/DBanno_detail.xls Diff_exprs.sig.anno.xls ; done
cd .. 



=================================================lncRNA 分析========================================================
=====ncRNA 注释, 与已知基因间位置关系，判断lncRNA种类，位置注释信息====
mkdir 14.ncRNA_genomic_anno && cd 14.ncRNA_genomic_anno 
python /home/rongzhengqin/pipe/RNAseq/lnc_anno.py ../6.lncRNA_predict/ref_novel_predicted.info.xls ../6.lncRNA_predict/refGTF.info.xls.gz
python /home/rongzhengqin/pipe/RNAseq/lnc_anno_summary.py lnc_gene_interact.resulut.xls ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls > intergenic.genes.list
cat intergenic.genes.list | python /home/rongzhengqin/scripts/select_tabfile.py - ../6.lncRNA_predict/ref_novel_predicted.info.xls 5 > intergenic.ncRNA.info.xls

===============================================================ncRNA 差异分析=========================================
mkdir 15.diff_noncoding && cd 15.diff_noncoding
awk -F"\t" '{print $1"|"$2"\t"$0}' ../norm/screen.isoforms.count_table.matrix.anno.noncoding.xls |  cut -f 1,4-100 > readcounts.xls
/data/tools/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix readcounts.xls --method edgeR --samples_file ../9.diff_coding/samples.list  --contrasts ../9.diff_coding/contrasts.list --output edgeR.diff.dir
#对无生物学重复差异分析
#/data/tools/trinityrnaseq-2.2.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix readcounts.xls --method edgeR --samples_file ../9.diff_coding/samples.list  --contrasts ../9.diff_coding/contrasts.list --output edgeR.diff.dir --dispersion 0.16
ls edgeR.diff.dir/readcounts.xls.*.edgeR.DE_results | while read i; do python /home/rongzhengqin/pipe/LNC/s2_fmt_DEresult.py $i 1; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do mkdir $i ;done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.Diff.sig.xls ${i}/Diff.sig.xls ;done 
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.Diff.total.xls ${i}/Diff.total.xls ;done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list | while read i; do cp edgeR.diff.dir/readcounts.xls.${i}.edgeR.DE_results.MA_n_Volcano.pdf ${i}/MA_n_Volcano.pdf ;done

cat  ../15.diff_noncoding/*/Diff.sig.xls | grep -v "^#" |  cut -f 2 | sort |  uniq | python /home/rongzhengqin/scripts/select_fa.py - ../3.merged_asm/exons.fa >diffexprs.transcripts_noncoding.fa
cd ..



===============================================================ncRNA -func =======================================

mkdir 16.diff_ncRNA_func && cd 16.diff_ncRNA_func
========1 position related ncRNA
cat  ../15.diff_noncoding/*_vs_*/Diff.sig.xls | grep -v "#" | cut -f 2 | python /home/rongzhengqin/scripts/select_tabfile.py - ../14.ncRNA_genomic_anno/lnc_gene_interact.resulut.xls 1 > pos_related_diffexprs_nc.interact.info.xls
========2 get diff exprs intergenic ncRNA result
cat  ../15.diff_noncoding/*_vs_*/Diff.sig.xls | grep -v "#"  | cut -f 2  | python /home/rongzhengqin/scripts/select_tabfile.py - ../14.ncRNA_genomic_anno/intergenic.ncRNA.info.xls 5 > intergenic_diffexprs_ncRNA.info.xls

========3 计算位置相关的ncRNA 和 mRNA 间的表达相关性
python /home/rongzhengqin/pipe/RNAseq/co_exprs/filter_interaction_lnc_mRNA_coexprs.py ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls ../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls noncoding,coding  pos_related_diffexprs_nc.interact.info.xls 1_6 -1  > coexprs_pos_related_diffexprs_nc.interact.result.xls

========4 计算 intergenic ncRNA 和 其他 mRNA 间的表达关系
grep -v "^#" intergenic_diffexprs_ncRNA.info.xls | cut -f 6 | sort | uniq > diffexprsNC.list
cat diffexprsNC.list | python /home/rongzhengqin/scripts/select_tabfile.py - ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls 1 > diffexprsNC.exprs.matrix.list
grep -v "#" ../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls | cut -f 2  >  coding.list
python /home/rongzhengqin/pipe/RNAseq/co_exprs/ncRNAvsmRNA_list.py diffexprsNC.list coding.list intergenic > intergenicNC_C.interact.list
python /home/rongzhengqin/pipe/RNAseq/co_exprs/filter_interaction_lnc_mRNA_coexprs.py ../norm/screen.isoforms.tpm_table.matrix.anno.noncoding.xls ../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls noncoding,coding intergenicNC_C.interact.list 0_1 2 > coexprs_diffexprs_intergenic_nc.interact.resulut.list
python /home/rongzhengqin/pipe/RNAseq/co_exprs/screenPvalue.py coexprs_diffexprs_intergenic_nc.interact.resulut.list 0.00001 > coexprs_diffexprs_intergenic_nc.interact.0.00001.xls 

mkdir plot && cd plot
python /home/rongzhengqin/bin/plot_gene/plot_lncVSgene.py ../coexprs_pos_related_diffexprs_nc.interact.result.xls ../../norm/screen.isoforms.tpm_table.matrix.anno.xls
python /home/rongzhengqin/bin/plot_gene/get_to_plot.py ../coexprs_pos_related_diffexprs_nc.interact.result.xls ../../6.lncRNA_predict/ref_novel_predicted.info.xls > tmp.xls
python /home/rongzhengqin/bin/plot_gene/plot_lncStruct.py ../coexprs_pos_related_diffexprs_nc.interact.result.xls tmp.xls > run.sh
doqsub -t 2 run.sh
cd ..

grep -v "#" coexprs_pos_related_diffexprs_nc.interact.result.xls | cut -f 2,7 > interact.list
grep -v "#" coexprs_diffexprs_intergenic_nc.interact.0.00001.xls | cut -f 1,2 >> interact.list
#cd ..

lncRNA和MRNA共表达的散点图
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do mkdir $i ;done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i;do cd ${i} && cat ../../15.diff_noncoding/$i/Diff.total.xls ../../9.diff_coding/$i/Diff.total.xls > total.diffgene.xls && cd ..; done
grep -v "#" coexprs_pos_related_diffexprs_nc.interact.result.xls|cut -f 2,7|sort -u > coexprs.list
？？awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i;do cd  ${i} &&  python /data/BIN/donglh/bin/coexprs-plot.py total.diffgene.xls ../coexprs.list ; done
如果报错说xy轴数目不同，则表明有一部分的lncRNA的转录本或者mRNA的转录本是在total中不显著（尤其是多个样本），所以需要对样本进行过滤后再构图
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/donglh/bin/matplot/coexprs/select-mRNAco.py total.diffgene.xls ../coexprs.list > mRNA-coexpr.xls && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/donglh/bin/matplot/coexprs/select-lncRNAco.py total.diffgene.xls ../coexprs.list > lncRNA-coexpr.xls && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  |while read i; do cd ${i} && python /data/BIN/donglh/bin/matplot/coexprs/select-mRNA-lncRNA.py mRNA-coexpr.xls lncRNA-coexpr.xls > new-coexprs.list && cd .. ; done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i;do cd ${i} &&  python /data/BIN/donglh/bin/matplot/coexprs/coexprs-plot.py total.diffgene.xls new-coexprs.list && cd .. ; done
cd ..

=================================================================enrich ==========================================================================
mkdir 17.diff_ncRNA_enrich && cd 17.diff_ncRNA_enrich
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do mkdir $i ;done
awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do python /home/rongzhengqin/pipe/RNAseq/get_lncRNAtarget_sigfile.py ../16.diff_ncRNA_func/interact.list ../15.diff_noncoding/${i}/Diff.sig.xls > ${i}/Diff.sig.xls ;done 

awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/lihuamei/bin/go_Annotation.py --nohuman --gene2ipr=../../8.seq_anno/ipr/known.novel.pep.fa.tsv.topo.xls  --genesigxls=../../17.diff_ncRNA_enrich/${i}/Diff.sig.xls --bg=../../9.diff_coding/${i}/Diff.total.xls --sampleinfo=../../0.info/sampleinfo.xls --annodata=../../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls --maxnum=30 && cd .. ; done

awk -F"\t" '{print $1"_vs_"$2}'  ../9.diff_coding/contrasts.list  | while read i; do cd ${i} && python /data/BIN/lihuamei/bin/kegg_ipr_Anotation.py --nohuman --gene2ipr=../../8.seq_anno/ipr/known.novel.pep.fa.tsv.topo.xls --ko=../../8.seq_anno/kegg/blast.ko.xls.KO.xls --genesigxls=../../17.diff_ncRNA_enrich/${i}/Diff.sig.xls --bg=../../9.diff_coding/${i}/Diff.total.xls --sampleinfo=../../0.info/sampleinfo.xls --annodata=../../norm/screen.isoforms.tpm_table.matrix.anno.coding.xls && cd .. ; done 


=============================================================== Call SNPs，Indels =====================================================
# call variants
#python /home/rongzhengqin/bin/variant_call/make_calling_sh.py 0.info/sampleinfo.xls  /hisat2.sort.bam  ../database/Cyprinus_carpio/GCA_000951615.1_common_carp_genome_genomic.fna  > total.variant.sh
python /home/rongzhengqin/bin/variant_call/make_calling_sh.py 0.info/sampleinfo.xls /Aligned.sortedByCoord.out.bam /data/database/ftp.ensembl.org/release-81/Rattus/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fmt.fa > total.variant.sh
doqsub -t 1 -m 40G total.variant.sh

mkdir 19.variants && cd 19.variants
grep -v "#" ../0.info/sampleinfo.xls  | cut -f  1 | while read i; do echo "/opt/software/bin/samtools mpileup -Q13 -q20 -f ../../database/GCA_000224295.1_LinUsi_v1.1_genomic.fna  ../${i}_tophat2/accepted_hits.bam -u -t DP | bcftools call -vm -O v -o ${i}.snp.vcf" > ${i}".samtools.snp.sh";done
ls *.samtools.snp.sh | while read i; do qsub -l ncpus=1 $i; done

ls *.vcf  | while read i ; do python /home/rongzhengqin/bin/variant_call/vcf2annovar.py $i ;done

python /home/rongzhengqin/bin/variant_call/vcf2annovar.py ../total.snp.vcf 

awk -F"\t" '{if($11>20) print $0 }' ../total.snp.vcf.fmt | sort -t $'\t' -k1,1 -k2,2n -k3,3n > ../qual20_total.snp.fmt


python /home/rongzhengqin/bin/variant_call/mergepop.py  ../qual20_total.snp.fmt ../0.info/sampleinfo.xls 0 ./  ./

#awk -F"\t" '{if($6>20) print $0 }' total.snp.fmt | sort -t $'\t' -k1,1 -k2,2n -k3,3n > qual20_total.snp.fmt
#awk -F"\t" '{if($6>20) print $0 }' total.indel.fmt | sort -t $'\t' -k1,1 -k2,2n -k3,3n > qual20_total.indel.fmt


==============================  Annotation SNPs and INDELs ========================================================
python  /home/rongzhengqin/bin/variant_call/gtf_gff_varanno.py  gtf Homo_sapiens.GRCh38.77.filter.gtf ../../database/hg38.fa > build.sh (若使用gff3 格式的注释文件，则采用gtf->gff   注释文件改为Homo_sapiens.GRCh38.77.filter.gff3)
sh build.sh
#/opt/software/annovar/table_annovar.pl total.snp.fmt test/ -buildver test -protocol  ensGene -operation g

ls *.snp.fmt | while read i; do /opt/software/annovar/table_annovar.pl $i test/ -buildver test -protocol  ensGene -operation g  ; done
ls *.indel.fmt | while read i; do /opt/software/annovar/table_annovar.pl $i test/ -buildver test -protocol  ensGene -operation g  ; done

# ls *.ensGene.variant_function   | while read i; do sed s/ncRNA/RNA/g $i tmp && mv tmp $i  ; done (#若 gtf 使用的是merged gtf)

python /home/rongzhengqin/bin/variant_call/mut_sub_pattern.py *.snp.fmt

python /home/rongzhengqin/bin/variant_call/varanno.stat.py  *.snp.fmt ensGene > SNP.anno.stat.xls 
python /home/rongzhengqin/bin/variant_call/plot_varanno.stat.py SNP.anno.stat.xls 
ls Annotation_Type_stat.* | while read i; do mv ${i} SNP_${i} ;done 

python /home/rongzhengqin/bin/variant_call/varanno.stat.py  *.indel.fmt ensGene > INDEL.anno.stat.xls
python /home/rongzhengqin/bin/variant_call/plot_varanno.stat.py INDEL.anno.stat.xls
ls Annotation_Type_stat.* | while read i; do mv ${i} InDel_${i} ;done

ls total.*.fmt  | while read i; do mv ${i} ${i}.xls ;done
ls total.*.test_multianno.txt |  while read i; do mv ${i} ${i}.xls ;done
rm -rf  *.log *.exonic_variant_function *.variant_function *.sh *.sh.e* *sh.o* ref* test *.vcf *.fmt *.txt *.invalid_input 
cd ..

# to do asprofile 最好分成 ncRNA 剪切模式 和 mRNA 剪切模式
=============================================可变剪接，仅模式生物，且放在最后跑此任务================================
mkdir asprofile && cd asprofile
python /home/rongzhengqin/bin/asprofile/asprofile_hdrs.py /data/database/ftp.ensemblgenomes.org/oryza_sativa/Oryza_sativa.IRGSP-1.0.25.dna_sm.toplevel.fa.gz genome
ls ../*/*_STR.gtf | while read i; do ln -s $i ./ ;done 
ls *_STR.gtf | while read i; do cuffcompare -r ../../database/Oryctolagus_cuniculus.OryCun2.0.83.filter.gtf ${i} ;done 

grep -v "#" ../0.info/sampleinfo.xls  | cut -f 1 | while read i; do /opt/software/ASprofile.b-1.0.4/extract-as ${i}_STR.gtf genome.hdrs -r cuffcmp.${i}_STR.gtf.tmap  ../../database/Oryctolagus_cuniculus.OryCun2.0.83.filter.gtf > ${i}.tmap.as ; done

grep -v "#" ../0.info/sampleinfo.xls  | cut -f 1 | while read i; do perl /opt/software/ASprofile.b-1.0.4/summarize_as.pl ${i}_STR.gtf ${i}.tmap.as -p ${i}_STR.gtf ;done  

ls *.as.nr | while read i; do python /home/rongzhengqin/bin/asprofile/result_fmt.py ${i} None > ${i}.splicing.events.stat.xls ;done
python /home/rongzhengqin/bin/asprofile/plot4samples.py *.as.nr.splicing.events.stat.xls
#python ~/bin/asprofile/result_fmt.py merged.as.nr ../norm/genes.attr_table   > splicing.events.stat.xls
rm *.gtf *.tmap *.refmap cuffcmp.* genome.hdrs *.tmap.as *.as.summary *.as.nr 
cd .. && mv asprofile 18.splicing_events



============================================== 生成报告：=======================================================
#mv merged_asm 3.merged_asm
mv norm 4.norm 
rm 4.norm/*_table ./*/*/*/*.lst */*.list ./*/*.sh ./*/*.sh.e* ./*/*.sh.o* 
cd 3.merged_asm && rm exons.fa.bal merged.fmt.gtf blastn.rRNA.result.xls  && cd ..
cd 4.norm && rm run.info samples.table coding.readcounts.xls noncoding.readcounts.xls *.anno && cd ..
cd 6.lncRNA_predict && rm  *tmp* class_code_stat.*  stringtie.info.xls sorted.ref.gtf  && cd ..
cd 8.seq_anno && rm coding.fa imm.fa lnc.fa predict.coding.fa predict.noncoding.fa pseudo.fa srna.fa  && cd ..

#把 result 目录中，以数字编号开头的文件夹《拷贝》到report 文件夹中
ls | grep -P '^[0-9]+' | while read i; do cp -R $i ./../report ; done 
cd ../report
cp ../result/*.html ./

#然后 拷贝配置文件：
cp -R /data/project/NCRNA/CSS .
cp -R /data/project/NCRNA/HELP .
cp -R /data/project/NCRNA/dist .

# 删除一些中间文件

python /data/project/NCRNA/genome_guild_ncRNA_v5.py 0.info/sampleinfo.xls 

================================================ 上传报告=======================================================
python /home/lizhenzhong/bin/ftp_client/ftp_client.py -I report.tar.gz --project HN16025000 --manager xuzhenhua


