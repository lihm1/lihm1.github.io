---
title: danseq-1
date: 2022-04-06 10:59:49
index_img: https://s2.loli.net/2022/04/06/JUSDs9f8GKChMH6.jpg
banner_img: /img/pic4.jpg 
categories:
- bioinformation
---
原始数据的处理流程：
主要用来记录，代码偏多，解释很少，方便自己后面使用

### STEP_1.数据清洗
下机数据的格式为fastq的格式，我们可以先查看数据排序是否一致：
~~~
gunzip -c test.fq.gz | paste - - - - | cut -f 1| head
~~~
进行数据清洗的工具：Soapnuke
~~~
SOAPnuke2.0.7 filter -1 test_1.fq.gz -2 test_2.fq.gz -l 5 -q 0.5 -n 0.1 -f adapter1 -r adapter2 -Q 2  -G 2 --seqType 0 -o 01-clean  -C test_1.clean.fq.gz -D test_2.clean.fq.gz
~~~
工具参数：
~~~
-f, –adapter1  : <s> 3′ adapter sequence of fq1 file

-r, –adapter2  : <s> 5′ adapter sequence of fq2 file [only for PE reads]

-1, –fq1       : <s> fq1 file

-2, –fq2       : <s> fq2 file, used to pe

–tile          : <s> tile number to ignore reads , such as [1101-1104,1205]
-M, –misMatch  : <i> the max mismatch number when match the adapter (default: [1])

-A, –matchRatio: <f> adapter’s shortest match ratio(default: [0.5])

-l, –lowQual   : <i> low quality threshold (default: [5])

-q, –qualRate  : <f> low quality rate (default: [0.5])

-n, –nRate     : <f> N rate threshold (default: [0.05])

-m, –mean      : <f> filter reads with low average quality, (<)

-p, –polyA     : <f> filter poly A, percent of A, 0 means do not filter (default: [ 0 ])

-d, –rmdup     : <b> remove PCR duplications

-i, –index     : <b> remove index

-c, –cut       : <f> the read number you want to keep in each clean fq file (unit:1024*1024, 0 means not cut reads)

-t, –trim      : <s> trim some bp of the read’s head and tail, they means:

read1’s head and tail and read2’s head and tail(default: [ 0,0,0,0 ])

-S, –small     : <b> filter the small insert size
-O, –overlap   : <i> minimun match length (default: [ 10 ])
-P, –mis       : <f> the maximum miss match ratio (default: [ 0.1 ])
-Q, –qualSys   : <i> quality system 1:illumina, 2:sanger (default: [ 1 ])

-L, –read1Len  : <i> read1 max length (default: all read1’s length are equal, and auto acquire)

-I, –read2Len  : <i> read2 max length (default: all read2’s length are equal, and auto acquire)

-G, –sanger    : <b> set clean data qualtiy system to sanger (default: illumina)

-a, –append    : <s> the log’s output place : console or file (default: [console])

-o, –outDir    : <s> output directory, directory must exists (default: current directory)

-C, –cleanFq1  : <s> clean fq1 file name

-D, –cleanFq2  : <s> clean fq2 file name

-E, –cutAdaptor: <i> cut sequence from adaptor index,unless performed -f/-r also in use,discard the read when the adaptor index of the read is less than INT

-b, –BaseNum   : <i> the base number you want to keep in each clean fq file,unless performed -E also in use

-R, –rawFq1    : <s> raw fq1 file name

-W, –rawFq2    : <s> raw fq2 file name

-5, –seqType   : <i> Sequence fq name type, 0->old fastq name, 1->new fastq name[default: 0]

old fastq name: [@FCD1PB1ACXX](https://weibo.com/n/FCD1PB1ACXX):4:1101:1799:2201#GAAGCACG/2

new fastq name: [@HISEQ](https://weibo.com/n/HISEQ):310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC

-6, –polyAType : <i> filter poly A type, 0->both two reads are poly a, 1->at least one reads is poly a, then filter, [default: 0]

-7, –outType: <i> Add /1, /2 at the end of fastq name, 0:not add, 1:add [default: 0]

-h, –help      : <b> help

-v, –version   : <b> show version
~~~
主要功能： 
● 主要功能：去adapter，index，read头尾剪切，低quality过滤，高比例N过滤。
● 使用平台：HiSeq2000, HiSeq 2500, HiSeq4000,HiSeqXTen,Zebra
● 测序策略：PE/SE 

##### 结果
我们可以看到得到的文件如下：
![](https://s2.loli.net/2022/04/06/fNqRSjle9s28wpK.png)

### STEP_2:数据质控
fqcheck也是一款常用的质控工具，但是它的质控结果并没有做可视化，而是直接输出一个文本文件，记录质控项结果，通常需要使用R或python或perl进行可视化。相比于fastqc，fqcheck可以接受pair-end的数据进行质控
软件位置：$fqcheck ||= "$BIN/fqcheck33"
~~~
fqcheck  -r 01-clean/test_1.clean.fq.gz -c 02-fqcheck/test-1.clean.fqcheck
~~~
输出结果：
包含有：总reads数，总碱基数，reads最大长度，reads平均长度，各个位置的碱基质量值分布，ATCG四种碱基的比例，GC含量，Q20，Q30等。

### STEP_3:sentieon
后面的步骤用的华大开发的sentieon，参考文档：https://support.sentieon.com/manual/DNAseq_usage/dnaseq/
#### 3.1 BWA比对
~~~
Sentieon bwa mem -M -t 40 -R "@RG\\tID:sampleID\\tLB:$KA\\tPU:sample\\tSM:135T-KAPA\\tPL:COMPLETE_GENOMICS\\tCN:BGI" Homo_sapiens_assembly38.fasta 01-clean/test_1.clean.fq.gz 01-clean/test_2.clean.fq.gz -o 03-sentieon/test-lane_barcode.sorted.bam
~~~

#### 3.2 Calculate data metrics （WGS)
A single command is run to generate 5 statistical summaries of the data quality and the pipeline data analysis quality results.
~~~
SentieonAuto driver -t 4 
-r Homo_sapiens_assembly38.fasta
-i 03-sentieon/test-lane_barcode.sorted.bam 
--algo GCBias 
--summary 04-sentieonauto/test.sorted.total.GC_SUMMARY.txt 04-sentieonauto/test.sorted.total.GC_METRIC.txt 
--algo MeanQualityByCycle 04-sentieonauto/test.sorted.total.MQ_METRIC.txt 
--algo QualDistribution 04-sentieonauto/test.sorted.total.QD_METRIC.txt 
--algo InsertSizeMetricAlgo 04-sentieonauto/test.sorted.total.IS_METRIC.txt
--algo AlignmentStat 04-sentieonauto/test.sorted.total.ALN_METRIC.txt 
--algo WgsMetricsAlgo 04-sentieonauto/test.sorted.total.WGS_METRIC.txt 
--algo SequenceArtifactMetricsAlgo 
--dbsnp dbsnp_151.hg38.vcf.gz 
04-sentieonauto/test.sorted.total.ARTIFACT_METRICS_BASE.txt
~~~

#### 3.3 Remove or mark duplicates 
分两步
~~~
Sentieon driver -t 8 -i 03-sentieon/test-lane_barcode.sorted.bam --algo LocusCollector --fun score_info  05-makdup/test-lane_barcode.mkdup.score.txt.gz

Sentieon driver -t 8 -i 03-sentieon/test-lane_barcode.sorted.bam --algo Dedup --score_info 05-makdup/test-lane_barcode.mkdup.score.txt.gz --metrics 05-makdup/test-lane_barcode.mkdup.metrics.txt  05-makdup/test-lane_barcode.mkdup.bam
~~~

#### 3.4 Indel realignment (optional)
~~~
SentieonAuto driver -t 8 -r Homo_sapiens_assembly38.fasta -i 05-makdup/test-lane_barcode.mkdup.bam --algo Realigner -k 1000G_phase1.snps.high_confidence.hg38.vcf.gz -k Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  06-realign/test-lane_barcode.mkdup.realn.bam
~~~
#### 3.5 Base quality score recalibration（WGS)
bqsr & snp, indel calling
BQSR 全称叫做 Base Quality Score Recalibration， 可以理解为碱基质量校正。对于变异位点的鉴定，碱基质量是非常重要的。比如测序识别到的一个位点，其碱基和参考基因组上的碱基不同，但是其质量值特别低，此时可以认为是一个测序错误，而不是一个SNP位点。
在测序的原始数据中，本身就提供了每个碱基对应的质量值，但是GATK官方认为测序仪提供的碱基质量值，是不准确的，存在误差的。
某个位点前后的碱基的种类，称之为上下文环境，会对这个碱基的质量值产生影响。对于A,T,C,G 4种碱基，共有4 x 4 =16 种上下文环境，左侧的图是利用fastq文件中测序仪给出的碱基质量值做的图，可以看到，对于不同的上下文环境，碱基质量值分布不同；右图为经过BQSR校正之后，不同上下文环境中碱基质量的分布。可以看到，校正之后，不同的上下文环境的碱基质量分布基本相同。也就是说，BQSR消除了上下文环境对碱基质量的影响。
https://cloud.tencent.com/developer/article/1626269

~~~
Sentieon driver -t 8 -i 06-realign/test-lane_barcode.mkdup.realn.bam
-r Homo_sapiens_assembly38.fasta 
--algo QualCal 
-k dbsnp_151.hg38.vcf.gz  
-k 1000G_phase1.snps.high_confidence.hg38.vcf.gz  
-k Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  
07-BQSR/test.bqsr.RECAL_PRE.txt

Sentieon driver -t 8 
-r Homo_sapiens_assembly38.fasta 
-i 06-realign/test-lane_barcode.mkdup.realn.bam 
-q 07-BQSR/test.bqsr.RECAL_PRE.txt --algo QualCal 
-k dbsnp_151.hg38.vcf.gz  
-k 1000G_phase1.snps.high_confidence.hg38.vcf.gz  
-k Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  
07-BQSR/test.bqsr.RECAL_POST.txt 
--algo ReadWriter 07-BQSR/test.bqsr.bam 
--algo Haplotyper 
-d test.Haplotyper.vcf.gz     ###Genotyper, Haplotyper: used for germline variant calling analysis.

Sentieon driver -t 8  --algo QualCal --plot 
--before 07-BQSR/test.bqsr.RECAL_PRE.txt --after 07-BQSR/135T-KAPA.bqsr.RECAL_POST.txt 
07-BQSR/test.bqsr.RECAL_RESULT.csv

Sentieon plot QualCal -o 07-BQSR/test.bqsr.RECAL_RESULT.pdf  
07-BQSR/test.bqsr.RECAL_RESULT.csv
~~~
##### 结果
![](https://s2.loli.net/2022/04/06/vOYlDhAC4qa859u.png)

#### germline mutation calling
~~~



####VARcal SNP



/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/SentieonPre10.02 driver \
-r Homo_sapiens_assembly38.fasta \
--algo VarCal \
-v $path/Haplotyper/<sample>.Haplotyper.vcf.gz \
--resource hapmap_3.3.hg38.vcf.gz \
--resource_param hapmap,known=false,training=true,truth=true,prior=15.0 \
--resource 1000G_omni2.5.hg38.vcf.gz \
--resource_param omni,known=false,training=true,truth=true,prior=12.0 \
--resource 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource_param 1000G,known=false,training=true,truth=false,prior=10.0 \
--resource dbsnp_b151_GRCh38p7_GATK.vcf.gz \
--resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 --annotation DP \
--annotation QD --annotation FS --annotation SOR --annotation MQ --annotation MQRankSum \
--annotation ReadPosRankSum \
--var_type SNP \
--plot_file $path/01-Germutation/<sample>.Haplotyper.snp.plotfile \
--tranches_file $path/01-Germutation/<sample>.Haplotyper.snp.tranches \
$path/01-Germutation/<sample>.Haplotyper.snp.recal && \




####applyVARcal SNP
SentieonPre10.02  driver \
-r Homo_sapiens_assembly38.fasta \
--algo ApplyVarCal -v $path/Haplotyper/<sample>.Haplotyper.vcf.gz \
--var_type SNP \
--tranches_file $path/01-Germutation/<sample>.Haplotyper.snp.tranches \
--sensitivity 99.0 \
--recal $path/01-Germutation/<sample>.Haplotyper.snp.recal \
$path/01-Germutation/01-snpvcf/<sample>.Haplotyper.snp.recaled.vcf.gz && \


Sentieon plot \
vqsr -o $path/01-Germutation/<sample>.vqsr.SNP.pdf \
$path/01-Germutation/<sample>.Haplotyper.snp.plotfile &&\

####VARcal INDEL

/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/SentieonPre10.02  driver \
-r Homo_sapiens_assembly38.fasta \
--algo VarCal \
-v $path/Haplotyper/<sample>.Haplotyper.vcf.gz \
--resource Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--resource_param Mills,known=false,training=true,truth=true,prior=12.0 \
--resource Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
--resource_param axiomPoly,known=false,training=true,truth=false,prior=10.0 \
--resource dbsnp_b151_GRCh38p7_GATK.vcf.gz \
--resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \
--annotation QD --annotation DP --annotation FS \
--annotation SOR --annotation MQ --annotation ReadPosRankSum \
--var_type INDEL \
--max_gaussians 4 \
--plot_file $path/01-Germutation/<sample>.Haplotyper.indel.plotfile \
--tranches_file $path/01-Germutation/<sample>.Haplotyper.indel.tranches \
$path/01-Germutation/<sample>.Haplotyper.indel.recal && \


####applyVARcal INDEL
SentieonPre10.02  driver \
-r Homo_sapiens_assembly38.fasta \
--algo ApplyVarCal -v $path/Haplotyper/<sample>.Haplotyper.vcf.gz \
--var_type INDEL \
--tranches_file $path/01-Germutation/<sample>.Haplotyper.indel.tranches \
--sensitivity 99.0 \
--recal $path/01-Germutation/<sample>.Haplotyper.indel.recal \
$path/01-Germutation/02-indel/<sample>.Haplotyper.indel.recaled.vcf.gz &&\


##
Sentieon/SentieonPre10.02  plot \
vqsr -o $path/01-Germutation/<sample>.vqsr.indel.pdf \
$path/01-Germutation/<sample>.Haplotyper.indel.plotfile &&\

#### merge
GATK/gatk-4.1.0.0/gatk \
MergeVcfs \
-I $path/01-Germutation/02-indel/<sample>.Haplotyper.indel.recaled.vcf.gz  \
-I $path/01-Germutation/01-snpvcf/<sample>.Haplotyper.snp.recaled.vcf.gz  \
-O $path/01-Germutation/03-merge/<sample>.Haplotyper.snp_indel.recaled.vcf.gz &&\

#### Filter

zcat $path/01-Germutation/03-merge/<sample>.Haplotyper.snp_indel.recaled.vcf.gz\
|awk '$1~/^#/ || $7=="PASS" {print}' |\
bcftools-1.9 norm \
-m-both \
-f Homo_sapiens_assembly38.fasta \
-o $path/01-Germutation/04-passvcf/<sample>.vqsr_pass.vcf.gz -O z - &&\

#####VEP annotation




export PATH="/share/app/perl-5.22.0/bin:$PATH"
export PERL5LIB="/share/app/vcftools-0.1.13/lib/perl5/site_perl:$PERL5LIB"
VEP/v102/bin/vep \
--cache --dir_cache VEP/Caches \
--fasta Homo_sapiens_assembly38.fasta \
--assembly GRCh38 --offline --everything --af_exac \
--format vcf --vcf \
--input_file $path/test.vqsr_pass.vcf.gz \
--output_file $path/test.vqsr_pass.vep.vcf \
--stats_file $path/test.vqsr_pass._summary
~~~