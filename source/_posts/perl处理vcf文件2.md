---
title: perl处理vcf文件2
date: 2022-04-08 17:18:53
index_img: /img/pic3.jpg  
categories:
- perl
tags:
- study
---
#### 方法一
~~~
vcftools
如何用vcftools从VCF文件中提取某条染色体信息
vcftools --gzvcf input.vcf --chr n --recode – recode-INFO-all --stdout | gzip -c > output.vcf.gz
说明：
–gzvcf：处理压缩格式的vcf文件（可替换为–vcf）
–chr n：选择染色体n，例：–chr 1
–recode：重新编码为vcf文件
–recode-INFO-all：将输出的文件保存所有INFO信息
–stdout：标准输出，后接管道命令
–gzip -c：压缩
> output.vcf.gz：将结果输出到output.vcf.gz
~~~
#### 方法二
按照染色体分离并保留header，且生成文件夹
~~~
use strict; 
my $bcftools = "/ldfssz1/ST_CANCER/POL/SHARE/CancerPipeline/Sentieon_V1/Bin_CentOS6/bcftools-1.9";

my $path = $ARGV[0];
my $varFile = $ARGV[1];
my $outdir = $ARGV[2];

#my $path = /hwfssz5/ST_HEALTH/P18Z10200N0098/CGR/PMO/CRC/Analysis2022/02.Germline/z6_Haplotyper/01-Germutation/05-vep/;
##my $varFile = 15TSLC0144TR3.vqsr_pass.vep.vcf.gz ;
#### $ARGV[0] : input the sample file path;
#### $ARGV[1] : input the sample name;


open FH, "$bcftools view -h $path/$varFile |" || die $!;   
 my @filname = (split /\./,$varFile,2);
  my $samplename = @filname[0];
   if ($samplename){
        `mkdir -p $outdir/$samplename`;
         }else{
                  print "please check your file";
              }
my (%normal_chr,%fh);
my @chrorder = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y);
foreach my $chr (@chrorder){
        open  $fh{"chr".$chr},">$outdir/$samplename/chr$chr.txt" or die $!;}


while (<FH>){
         chomp();
             if (/##fileformat/ || /##FORMAT=<ID=GT/ || /##FORMAT=<ID=AD/ || /##FORMAT=<ID=DP/ || /##VEP/ || /##INFO=<ID=CSQ/ || /#CHROM/ ){
                         foreach my $chr (@chrorder){ 
                                     $fh{"chr".$chr} -> print("$_\n");
                                 }
                             }
                         }
open IH, "$bcftools view -H $path/$varFile |" || die $!;   
open other,">$outdir/$samplename/other.vcf" || die $!;


foreach my $chr (@chrorder){
        $normal_chr{"chr$chr"}=1 ; 
        }
while(<IH>){
        chomp();
            my  @seq = split /\s+/, $_ ;
            if (exists $normal_chr{$seq[0]}){
                print $seq[0];
                $fh{$seq[0]}->print( "$_\n");
            }else{
                print other $_;}
        }
foreach my $chr (@chrorder){
    close $fh{"chr".$chr};
    }
close other;

close IH;

close FH;
~~~
~~~
use strict; 
my $path = $ARGV[0];
my $varFile = $ARGV[1];
my $outdir = $ARGV[2];

#my $path = /hwfssz5/ST_HEALTH/P18Z10200N0098/CGR/PMO/CRC/Analysis2022/02.Germline/z6_Haplotyper/01-Germutation/05-vep/;
###my $varFile =16TSRC0715TR3.vqsr_pass.vep.vcf.gz

open FH, "gzip -dc $path/$varFile|" || die $!;
      
 my @filname = (split /\./,$varFile,2);
   my $samplename = @filname[0];
      if ($samplename){
                  `mkdir -p $outdir/$samplename`;
                           }else{
                               print "please check your file";}
                                                               
my (%normal_chr,%fh);
my @chrorder = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y);
foreach my $chr (@chrorder){
open  $fh{"chr".$chr},">$outdir/$samplename/chr$chr.txt" or die $!;}
while (<FH>){
    chomp();
    if (/##fileformat/ || /##FORMAT=<ID=GT/ || /##FORMAT=<ID=AD/ || /##FORMAT=<ID=DP/ || /##VEP/ || /##INFO=<ID=CSQ/ || /#CHROM/ ){
        foreach my $chr (@chrorder){
            $fh{"chr".$chr} -> print("$_\n");
        }}
        }
 open IH, "gzip -dc $path/$varFile |" || die $!;   
 open other,">$outdir/$samplename/other.vcf" || die $!;
 foreach my $chr (@chrorder){
     $normal_chr{"chr$chr"}=1 ; }
 while(<IH>){
     chomp();
     my @seq = split /\s+/, $_ ;
      if (exists $normal_chr{$seq[0]}){
         print $seq[0];
         $fh{$seq[0]}->print( "$_\n");}else{
          print other $_;}
  }
  
  
  foreach my $chr (@chrorder){
      close $fh{"chr".$chr};
    }

    
    close other;
    close IH;
    close FH;
~~~