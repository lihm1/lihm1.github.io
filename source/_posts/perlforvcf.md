---
title: perl处理vcf文件3
date: 2022-04-18 11:17:25
index_img: /img/pic3.jpg  
categories:
- perl
tags:
- study
---
### 方法一
conding 区域的bed文件，按照染色体查找在区间内的区域
~~~
#!user/bin/perl


use strict;
use warnings;
use Data::Dumper;

my $path = $ARGV[0];
my $varFile = $ARGV[1];
my $outdir = $ARGV[2];
my $bed = $ARGV[3];

open FH, "gzip -dc $path/$varFile|" || die $!;
      
 my @filname = (split /\./,$varFile,2);
 my $samplename = $filname[0];
 if ($samplename){
                  `mkdir -p $outdir/$samplename`;
                        }else{
          print "please check your file";}
                                                                                                                                                         

 open CODE,">$outdir/$samplename/$samplename.vcf" or die $!;

 while (<FH>){
             chomp();                                                                                             if (/##fileformat/ || /##FORMAT=<ID=GT/ || /##FORMAT=<ID=AD/ || /##FORMAT=<ID=DP/ || /##VEP/ || /##INFO=<ID=CSQ/ || /#CHROM/ ){
   print CODE "$_\n";}
        }

open IH, "gzip -dc $path/$varFile |" || die $!;   

open BED,"$bed" || die $!;
                                                                                                             my (%hashBed,%hashVar);
                                                                                                             while (<BED>){
                                                                                                             chomp;
my @seq = split /\s+/, $_;
 my ($chr, $start,$end) = @seq[0, 1, 2];
 $hashBed{$chr}{$start}=$end;
 }


while (<IH>){
    chomp;
if (/^#/){
    next;}else{
    my @seq = split /\s+/, $_,3;
    my ($chr, $start,$other) = @seq[0, 1, 2];
   $hashVar{$chr}{$start}=$other;
}
}

my @chrorder =qw (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chrX chrY);


my (@test);
foreach my $chr(@chrorder){
    if(exists $hashBed{$chr}){
        print "$chr\n";
        while (my ($start1) = each $hashVar{$chr}){
            print "$start1\n";
            while ( my ($start2,$end) = each $hashBed{$chr}){ 
            if ($start1 >= $start2 && $start1<=$hashBed{$chr}{$start2}){
                 my $mutation = $hashVar{$chr}{$start1};
                 print CODE join ("\t",$chr, $start1,$mutation ), "\n";
            }else{
                print "no";}
         
}}}}



#print Dumper (\%hashVar);

#print Dumper (\%hashBed);
~~~
### 方法二
bcftools
调整 bed文件的格式 9:4700000-4800000，一次读入进bcftools

~~~
bcftools filter 1000Genomes.vcf.gz --regions 9:4700000-4800000 > 4700000-4800000.vcf
~~~

