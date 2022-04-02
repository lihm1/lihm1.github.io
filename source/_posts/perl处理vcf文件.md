---
title: perl处理vcf文件
date: 2022-04-02 14:37:03
index_img: /img/pic3.jpg  
categories:
- perl
tags:
- study
---
最近尝试着用perl处理vcf文件得到自己想要的内容，在这里做个总结。

### 1.向vcf文件添加一列：
提取vep注释后的vcf中的样本名并将样本名添加到文件的每条染色体后面
解决思路：
读入文件，join ;push ;
~~~
#!user/bin/perl

use strict;
my $bcftools = "bcftools-1.9";
my $varFile = "vqsr_pass.vep.vcf";

my(@name,$samplename);
    open IH, "$bcftools view -h $varFile |" || die $!;
    open VCF, ">vctout";
    while(<IH>){
        chomp();
        if (/^#CHROM/){
            @name = split /\s+/,$_;
            $samplename = @name[9]
        }
    }
print $samplename;

open FA, "$bcftools view -H $varFile |" || die $!;
while (<FA>){
    chomp();
    if (/^#/){
        next;
    }else{
        my @seq = split /\s+/, $_;
        my ($chr, $pos, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @seq[0, 1, 2, 3, 4, 5, 6,7, 8];
        my $FORMAT = join (",",@seq[8],@seq[9]);
        print VCF join ("\t",$chr, $pos, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT,$samplename), "\n";    }
}  
close FA;
close VCF;
close IH;
~~~

### 2.读取vcf文件
保留#CHROM 这一行
~~~
#!user/bin/perl

use strict;

open VCF,"<$ARGV[0]";
open OUT,">mmary.txt";

my @linelist=<VCF>;
foreach my $eachline(@linelist){
if ($eachline=~/^##/){
next;}
elsif($eachline=~/^#CHROM/){
  $eachline =~ s/#//;
  print OUT $eachline;
}else{
print OUT $eachline;
}
}
~~~
### 3.将多个vcf文件合并为一个并且按照染色体进行排序
排序的解决思路：哈希表的排序
~~~
#!user/bin/perl

use warnings;
use Data::Dumper;

my $bcftools = "bcftools-1.9";

my $indir = $ARGV[0];

opendir  DIR, $indir or die "Connot open $indir: $!";

@allfile = grep { /^[^\.]/ } readdir DIR;


open mergeVCF, ">mergevcf" || die $!;

my (%hashVar);

foreach my $file (@allfile){
open FA,"$indir/$file"|| die $!;
while(<FA>){
chomp();
my @seq = split /\s+/, $_;
my ($chr, $pos, $ID, $REF, $ALT, $QUAL, $FILTER) = @seq[0, 1, 2, 3, 4, 5, 6];
$hashVar{$chr}{$pos}{$ID}{$REF}{$ALT}{$QUAL}{$FILTER}="sample";
}
}


my @chrorder =qw (chr1 chr2);

 
for  my $chr (@chrorder){
next unless (exists $hashVar{$chr});
for  my $pos (sort {$a<=>$b} keys %{$hashVar{$chr}}){
print "$pos\n";
for  my $ID (keys %{$hashVar{$chr}{$pos}}){
for  my $REF (keys %{$hashVar{$chr}{$pos}{$ID}}){
for  my $ALT (keys %{$hashVar{$chr}{$pos}{$ID}{$REF}}){
for  my $QUAL (keys %{$hashVar{$chr}{$pos}{$ID}{$REF}{$ALT}}){
for  my $FILTER (keys %{$hashVar{$chr}{$pos}{$ID}{$REF}{$ALT}{$QUAL}}){
print mergeVCF join ("\t",$chr, $pos, $ID, $REF, $ALT,$QUAL, $FILTER), "\n";

}}}}}}}


print Dumper (\%hashVar);
close mergeVCF;
~~~

### 4.按照染色体提取分别保存在不同的文件里

~~~use strict;


my $bcftools = "bcftools-1.9";
my $varFile = "15TSLC.vqsr_pass.vep.vcf.gz ";

open IH, "$bcftools view -H $varFile |" || die $!;   
open other,">other.vcf" || die $!;

my (%normal_chr,%fh);
my @chrorder = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y);

foreach my $chr (@chrorder){
    $normal_chr{"chr$chr"}=1 ; 
    open  $fh{"chr".$chr},">chr$chr.txt" or die $!;
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
 ~~~


