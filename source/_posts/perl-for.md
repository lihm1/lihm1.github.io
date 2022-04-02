---
title: perl_for
date: 2022-03-30 10:58:22
index_img: /img/pic3.jpg  
categories:
- perl
tags:
- study
---
### perl 循环用法
#### while
~~~
my $i = 0;  
while($i < scalar @array) {  
    print $i, ": ", $array[$i];  
    $i++;  
} 
~~~
举例：
~~~
@tt = 5..9 ;
my $i = 0;  
while($i < scalar @tt) {  
    print $i, " : ", "$tt[$i]\n";  
    $i++;  
} 
~~~
#### until
~~~
my $i = 0;  
until($i >= scalar @array) {  
    print $i, ": ", $array[$i];  
    $i++;  
} 
#!/usr/bin/perl
@tt = 5..9 ;

my $i = 0;  
until($i >= scalar @tt) {  
    print $i, " : ", "$tt[$i]\n";  
    $i++;  
} 
~~~
#### Do
将until 或者while 加到后面
~~~
#!/usr/bin/perl
@tt = 5..9 ;

my $i = 0;  
do {  
    print $i, " : ", "$tt[$i]\n";  
    $i++;  
}until ($i >= scalar @tt)
~~~

#### 迭代哈希键值
你不能去迭代一个哈希变量。但是，你可以迭代它的键值。使用keys内置函数，获取哈希变量的所有键值数组。然后使用foreach方法，就像数组一样：

foreach my $key (keys %scientists) {  
    print $key, ": ", $scientists{$key};  
} 

因为哈希变量没有顺序，键值可能以任何顺序被返回。使用sort内置函数对键值数组排序，按照字母表从小到大的方式：
foreach my $key (sort keys %scientists) {  
    print $key, ": ", $scientists{$key};  
} 


#### 跳出循环
Next last 允许在循环中改变程序
next跳至本次循环的结束，跳出迭代
last 退出循环语句块，从而结束循环，last语句之后的语句不再执行，continue语句块也不再执行


#### for 循环
for循环不仅仅可以支持数字递增递减的循环方式，还可以支持其他类型的循环，只要能进行判断就可。
数字的递增递减
for (my $i=1;$i<=10;$i++){
print $i,"/n";}

for括号里面的三个表达式都可以省略，但是两个分号不能省略。
~~~
for(my $str="malongshuai";$str =~ s/(.)//;){
    print $str,"\n";
}
~~~

像shell的for i in ..的操作也可以使用，遍历整个列表；
~~~
my @arr=qw(shel cuabc bjab)
for $i (@arr){
print "$i/n";
}
for @arr{
print "$_/n";
}
~~~

#### perl 控制语句
#### perl 三目运算符 ？

如果表达式为真，表达式返回if_true，否则返回if_false,
expression ? if_true:if_false
~~~
1.$avg = $n ? $sum/$n : "------";
2.if($n){
    $avg = $sum / $n;
}else{
    $avg = "------";
}
~~~
上述代码表达内容相等
三目还可以更复杂：
~~~
use 5.010;

$score = $ARGV[0];
$mark = ($score < 60) ? "a" : 
        ($score < 80) ? "b" : 
        ($score < 90) ? "c" : 
        "d";          # 默认值
say $mark;
~~~
输出的结果：
~~~
$ perl test.plx
a
$ perl test.plx 33
a
$ perl test.plx 63
b
$ perl test.plx 83
c
$ perl test.plx 93
d
~~~
#### and 操作符：
用于连接两个行为，左边为真，就执行右边的操作。
$m < $n and $m = $n;   # 以$m的姿态取出$m和$n之间较大值
下面代码所表达的意思是一样的。
~~~
if ($m<$n){$m=$n}
$m=$n if $m<$n;
$m=($m < $n) ? $n : $m;
~~~
