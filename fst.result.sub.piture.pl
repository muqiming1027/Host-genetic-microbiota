#!/bin/perl -w
use strict;
die "usage : perl $0 <input fst result>\n" unless @ARGV == 1;
open IN,"$ARGV[0]" or die "fst result file not found !\n";
my %hash;
open INDEX,"/home/qmmou/code/index_chr" or die "can not open chr index file !\n";
open OUT,">","$ARGV[0]\.picture" or die "permission denied !\n";
while (<INDEX>){
	chomp;
	my @a = split;
	$hash{$a[0]} = $a[1];
}
while (<IN>){
	next if /CHR/;
	my $line =$_ ;
	foreach (keys %hash){
		if ($line =~ /$_/){
			$line =~ s/$_/$hash{$_}/g;
			last;
		}
	}
	chomp $line ;
	my @a = split /\t/,$line;
	print OUT "$a[0]\t$a[0]\t$a[1]\t$a[-1]\n";
}
