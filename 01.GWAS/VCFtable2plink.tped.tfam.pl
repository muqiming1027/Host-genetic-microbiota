#!/usr/bin/perl -w 

# Read a vcf file and extract useful inforamtion into a new file
###################################
# write by Zhengkui Zhou 8-7-2012 #
###################################
use warnings;
use strict;

# Declare and initialize variables 
my @file1_data1 = (  ); 
my $file1_data1;
my @file1_data2 = (  ); 
my $file1_data2;

#INPUT: IGDB-hmp format; OUTPUT: PLINK tped & tfam
die  "Version 1.0\t2013-07-08;\nUsage: $0 <InPut><Ind-number><OutDir1:tped><OutDir2:tfam>\n" unless (@ARGV ==4);

open     INFile,"$ARGV[0]"  || die "input file can't open $!" ;
open     OutFile,">$ARGV[2]" || die "output file can't open $!" ;
open     OutFile2,">$ARGV[3]" || die "output file can't open $!" ;

@file1_data1 = <INFile>; 
my $query1 = @file1_data1;
my $i = 0;
my $SNP;
my $allele;
close  INFile;

my $m;
my $k;
for ($k=6; $k<($ARGV[1]+6); $k++)   # The number of individual $k < 303 +6 ;
{
    for ($m=0; $m<1; $m++)   # The individual name;
    {
    my @temp3 = split (/\s+/, $file1_data1[$m]);
    my $temp3_len = @temp3 ;
    my @temp4 = split (/\./, $temp3[$k]);
    print OutFile2 "$temp4[0] $temp4[0] 0 0 1 -9\n";     
    }
}
close (OutFile2) or die( "Cannot close file : $!");

my $n = 1;
for ($i=1; $i<$query1; $i++)   # care the start line number;
{
    my @temp1 = split (/\s+/, $file1_data1[$i]);
    my $temp1_len = @temp1 ;
                $temp1[0] =~ s/chr//;
                #PLINK tped format;
                #chrom rs genetics-pos physics-pos
                print OutFile "$temp1[0] snp$n 0 $temp1[1] ";
                $n ++;            
           
        my $m ;
        for ($m=0; $m<( $temp1_len-6 ); $m++)  ####
        {
            my @snp = split (/\s+/, $temp1[$m+6]);
            my $allele1=substr($snp[0],0,1);
            my $allele2=substr($snp[0],2,1);
                if ($snp[0] eq "./.")
                {
                    $SNP="0 0"   
                }
                else
                {
                    $SNP = "$allele1 $allele2"
                }
                print OutFile "$SNP ";
        }
       print OutFile "\n";
}     
close (OutFile) or die( "Cannot close file : $!");
