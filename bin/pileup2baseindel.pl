#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Fri 03 Aug 2012 10:34:54 PM CDT  
#Vanderbilt Center for Quantitative Sciences
#############################################
#Modified by Matthieu Foll (follm@iarc.fr) and Tiffany Delhomme 2015
#from https://github.com/riverlee/pileup2base (see https://github.com/mfoll/pileup2base)
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::File;

######### Update from pile2base.pl ####################
# 1).it will parse the insert and deletion too (and will not if option -no-indels), 
# 2).it can parse pileups from multiple sample,
# 3).you don't need to define the outputfile,instead, 
#    you will need to provide 'prefix' which will name result from each sample
#    as prefix1.txt, prefix2.txt,etc. Default, prefix is sample.
# 4).Provide option to define the base quality score offset, default is 33 (Sanger standard).
# 5).Read parameters from command line with options

#Usage: perl pileup2baseindel.pl -i <pileupfile> -bq [BQcutoff] -prefix [sample] -offset [33] (-no-indels)
my $usage = <<USAGE;
Usage: perl pileup2base.pl -i <pileupfile> -bq [BQcutoff] -prefix [sample] -offset [33]
        -i        input pileup file, could be from 1 sample or multiple samples
        -bq       base quality score cutoff for each mapped/unmapped base,
                  only those larger than cutoff will be output in the result, default is -5, means no filter
        -prefix   output file prefix, default is sample, the output will be named as prefix1.txt, prefix2.txt, etc.
        -offset   Offset to change ASCII character to base quality score, default is 33 (sanger format).
        -h        print out this
USAGE

my ($input,$BQcut,$offset,$prefix,$help,$noindel) = (undef,-5,33,"sample",undef,undef);
GetOptions(
	"i=s"=>\$input,
	"bq=i"=>\$BQcut,
	"offset=i"=>\$offset,
	"prefix=s"=>\$prefix,
	"h"=>\$help,
	"no-indels"=>\$noindel
);
if($help){
	print $usage;
	exit(0);
}
unless ($input){
	print "Input file does not provide yet\n";
    print "\n$usage\n";
    exit(1);
}
if(! -e $input){
    print "Input file '$input' does not exists\n";
    print "\n$usage\n";
    exit(1);
}

#Do the parsing
open FILE, $input or die "error, can not open $input";
print "[",scalar(localtime),"] Begin parsing...\n";

my $line = <FILE>;
$line=~s/\r|\n//g;
my ($chr,$loc,$ref,@dp_bases_bq) = split /\s+/, $line;
my $n = int (scalar(@dp_bases_bq)/3); #determine how many samples, use int just in safe
my %files;
foreach my $i (1..$n){
	my $fh = new IO::File;
	$fh->open("> ${prefix}${i}.txt");
	print $fh "chr\t"."loc\t"."ref\t"."A\t"."T\t"."C\t"."G\t"."a\t"."t\t"."c\t"."g\t"."Insertion\t"."Deletion\n";
	
	my @region=(3*($i-1),3*($i-1)+1,3*($i-1)+2);
	my ($dp,$bases,$bq) = @dp_bases_bq[@region];
	my $str = parsePileup($bases,$bq,$BQcut,$offset,$ref);
	
	if($str ne "*"){
		print $fh join "\t",($chr,$loc,$ref,$str);
	}
	$files{$i}=$fh;
}

while(<FILE>){
	s/\r|\n//g;
	my ($chr,$loc,$ref,@dp_bases_bq) = split /\s+/;
	my $n = int (scalar(@dp_bases_bq)/3); #determine how many samples, use int just in safe
	foreach my $i (1..$n){
		my $fh = $files{$i};
		my @region=(3*($i-1),3*($i-1)+1,3*($i-1)+2);
		my ($dp,$bases,$bq) = @dp_bases_bq[@region];
		my $str = parsePileup($bases,$bq,$BQcut,$offset,$ref);
		#if($str ne "*"){
			print $fh join "\t",($chr,$loc,$ref,$str);
		#}
	}
}

close FILE;
foreach my $k (keys %files){
	$files{$k}->close;
}
print "[",scalar(localtime),"] Finished\n";


sub parsePileup{
	my ($bases,$bq,$BQcut,$offset,$ref) = @_;
	
	if($bases eq "*"){
		my $str="0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\n";
		return $str;
		#return "*";
	}
	
	#do some modificaton on $base to remove additional characters
	#1,remove the ^. pattern
	$bases=~s/\^.//g;
	#2,remove the $ pattern
	$bases=~s/\$//g;
	#3,remove -[0-9]+[ACGTNacgtn]+ pattern
	my %hash=();
	my %deletion=();
	my %insertion=();
	while($bases=~/-(\d+)/g){
		$hash{$1}=1;
	}

	#if !noindels get the deletion and insertion sequences and delete them
	foreach my $k (keys %hash){
		if(!$noindel){
			while($bases=~/-$k([ACGTNacgtn]{$k})/g){
	            		$deletion{$1}++;
	        	}
		}
	        $bases=~s/-$k[ACGTNacgtn]{$k}//g;
	}
	%hash=();
	while($bases=~/\+(\d+)/g){
		$hash{$1}=1;
	}
	foreach my $k (keys %hash){
		if(!$noindel){
		        while($bases=~/\+$k([ACGTNacgtn]{$k})/g){
			    $insertion{$1}++;
		        }
		}
		$bases=~s/\+$k[ACGTNacgtn]{$k}//g;
	}

	#Now @base and @bq have the same length
	my @base=split (//,$bases);
	my @bq=split(//,$bq);
	#I have check it
	#if(scalar(@base) ne scalar(@bq)){
	#	print $_,"\n";
	#}
	#foreach my $c (@base){
	#	$check{$c}++;
	#}
	my $forward_A=0;
	my $forward_T=0;
	my $forward_C=0;
	my $forward_G=0;
	my $reverse_A=0;
	my $reverse_T=0;
	my $reverse_C=0;
	my $reverse_G=0;

	#start the loop
	for(my $i=0;$i<@base;$i++){
		my $ch=$base[$i];
		my $score=ord($bq[$i])-$offset;
		if($score>=$BQcut){
		    if($ch eq "A"){
		        $forward_A++;
		    }elsif($ch eq "T"){
		        $forward_T++;
		    }elsif($ch eq "C"){
		        $forward_C++;
		    }elsif($ch eq "G"){
		        $forward_G++;
		    }elsif($ch eq "a"){
		        $reverse_A++;
		    }elsif($ch eq "t"){
		        $reverse_T++;
		    }elsif($ch eq "c"){
		        $reverse_C++;
		    }elsif($ch eq "g"){
		        $reverse_G++;
		    }elsif($ch eq "."){
		        if($ref eq "A"){
		            $forward_A++;
		        }elsif($ref eq "T"){
		            $forward_T++;
		        }elsif($ref eq "C"){
		            $forward_C++;
		        }elsif($ref eq "G"){
		            $forward_G++;
		        }
		    }elsif($ch eq ","){
		        if($ref eq "A"){
		            $reverse_A++;
		        }elsif($ref eq "T"){
		            $reverse_T++
		        }elsif($ref eq "C"){
		            $reverse_C++;
		        }elsif($ref eq "G"){
		            $reverse_G++;
		        }
		    }
		}#end the condition  $score>=$BQcut
	}#end the loop
	
    #my $str="$chr\t$loc"."\t".$ref."\t".$forward_A."\t".$forward_T."\t".$forward_C."\t".$forward_G."\t".$reverse_A."\t".$reverse_T."\t".$reverse_C."\t".$reverse_G."\t";
	my $str=$forward_A."\t".$forward_T."\t".$forward_C."\t".$forward_G."\t".$reverse_A."\t".$reverse_T."\t".$reverse_C."\t".$reverse_G."\t";
	my $insertion="NA";
	my $deletion="NA";
	if(!$noindel){
		if(scalar(keys %insertion)){
			$insertion="";
			foreach my $k (sort {$insertion{$b}<=>$insertion{$a}} keys %insertion){
				$insertion.=$insertion{$k}.":".$k."|";
			}
			chop($insertion);
		}
	
		if(scalar(keys %deletion)){
			$deletion="";
			foreach my $k (sort {$deletion{$b}<=>$deletion{$a}} keys %deletion){
				$deletion.=$deletion{$k}.":".$k."|";
			}
			chop($deletion);
		}
	}
	$str.=$insertion."\t".$deletion."\n";
	return $str;
}
