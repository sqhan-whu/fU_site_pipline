#!/usr/bin/perl


#  Author     : Han Shaoqing
#  Email      ：sqhan1@whu.edu.cn
#  Last Edited: 20201021


use strict;
use warnings;
use Getopt::Long;
my %opt;


GetOptions (\%opt,"bg:s","in:s","help");

my $help=<<USAGE;
This script could correct TFEA mRNA by background mRNA
Usage: perl $0 -bg bg.tsv -in in.tsv
USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help\n";
		exit();
		}

		my $in1 = $opt{bg};
		my $in2 = $opt{in};

		my %hash;

		open (IN1,$in1) or die $!;
		while (my $line=<IN1>) {
			chomp $line;
				my @array = split "\t",$line;
				        my $loc = join "_",($array[2],$array[6]);
				                $hash{$loc} = "exist";
				                }
				                close IN1;


				                open (IN2,$in2) or die $!;
				                open (OUT,">${in2}_corrected.tsv") or die $!;
				                open (OUT2,">${in2}_discard.tsv") or die $!;
				                while (my $line=<IN2>) {
				                       chomp $line;
				                              my @array = split "\t",$line;
				                                     my $loc = join "_",($array[2],$array[6]);
				                                            if (! $hash{$loc}) {
				                                                      print OUT $line,"\n";
				                                                             } else {
				                                                                       print OUT2 $line,"\n";
				                                                                              }
				                                                                              }
				                                                                              close OUT;
				                                                                              close OUT2;
				                                                                              close IN2;


