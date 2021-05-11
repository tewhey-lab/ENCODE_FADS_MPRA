#!/usr/bin/perl

use strict;
use warnings;
use Text::LevenshteinXS qw(distance);
use Getopt::Std;

my %options=();
getopts('B:', \%options);
 
my $fasta = $ARGV[0];
my $out = $ARGV[1];

open (FASTA, "$fasta") or die("ERROR: can not read file ($fasta): $!\n");
open (MATCH, ">$out".".match") or die("ERROR: can not create $out .matched: $!\n");
open (REJECT, ">$out".".reject") or die("ERROR: can not create $out .rejected: $!\n");

my $L2_SEQ_BARCODE_LEN = $options{B} || 20;
print STDERR "Barcode length: $L2_SEQ_BARCODE_LEN\n";


my $no_L1_match = 1;
my $no_L2_match = 1;

my $L1_SEQ = "AACTGGCCGCTTGACG";
my $L1_SEQ_ANCHOR = "CG";
my $L1_SEQ_ANCHOR_POS = 14;
my $L2_SEQ = "CACTGCGGCTCCTGCGATCGCGTCGACGAACCTCTAGA";  #Change substr length if this is not 26bp
my $L2_SEQ_ANCHOR = "GA";
my $L2_SEQ_LEFT_ANCHOR = "CA";
my $L2_SEQ_POS = 216;

my $MIN_SEQ_SIZE = 100;

my $MIN_ENH_SIZE = 50;
my $MAX_ENH_SIZE = 210;

my $r1;
my $id;
my $L1_read_seq;
my $L1_index;
my $L2_index;
my $L2_read_seq;
my $dist;

my $enh_seq;
my $barcode_seq;

my $match1 = -9;
my $match2 = -9;
my $L2_match_len = 0;
my $L2_match_start = 0;

my $match_dist = 3; #max levenshtein distance for match
my $anchor_wiggle = 2; #Distance in bp to move around anchor looking for a match

my $i;

while (<FASTA>)
{
	chomp;
	$id=$_;
	$id =~ s/^>//;
	$id =~ s/\/1$//;

	$r1=<FASTA>;
	chomp $r1;

	next if(length($r1) < $MIN_SEQ_SIZE);

##Left hand adapter match
	$match1 = -9;
	if(substr($r1,$L1_SEQ_ANCHOR_POS,length($L1_SEQ_ANCHOR)) eq $L1_SEQ_ANCHOR) #Anchor match
		{
		$match1 = 0 if($no_L1_match ==  1);		
		$match1 = 0 if(abs(distance($L1_SEQ, substr($r1,0,length($L1_SEQ)))) <= $match_dist);
		}
	else
		{
MATCH1:
		for($i=1;$i<=$anchor_wiggle;$i++)
			{
			if(substr($r1,($L1_SEQ_ANCHOR_POS-$i),length($L1_SEQ_ANCHOR)) eq $L1_SEQ_ANCHOR) #Anchor match
				{
				if(abs(distance($L1_SEQ, substr($r1,0,(length($L1_SEQ)-$i)))) <= $match_dist)
					{
					$match1 = -$i;
					last MATCH1;
					}
				}
			elsif(substr($r1,($L1_SEQ_ANCHOR_POS+$i),length($L1_SEQ_ANCHOR)) eq $L1_SEQ_ANCHOR) #Anchor match
				{
					if(abs(distance($L1_SEQ, substr($r1,0,(length($L1_SEQ)+$i)))) <= $match_dist)
					{
					$match1 = $i;
					last MATCH1;
					}	
				}
			elsif($no_L1_match == 1)
				{
					$match1 = 0.0;
					last MATCH1;				
				}
			}
		}
		
##Middle adapter match
	$match2 = -9;
	$L2_match_len = -9;
	$L2_match_start = -9;
	$dist = -9;
	if(substr($r1,-$L2_SEQ_BARCODE_LEN-length($L2_SEQ_ANCHOR),length($L2_SEQ_ANCHOR)) eq $L2_SEQ_ANCHOR) #Anchor match
		{
		$match2 = -99;
		$match2 = 0.0 if($no_L2_match == 1);		
MATCH2:
		for($i=0;$i<=$anchor_wiggle;$i++)
			{
			if(substr($r1,-($L2_SEQ_BARCODE_LEN+length($L2_SEQ)+$i),length($L2_SEQ_LEFT_ANCHOR)) eq $L2_SEQ_LEFT_ANCHOR) ##Left anchor match
				{
					$match2 = $i;
					last MATCH2;
				}
			if(substr($r1,-($L2_SEQ_BARCODE_LEN+length($L2_SEQ)-$i),length($L2_SEQ_LEFT_ANCHOR)) eq $L2_SEQ_LEFT_ANCHOR)  ##Left anchor match
				{
					$match2 = -$i;
					last MATCH2;
				}		
			}
		if($match2 > -9)
			{
			if(abs(distance($L2_SEQ, substr($r1,length($r1)-$L2_SEQ_BARCODE_LEN-length($L2_SEQ)-$match2,length($L2_SEQ)+$match2))) <= $match_dist)  ##Adapter match
				{
				$dist = distance($L2_SEQ, substr($r1,length($r1)-$L2_SEQ_BARCODE_LEN-length($L2_SEQ)-$match2,length($L2_SEQ)+$match2));
				#Matched internal adapter
				}
			elsif($no_L2_match == 1)
				{
				$dist = distance($L2_SEQ, substr($r1,length($r1)-$L2_SEQ_BARCODE_LEN-length($L2_SEQ)-$match2,length($L2_SEQ)+$match2));
				$match2 = 0.0;
				}	
			else
				{
				$match2 = -999
				}
			}
		}
	  elsif($no_L2_match == 1)
		{
			$match2 = 0.0;
		}
	
	
	
	if($match1 != -9 && $match2 > -9)
		{
		$L2_match_len = length($L2_SEQ)+$match2;
		$L2_match_start =  length($r1)-$L2_SEQ_BARCODE_LEN-length($L2_SEQ)-$match2;
		
		$enh_seq = substr($r1,length($L1_SEQ)+$match1,$L2_match_start-length($L1_SEQ)+$match1);
		$barcode_seq = substr($r1,-$L2_SEQ_BARCODE_LEN);

		$L1_read_seq = substr($r1,0,length($L1_SEQ)+$match1);
		$L2_read_seq = substr($r1,$L2_match_start,$L2_match_len);
		if(length($enh_seq) >= $MIN_ENH_SIZE && length($enh_seq) <= $MAX_ENH_SIZE)
			{
			print MATCH "$id\t$match1\t$match2\t$enh_seq\t$barcode_seq\t$L1_read_seq\t$L2_read_seq\n";
			}
		else
			{
			print REJECT join("\t",$id,$match1,$match2,substr($r1,0,length($r1)-$L2_SEQ_BARCODE_LEN),substr($r1,-$L2_SEQ_BARCODE_LEN),$dist,length($L2_SEQ),0)."\n";		
			}
		}
	else
		{
		print REJECT join("\t",$id,$match1,$match2,substr($r1,0,length($r1)-$L2_SEQ_BARCODE_LEN),substr($r1,-$L2_SEQ_BARCODE_LEN),$dist,length($L2_SEQ),1)."\n";
		}
}
close FASTA;