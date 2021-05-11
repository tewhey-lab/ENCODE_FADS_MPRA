#!/usr/bin/perl

use strict;
use warnings;

my $fastq_r1_in = $ARGV[0];

open (FASTQ_R1, "$fastq_r1_in") or die("ERROR: can not read file ($fastq_r1_in): $!\n");

my $r1;
my $tmp;
my @id;

while (<FASTQ_R1>)
{
	$r1 = $_;
	$r1 =~ s/^@/>/;
	chomp $r1;
	@id=split(/\s/,$r1);
	print "$id[0]\n";
	#print $r1;
	$r1 = <FASTQ_R1>;
	print $r1;
	$tmp = <FASTQ_R1>;
	$tmp = <FASTQ_R1>;
}
close FASTQ_R1;
