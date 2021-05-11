#!/usr/bin/perl

use strict;
use warnings;

my $fastq = $ARGV[0];

open (FASTQ, "$fastq") or die("ERROR: can not read file ($fastq): $!\n");

my $line;
my $revcomp;
while (<FASTQ>)
{
print "$_";

$line = <FASTQ>;
chomp $line;
$revcomp = reverse($line);
$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
print $revcomp."\n";

$line = <FASTQ>;
print "$line";

$line = <FASTQ>;
chomp $line;
$revcomp = reverse($line);
print $revcomp."\n";
}
close FASTQ;
