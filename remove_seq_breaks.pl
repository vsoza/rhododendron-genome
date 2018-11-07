#!/usr/bin/perl
#written by Adam Waalkes to remove line breaks (newline characters), which occur every 80bp, in genome fasta file

use warnings;
use strict;

my $vcf_file = $ARGV[0];
my $first_line="";

open(IP, "$vcf_file");

$first_line=<IP>;
print $first_line;

while(my $line = <IP>) {
   if((substr($line, 0,1)) eq ">") {
      print "\n$line";
   }
   else {
      chomp $line;
      print $line;
   }
}
print "\n";
close IP;
exit;
