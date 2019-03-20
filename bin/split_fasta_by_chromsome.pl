#!/usr/bin/perl
#usage: perl split_fasta_by_chrom.pl <genome.fa> <outDIR>
$f = $ARGV[0]; #get the file name
$o = $ARGV[1]; #get the DIR name with "/"

open (INFILE, "<$f") or die "Can't open: $f $!";

while (<INFILE>) {
$line = $_;
chomp $line;
if ($line =~ /\>/) { #if has fasta >
	#split (/\n/,$_)
	@tmp = split (/\s/, substr($line,1)) ;
	$new_file = $tmp[0];
	#print $new_file."\n";
	$new_file .= ".fa";
	open (OUTFILE, ">$o/$new_file") or die "Can't open: $new_file $!";
	}
	print OUTFILE "$line\n";
}
	close OUTFILE;