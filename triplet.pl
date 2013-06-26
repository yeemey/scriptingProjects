#!/usr/bin/perl
my $filename = $ARGV[0];

chomp $filename;

open MYFILE1, $filename || die "Sequence file not found!\n";
<MYFILE1>;
my @SEQFILE1 = <MYFILE1>;
close MYFILE1;
#Open first sequence FASTA file and read into array.

my $oneseq = join(' ', @SEQFILE1);
$oneseq =~s/\s//g;
print "\nYour sequence is $oneseq\n\n";
#Print sequence to screen

for (my $position = 3; $position < length($oneseq); $position=$position+3){
 my $codon = substr($oneseq, $position, 3);
 my $realpos = $position+1;
 if ($codon eq "TCT"){
  print "Alternative start codon TCT at position $realpos\n\n";
 }
 elsif ($codon eq "GTG"){
   print "Alternative start codon GTG at position $realpos\n\n";
 }
  elsif ($codon eq "ATA"){
   print "Alternative start codon ATA at position $realpos\n\n";
 }
 elsif ($codon eq "ATG"){
   print "Alternative start codon ATG at position $realpos\n\n";
 }
 else {next;}
}


#my @SEQTRIPLET = unpack('a3' x (length($oneseq)-2), $oneseq);
#print "You have split your sequence into the following triplets:\n";
#print "@SEQTRIPLET \n\n";
#Split sequence into triplets
