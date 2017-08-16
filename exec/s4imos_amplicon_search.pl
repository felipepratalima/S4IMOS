#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::Tools::AmpliconSearch;
use Bio::SeqIO;
use Bio::PrimarySeq;

## Helper messages
my $DIE_MESSAGE = "Incorrect arguments... USAGE: s4imos_amplicon_search.pl --fasta FASTA_FILE --primers PRIMERS_FILE > OUTPUT_FILE\n";

## Get command line options
GetOptions(\my %options,
  "fasta=s",
  "primers=s",
)
or die ($DIE_MESSAGE);

## Check options
die ($DIE_MESSAGE) unless exists $options{fasta};
die ($DIE_MESSAGE) unless exists $options{primers};

## Get file names
my $fasta = $options{fasta};
my $primers = $options{primers};

## Test if files exist
die ("Fasta file not found: $fasta\n") unless -e $fasta;
die ("Primers file not found: $primers\n") unless -e $primers;

## Object to access sequences in the fasta
my $sequences = Bio::SeqIO->new(
  -file => $fasta,
  "-format" => "Fasta"
) or die ("Fasta file \n");

## For each sequence in the fasta
while(my $sequence = $sequences->next_seq) {
  ## Search the amplicon
  my $genome_sequence = Bio::PrimarySeq->new(
    -seq => $sequence->seq,
  );
  my $amplicon_search = Bio::Tools::AmpliconSearch->new(
    -template => $genome_sequence,
    -primer_file => $primers,
  );
  
  ## If amplicon(s) found, print to output
  if (my $amplicon = $amplicon_search->next_amplicon) {
    print $sequence->id . "\n";
  }
}
