#!/usr/bin/perl -w

# Script use to collect each miRNA information from gtf annotation file such as start_position, end_position, strand, chromosome_name and gene_id

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use List::Util 'first';

# Define the different input/output files
my $fasta; # Input fasta file containing miRNA sequence
my $gff;    # Input gff file
my $output; # Output file which will contain all the different genes info

# Define the parameter in order to submit input files to this script
&GetOptions (
    'fasta=s' => \$fasta,
    'gff=s' => \$gff,
    'output=s' => \$output,
);

my $start_date = localtime;
print STDERR "\n################################\nSTART = $start_date\n################################\n\n";

# Open the fasta input file
unless ($fasta) {
    die "Please specify the fasta file containing the miRNA sequence via -fasta parameter!\n";
}
open (FASTA, "<$fasta") || die "Cannot open $fasta: $!\n"; $_="1";

# Open the gtf input file
unless ($gff) {
    die "Please specify the gff file containing the gene annotation via -gff parameter!\n";
}
open (GFF, "<$gff") || die "Cannot open $gff: $!\n"; $_="1";

# Open the output file
unless ($output) {
    die "Please specify the output file via -output parameter!\n";
}
if (-e $output) {
    die "This file: $output already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Define variables required for reading the fasta input file
my %mirna_seq;   # Hash containing each sequence of mature miRNA
my $total = 0;  # Variable to count line in input file
my $ig_count = 0;  # Variable to count line ignored from input file
my $keep_count = 0;  # Variable to count line kept from input file

# Read in the fasta sequence file
while(1) {
    my $line;        # Scalar containing line from the fasta file
    chomp($line = <FASTA>);    # Read lines one by one from the fasta file
    if ($line =~ /^>bta-/){   # Try to match the value to collect
        my ($name, $id, $full_name) = (split(/\s/, $line));
        $name =~ s/^>//;
        my $sequence;
        chomp($sequence = <FASTA>);    # Read following line from the fasta file which contains associated sequence
        unless (exists $mirna_seq{$name}) {
            $mirna_seq{$name} = $sequence;
            $total += 2;
            $keep_count += 2;
        }
        else {
            die "The miRNA name $name already has a sequence value!\n";
        }
    }
    else {
        $total ++;
        $ig_count ++;
    }
    last if eof (FASTA);  # If the sequence fasta file was fully read, then exit reading the fastq file
}
close (FASTA);    # Close sequence fasta file

print STDERR "The total line count from $fasta file is: $total; ignored line is: $ig_count and kept line is: $keep_count!\n\n";

# Define variables required for reading gff input file
my %premirna;   # Hash containing each precursor miRNA information
$total = 0; # Variable to count line in input file
my $out = 0;   # Variable to count line in output file

# Read and split the gff file for information collection
print OUTPUT "gene_id\tgene_name\tchromosome\tstart_position\tend_position\tstrand\tsequence\tprecursor_id\tprecursor_name\tprecursor_start\tprecursor_end\n";
while (1) {
    chomp (my $line = <GFF>);
    $total ++;
    unless ($line =~ /^#/) {
        my ($chromosome, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = (split(/\t/, $line));
        my ($id, $alias, $name, $derives) = (split(/;/, $attributes));
        $id =~ s/ID=//;
        $name =~ s/Name=//;
        if ($feature eq "miRNA_primary_transcript") {
            $premirna{$id}{name} = $name;
            $premirna{$id}{start} = $start;
            $premirna{$id}{end} = $end;
        }
        elsif ($feature eq "miRNA") {
            $derives =~ s/Derives_from=//;
            unless (exists $mirna_seq{$name}) {
                $mirna_seq{$name} = "Undefined";
            }
            print OUTPUT "$id\t$name\t$chromosome\t$start\t$end\t$strand\t$mirna_seq{$name}\t$derives\t$premirna{$derives}{name}\t$premirna{$derives}{start}\t$premirna{$derives}{end}\n";
            $out ++;
        }
        else {
            die "Biotype value not recognised at gff file line: $total!\n";
        }
    }
    last if eof GFF;
}
close (GFF);
close (OUTPUT);

print STDERR "There were $total lines from $gff input file and $out lines in the $output output file!\n";

my $finish_date = localtime;
print STDERR "\n##############################\nEND = $finish_date\n##############################\n\n";

__END__