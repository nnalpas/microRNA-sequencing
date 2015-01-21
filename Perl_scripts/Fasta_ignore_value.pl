#!/usr/bin/perl -w

# Script used on sequence fasta file for parsing into smaller fasta files

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the different input/output files such as fastq sequence file, indices list file, each index output fastq file
my $fasta;    # Input file listing all indices used in the pool RNA-seq library
my $ignore;         # Header value (and subsequent sequence) to ignore
my $output;     # Ouput file containing fasta sequences

# Define the parameter in order to submit input files to this script
&GetOptions (
    'fasta=s' => \$fasta,
    'ignore=s' => \$ignore,
    'output=s' => \$output,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input fasta file
unless ($fasta) {
    die "Please specify the fasta file via -fasta parameter!\n";
}
open (FASTA, "<$fasta") || die "Cannot open $fasta: $!\n"; $_="1";

# Obtain the value to ignore
unless ($ignore) {
    die "Please specify the value to ignore in the fasta header via -ignore parameter!\n";
}

# Open the ouput file which will be in gtf format
unless ($output) {
    die "Please specify the output fasta file via -output parameter!\n";
}
if (-e $output) {
    die "This file: $output already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Define variables required for reading the input file
my $action;
my $total = 0;  # Variable to count line in input file
my $ig_count = 0;  # Variable to count line ignored from input file
my $keep_count = 0;  # Variable to count line kept from input file

# Read in the fastq sequence file
while(1) {
    my $line;        # Scalar containing line from the fasta file
    chomp($line = <FASTA>);    # Read lines one by one from the fasta file
    if ($line =~ /${ignore}/){   # Try to match the value to ignore
        $action = "remove";
        $total ++;
        $ig_count++;
    }
    elsif ($line =~ /^>/) {
        $action = "keep";
        $total ++;
        $keep_count++;
        print OUTPUT ("$line\n");
    }
    elsif ($action eq "remove") {
        $total ++;
        $ig_count++;
    }
    else {
        $total ++;
        $keep_count++;
        print OUTPUT ("$line\n");
    }
    last if eof (FASTA);  # If the sequence fastq file was fully read, then exit reading the fastq file
}
close (FASTA);    # Close sequence fasta file
close (OUTPUT);     # Close the output fasta file

print STDERR "The total line count from input file is: $total; ignored line is: $ig_count and kept line is: $keep_count!\n\n";

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__