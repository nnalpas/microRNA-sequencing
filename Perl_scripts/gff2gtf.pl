#!/usr/bin/perl -w

# Script used to convert gff format to gtf format

# Define all modules to be used in script
use strict;
use warnings;
use Getopt::Long;
use IO::File;

# Define the input file
my $input;    # Input file containing gff format annotation
my $output;  # Output file containing gtf format annotation

# Define the parameter in order to submit input files to this script
&GetOptions (
    'i=s' => \$input,
    'o=s' => \$output,
);

my $start_date = localtime;
print STDERR "START = $start_date\n\n";

# Open the input gtf file
unless ($input) {
    die "Please specify the gff file to convert via -i parameter!\n";
}
open (INPUT, "<$input") || die "Cannot open $input: $!\n"; $_="1";

# Define output file
unless ($output) {
    $output = $input;
    chomp ($output);
    $output =~ s/(.*)\.gff3$/$1\.gtf/;
    print STDERR "File name for gtf output file not provided (-o parameter), name generated from input gff file: $output\n";
}

# Open the ouput file which will be in gtf format
if (-e $output) {
    die "This file: $output already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Define variables required for reading input file
my $gff_line = 0;

# Read and split the gff file for converting to gtf format
while (1) {
    chomp (my $line = <INPUT>);
    $gff_line ++;
    my $biotype;
    unless ($line =~ /^#/) {
        my ($chromosome, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = (split(/\t/, $line));
        my ($id, $alias, $name, $derives) = (split(/;/, $attributes));
        if ($feature eq "miRNA_primary_transcript") {
            $biotype = "pre-miRNA";
        }
        elsif ($feature eq "miRNA") {
            $biotype = "miRNA";
        }
        else {
            die "Biotype value not recognised at gff file line: $gff_line!\n";
        }
        $id =~ s/ID\=(.*)/gene_id \"$1\"/;
        $name =~ s/Name\=(.*)/gene_name \"$1\"/;
        print OUTPUT "$chromosome\t$biotype\texon\t$start\t$end\t$score\t$strand\t$frame\t$id\; $name\;\n";
    }
    last if eof INPUT;
}

# Close the files
close (INPUT);
close (OUTPUT);

print STDERR "Converting from gff format to gtf format completed!\n\n";

my $finish_date = localtime;
print STDERR "Finish = $finish_date\n\n";

__END__