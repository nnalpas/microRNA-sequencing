#!/usr/bin/perl -w

########################################################################################
########################################################################################
##                                                                                    ##
## Script to collect the isomiR counts and information from miRDeep2 software outputs ##
##                                                                                    ##
########################################################################################
########################################################################################

######################
# Load Perl packages #
######################

# Define the various packages to use within the script
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Array::Utils qw(:all);

# Print the local time at the start
my $start_date = localtime;
print STDERR ("Start = $start_date\n\n");

####################################################
# Define and open the different input/output files #
####################################################

# Input MRD file containing pri-miRNA counts and isomiR details
my $mrd;

# Output TXT file containing isomiR counts
my $output;

# Define the parameters in order to submit file names to this script
&GetOptions (
    'mrd=s' => \$mrd,
    'output=s' => \$output,
);

# Check for submission of input file
unless ($mrd) {
    die "Please specify the MRD file from miRDeep2 via -mrd parameter!\n";
}
# Open the input file
open (MRD, "<$mrd") || die "Cannot open $mrd: $!\n"; $_="1";

# Check for submission of output file
unless ($output) {
    die "Please specify the output TXT file via -output parameter!\n";
}

###############################################
# Read in the information from the input file #
###############################################

# Define variables required for reading the input file
# Scalar containing line from the MRD input file
my $line;

# Number of line already read from the input file
my $line_count = 0;

# Number of mature-miRNA and isomiR in MRD file for the current pri-miRNA (sanity scalars)
my $isomir_total_count = 0;
my $isomir_match_count = 0;
my $isomir_remain_count = 0;
my $bta_counter = 0;

# Define scalar that will contain the total isomiR and mature-miRNA counts for the current pri-miRNA (sanity scalars)
my $total_iso_count = 0;
my $total_mat_count = 0;
my $total_remain_count = 0;

# Hash containing all mature-miRNA information
my %mirna;

# Hash containing all isomiR information
my %isomir;

# Define the scalars that will contain pri-miRNA information
my ($primirna_name, $primirna_id, $primirna_count, $primirna_seq) = "";

# Define the array that will contain all mature-miRNA ID for the current pri-miRNA
my @maturemirna;

# Process input files
PROCESSING: while(1) {
    # Read lines one by one from the MRD input file
    chomp ($line = <MRD>);
    # Add one to the processed line counter
    $line_count ++;
    # Check if this line contains a non-required value
    if ($line =~ /^pri_struct/) {
        # Go to next line in file
        next PROCESSING;
    }
    # Check if this line is an empty line
    elsif ($line eq "") {
        # Go to next line in file
        goto FILE_END;
    }
    # Check if this is the first line of the record for a specific pri-miRNA
    elsif ($line =~ /^>/) {
        # Check if this is the first line read from the file
        unless ($line_count == 1) {
            # Check whether the current pri-miRNA count equals the current mature-miRNAs counts
            if ($primirna_count != $total_mat_count) {
                # Stop the script with the following message
                die "Error! The overall check of pri-miRNA count versus mature-miRNA count does not match for $primirna_id within file $mrd, with pri count: $primirna_count and mature count: $total_mat_count!\n";    
            }
            # Check whether the current pri-miRNA contains only regular mature-miRNA name and whether all counts match (pri-miRNA vs mature-miRNA vs isomiR)
            elsif (($primirna_count != $total_mat_count) | ($isomir_total_count != ($isomir_match_count + $isomir_remain_count))) {
                # Stop the script with the following message
                die "Error! The overall check of pri-miRNA count versus isomiR count does not match for $primirna_id within file $mrd!\n";
            }
            # Re-set to zero the number of matching and non-matching isomiR and mature-miRNA for the sanity scalars
            $isomir_total_count = 0;
            $isomir_match_count = 0;
            $isomir_remain_count = 0;
            # Re-set to zero the counts of isomiR and mature-miRNA for the sanity scalars
            $total_iso_count = 0;
            $total_mat_count = 0;
            $total_remain_count = 0;
        }
        # Split line into scalars for the current pri-miRNA name and id
        ($primirna_name, $primirna_id) = split (/_/, $line);
        # Remove the special character at the start of the pri-miRNA name
        $primirna_name =~ s/^>//;
        # Re-set to empty the array containing mature-miRNA for the current pri-miRNA
        undef @maturemirna;
    }
    # Check if this line contains mature-miRNA information
    elsif ($line =~ /^bta-/) {
        # Add one to the counter of bta matching line that correspond to mature-miRNA read
        $bta_counter ++;
        # Define the scalars that will contain the current mature-miRNA information
        my ($maturemirna_name, $maturemirna_id, $word1, $word2, $count);
        # Check whether the word count is separated by space from the actual count value for the current mature-miRNA
        if ($line =~ /count\d+$/) {
            # Split line into scalars for the current mature-miRNA name, id and count (two words are ignored)
            ($maturemirna_name, $maturemirna_id, $word1, $word2) = split (/\s+|_/, $line);
            # Due to the lack of space between the word count and the actual count value, the count value is defined separately
            ($count = $word2) =~ s/^count(\d+)$/$1/;
        }
        else {
            # Split line into scalars for the current mature-miRNA name, id and count (two words are ignored)
            ($maturemirna_name, $maturemirna_id, $word1, $word2, $count) = split (/\s+|_/, $line);
        }
        # Add the current mature-miRNA count value to the total mature-miRNA count for the current pri-miRNA
        $total_mat_count += $count;
        # Check that the mature-miRNA ID does not already exists in the corresponding hash
        if (exists $mirna{$maturemirna_id}) {
            # Find how many mature-miRNA ID match the current value in the corresponding hash
            my @previous_mat = grep {/$maturemirna_id/} keys %mirna;
            # Collate the number of matching ID plus one to the current mature-miRNA ID
            $maturemirna_id .= join ("", '_', (scalar @previous_mat + 1));
        }
        # Include all mature-miRNA ID for the current pri-miRNA into an array
        push (@maturemirna, $maturemirna_id);
        # Include this new mature-miRNA information in the corresponding hash
        $mirna{$maturemirna_id}{name} = $maturemirna_name;
        $mirna{$maturemirna_id}{count} = $count;
        $mirna{$maturemirna_id}{pri_id} = $primirna_id;
        $mirna{$maturemirna_id}{pri_name} = $primirna_name;
        $mirna{$maturemirna_id}{pri_count} = $primirna_count;
        $mirna{$maturemirna_id}{pri_seq} = $primirna_seq;
        # Define scalar that contains the size of the %mirna hash
        my $size = keys %mirna;
        # Check if the number of mature-miRNA read matches the %mirna size
        if ($bta_counter != $size) {
            # Stop the script with the following message
            die "Error! The number of mature-miRNA read does not correspond to the number of mature-miRNA recoreded in the hash at line $line_count within file $mrd, please seek advice!\n";
        }
    }
    # Check if this line contains total read count information
    elsif ($line =~ /^total read count/) {
        # Check whether there is an actual count value for the total read count line
        if ($line !~ /\d+$/) {
            # If not define the total read count value as zero
            $primirna_count = 0;
        }
        else {
            # Keep only the count value for the pri-miRNA
            ($primirna_count = $line) =~ s/^.*\s+(\d+)$/$1/;
        }
    }
    # Check if this line contains remaining read count information
    elsif ($line =~ /^remaining read count/) {
        # Keep only the count value for the remaining isomiR count
        ($total_remain_count = $line) =~ s/^.*\s+(\d+)$/$1/;
        # Add the count value for the remaining reads to the total mature-miRNA count for sanity check
        $total_mat_count += $total_remain_count;
    }
    # Check if this line contains pri-miRNA sequence information
    elsif ($line =~ /^pri_seq/) {
        # Keep only the sequence value for the pri-miRNA
        ($primirna_seq = $line) =~ s/^.*\s+(\w+)$/$1/;
    }
    # Check if this line contains mature-miRNA position information relative to the pri-miRNA sequence
    elsif ($line =~ /^exp/) {
        # Keep only the position of mature-miRNA relative to start of pri-miRNA
        (my $mirna_pos = $line) =~ s/^exp\s+(\w+)$/$1/;
        # Define the array that will contain the mature-miRNA ID that have the same position marking
        my @mirna_5;
        my @mirna_3;
        my @mirna_m;
        # Go successively through all the mature miRNA ID for the current pri-miRNA
        foreach my $mat (@maturemirna) {
            # Check whether the mature-miRNA name contains a 5p
            if ($mirna{$mat}{name} =~ /^.*-(5)p$/) {
                # Define the position value to look for in pri-miRNA sequence based on mature-miRNA name
                $mirna{$mat}{look_pos} = $1;
                # Add the current mature-miRNA ID to the array of 5' located mature-miRNA
                push (@mirna_5, $mat);
            }
            # Check whether the mature-miRNA name contains a 3p
            elsif ($mirna{$mat}{name} =~ /^.*-(3)p$/) {
                # Define the position value to look for in pri-miRNA sequence based on mature-miRNA name
                $mirna{$mat}{look_pos} = $1;
                # Add the current mature-miRNA ID to the array of 3' located mature-miRNA
                push (@mirna_3, $mat);
            }
            # Check whether the mature-miRNA name finishes by a alphanumeric value
            elsif ($mirna{$mat}{name} =~ /^.*\w$/) {
                # Define the position value to look for in pri-miRNA sequence based on mature-miRNA name
                $mirna{$mat}{look_pos} = "M";
                # Add the current mature-miRNA ID to the array of undetermined located mature-miRNA
                push (@mirna_m, $mat);
            }
            # Stop the script if the current mature-miRNA name is not recognized
            else {
                # Stop the script with the following message
                die "Error! The name of the mature-miRNA $mat is not recognized for $primirna_id within file $mrd, please seek advice!\n";
            }
        }
        # Define the array that will contain the mature-miRNA that require a zero as start and end position value
        my @zero_pos;
        # Check if several mature-miRNA for the current pri-miRNA are 5' located or if the pri-miRNA position does not contain 5 markings or there are several 5 markings or overlapping markings
        if ((scalar @mirna_5 > 1) | $mirna_pos !~ /5+/ | $mirna_pos =~ /5+[^5]+5+/ | $mirna_pos =~ /5+[M3]+/) {
            # Add the 5' located mature-miRNA ID to the zero position array
            push (@zero_pos, @mirna_5);
        }
        # Check if several mature-miRNA for the current pri-miRNA are 3' located or if the pri-miRNA position does not contain 3 markings or there are several 3 markings or overlapping markings
        if ((scalar @mirna_3 > 1) | $mirna_pos !~ /3+/ | $mirna_pos =~ /3+[^3]+3+/ | $mirna_pos =~ /[5M]+3+/) {
            # Add the 3' located mature-miRNA ID to the zero position array
            push (@zero_pos, @mirna_3);
        }
        # Check if several mature-miRNA for the current pri-miRNA are undetermined located or if the pri-miRNA position does not contain M markings or there are several M markings or overlapping markings
        if ((scalar @mirna_m > 1) | $mirna_pos !~ /M+/ | $mirna_pos =~ /M+[^M]+M+/ | $mirna_pos =~ /5+M+|M+3+/) {
            # Add the undetermined located mature-miRNA ID to the zero position array
            push (@zero_pos, @mirna_m);
        }
        # Get the mature-miRNA ID from array @maturemirna that are not in array @zero_pos
        my @evaluate_pos = array_minus(@maturemirna, @zero_pos);
        # Check if the array @evaluate_pos is non-empty
        if (@evaluate_pos) {
            # Go successively through all the mature miRNA ID present in array @evaluate_pos
            foreach my $mat (@evaluate_pos) {
                # Define for the mature-miRNA the position value to look for in the pri-miRNA position
                my $pos_val = $mirna{$mat}{look_pos};
                # Check whether the position value to look for in the pri-miRNA position is present or not
                if ($mirna_pos =~ /f*($pos_val+)f*/) {
                    # Include the start and end position information for mature-miRNA in the corresponding hash
                    $mirna{$mat}{start} = ($-[1]);
                    $mirna{$mat}{end} = ($+[1]);
                }
                else {
                    # Stop the script with the following message
                    die "Error! The line $line contains an unrecognized pattern of position for $primirna_id within file $mrd, please seek advice!\n";
                }
            }
        }
        # Check if the array @zero_pos is non-empty
        if (@zero_pos) {
            # Go successively through all the mature miRNA ID present in array @zero_pos
            foreach my $mat (@zero_pos) {
                # Define and include as zero the start and end position information for mature-miRNA in the corresponding hash
                $mirna{$mat}{start} = 0;
                $mirna{$mat}{end} = 0;
            }
        }    
    }
    # Check if this line contains isomiR information
    elsif ($line =~ /^seq_/) {
        # Add one to the number of isomiR processed for the current pri-miRNA
        $isomir_total_count ++;
        # Split line into scalars for the current isomiR id, sequence, number of mismatches and count
        my ($isomir_id, $isomir_count, $isomir_seq, $isomir_mm) = split (/_x|\s+/, $line);
        # Add the count value for the current isomiR to the total count for current pri-miRNA
        $total_iso_count += $isomir_count;
        # Keep only the actual isomiR sequence
        $isomir_seq =~ s/^\.*(\w+)\.*$/$1/;
        # Define the scalar containing the start and end of the isomiR
        my $start = ($-[1]);
        my $end = ($+[1]);
        # Calculate the middle position value for the isomiR
        my $middle = ($start + ($end - $start) / 2);
        # Determine to which mature-miRNA the isomiR matches based on position
        my $check_count = 0;
        # Define the scalar that will contain the list of mature-miRNA that the current isomiR matches to
        my $list_mature = "";
        # Go successively through all the mature miRNA ID for the current pri-miRNA
        foreach my $val (@maturemirna) {
            # Check that the isomiR position matches a mature-miRNA position
            if ($middle >= $mirna{$val}{start} && $middle <= $mirna{$val}{end}) {
                # Add one to the counter for checking current match per pri-miRNA
                $check_count ++;
                # Attribute the matching mature-miRNA ID to the matching list
                $list_mature = $val;
            }
        }
        # Define scalar that will contain the mature-miRNA and pri-miRNA fields for the current isomiR
        my $isomir_pri;
        # Check to how many of the current mature-miRNA the isomiR matched
        if ($check_count == 1) {
            # Add one to the counter for matching isomiR
            $isomir_match_count ++;
            # Define the values for the pri field in the %isomir hash for the current isomiR to contain information about pri-miRNA and mature-miRNA
            $isomir_pri = join ("", $primirna_id, ' (', $start, '-', $end, ', #', $mirna{$list_mature}{pri_count}, ', mm', $isomir_mm, ' [', $list_mature, ', ', $mirna{$list_mature}{start}, '-', $mirna{$list_mature}{end}, ', #', $mirna{$list_mature}{count}, '])');
        }
        elsif ($check_count == 0) {
            # Add one to the counter for remaining (unmatching) isomiR
            $isomir_remain_count ++;
            # Define the values for the pri field in the %isomir hash for the current isomiR to contain information about pri-miRNA and mature-miRNA
            $isomir_pri = join ("", $primirna_id, ' (', $start, '-', $end, ', #', $primirna_count, ', mm', $isomir_mm, ' [NA, NA, #NA])');
        }
        else {
            # Stop the script with the following message
            die "Error! The isomiR $isomir_id is matching several mature-miRNA for the pri-miRNA $primirna_id within file $mrd, please seek advice!\n";
        }
        # Check whether the isomiR ID already exists in the corresponding hash
        if (exists $isomir{$isomir_id}) {
            # Check if the values previously registered for the current isomiR are still correct
            if (($isomir_count != $isomir{$isomir_id}{count}) | ($isomir{$isomir_id}{seq} !~ /$isomir_seq/i)) {
                # Stop the script with the following message
                die "Error! The isomiR $isomir_id within file $mrd was found in the isomiR hash with either different count or sequence, please seek advice!\n";
            }
            # Check if the pri-miRNA has already been recorded for the current isomiR
            elsif ($isomir{$isomir_id}{pri} =~ $primirna_id) {
                # Replace the pri-miRNA field values for the current isomiR by stating that it is a multihit for the current pri-miRNA
                $isomir{$isomir_id}{pri} =~ s/($primirna_id) \(.+?, (#\d+), .+? \[.+?, .+?, #.+?\]\)/$1 \(Multihit, $2, mmNA \[NA, NA, #NA\]\)/;
            }
            # If none of the above are true then concatenate the new value
            else {
                # Include these new information for the existing isomiR in the corresponding hash
                $isomir{$isomir_id}{pri} .= join ("", ', ', $isomir_pri);
            }
        }
        else {
            # Include this new isomiR information in the corresponding hash
            $isomir{$isomir_id}{count} = $isomir_count;
            $isomir{$isomir_id}{seq} = $isomir_seq;
            $isomir{$isomir_id}{pri} = $isomir_pri;
        }
    }
    # If none of the above checks are true then stop the script
    else {
        # Stop the script with the following message
        die "Error! The line $line contains an unrecognized value within file $mrd, please seek advice!\n";
    }
    # If the MRD input file was fully read, then exit reading
    FILE_END:if (eof (MRD)) {
        # Check whether the current pri-miRNA count equals the current mature-miRNAs counts
        if ($primirna_count != $total_mat_count) {
            # Stop the script with the following message
            die "Error! The overall check of pri-miRNA count versus mature-miRNA count does not match for $primirna_id within file $mrd, with pri count: $primirna_count and mature count: $total_mat_count!\n";    
        }
        # Check whether the current pri-miRNA contains only regular mature-miRNA name and whether all counts match (pri-miRNA vs mature-miRNA vs isomiR)
        elsif (($primirna_count != $total_mat_count) | ($isomir_total_count != ($isomir_match_count + $isomir_remain_count))) {
            # Stop the script with the following message
            die "Error! The overall check of pri-miRNA count versus isomiR count does not match for $primirna_id within file $mrd!\n";
        }
        last;
    }
}
# Close the input MRD file
close (MRD);

#####################################
# Open and write on the output file #
#####################################

# Open the output file only if the outpout file does not already exists
if (-e $output) {
    die "This file: $output already exists, cannot overwrite it: $!\n"; $_="1";
}
else {
    open (OUTPUT, ">$output") || die "Cannot open $output: $!\n"; $_="1";
}

# Print the columns header to the output file
print OUTPUT ("IsomiR_ID\tSequence\tMatching_primiRNA_(isomiR_position,pri_count,mismatch[Matching_miRNA,miRNA_position,count])\tCount\n");

# Go successively through all the isomiR ID recorded for the input file
foreach my $key (keys %isomir) {
    # Print to the output file all the recorded isomiR information
    print OUTPUT ("$key\t$isomir{$key}{seq}\t\"$isomir{$key}{pri}\"\t$isomir{$key}{count}\n");
}

# Close the ouput TXT file
close (OUTPUT);

###########################################
# Finalisation and script summary results #
###########################################

# Define the scalar containing the number of mature miRNAs and nnumber of isomiRs
my $length_iso = keys %isomir;
my $length_mature = keys %mirna;

# Print final message in STDERR to indicate script run completing
print STDERR ("There are $length_iso isomiRs for $length_mature mature-miRNAs in the MRD file $mrd!\n\n");

# Print the local time at the end
my $finish_date = localtime;
print STDERR ("Finish = $finish_date\n\n");

__END__