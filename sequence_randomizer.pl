#!/usr/local/bin/perl

use warnings;
use strict;

##################################################################
=head1 # SEQUENCE RANDOMIZER #####################################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program accepts a tab-delineated input file filled with nucleotide
or protein sequences (it can actually accept any string of characters).
It then randomizes the sequences and outputs them to a tab-delineated file. 

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

After the executable, one should type the path and filename of the
input tab-delineated file. The first column in the tab-delineated
file should be an ID for the gene/protein and the second column
should be the sequence. 

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($input_file_line,
    $original_sequence,
    $sequence_name,    
    $randomized_sequence,
    $output_file_name,
    );

$output_file_name = 'randomized_sequences.tab';
  
open(RANDOM_SEQ_OUTPUT, ">", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

# We will now randomize each sequence.

while (defined( $input_file_line = <> ))
    {
    chomp $input_file_line;
    
    $input_file_line =~ m/^(\d+)\t([a-z]+)$/gi;
    
    $sequence_name = $1;
    $original_sequence = $2;    
    
    $randomized_sequence = &randomizer($original_sequence);    
    
    # Now that the sequence is randomized, we can now print it
    # to the output file.
    
    print RANDOM_SEQ_OUTPUT "$sequence_name\t$randomized_sequence\n";
    };

close RANDOM_SEQ_OUTPUT;

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# RANDOMIZER

# The following subroutine will randomize the input sequence. That is,
# it will randomly order the characters that are present in the
# input string. The randomizing engine is the perl "rand" command. 

sub randomizer
    {
    
    # Input Variables
    
    my  ($original_sequence) = @_;
    
    # Subroutine Variables
    
    my  (@original_seq_array,
        $remaining_aa_number,
        $index,
        @randomized_seq_array, 
        $randomized_sequence,
        );
    
    @original_seq_array = split('', $original_sequence);
    
    # We can now begin randomizing the sequence.
    
    until ( scalar(@original_seq_array) == 0 )
        {    
        $remaining_aa_number = scalar(@original_seq_array);
        
        # The value that is output from rand() ranges from 0 to
        # LESS THAN the number input into rand() (in this case
        # $remaining_aa_number). Therefore, although we are using
        # the "scalar" command instead of "$#", we are still
        # giving rand() the correct range.
        
        $index = int( rand($remaining_aa_number) );
        
        # We will now place the amino acid at $index into the
        # new randomized array
        
        push( @randomized_seq_array, $original_seq_array[$index] );
        
        splice(@original_seq_array, $index, 1);
        };
    
    undef @original_seq_array;
    
    $randomized_sequence = join('', @randomized_seq_array);
    
    undef @randomized_seq_array;

    return $randomized_sequence;
    
    };

########################################################################################  
