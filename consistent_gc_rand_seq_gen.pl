#!/usr/local/bin/perl

##################################################################
=head1 # CONSISTENT GC RAND SEQ GEN ##############################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will generate random nucleotide sequences that have
the same length and GC content as the input nucleotide sequences.
These random sequences are output to STDOUT in a tab-delineated
format. The first column is the Unique ID from the input database
and the second column is the random nucleotide sequence.

Additionally, this script should be executed on the computer running
the bioinformatics program; the MySQL database does not need to
be located on this same computer.

=head2 REQUIREMENTS:

This program requires the DBD::mysql module (and the
library [DBD::mysql directory]). Although DBD and DBI are both
part of the core modules, the DBD:mysql driver is not and must be
downloaded from CPAN.

Finally, no protein unique ID may contain a space or tab.

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

The hardcoded values below are the default values for each variable.
Additionally, below these are descriptions of each variable.

=head2 ARGUMENT VALUES

The following input variables are package variables and can be
utilized thoughout the entire CONSISTENT_GC_RAND_SEQ_GEN program.

=cut

use warnings;
use strict;

package Input_Variables;   

use vars qw
    ($stop_codon_present
    
    $mysql_extraction_query

    $database_host_ip_address
    $mysql_account_username
    $mysql_account_password
    $mysql_database_name
    );

# Default Program Execution Options:
# _________________________________

$stop_codon_present = 'yes';
    
# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT PrimaryKey, NucleotideSequence FROM TableName'; 

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$stop_codon_present = This option can be either "yes" or "no". If this option
is "yes", then that means that the sequences end (or are supposed to end)
with a stop codon. If "yes", then stop codons are removed when When calculating
the GC content. The stop codons are later added back to the output nucleotide
sequences. This prevents the stop codon from influencing the GC content
calculations. It also prevents having to adjust the GC content in order to
generate a new stop codon at the end of the output nucleotide sequence. If
this option is "no", then the program assumes that the input nucleotide
sequences do not contain a stop codon, and thus it does not need to be removed.

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
aggregation propensity of only certain protein sequences, the user may add a
"WHERE" statement to the MySQL command. 

$database_host_ip_address = This is the IP address of the computer hosting
the MySQL database. If the MySQL database is located on the local computer
that is executing this program, the IP address should be set as
"localhost" or "127.0.0.1". 

$mysql_account_username = This is the username required to access the MySQL
database. Importantly, this is NOT the username required to access the
remote computer. That is, this is the MySQL username, NOT the computer login
username.

$mysql_account_password = This is the password required to access the MySQL
database. Importantly, this is NOT the password required to access the
remote computer. That is, this is the MySQL password, NOT the computer login
password.

$mysql_database_name = This is the name of the MySQL database containing the
information to be extracted. This is also the database into which the newly
generated aggregation information will be placed (i.e., the one containing the
output information from the bioinformatics program). 

=cut

package main;

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $row_array,
    $output_file_name,
    
    $unique_id,
    $gene_sequence,
    $stop_codon,
    
    $gc_content,
    
    $randomized_sequence,
    $random_binary,
    
    $three_prime_seq,
    $five_prime_seq,
    $codon,    
    $stop_codon_position,
    $random_nucleotide_position,
    );

# We must first access the MySQL database.

use DBI;

print 'Connecting to host database = ' . "$Input_Variables::database_host_ip_address\n";

$database_connection_object = DBI->connect("dbi:mysql:$Input_Variables::mysql_database_name;$Input_Variables::database_host_ip_address",
                                           $Input_Variables::mysql_account_username, $Input_Variables::mysql_account_password);

print "Database connected\n";

# We now pre-prepare our MySQL commands that will be used to extract data.

$database_extraction_object = $database_connection_object->prepare($Input_Variables::mysql_extraction_query);

$database_extraction_object->execute();

print "Database extraction command executed\n";

$output_file_name = './gc_rand_seq_output.tab';

open (GC_RAND_OUTPUT, ">", "$output_file_name") || die("Cannot open the $output_file_name file: $!\n"); 

while ($row_array = $database_extraction_object->fetchrow_arrayref)
    {
    chomp @$row_array;
    
    # The $row_array is a reference to an array, each element of which
    # is a value from a row of the database. Specifically, each element
    # contains the value returned by the "SELECT" statement in the
    # $mysql_extraction_query.In this case, the zeroth
    # element is the PrimaryKey and the first element is the
    # NucleotideSequence.
    
    $unique_id = $row_array->[0];    
    $gene_sequence = $row_array->[1];    
    
    undef $row_array;
    
    # We must first check the integrity of the gene sequence.

    $gene_sequence =~ s/[^atgcn]//ig;
    
    # If the $stop_codon_present is set to "yes", we must subtract the
    # stop codon from the sequence (we will add it back later). If the
    # nucleotide sequence does not end in a stop codon, then nothing
    # will be removed from the 5' end of the DNA sequence.
    
    $stop_codon = '';
    
    if ($Input_Variables::stop_codon_present =~ m/^yes$/i)
        {
        $stop_codon = substr($gene_sequence, -3, 3, ''); 
    
        unless ($stop_codon =~ m/^((taa)|(tag)|(tga))$/i)
            {
            $gene_sequence = $gene_sequence . $stop_codon;
            $stop_codon = '';
            };    
        };

    # We can now calculate the GC content percentage.
    # This calculation is performed AFTER the stop codon
    # has been removed (assuming $stop_codon_present = 'yes').
    # That is, we do not want to include the stop codon in
    # the GC content percentage.
    
    $gc_content = &gc_content_percent($gene_sequence);
    
    # The following "foreach" loop will generate A relative to T and G
    # relative to C with 50% probability. For example, there is a 50%
    # probability of each "A" being converted into a "T" and a 50%
    # probability of each one remaining an "A". Note that removing the
    # following "foreach" loop will scramble the sequences while
    # keeping the nucleotide frequency intact. 
    
    $randomized_sequence = '';
    
    foreach ( split('', $gene_sequence) )
        {
        $random_binary = int(rand(2));
        
        if ( ($_ =~ m/(a|t)/i) and ($random_binary == 0) ) {$_ = "A"}            
            
        elsif ( ($_ =~ m/(c|g)/i) and ($random_binary == 0) ) {$_ = "C"}            
        
        elsif ( ($_ =~ m/(a|t)/i) and ($random_binary == 1) ) {$_ = "T"}
        
        elsif ( ($_ =~ m/(c|g)/i) and ($random_binary == 1) ) {$_ = "G"};
        
        $randomized_sequence = $randomized_sequence . $_;
        };
    
    # We will now scramble the sequence. 
    
    $randomized_sequence = &randomizer($randomized_sequence); 
    
    # We must now check for stop codons. The program will remove
    # a random nucleotide from any detected stop codon and place it 
    # into a randomly chosen position in the sequence. Importantly,
    # the following variable names are based upon the nucleotide
    # sequence being a DNA sequence, not a transcript sequence. 
    
    $five_prime_seq = $randomized_sequence;
    
    undef $randomized_sequence;
    
    $three_prime_seq = '';
 
    # The following "until" loop must be named so that the
    # variables (i.e., $five_prime_seq) can be reset each
    # time the loop is reinitiated.
 
    STOP_CODON_REMOVE: until ( length($five_prime_seq) == 0 )
        {
        $codon = substr($five_prime_seq, 0, 3, '');
        
        if ( $codon =~ m/^((tag)|(taa)|(tga))$/i )
            {
            # Because a stop codon was found, the algorithm must reexamine
            # the entire sequence. Therefore, the $five_prime_seq will
            # revert to again being the entire sequence.
            
            $five_prime_seq = $three_prime_seq . $codon . $five_prime_seq;
            
            # We must first determine the positions of the two nucleotides
            # that need to be switched.
            
            $random_nucleotide_position = int( rand( length($five_prime_seq) ) );
            
            # The following "rand" command will produce a number that ranges from
            # 0 to 2.9999. The stop codon has positions "0", "1", & "2".
            
            $stop_codon_position = length($three_prime_seq) + int( rand(3) );
            
            # The positions entered into the following subroutine assume the
            # input sequence begins at position zero (0), not one (1). This is
            # why we only add 0, 1, or 2 to the $three_prime_seq length when
            # determining the $stop_codon_position instead of adding 1, 2, or 3.
            
            $five_prime_seq = switch_characters($five_prime_seq, $stop_codon_position, $random_nucleotide_position);
    
            $three_prime_seq = '';
            
            redo STOP_CODON_REMOVE;    
            };
        
        $three_prime_seq = $three_prime_seq . $codon;
        };
    
    $randomized_sequence = $three_prime_seq . $stop_codon;
    
    print GC_RAND_OUTPUT "$unique_id\t$gc_content\t$randomized_sequence\n";
    
    undef $stop_codon;
    undef $three_prime_seq;
    undef $five_prime_seq;
    undef $randomized_sequence;
    undef $gc_content;
    };
    
close GC_RAND_OUTPUT;

$database_connection_object->disconnect;

no DBI;  

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# RANDOMIZER

# The following subroutine will randomly order the characters that
# are present in the input string of characters.

sub randomizer
    {
    
    # Input Variables
    
    my  ($original_string) = @_;
    
    # Subroutine Variables
    
    my  (@original_string_array,
        $remaining_character_number,
        $index,
        @randomized_string_array, 
        $randomized_string,
        );
    
    @original_string_array = split('', $original_string);
    
    # We can now begin randomizing the sequence.
    
    until ( scalar(@original_string_array) == 0 )
        {    
        $remaining_character_number = scalar(@original_string_array);
        
        # The value that is output from rand() ranges from 0 to
        # LESS THAN the number input into rand() (in this case
        # $remaining_character_number). Therefore, although we are
        # using the "scalar" command instead of "$#", we are still
        # giving rand() the correct range.
        
        $index = int( rand($remaining_character_number) );
        
        # We will now place the character at $index into the
        # new randomized array
        
        push( @randomized_string_array, $original_string_array[$index] );
        
        splice(@original_string_array, $index, 1);
        };
    
    undef @original_string_array;
    
    $randomized_sequence = join('', @randomized_string_array);
    
    undef @randomized_string_array;

    return $randomized_string;
    
    };

######################################################################################## 

# SWITCH_CHARACTERS

# The following subroutine will take two characters in a string and
# switch them. The characters are identified by position, not by
# the actual identity of the character. Importantly, the string begins
# at position zero (0), not one (1).

sub switch_characters
    {
    
    # Input Variables
    
    my  ($string, $first_switch_position, $second_switch_position) = @_;
    
    # Subroutine Variables
    
    my  ($first_character_identity,
        $second_character_identity,
        );
    
    if ($first_switch_position == $second_switch_position) {return $string};
    
    $string =~ m/^[a-z]{$first_switch_position}([a-z])/i;
    $first_character_identity = $1;
    
    $string =~ m/^[a-z]{$second_switch_position}([a-z])/i;
    $second_character_identity = $1;
    
    substr($string, $first_switch_position, 1, "$second_character_identity");
    substr($string, $second_switch_position, 1, "$first_character_identity");
    
    return $string;
    
    };
    
########################################################################################

# GC_CONTENT_PERCENT

# The following subroutine will calculate the GC content percent
# given an input nucleotide sequence. Any unknown nucleotides will
# not be counted as G or C. Additionally, any abbreviation representing
# multiple nucleotides will not be counted as G or C except S (which
# means "G or C").

sub gc_content_percent
    {
    
    # Input Variables
    
    my  ($gene_sequence) = @_;
    
    # Subroutine Variables
    
    my  ($total_length,
        $gc_content,
        );
    
    $total_length = length($gene_sequence);
    
    $gene_sequence =~ s/g|c|s//gi;
    
    $gc_content = ( $total_length - length($gene_sequence) ) / $total_length;
    
    return $gc_content;
    
    };
    
########################################################################################




