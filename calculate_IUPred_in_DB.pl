#!/usr/local/bin/perl

##################################################################
=head1 # CALCULATE IUPRED ISD IN DB ##############################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will calculate the intrinsic structural disorder (ISD)
of polypeptide sequences located inside a MySQL database. This
program is similar to a shell script in that it performs several
small, independent tasks. The first task extracts data from the
MySQL database one row at a time. For each row, the IUPred program
is executed on the row's polypeptide sequence in order to generate
the ISD of that sequence. Note that we must extract only individual
rows from MySQL because IUPred only accepts one sequence per execution.
Finally, the generated ISD score will be placed back into the database.

Finally, this script should be executed on the computer running
the bioinformatics program; the MySQL database does not need to
be located on this same computer.

=head2 REQUIREMENTS:

The "IUPred" program. 

Additionally, this program requires the DBD::mysql module (and the
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

The following input variables are lexical and can be utilized
thoughout the entire "Calculate IUPred ISD in DB" program.

=cut

use warnings;
use strict;

print "Time of program execution:\n";
main::print_time();
print "\n";

my  ($execute_iupred_command,

    $mysql_extraction_query,
    $mysql_insertion_statement,
    
    $database_host_ip_address,
    $mysql_account_username,
    $mysql_account_password,
    $mysql_database_name,
    );

# Default ISD Program:
# ____________________________

$execute_iupred_command = '/home/USER_NAME/bin/iupred/iupred ./temp_input_file.fasta long';

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT GeneUniqueID, NoCysteineProteinSequence FROM TableName'; 
$mysql_insertion_statement = 'UPDATE TableName SET NoCysOrigDataIUPredAminoAcidISDScores = ?, NoCysIUPredMeanISD = ? WHERE MouseGenesPrimaryKey = ?';

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$execute_iupred_command = This is the Linux command that executes the
IUPred program. According to this command, IUPred must receive
input from a file called "temp_input_file.fasta"; the output for IUPred
is printed to STDOUT and, therefore, is input directly to an array (using
backticks[`]). The $execute_iupred_command must end with the word "long"
or "short". These are arguments that are required for IUPred to execute.
See the IUPred manual for more information on these arguments. 

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
ISD of only certain protein sequences, the user may add a "WHERE" statement
to the MySQL command. 

$mysql_insertion_statement = This is the MySQL command that inserts the
various ISD scores back into the MySQL database. For the most part, this
command should be input exactly as it would be input into a MySQL command
line interface. However, the command requires 3 placeholders (question marks).
The first placeholder is for the string of ISD values per amino acid, the
second is the mean ISD score for the entire polypeptide. Thie latter score is
derived calculating the mean of the former amino acid scores. The third
placeholder is the unique ID so we know where to insert the ISD information.

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
generated ISD information will be placed (i.e., the one containing the
output information from IUPred). 

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $database_insertion_object,
    $row_array,
    
    $unique_id,
    $protein_sequence,
    $print_fasta_command,
    
    @iupred_output_lines,
    $mean_isd_score,
    $amino_acid_isd_scores_string,
    );

# We must first access the MySQL database.

use DBI;

print 'Connecting to host database = ' . "$database_host_ip_address\n";

$database_connection_object = DBI->connect("dbi:mysql:$mysql_database_name;$database_host_ip_address",
                                           $mysql_account_username, $mysql_account_password);

print "Database connected\n";

# We now pre-prepare our MySQL commands that will be used to both
# extract and insert data.

$database_extraction_object = $database_connection_object->prepare($mysql_extraction_query);
$database_insertion_object = $database_connection_object->prepare($mysql_insertion_statement);

$database_extraction_object->execute();

print "Database extraction command executed\n";

while ($row_array = $database_extraction_object->fetchrow_arrayref)
    {
    # The $row_array is a reference to an array, each element of which
    # is a value from a row of the database. Specifically, each element
    # contains the value returned by the "SELECT" statement in the
    # $mysql_extraction_query.In this case, the zeroth
    # element is the PrimaryKey and the first element is the
    # ProteinSequence.
    
    # IUPred only accepts a single polypeptide sequence per execution.
    # We must now generate a FASTA containing the sequence from our current
    # database row.

    $unique_id = $row_array->[0];    
    $protein_sequence = $row_array->[1];

    undef $row_array;   

    # A FASTA file is required to execute IUPred.
    
    &print_fasta($unique_id, $protein_sequence);

    # We will now execute the IUPred command.
    
    @iupred_output_lines = `$execute_iupred_command` or die "Cannot input the temporary FASTA file into IUPred";

    !system 'rm ./temp_input_file.fasta' or die ("Cannot delete temp_input_file.fasta: $!\n");

    chomp @iupred_output_lines;
    
    ($amino_acid_isd_scores_string, $mean_isd_score) = &filter_iupred_output_for_isd_score(\@iupred_output_lines);

    undef @iupred_output_lines;

    # We must now insert the $mean_isd_score and the $amino_acid_isd_scores_string
    # into the database for this specific row.

    $database_insertion_object->execute($amino_acid_isd_scores_string, $mean_isd_score, $unique_id);

    undef $mean_isd_score;
    undef $unique_id;
    undef $protein_sequence;
    undef $amino_acid_isd_scores_string;
    };

$database_connection_object->disconnect;

no DBI;  

print "\nTime of program completion\n";
main::print_time();
print "\nProgram finished";

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# PRINT_TIME

# The following subroutine will use the "localtime" builtin perl function
# to print the date and time. Note that this is the local time on the machine.
# If the machine is a remote login server in a different time zone, then the
# time will be incorrect according to the time of the user's local computer.

sub print_time
    {
    
    # There are no input variables	
    
    # Subroutine Variables
    
    my (@master_time_variables,
        @month_array,
        
        $seconds,
        $minutes,
        $hours,
        $month,
        $day_number,
        $year,
        
        $time,
        $date,
        ); 
    
    @master_time_variables = localtime(time);

    # The @master_time_variables array contains all the information one needs
    # to print the local time. The array indexes contain information as follows:
    
    # [0] => Seconds past the minute
    # [1] => Minutes past the hour
    # [2] => Hours past midnight
    # [3] => Day of the month
    # [4] => Months past the start of the year
    # [5] => Number of years since 1900
    # [6] => Number of days since the start of the week (Sunday)
    # [7] => Number of days since the start of the year
    # [8] => Whether or not daylight savings is active
    
    @month_array = qw/January February March April May June July August September October November December/;

    $seconds = $master_time_variables[0];
    $minutes = $master_time_variables[1];
    $hours = $master_time_variables[2];
    $month = $month_array[ $master_time_variables[4] ];
    $day_number = $master_time_variables[3];
    $year = 1900 + $master_time_variables[5];
    
    $time = "$hours". ':' . "$minutes" . ':' . "$seconds";
    $date = "$year $month $day_number";
    
    print "\nDATE: $date, TIME: $time\n\n";
    
    };

########################################################################################

# FILTER_IUPRED_OUTPUT_FOR_ISD_SCORE

# The following subroutine will utilize the output from the IUPred
# program and calculate an overall ISD score for the polypeptide
# sequence. Because IUPred outputs an ISD score for each amino acid,
# this subroutine will calculate the mean ISD score for the entire
# polypeptide sequence. 

sub filter_iupred_output_for_isd_score
    {
    
    # Input Variables
    
    my 	($iupred_output_lines) = @_;
    
    # Subroutine Variables
    
    my  (@isd_scores,
        @unfinished_isd_scores,
        $summation,
        $mean_isd_score,
        $amino_acid_isd_scores_string,
        );
    
    $summation = 0;
    
    @unfinished_isd_scores = grep
        {
        ( $_ !~ m/^#/ ) and ( $_ =~ m/\d$/g );
        } @$iupred_output_lines;
    
    foreach (@unfinished_isd_scores)
        {
        $_ =~ m/\w\s+(\d+\.\d+)$/g;
            
        $summation += $1;
        
        push(@isd_scores, $1)
        };
    
    $mean_isd_score = $summation / ( scalar(@isd_scores) ); 

    # The following "sprintf" command will limit the $isd_score
    # variable to only four digits after the decimal point.

    $mean_isd_score = sprintf('%.4f', "$mean_isd_score");
    
    # We will now calculate the individual ISD scores for the amino acids.
    
    $amino_acid_isd_scores_string = join(',', @isd_scores);
    
    return ($amino_acid_isd_scores_string, $mean_isd_score);
    
    };

########################################################################################

# PRINT_FASTA

# The following program accepts a title line and sequence and prints
# out a FASTA file. This program will print the results to a file
# called "temp_input_file.fasta" in the current directory.

sub print_fasta
    {
    
    # Input Variables
    
    my 	($title_line,
        $full_sequence,
        ) = @_;
    
    # Subroutine Variables
    
    my  ($iupred_input_file,
        $small_sequence,
        );
    
    $summation = 0;

    $iupred_input_file = './temp_input_file.fasta';
    
    open(FASTA_INPUT_FOR_IUPRED, ">", $iupred_input_file) || die("Cannot open the $iupred_input_file file: $!\n");
    
    print FASTA_INPUT_FOR_IUPRED '>' . "$title_line\n";    
    
    # The following subroutine will separate the $full_sequence into
    # several subsequences, each consisting of only 80 characters.
    
    while ( length($full_sequence) > 0 )
        {
        # The following commands means:
        # Read the first (zeroth position) substring of 80 characters and
        # replace it with an empty string.
        
        $small_sequence = substr($full_sequence, 0, 80, ""); 
        
        print FASTA_INPUT_FOR_IUPRED "$small_sequence\n";
        };
    
    print FASTA_INPUT_FOR_IUPRED "\n";
    
    close FASTA_INPUT_FOR_IUPRED;
    
    };

########################################################################################

