#!/usr/local/bin/perl

##################################################################
=head1 # DETERMINE GENE FAMILIES #################################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program acts as a filter for gene families. It prevents
one older gene from assigning the phylostratum (PS) of the entire
gene family filled with newer genes. For example, traditionally,
if I have 20 genes in PS18 and one gene (in the same family) in PS7,
then the entire family will be assigned to PS7. Instead, this program
will calculate a series of decisions to determine if this gene
family should be split into two families: one family for the 20
newer genes and the other family for the one older gene. For genes
to be in the same family, the following conditions must be met:

1) Oldest PS must must contain two or more genes.
2) More than 50% of the newer genes must be homologous to at least
    one of the older genes.
3) There must be more than one newest gene.

If the aforementioned conditions are not met, the genes are split
into two (or possibly more) families.

This program requires input from two sources: an input text file
and a MySQL database. The text file is tab-delineated and lists
the homologous pairs of genes (e.g., mouse_gene_1 in the first
column and a homologous mouse_gene_2 in the second column). The
first variable obtained from the MySQL database is the PS number of
the gene, then the family of that gene. Family information can be
generated using the script entitled "pairwise_homolog_phylostrata_clustering.pl".
This script will return the modified family numbers back into
the database. Importantly, most likely one will only want to run this
script on coding genes in the MySQL table. If the MySQL table also
contains control polypeptide sequences, use a "WHERE" command in the
$mysql_extraction_query statement to remove any control sequences.  

The program will output to a tab-delineated file that contains the
gene's unique ID, the gene family number, and the PS number of
that gene family. It outputs to a file instead of back into the
database to allow for maximum flexibility. 

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
utilized thoughout the entire DETERMINE_GENE_FAMILIES program.

=cut

use warnings;
use strict;

package Input_Variables;   

use vars qw
    ($new_gene_cutoff
    $output_file_name
    $initial_family_number
    
    $mysql_extraction_query

    $database_host_ip_address
    $mysql_account_username
    $mysql_account_password
    $mysql_database_name
    );

# Default Program Variables:
# _________________________
    
$new_gene_cutoff = 2;
$output_file_name = 'gene_names_families_and_ages.tab';
$initial_family_number = 1;

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT GeneUniqueID, Phylostratum, GeneFamilyNumber FROM TableName'; 

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$new_gene_cutoff = This is the minimum number of new genes are that required
for a group of new genes to be considered their own family. That is, a
single homologous old gene will be assigned the phylostratum of these new
genes.

$output_file_name = The name of the output file that contains the 
gene unique ID, the family number, and the family phylostratum.

$initial_family_number = This number will be the lowest possible family
number that the user wants to insert into the database. If this variable
did not exist, the gene families would begin at 0 or 1 everytime this
program is run. This is a problem if someone wishes to conglomerate two
sets of genes together. The combined set would then have two families
labelled as "1", two labelled as "2", etc, even though they are really two
distinct families. 

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
information to be extracted. 

=cut

package main;

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

# Lexical Variables:

my  ($input_file_line,
    $gene_uid,
    $paralog_list_ds, # "ds" stands for "data strcture"

    $database_connection_object,
    $database_extraction_object,
    $database_insertion_object,
    $row_array,
    
    $db_primary_key,
    $gene_phylostratum,
    $gene_family_number,
    $genes_in_family_temp_hash,
    %phylostratum_per_gene_ds,
    $genes_in_family_ds,
    
    %temp_phylostrata_in_family,
    @phylostrata_in_family,
    $completed_families_ds,
    $old_phylostratum,
    $new_phylostratum,
    @old_gene_names,
    @new_gene_names,
    $new_gene_id,
    $old_gene_id,
    %new_gene_homologies,
    $homology_fraction_numerator,
    $new_gene_homology_fraction,
    $index_to_add,
    
    $family_number_for_db,
    );

# Statistical Counters:

my  ($transferred_to_complete_counter,
    $not_transferred_to_complete_counter, 
    $new_gene_cutoff_counter,
    $old_gene_is_one_counter,
    $genes_remain_as_single_family_counter,
    $genes_split_into_two_families_counter,    
    %status_hash,
    );

# The %status_hash contains gene names as keys and
# abbreviations representing the status of how that gene
# was filed. The abbreviations are as follows:

# A = Family has only one phylostratum (all genes have this; this is for printer purposes)
# N = New phylostratum only had one gene
# O = Old phylostratum only had one gene
# R = Family remained as a single family
# S = Family split into two families

# Families with more than two total phylostrata will have
# multiple abbreviations. Additionally, only gene families
# that begin having multiple phylostrata are included in
# this hash.

$transferred_to_complete_counter = 0;
$not_transferred_to_complete_counter = 0;
$new_gene_cutoff_counter = 0;
$old_gene_is_one_counter = 0;
$genes_remain_as_single_family_counter = 0;
$genes_split_into_two_families_counter = 0;

#####################################################################
# INPUT #############################################################
#####################################################################

# We will first read the input text file and generate a
# data structure for the homologous genes.

print 'Inputting text file...';

while (defined( $input_file_line = <> ))
    {
    chomp $input_file_line;
    
    # In the following regular expression, the "?:" after the opening parenthesis
    # tells the perl interpreter not to count this set of parentheses when saving
    # to the memory variables (i.e., $1, $2, etc.). Furthermore, having the
    # "?<LABEL>" after the opening parenthesis adds a label to the expression.
    # The label is stored as a hash key in the "$+" hash. Each key's value is
    # the matching pattern inside the parentheses of the label.

    $input_file_line =~ m/^(?<ref_gene_uid>\S+)\t(?<paralog_uid>\S+)\t/gi;
    
    $gene_uid = $+{ref_gene_uid};    

    push( @{ $paralog_list_ds->{$gene_uid} }, $+{paralog_uid});
    };

print "Complete\n";

# Next we must access the MySQL database, extract data from it,
# and use this data to generate two data structures. One data
# structure will be a hash of GeneUniqueIDs and their phylostrata,
# the other will be a hash containing the gene family number and
# all the genes (labelled by GeneUniqueID) currently within that
# family. 

use DBI;

print 'Connecting to host database = ' . "$Input_Variables::database_host_ip_address\n";

$database_connection_object = DBI->connect("dbi:mysql:$Input_Variables::mysql_database_name;$Input_Variables::database_host_ip_address",
                                           $Input_Variables::mysql_account_username, $Input_Variables::mysql_account_password);

print "Database connected\n";

# We now pre-prepare our MySQL commands that will be used to both
# extract and insert data.

$database_extraction_object = $database_connection_object->prepare($Input_Variables::mysql_extraction_query);

$database_extraction_object->execute();

print "Database extraction command executed\n";

while ($row_array = $database_extraction_object->fetchrow_arrayref)
    {    
    unless ( defined( $row_array->[4] ) ) {next};
    
    chomp @$row_array;
    
    # The $row_array is a reference to an array, each element of which
    # is a value from a row of the database. In this case, the zeroth
    # element is the GeneUniqueID, the first element is the
    # Phylostratum, and the second element is the GeneFamilyNumber.

    $gene_uid = $row_array->[0];
    $gene_phylostratum = $row_array->[1];
    $gene_family_number = $row_array->[2];
    
    undef $row_array;
    
    print "Input: Gene: $gene_uid\n";

    push( @{$genes_in_family_temp_hash->{$gene_family_number}}, $gene_uid );
    $phylostratum_per_gene_ds{$gene_uid} = $gene_phylostratum;
    $status_hash{$gene_uid} = 'A';
    };

print "Database data extraction complete\n";

$database_connection_object->disconnect;

no DBI;  

print "Database disconnected\n";

# We must now turn the $genes_in_family_ds into an array-of-arrays.
# The index of the initial array will be the new family number. This
# will allow us to consolidate the family numbers and easily add and
# remove families from the array.

foreach (values %$genes_in_family_temp_hash)
    {
    push(@$genes_in_family_ds, $_); 
    };

undef %$genes_in_family_temp_hash;

print 'Gene families extracted: ' . scalar(@$genes_in_family_ds) . "\n";

#####################################################################
# ALGORITHM #########################################################
#####################################################################
    
# Now that the input is finished, we can modify the appropriate
# gene family numbers.

print "Calculating the family number and phylostratum\n\n";

while ( scalar(@$genes_in_family_ds) > 0 )
    {    
    # We will begin with the last family in the array (i.e., the
    # last element of the array). When each family
    # is assigned to only one phylostratum, the family will be moved
    # from the $genes_in_family_ds to the $completed_families_ds. 
    
    # We must order all the phylostrata in this gene family.
    
    foreach ( @{$genes_in_family_ds->[-1]} )
        {
        # The intermediate step of creating the %temp_phylostrata_in_family
        # hash will remove duplicate phylostrata.
        
        $temp_phylostrata_in_family{ $phylostratum_per_gene_ds{$_} } = 0;     
        };
    
    undef @phylostrata_in_family;
    
    @phylostrata_in_family = (keys %temp_phylostrata_in_family);    
    
    undef %temp_phylostrata_in_family;    
    
    # If all genes in the family are the same age, the family number does not change.
    
    if ( scalar(@phylostrata_in_family) == 1 )
        {
        # We will move this gene family into the $completed_families_ds
        # data structure.
        
        push( @$completed_families_ds, pop(@$genes_in_family_ds) );        
        
        $transferred_to_complete_counter++;
        
        next;
        };    
    
    $not_transferred_to_complete_counter++;
    
    # It is possible to have more than two phylostrata representing each
    # gene family. If this occurs, we will assess the two oldest phylostrata
    # first, condense them, and then iterate to the newer phylostrata.
    # Eventually, we will condense the @phylostrata_in_family array
    # down to a single number, and this will be the phylostrata for the
    # entire family. If a split occurs, we will add the family 
    
    # The following "sort" command will sort the phylostrata in ascending order.

    @phylostrata_in_family = sort { $a <=> $b } @phylostrata_in_family;
    
    # We must now determine which genes are in the oldest two phylostrata.
    # Those in the oldest phylostratum will be called "old" and those in
    # the second oldest phylostratum will be called "new".
    
    $old_phylostratum = $phylostrata_in_family[0];
    $new_phylostratum = $phylostrata_in_family[1];
    
    undef @old_gene_names;
    undef @new_gene_names;
    
    foreach ( @{$genes_in_family_ds->[-1]} )
        {
        if ( $phylostratum_per_gene_ds{$_} == $old_phylostratum )
            {
            # Old genes
            
            push(@old_gene_names, $_);    
            }
            
        elsif ( $phylostratum_per_gene_ds{$_} == $new_phylostratum )
            {
            # New genes
            
            push(@new_gene_names, $_);    
            };
        };    

    # We will now perform the algorithms that determine to which
    # phylostrata each gene belongs.
    
    # If less than X number of new genes exist (where X is the
    # $new_gene_cutoff), they are considered homologous with the
    # old gene(s) and they will obtain the phylostratum of the old
    # gene, and the family will be reassessed with these updated
    # phylostrata.  
    
    if ( scalar(@new_gene_names) < $Input_Variables::new_gene_cutoff )
        {
        foreach (@old_gene_names) { $status_hash{$_} = $status_hash{$_} . 'N' };
        foreach (@new_gene_names) { $status_hash{$_} = $status_hash{$_} . 'N' };
        
        foreach (@new_gene_names)
            {
            $phylostratum_per_gene_ds{$_} = $old_phylostratum;  
            };
        
        $new_gene_cutoff_counter++;
        
        next;    
        };
    
    # To prevent possible errors, the old genes must have
    # more than one gene. A single gene in the old genes
    # group is considered to be an error and is therefore
    # given the phylostratum of the new genes. The gene
    # family is then reassessed using the updated phylostrata.
    
    if ( scalar(@old_gene_names) < 2 )
        {
        foreach (@old_gene_names) { $status_hash{$_} = $status_hash{$_} . 'O' };
        foreach (@new_gene_names) { $status_hash{$_} = $status_hash{$_} . 'O' };
        
        foreach (@old_gene_names)
            {
            $phylostratum_per_gene_ds{$_} = $new_phylostratum;  
            };
        
        $old_gene_is_one_counter++;
        
        next;    
        };
    
    # If the new and old genes both meet their quantity cutoffs,
    # then at least 50% of the new genes must be homologous with
    # at least one of the old genes in order for the family to
    # remain intact. If this goal is met, the new genes will be
    # assigned the phylostratum of the old genes. If fewer than
    # 50% of the new genes are homologous to at least one of the
    # old genes, then the gene family is split into two families.
    # Both of these families will be assessed.
    
    undef %new_gene_homologies;
    
    foreach $new_gene_id (@new_gene_names)
        {
        # In the %new_gene_homologies hash, a "0" means that this
        # gene is not homologous to any of the old genes. A "1"
        # means that a homology is found.  
        
        $new_gene_homologies{$new_gene_id} = 0;        
        
        foreach $old_gene_id (@old_gene_names)
            {
            foreach ( @{$paralog_list_ds->{$new_gene_id}} )
                {
                if ($old_gene_id eq $_)
                    {
                    # The new gene and the old gene are paralog matches.
                    
                    $new_gene_homologies{$new_gene_id} = 1;
                    
                    last; # For speed; only one match will exist in the $paralog_list_ds.
                    };
                };
            };
        };
    
    $homology_fraction_numerator = 0;
    
    foreach (values %new_gene_homologies) { $homology_fraction_numerator += $_ };
    
    $new_gene_homology_fraction = $homology_fraction_numerator / scalar(@new_gene_names);
    
    if ( $new_gene_homology_fraction >= 0.5 )
        {
        # The genes remain as a single family.    
        
        $genes_remain_as_single_family_counter++;
        
        foreach (@old_gene_names) { $status_hash{$_} = $status_hash{$_} . 'R' };
        foreach (@new_gene_names) { $status_hash{$_} = $status_hash{$_} . 'R' };
        
        foreach (@new_gene_names)
            {
            $phylostratum_per_gene_ds{$_} = $old_phylostratum;  
            };
        }
        
    elsif ( $new_gene_homology_fraction < 0.5 )
        {
        # The genes split into two families.    
        
        $genes_split_into_two_families_counter++;
        
        foreach (@old_gene_names) { $status_hash{$_} = $status_hash{$_} . 'S' };
        foreach (@new_gene_names) { $status_hash{$_} = $status_hash{$_} . 'S' };

        # We cannot just use the @new_gene_names because this array does
        # not contain any genes in the gene family that are newer than
        # those in the @new_gene_names.  
        
        undef @new_gene_names; 
        
        foreach ( @{$genes_in_family_ds->[-1]} )
            {
            unless ( $phylostratum_per_gene_ds{$_} == $old_phylostratum )
                {
                push(@new_gene_names, $_);      
                };
            };
            
        pop(@$genes_in_family_ds);
        
        # We cannot just add the \@old_gene_names reference to the
        # @$genes_in_family_ds data structure because anything
        # that influences the \@old_gene_names reference (such as
        # an "undo" command) will affect what is in the
        # $genes_in_family_ds data structure. Therefore, we must add
        # each gene name individually.
        
        undef $index_to_add;
        
        $index_to_add = scalar(@$genes_in_family_ds);
        
        foreach (@old_gene_names)
            {
            push(@{$genes_in_family_ds->[$index_to_add]}, $_);    
            };
            
        $index_to_add = scalar(@$genes_in_family_ds);
        
        foreach (@new_gene_names)
            {
            push(@{$genes_in_family_ds->[$index_to_add]}, $_);    
            };
        };    
    };
    
#####################################################################
# OUTPUT ############################################################
#####################################################################    

# First we will print the various statistics to STDOUT.

print "\nTransferred to complete: $transferred_to_complete_counter\n";
print "Not transferred to complete: $not_transferred_to_complete_counter\n\n";
print "New gene cutoff: $new_gene_cutoff_counter\n";
print "Old genes are one: $old_gene_is_one_counter\n\n";
print "Genes remain as single family: $genes_remain_as_single_family_counter\n";
print "Genes split into two families: $genes_split_into_two_families_counter\n\n";
print 'Size of $completed_families_ds: ' . scalar(@$completed_families_ds) . "\n";    
    
# We will now output the GeneUniqueID, the GeneFamilyNumber, and
# the GeneFamilyPhylostratum into a tab-delineated file. 

print "Printing to output file\n";
    
open(GENES_FAMILIES_AGES_OUTPUT, ">", $Input_Variables::output_file_name) || die("Cannot open the $Input_Variables::output_file_name file: $!\n");

foreach $gene_family_number ( 0..( scalar(@$completed_families_ds) - 1  ) )
    {    
    # To avoid having a family number of zero, we will add one to
    # the $gene_family_number.
    
    $family_number_for_db = $gene_family_number + $Input_Variables::initial_family_number;
    
    foreach $gene_uid ( @{$completed_families_ds->[$gene_family_number]} )
        {    
        $gene_phylostratum = $phylostratum_per_gene_ds{$gene_uid};
        
        print GENES_FAMILIES_AGES_OUTPUT "$gene_uid\t$family_number_for_db\t$gene_phylostratum\t" . $status_hash{$gene_uid} . "\n";
        };
    };

close (GENES_FAMILIES_AGES_OUTPUT);

print "\nProgram is finished\n";

exit; 
