#!/usr/local/bin/perl

use warnings;
use strict;

##################################################################
=head1 # PAIRWISE HOMOLOG PHYLOSTRATA CLUSTERING #################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will take a pairwise list of homologous genes from various
phylostrata and place them into paralogous families (i.e, clusters).
The output from this program is a tab-delineated file. The first column
is the gene name, the second is the original phylostratum of the gene,
the third is the phylostratum of the group, and the fourth is the group
number. Group numbers are unordered but cumulative. They are not repeated
over differing phylostrata.  

=head2 INPUT:

Bash-style input. The name of the input file will be listed after
the executable. No argument prefix required. The input file should
be tab-delineated with the first column being a UID of the reference
gene, the second column being the UID of a gene that is a paralog of
the reference gene, and the third column should be the phylostratum of
the reference gene. Matches are one-way; each paralog gene must also
be listed as a reference gene somewhere in the input file.

=head2 REQUIREMENTS:

None.

=head2 ARGUMENTS:

None.

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($input_file_line,
    %reference_phylostrata,
    $paralogs,    
    $reference_gene,
    $paralog_gene,
    $ref_gene_phylostratum,
    
    $cluster_counter,
    %gene_clusters,
    
    %original_cluster_numbers,
    %new_cluster_numbers,
    %cluster_phylostrata,

    $cluster,
    
    $output_file_name,        
    $cluster_phylostratum,
    );

# First we must read data from the input file.

while (defined( $input_file_line = <> ))
    {
    chomp $input_file_line;
    
    # In the following regular expression, the "?:" after the opening parenthesis
    # tells the perl interpreter not to count this set of parentheses when saving
    # to the memory variables (i.e., $1, $2, etc.). Furthermore, having the
    # "?<LABEL>" after the opening parenthesis adds a label to the expression.
    # The label is stored as a hash key in the "$+" hash. Each key's value is
    # the matching pattern inside the parentheses of the label.

    $input_file_line =~ m/^(?<ref_gene_uid>\S+)\t(?<paralog_uid>\S+)\t(?<ref_gene_phylostrat>\d+)$/gi;
    
    $reference_gene = $+{ref_gene_uid};
    $paralog_gene = $+{paralog_uid};
    $ref_gene_phylostratum = $+{ref_gene_phylostrat};
    
    # We will now make sure all data is present. To do this
    # we must ensure the input scalars are defined and that
    # they approximately fit the characters they are supposed to.
    
    unless ($reference_gene) {next};
    unless ($paralog_gene) {next};
    unless ($ref_gene_phylostratum) {next}; 
    
    unless ( $reference_gene =~ m/(\w|\d)/g ) {next};
    unless ( $paralog_gene =~ m/(\w|\d)/g ) {next};
    unless ( $ref_gene_phylostratum =~ m/^\d+$/ ) {next};    
    
    $paralogs->{ $reference_gene }->{ $paralog_gene } = 0;
    $reference_phylostrata{ $reference_gene } = $ref_gene_phylostratum;
    };    

# As a final check, we must make sure that each paralogous gene
# has a known phylostratum.

foreach $reference_gene (keys %$paralogs)
    {
    foreach $paralog_gene ( keys %{$paralogs->{$reference_gene}} )
        {
        unless ( $reference_phylostrata{$paralog_gene} )
            {
            delete $paralogs->{$reference_gene}->{$paralog_gene};
            };
        };
    };

# We must now fill-in the cluster names.

$cluster_counter = 0;

foreach (keys %$paralogs)
    {
    $cluster_counter++;
    
    $gene_clusters{$_} = $cluster_counter;
    };

# We can now generate clusters for each group of paralogous
# genes.

foreach $reference_gene (keys %$paralogs)
    {
    foreach $paralog_gene ( keys %{$paralogs->{$reference_gene}} )
        {
        $gene_clusters{$paralog_gene} = $gene_clusters{$reference_gene};
        };
    };

# Although each cluster now has a number, these numbers are
# not consecutive. We will now make each cluster number
# consecutive.

foreach (values %gene_clusters)
    {
    $original_cluster_numbers{$_} = 0;    
    };

$cluster_counter = 0;

foreach (keys %original_cluster_numbers)
    {
    $cluster_counter++;
    
    $new_cluster_numbers{$_} = $cluster_counter;     
    
    # The %cluster_phylostrata hash is used in a
    # subsequent section.
    
    $cluster_phylostrata{$cluster_counter} = 999999;
    };

undef %original_cluster_numbers;

foreach (keys %gene_clusters)
    {
    $gene_clusters{$_} = $new_cluster_numbers{ $gene_clusters{$_} };
    };

undef %new_cluster_numbers;

# We will now calculate the minimum phylostratum for
# all the genes in a given cluster. This is the
# phylostratum number for for each cluster.

while (($reference_gene, $ref_gene_phylostratum) = each %reference_phylostrata)
    {
    $cluster = $gene_clusters{$reference_gene};
    
    if ( $ref_gene_phylostratum < $cluster_phylostrata{$cluster} )
        {
        $cluster_phylostrata{$cluster} = $ref_gene_phylostratum; 
        };
    };

# We will now print output for the information we have calculated.

$output_file_name = 'gene_clusters_output.tab';

open(GENE_OUTPUT, ">", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

while (($reference_gene, $ref_gene_phylostratum) = each %reference_phylostrata)
    {
    $cluster = $gene_clusters{$reference_gene};
    $cluster_phylostratum = $cluster_phylostrata{$cluster};
    
    print GENE_OUTPUT "$reference_gene\t$ref_gene_phylostratum\t$cluster_phylostratum\t$cluster\n";
    };

close GENE_OUTPUT;

exit;