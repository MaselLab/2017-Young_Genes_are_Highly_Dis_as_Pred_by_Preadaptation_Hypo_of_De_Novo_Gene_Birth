#!/usr/local/bin/perl

##################################################################
=head1 # INTERGENIC CONTROL SEQUENCE GENERATOR ###################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will extract intergenic sequences from a genome using a
list of proteins associated with the genome. Each extracted intergenic
sequence is the same length as an associated protein. Multiple intergenic
sequences can be extracted per protein. A buffer of nucleotides separating
the genes and the extracted intergenic regions is also included to
prevent the extraction of a promoter. Each extracted intergenic sequence
is extracted as close as possible to its associated protein's gene.
This will limit the effect of chromosomal location differences between
the protein and the intergenic control sequence. 

Although the length of the extracted intergenic sequences is based upon
an input list of proteins, a separate list of gene loci is input into
the program to determine the location of the genes on the genome, thus
giving the program the location of the intergenic regions.

As a sequence begins to be extracted from the genome, this sequence has a
frame. Any stop codon (given the frame of the sequence) within the extracted
intergenic sequence is skipped. The length of the extracted sequence is
extended appropriately until the final sequence that is extracted is the
appropriate length given the associated query protein sequence. In addition
to skipping stop codons, this program also skips any "N"s or "X"s in the
intergenic regions. Although "X" is not a standard nucleotide designation,
it can be used as an arbitrary nucleotide designation that the user wishes
to skip. For example, when using RepeatMasker on a genome before inputting
it into this program, the masked nucleotides can be given an "X" designation.

=head2 INPUT:

This program requires three input files:

1)  FASTA file for the genome. The title line of the FASTA file must contain
    the chromosome. The regular expression used to extract the chromosome
    from the title line can be found in the "extract_chromosome_sequence"
    subroutine.

2)  File containing the proteins for which we want the intergenic controls.
    This is a tab-delineated file containing the protein name in the first
    column, the chromosome of the protein in the second column, and the
    protein's sequence in the final column.
    
3)  File containing the genes, enabling the program to mark the genic and
    intergenic regions of the genome. This is a tab-delineated file with
    four columns: the first column is the gene name, the second column is
    the chromosome, and the third and fourth columns contain the start and
    stop locus respectively. Importantly, the loci contained in this column
    must match to the genome FASTA file.
    
Importantly, all gene/protein names and chromosomes must match between the
three aforementioned input files. For example, the names of the proteins
from the protein input file must match the gene names in the gene input file.
The chromosomes between all three files must also match.  

=head2 REQUIREMENTS:

None.

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

The hardcoded values below are the default values for each variable.
Additionally, below these are descriptions of each variable.

=head2 ARGUMENT VALUES

The following input variables are lexical and can be utilized
thoughout the entire "Intergenic Control Sequence Generator" program.

=cut

use warnings;
use strict;

print "Time of program execution:\n";
main::print_time();
print "\n";

my  ($protein_input_file_name,
    $gene_input_file_name,
    $genome_input_file_name,
    $output_file_name,
    
    $control_seq_per_gene,
    $buffer,
    );

$genome_input_file_name = '/PATH/GENOME_FILE_NAME.fasta';
$protein_input_file_name = '/PATH/PROTEIN_FILE_NAME.fasta.tab';
$gene_input_file_name = '/PATH/GENE_FILE_NAME.fasta..tab';

$output_file_name = '/PATH/OUTPUT_FILE_NAME.tab';
$control_seq_per_gene = 2;
$buffer = 1000; # Number of nucleotides.

=head2 ARGUMENT DESCRIPTIONS

$genome_input_file_name = File name of the genome FASTA file. Input #1
from the INPUT section. 

$protein_input_file_name = File name of the tab-delineated file containing
the protein information. Input #2 from the INPUT section. 

$gene_input_file_name = File name of the tab-delineated file containing
the gene information. Input #3 from the INPUT section. 

$output_file_name = Name and path of the output file containing the intergenic
sequences. This is a tab-delineated file where the first column is the name of
the protein associated with the intergenic sequence, the second column is
the chromosome from where the intergenic sequence was extracted, and the
third column is the nucleotide sequence of the extracted intergenic sequence.

$control_seq_per_gene = The number of control intergenic sequences that are
extracted per protein.

$buffer = The number of nucleotides separating the genic and intergenic regions
in the genome. Intergenic sequences will not be extracted from the buffer regions.

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($query_proteins_ds,
    $all_genes_ds,
    %chromosome_for_each_gene,
    %chromosomes_hash,
    @chromosomes,

    $chromosome,
    $chromosome_sequence,
    $intergenic_sequences,
    
    $gene_uid,
    $orf_length,
    $single_intergenic_seq,
    $codon,
    $control_sequence,
    );

# First we must read data from the protein and gene
# input files.

$query_proteins_ds = &protein_input($protein_input_file_name);
$all_genes_ds = &gene_input($gene_input_file_name);

# Now we must order the query proteins by the length of
# their associated nucleotides sequences.

# The following sort command will sort the lengths in ascending order.

@$query_proteins_ds = sort
    {
    if ($a->{'orf_length'} < $b->{'orf_length'}) {-1} 

    elsif ($a->{'orf_length'} > $b->{'orf_length'}) {1} 

    else {0};
    
    } @$query_proteins_ds;

# We must generate a list of chromosomes.

foreach (@$query_proteins_ds)
    {
    $chromosomes_hash{ $_->{'chromosome'} } = 0;    
    };

@chromosomes = keys %chromosomes_hash;

undef %chromosomes_hash;

# We can now look into each chromosome to generate control sequences.

open(CONTROL_OUTPUT, ">", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

foreach $chromosome (@chromosomes) 
    {
    # For each chromosome, we must first extract the
    # chromosome sequence from the genome file.
    
    print "\nExtracting sequence for Chromosome $chromosome...\n";
    
    $chromosome_sequence = &extract_chromosome_sequence($chromosome);
    
    print "Sequence extracted.\n";    

    # It is now time to generate an array of intergenic sequences.
    
    print "Generating intergenic sequences...\n";
    
    $intergenic_sequences = &extract_intergenic_sequences($chromosome, $chromosome_sequence, $all_genes_ds, $buffer);

    print "Intergenic sequences generated.\n";
    
    undef $chromosome_sequence;
    
    # We must now sort the intergenic sequences from smallest to largest.
    # This will allow us to maximize the usage of the intergenic sequences.
    # That is, it will allow us to get the most control sequences possible.

    # The following sort command will sort the lengths in ascending order.    
    
    @$intergenic_sequences = sort
        {
        if ( length($a) < length($b) ) {-1} 
    
        elsif ( length($a) > length($b) ) {1} 
    
        else {0};
        
        } @$intergenic_sequences;

    # We can now match the query protein sequences to the
    # intergenic sequences.
    
    foreach (@$query_proteins_ds)
        {
        if ( $_->{'chromosome'} =~ m/^$chromosome$/i )
            {
            $gene_uid = $_->{'gene_uid'};
            $orf_length = $_->{'orf_length'};            
            
            foreach (1..$control_seq_per_gene)
                {
                # For the following loop, if $_ is the
                # sequence, then we won't actually be modifying the
                # array. Therefore, we must make $_ equal to the index.
                
                INTERGENIC: foreach (0..$#{$intergenic_sequences})
                    {
                    $single_intergenic_seq = $intergenic_sequences->[$_];
                    
                    if ( $orf_length > length($single_intergenic_seq) ) {next};

                    # We must first do a check to make sure the ORF length
                    # can be divided by 3.
                    
                    unless( $orf_length % 3 == 0 )
                        {
                        die("The ORF sequence cannot be divided by 3");    
                        };                    
                    
                    # We can now extract the control sequence one codon at a
                    # time. We must do it this way because we must also remove
                    # any stop codons in the sequence. If a stop codon is found,
                    # the code will remove it and choose the next three
                    # nucleotides from the genome and add them on to the 5' end
                    # (i.e., it will extend the genome to keep the length
                    # consistent). Although it is possible for a control
                    # sequence to extend past the start locus of the next gene
                    # (after stop codons are removed), the extentions should be
                    # short enough that this will not be too much of a
                    # problem--especially with the nucleotide buffer.  
                    
                    $control_sequence = '';
                    
                    until ( length($control_sequence) == $orf_length )
                        {
                        # If stop codons are present, one might have to extend
                        # the codons past the end of the $single_intergenic_seq.
                        # To prevent this from becoming a problem, we will name
                        # the previous loop "INTERGENIC" and permit the next
                        # element in it using the following conditional statement.                        
                        
                        if ( length($single_intergenic_seq) < 3) {next INTERGENIC};                        
                        
                        $codon = substr($single_intergenic_seq, 0, 3, '');
                        
                        unless ($codon =~ m/TAG|TAA|TGA|X|N/gi)
                            {
                            $control_sequence = $control_sequence . $codon;
                            };
                        
                        undef $codon;
                        };
 
                    print CONTROL_OUTPUT "$gene_uid\t$chromosome\t$control_sequence\n";
                    
                    undef $control_sequence;
 
                    # We must now return the new $single_intergenic_seq
                    # to the @$intergenic_sequences array.
                    
                    $intergenic_sequences->[$_] = $single_intergenic_seq;
                    
                    last;
                    };
                };
            };
        };
    
    undef $intergenic_sequences;
    };

close CONTROL_OUTPUT;

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

# PROTEIN_INPUT

# The following subroutine will input the protein's UID and the
# amino acid protein sequence associated with it. 

sub protein_input
    {
    
    # Input Variables
    
    my  ($protein_input_file_name) = @_;
    
    # Subroutine Variables
    
    my  ($query_proteins_ds_ref,
        $query_proteins_ds,
        );

    open(PROTEIN_FILE_INPUT, "<", $protein_input_file_name) || die("Cannot open the $protein_input_file_name file: $!\n");
    
    foreach (<PROTEIN_FILE_INPUT>)
        {
        chomp $_;
        
        $_ =~ m/^(\S+)\t(\S+)\t(\w+)/;
        
        $query_proteins_ds_ref =
            {            
            'gene_uid' => $1,            
            'chromosome' => $2,
            'orf_length' => ( length($3) ) * 3,
            };        
        
        push(@$query_proteins_ds, $query_proteins_ds_ref);        
        
        undef $query_proteins_ds_ref;
        };
    
    close PROTEIN_FILE_INPUT;
    
    return $query_proteins_ds;

    };
    
########################################################################################    

# GENE_INPUT

# The following subroutine will input the gene's UID, the
# chromosome, and the start and stop loci of the genes on the genome.  

sub gene_input
    {
    
    # Input Variables
    
    my  ($protein_input_file_name) = @_;
    
    # Subroutine Variables
    
    my  ($chromosomes,
        $gene_start_loci,
        $gene_stop_loci,
        $all_genes_ds_ref,
        $all_genes_ds,
        );

    open(GENE_FILE_INPUT, "<", $gene_input_file_name) || die("Cannot open the $gene_input_file_name file: $!\n");
    
    foreach (<GENE_FILE_INPUT>)
        {
        chomp $_;
        
        $_ =~ m/^(\S+)\t(\S+)\t(\S+)\t(\S+)$/;
        
        $all_genes_ds_ref =
            {'gene_uid'     => $1,
            'chromosome'    => $2,
            'start_locus'   => $3,
            'stop_locus'    => $4,
            };
        
        push(@$all_genes_ds, $all_genes_ds_ref);
        
        undef $all_genes_ds_ref;
        };
    
    close GENE_FILE_INPUT;
    
    return $all_genes_ds;

    };

########################################################################################

# EXTRACT_CHROMOSOME_SEQUENCE

# The following subroutine will extract a single
# chromosomal sequence from the genome FASTA file.
# Importantly, this information cannot be put into
# a hash because it requires too much memory. Therefore,
# we will extract each chromosome's sequence directly out
# of the genome file as it is needed. 

sub extract_chromosome_sequence
    {
    
    # Input Variables
    
    my  ($query_chromosome) = @_;
    
    # Subroutine Variables
    
    my  ($chromosome_sequence,
        @chromosome_sequence_array,
        $current_chromosome,
        );
    
    $current_chromosome = '99999';
    
    open(GENOME_FILE_INPUT, "<", $genome_input_file_name) || die("Cannot open the $genome_input_file_name file: $!\n");
    
    foreach (<GENOME_FILE_INPUT>)
        {
        chomp $_;
        
        if ($_ =~ m/^\>(\S+)\sdna/)
            {
            # Title line.
            
            if ($current_chromosome =~ m/^$query_chromosome$/i)
                {
                # This loop will be entered on the chromosome
                # that follows our query chromosome.
                
                $chromosome_sequence = join('', @chromosome_sequence_array);
                
                undef @chromosome_sequence_array;                            
                
                close GENOME_FILE_INPUT;
    
                return $chromosome_sequence;
                };
            
            undef $current_chromosome;
            
            $current_chromosome = $1;
            }

        elsif ($_ =~ m/^[a-zA-Z]/)
            {
            # Sequence.
            
            if ($current_chromosome =~ m/^$query_chromosome$/i)
                {
                # This sequence should be saved.
                
                push(@chromosome_sequence_array, $_);
                };            
            };
        };
    
    undef @chromosome_sequence_array;                            

    close GENOME_FILE_INPUT;

    # If for some reason the query chromosome cannot be found
    # within the genome, we will use Chromosome 1 to obtain the
    # control sequences. Although these control sequences will
    # be duplicates (except for their lengths), this should not
    # occur often enough to influence any results obtained by
    # these controls. To obtain the sequence of Chromosome 1,
    # we will use the current subroutine (&extract_chromosome_sequence)
    # recursively.
    
    print "The chromosome: $query_chromosome was not found in the genome file.\n";
    
    $chromosome_sequence = &extract_chromosome_sequence('1');
    
    return $chromosome_sequence;    
    
    };

########################################################################################

# EXTRACT_INTERGENIC_SEQUENCES

# The following subroutine extracts intergenic sequences
# from each chromosome. These sequences are stored in an
# array called @intergenic_sequences. Importantly, the
# sequences range from the end of the most 3' gene to the
# beginning of the most 5' gene.

sub extract_intergenic_sequences
    {
    
    # Input Variables
    
    my  ($chromosome,
        $chromosome_sequence,
        $all_genes_ds,
        $buffer,
        ) = @_;
    
    # Subroutine Variables
    
    my  (@chromosome_nucleotides,
        $single_intergenic_sequence,
        @intergenic_sequences,
        );

    # We must first create an array containing all the nucleotides
    # in the chromosome sequences.
    
    print "Generating chromosome array\n";
    
    @chromosome_nucleotides = split('', $chromosome_sequence);
    
    # For ease, we will name the zeroth element "zero".
    # Therefore, index 1 will represent the 3' end of the
    # chromosome sequences.
    
    unshift(@chromosome_nucleotides, 'zero');

    # Now we can label all the gene positions in the array.

    foreach (@$all_genes_ds)
        {
        if ( $_->{'chromosome'} =~ m/^$chromosome$/i )
            {
            foreach ( ($_->{'start_locus'})..($_->{'stop_locus'}) )
                {
                $chromosome_nucleotides[$_] = 'gene';
                };
            };
        };

    # We will now label the 3' end of the chromosome "start" and
    # the 5' end of the chromosome "end".
    
    foreach (@chromosome_nucleotides)
        {
        if ($_ =~ m/gene/i) {last}

        elsif ($_ =~ m/^[a-zA-Z]$/i) { $_ = 'start' }
        
        elsif ($_ =~ m/zero/i) {next}
        };
    
    foreach (reverse @chromosome_nucleotides)
        {
        if ($_ =~ m/gene/i) {last}

        elsif ($_ =~ m/^[a-zA-Z]$/i) { $_ = 'end' };
        };
    
    reverse @chromosome_nucleotides;

    print "Chromosome array generated\n";
    
    # We can now begin to build the @intergenic_sequences array.
    
    undef $single_intergenic_sequence;
    
    foreach (@chromosome_nucleotides)
        {
        if ($_ =~ m/^[a-zA-Z]$/i)
            {
            $single_intergenic_sequence = $single_intergenic_sequence . $_;
            }
            
        else
            {
            push(@intergenic_sequences, $single_intergenic_sequence);    
            
            undef $single_intergenic_sequence;
            };
        };    
    
    undef @chromosome_nucleotides;
    
    # Our @intergenic_sequences array will be filled with empty
    # elements; we must now delete these elements.
    
    @intergenic_sequences = grep { defined } @intergenic_sequences;
    
    # We must now give a "buffer zone" of $buffer nucleotides
    # on both sides of each gene. Note that if $_ is the
    # sequence, then we won't actually be modifying the
    # array. Therefore, we must make $_ equal to the index.
    
    foreach (0..$#intergenic_sequences)
        {
        # First we must make sure the intergenic sequence is long
        # enough.
        
        if ( length($intergenic_sequences[$_]) <= (2 * $buffer) )
            {
            undef($intergenic_sequences[$_]);
            next;
            };
        
        substr( $intergenic_sequences[$_], 0, $buffer, '' ); # Truncates the $buffer length of 3' nucleotides.
        substr( $intergenic_sequences[$_], (-1 * $buffer), $buffer, '' ); # Truncates the $buffer length of 5' nucleotides.
        };

    @intergenic_sequences = grep { defined } @intergenic_sequences;

    return \@intergenic_sequences;
    
    };

########################################################################################





