#!/usr/bin/perl 
use strict;
use warnings;

my $dir=$ARGV[0];
my $uniprot_rep_path=$ARGV[1];
my $trinotate_rep_path=$ARGV[2];
my $node_protein_path=$ARGV[3];

# USAGE example:
# for DIR in *files
# do 
# 	perl $0 $DIR ${DIR}/${DIR/_files}_uniprot_ann.txt \
# 	     ggallus_transc_Trinotate_report.tsv \
#	     ${DIR}/${DIR/_files}_node_and_protID.txt
# done

# Open files
open my $tri_out, ">", "$dir/trinotate_report_table.tsv" or die "Cannot open output file: $!";
open my $ua_fh,   "<", "$uniprot_rep_path" or die "Cannot open uniprot annotation file: $!";
open my $tr_fh,   "<", "$trinotate_rep_path" or die "Cannot open trinotate annotation file: $!";
open my $n2p_fh,  "<", "$node_protein_path" or die "Cannot open node to protein file: $!";

### Fill hashes ### 
# uniprot annotations
my %uniprot;
while (my $line = <$ua_fh>){
        chomp $line;
        my @temp = split /\t/, $line;
        my $id = $temp[0];
        unless (exists $uniprot{$id}){
                $uniprot{$id} = $temp[1];
        }
}
close $ua_fh;

# Trinotate results
my %report;
<$tr_fh>;
while (my $line = <$tr_fh>){
        chomp $line;
        my @temp = split /\t/, $line;
        my $id = $temp[0];
        my @blast_hit = split '\^', $temp[2];
        my $prot_hit = $blast_hit[0];
        my @query_hit = split ',', $blast_hit[2];
        (my $query = $query_hit[0]) =~ s/Q://;
        my @queries = split '-', $query;
        my $query_length = abs($queries[1] - $queries[0]);
        (my $hit = $query_hit[1]) =~ s/H://;
        (my $identity = $blast_hit[3]) =~ s/%ID//;
        my $eValue   = $blast_hit[4];
        (my $hit_name = $blast_hit[5]) =~ s/RecName: Full=//;
        $hit_name =~ s/;//;
        my $kegg = $temp[11];
        my $go = $temp[12];

        if (($identity > 50) && ($query_length >= 100)){
          unless (exists $report{$id}){
                $report{$id}[0] = $prot_hit;
                $report{$id}[1] = $query;
                $report{$id}[2] = $hit;
                $report{$id}[3] = $identity;
                $report{$id}[4] = $eValue;
                $report{$id}[5] = $hit_name;
                $report{$id}[6] = $kegg;
                $report{$id}[7] = $go;
          }
        }
}
close $tr_fh;

# node to protein list
my %node2pro;
while (my $line = <$n2p_fh>){
        chomp $line;
        my @temp = split /\t/, $line;
        my $id = $temp[0];
        unless (exists $node2pro{$id}){
                $node2pro{$id} = $temp[1];
        }
}
close $n2p_fh;

sub read_fasta {
    my ($filename) = @_;
    my %sequences;
    my $header;

    open(my $fh, '<', $filename) or die "Cannot open file $filename: $!";
    while (my $line = <$fh>) {

        chomp $line;
        if ($line =~ /^>(.*)/) {
            ($header = $1) =~ s/>//;
            $header =~ s/ .*$//;
        } else {
            $sequences{$header} .= $line;
        }
    }
    close $fh;

    return %sequences;
}

# Read both FASTA files
my %nucleotide_seqs = read_fasta("$dir/transcripts.fasta.transdecoder.complete.cds");
my %protein_seqs = read_fasta("$dir/transcripts.fasta.transdecoder.complete.pep");

open my $nt_out, ">", "$dir/transcripts.fasta.transdecoder.complete_with_annotation.fasta";
open my $pro_out, ">", "$dir/transcripts.fasta.transdecoder.complete_with_annotation.faa";


# Print header of trinotate table
print $tri_out "Transc_ID\tNode\tProtein\tQuery\tHit\tIdentity\te-value\tname\tKeggAccession\tGeneOntology\n";

# Iterate over the nucleotide sequences
foreach my $id (keys %nucleotide_seqs) {
      if (exists $report{$id}){
        my $nucleotide_seq = $nucleotide_seqs{$id};
        my $protein_seq = $protein_seqs{$id};
        my $protein = $node2pro{$id};
	if (($protein) && ($uniprot{$protein})) {
                my $fasta_annotations = $uniprot{$protein};
                my $long_id = ">" . $id . " " . $fasta_annotations;
                print $nt_out $long_id, "\n", $nucleotide_seq, "\n";
                print $pro_out $long_id, "\n", $protein_seq, "\n";
                print $tri_out "$id\t$id\t$report{$id}[0]\t$report{$id}[1]\t$report{$id}[2]\t$report{$id}[3]\t$report{$id}[4]\t$report{$id}[5]\t$report{$id}[6]\t$report{$id}[7]\n";
          }
        }
}

close $tri_out;
close $nt_out;
close $pro_out;
exit;
