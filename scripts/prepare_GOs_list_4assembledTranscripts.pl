#!/usr/bin/perl

use strict;
use warnings;

my $GOs_table  = $ARGV[0];
my $pseudo_ids = $ARGV[1];

open my $gos_desc_fh, "<", $GOs_table;
open my $pseudo_ids_fh, "<", $pseudo_ids;

# Fill hashes
my %pseudo_ids;

while (my $line = <$pseudo_ids_fh>){
    chomp $line;
    my @temp = split /\t/, $line;
    $pseudo_ids{$temp[0]} = $temp[1];
}
close $pseudo_ids_fh;


while (my $line = <$gos_desc_fh>){
    chomp $line;

    my @temp = split /\t/, $line;
    my $prot_id = $temp[2];
    my $seudoID = $pseudo_ids{$prot_id};
    my $gos     = $temp[9];
    my @ind_GOs = split /(?<=`)(?=GO:)/, $gos;
    foreach my $term(@ind_GOs){
        my @GOs_parts = split /\^/, $term;
        my $GOterm    = $GOs_parts[0];
        my $GOdesc    = $GOs_parts[2];
        my $ontology  = $GOs_parts[1];
        if ($GOterm && $GOdesc && $ontology) {
            $GOdesc =~ s/`//;
            print "$seudoID\t$GOterm\t$GOdesc\t$ontology\n";
        }
    }
}
close $gos_desc_fh;

exit;

    

    
    
    

