#!/usr/bin/env perl                                                                                                                                                            
use strict;
use warnings;


my %refseqs;
my %geneids;

open(REFSEQ,$ARGV[0]);
while(<REFSEQ>){
    chomp;
    $refseqs{$_} = 1;
}
close(REFSEQ);

open(GENE2ACCESSION,"unpigz -c ".$ARGV[1]."|");
while(<GENE2ACCESSION>){
    chomp;
    my @cols = split(/\t/);

    if(exists $refseqs{$cols[5]}){
        $geneids{$cols[1]} = $cols[5];
    }
}
close(GENE2ACCESSION);


open(GENEINFO,"unpigz -c ".$ARGV[2]."|");
while(<GENEINFO>){
    chomp;
    my @cols = split(/\t/);

    if(exists $geneids{$cols[1]}){
        if($cols[5]=~/HGNC:(HGNC:\d+)/){
            print $geneids{$cols[1]}."\t".$cols[1]."\t".$1."\n";
        }
    }
}
close(GENEINFO);
