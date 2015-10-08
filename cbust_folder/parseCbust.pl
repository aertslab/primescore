#!/usr/bin/env perl
#
# Purpose: Parse ClusterBuster output and return a gff file.
#
# Usage:   cat ./cbust-output.txt |./parseCbust.pl <customclustername> > ./cbust-output.gff


my $line;
if($ARGV[0]){
    $name=$ARGV[0]; #e.g. dl-I
}
else{
    $name="cbust";
}
my $seq;
my $cnt=0;

while(<STDIN>){
    $line=$_;
    chomp $line;
    if($line=~/CLUSTER\s+(\d+)/){
        $cluster=$1;
    }
    if($line =~/^>(.*)\s+\(/){
	$seq=$1;
    }
    if($line=~/Location:\s+(\d+)\s+to\s+(\d+)/){
	$from=$1;
	$to=$2;
    }
    if($line=~/Score:\s+([\d,.]+)/){
	$score=$1;
	print $seq."\tCbust\tCRM\t".$from."\t".$to."\t".$score."\t+\t.\tcluster \"".$name."-cluster-".$cluster."\"\n";
    }

    if($line=~/^(.*\w)\s+(\d+)\s+(\d+)\s+([+,-])\s+([\d,\.]+)\s+(\w+)/){      
        $strand=$4;
        if($strand eq 'p'){$strand='.'}
        print $seq."\tCbust\tmisc_feature\t".$2."\t".$3."\t".$5."\t".$strand."\t.\tid \"".$1."\"; site \"".$6."\"; cluster \"".$name."-cluster-".$cluster."\"\n";
    }

}

