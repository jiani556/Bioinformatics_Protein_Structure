#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

#####################################################################################################################
#A Protein-Protein Interaction (PPI) network is a graph where nodes represent proteins and an edge between two      #
#               nodes represents interacting proteins (either physically or functionally).                          #
#This script will compute the following local properties of a PPI network and compare those in random graphs:       #
#                             degree distribution and clustering coefficient                                        #
#Input: the Protein-Protein Interaction (PPI) graph in sif format. Each line of the file represents an edge.        #
#Output:                                                                                                            #
#1.a table (histogram) of the degree frequencies, where each line consists of a value of k and the corresponding 𝑝" #
#2.the average clustering coefficient AVG_C of the network                                                          #
#3.for all k, the average clustering coefficient AVG_C(k) of all nodes of degree k                                  #
#####################################################################################################################

sub histogram;
sub cu;
sub random;
my $sif=$ARGV[0];
open(SIF,$sif)or die "Unable to open $sif: $!"; #read filename from command line and open file

my %ar;
my %cu;
while(<SIF>){
  my %degree;
  my @columns=split;
  if($columns[2] eq $columns[0]){next;}
  if(!grep(/^$columns[2]$/,@{$ar{$columns[0]}})){
    push(@{$ar{$columns[0]}},$columns[2]);
  }
  if(!grep(/^$columns[0]$/,@{$ar{$columns[2]}})){
    push(@{$ar{$columns[2]}},$columns[0]);
  }
  foreach my $key(keys %ar){
    $degree{$key}=scalar @{$ar{$key}};
    print $key."\t".$degree{$key}."\n";}
}#make a hash to store array named by each node which contains every node have a connection.
close SIF;

print "\nhistogram for input network\n";
histogram(%ar);
cu(%ar);#make histgoram count coefficient for the input network file

my $p1=5;
my %p1=random(\%ar,\$p1);
print "histogram for p=1/2\n";
histogram(%p1);
cu(%p1);#to generate a new network file by the same nodes, random sub make the possiblity of connection between two node 0.5,

my $p2=1;
my $p3=8;
my %p2=random(\%ar,\$p2);#to generate a new network file by the same nodes, random sub make the possiblity of connection between two node 0.1,
my %p3=random(\%ar,\$p3);#to generate a new network file by the same nodes, random sub make the possiblity of connection between two node 0.8,

print "histogram for p=1/10\n";
histogram(%p2);
cu(%p2);

print "histogram for p=8/10\n";
histogram(%p3);
cu(%p3);


sub histogram {
  my (%degree)=@_;
  my $width=2;
  my ($max,$min)=();
  my %histogram; #make a hash


  foreach my $keys(keys %degree){
    my $count=scalar @{$degree{$keys}};
    my $range=floor($count/$width);

    if (!defined($max) || $count > $max) {
    $max = $count;
  }
   if (!defined($min) || $count < $min) {
    $min = $count;
  }

    if(exists $histogram{$range}){$histogram{$range}+= 1;}
    else{$histogram{$range}= 1;}
  } #find the max and the min degree# range as key of hash; numbers of distance in certian range as vaule of hash

  for (my $i=floor($min/$width);$i<=floor($max/$width);$i++){
    my $frequency = $histogram{$i};
    if(!$frequency){next;}
    else{print $i,"\t",$frequency = "#" x $frequency,"\n";}
  }#print the range number and display the frequency of distance in this range with "#" counts
  print     "------------------------------------------------------\n";
  print     "Width:",$width,"\n";
  print     "Range:",$min,"-",$max,"\n\n";
}

sub cu {
  my (%degree)=@_;
  my %nu;
  my %cluster;
  my $sum=0;
  my $count=0;

  foreach my $key1(keys %degree){
    my $n=scalar @{$degree{$key1}};
    for(my $i=0;$i<$n-1;$i++){
      for(my $j=$i+1;$j<$n;$j++){
        my $c=@{$degree{$key1}}[$j];
        my $d=@{$degree{$key1}}[$i];
        if(grep(/^$c$/,@{$degree{$d}})){
           $nu{$key1}+=1;}#new hash to store every nodes' cluster number
        }}}
  print "======================================================\n\n";
  print "Node\tClustering coefficient\n";
  foreach my $keys(keys %degree){
          my $n=scalar @{$degree{$keys}};
          if($n>1 && exists $nu{$keys}){
          my $cu=2*$nu{$keys}/($n*($n-1));
          $sum=$sum+$cu;
          print $keys."\t".$cu."\n";}
          $count++;}
        my $ave=$sum/$count;
        print     "------------------------------------------------------\n";
        print "AVG-C=$ave\n";
        print "======================================================\n\n";
      }

sub random {
    my ($hash_ref,$p_ref)=@_;
    my %degree=%{$hash_ref};
    my $p=${$p_ref};
    my %random;
    my $filename="p_$p.sif";
    open(FW,">$filename") or die("could not open");#generate new sif file
        foreach my $a(sort keys %degree){
        foreach my $b(sort keys %degree){
        if($a ne $b && !grep(/^$b$/,@{$random{$a}})){
        my $i=int(rand(10));#$i equal to any int number between 0-9 randomly
        if($i<$p){push(@{$random{$a}},$b);
        if(!grep(/^$a$/,@{$random{$b}})){print FW "$a\t1\t$b\n"}}
      }
        else{next;}
         }
        }
        close FW;
        return %random;
      }



exit;
