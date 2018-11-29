#!/usr/bin/perl
use strict;
use warnings;

###Comput The top 10 nodes (sorted) with the highest closeness centrality value###

my $sif=$ARGV[0];
open(SIF,$sif)or die "Unable to open $sif: $!"; #read filename from command line and open file
my @data=<SIF>;
close SIF;


sub Create_Adjacency_Matrix;
sub Degree;
sub All_Pairs_Shortest_Paths;
sub Closeness_Centrality;
sub Get_shortest_path;

my($hash_ref,$matrix_ref)=Create_Adjacency_Matrix(@data);
Degree(@data);
my($path_ref,$p_ref)=All_Pairs_Shortest_Paths($hash_ref,$matrix_ref);
my(%cu)=Closeness_Centrality($hash_ref,$path_ref);
my ($max,$min,$max_note,$min_note);
foreach my $key(keys %cu){#find the max and min
  if(!defined $max|| $cu{$key}>$max){$max=$cu{$key};$max_note=$key;};
  if(!defined $min|| $cu{$key}<$min){$min=$cu{$key};$min_note=$key;};
}
print "-------------------------------------------------------------\n";
print "The lowest Closeness Centrality nodes:$min_note\nThe highest Closeness Centrality nodes:$max_note\n\n";

# Get_shortest_path(\$max_note,\$min_note,$p_ref);

sub Create_Adjacency_Matrix(){
  my (@input)=@_;
  my %hash;
  my @matrix;
  my $i=0;

  for(@input){
    my @columns=split;
    my $a=$columns[0];
    my $distance=$columns[1];
    my $b=$columns[2];
    if (!exists $hash{$a}){$hash{$a}=$i; $i++};
    if (!exists $hash{$b}){$hash{$b}=$i; $i++};

    $matrix[$hash{$a}][$hash{$b}]=$distance;
    $matrix[$hash{$b}][$hash{$a}]=$distance;
  }
  foreach my $i(sort keys %hash){
    $matrix[$hash{$i}][$hash{$i}]=0;
    foreach my $j(sort keys %hash){
      if(!defined $matrix[$hash{$i}][$hash{$j}]){
        $matrix[$hash{$i}][$hash{$j}]=9999;}
      if(!defined $matrix[$hash{$j}][$hash{$i}]){
        $matrix[$hash{$j}][$hash{$i}]=9999;}
}}#generate matrix and hash
  return(\%hash,\@matrix);
}

sub Degree(){
  my (@input)=@_;
  my %ar;
  my %degree;
  my $i=0;
  for(@input){
    my @columns=split;
    if($columns[2] eq $columns[0]){next;}
    if(!grep(/^$columns[2]$/,@{$ar{$columns[0]}})){
      push(@{$ar{$columns[0]}},$columns[2]);
    }
    if(!grep(/^$columns[0]$/,@{$ar{$columns[2]}})){
      push(@{$ar{$columns[2]}},$columns[0]);}
    }
  foreach my $key(keys %ar){
    $degree{$key}=scalar @{$ar{$key}};}
    print "-------------------------------------------------------------\n";
    print "The top 10 nodes (sorted) with the highest degree\n";
    print "Nodes\t\tDegree\n";
  foreach my $key(sort { $degree{$b} <=> $degree{$a} } keys %degree){
    last if $i==10;
    print $key."\t".$degree{$key}."\n";
    $i++;}
  }#count drgree for each nodes
#
sub All_Pairs_Shortest_Paths(){
  my ($refhash,$refmat)=@_;
  my %hash=%{$refhash};
  my @path=@{$refmat};
  my %p;
  foreach my $i(keys %hash){
    foreach my $j(keys %hash){
      foreach my $k(keys %hash){
        if($path[$hash{$i}][$hash{$j}] + $path[$hash{$i}][$hash{$k}] < $path[$hash{$j}][$hash{$k}] ){
          $path[$hash{$j}][$hash{$k}]=$path[$hash{$i}][$hash{$j}] + $path[$hash{$i}][$hash{$k}];
        }}}}

return(\@path,\%p);}


sub Closeness_Centrality(){
  my($refhash,$refmat)=@_;
  my %hash=%{$refhash};
  my @path=@{$refmat};
  my %cu;
  my $i=0;
  foreach my $i(keys %hash){
    my $sum=0;
    foreach my $j(keys %hash){
      $sum= $sum + $path[$hash{$i}][$hash{$j}];
    }
      $cu{$i}=1/$sum;
  }
  print "-------------------------------------------------------------\n";
  print "The top 10 nodes (sorted) with the highest closeness centrality value\n";
  print "Nodes\t\tCloseness Centrality Value\n";
  foreach my $key(sort { $cu{$b} <=> $cu{$a} } keys %cu){
    last if $i==10;
    print $key."\t".$cu{$key}."\n";
    $i++;
  }
  return %cu;
}
#
# sub Get_shortest_path(){
#   my($max_ref,$min_ref,$p_ref)=@_;
#   my $highest=${$max_ref};
#   my $lowest=${$min_ref};
#   my %p=%{$p_ref};
#   print "$lowest\t$highest\n";
#   for(@{$p{$highest}{$lowest}}){print $_;}
#}

exit;
