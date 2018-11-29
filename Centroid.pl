#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

##
##Input: It takes in input the name of a protein with extension .pdb
##It reads the corresponding PDB file, extracts from it the coordinates of the C_alpha atoms
##and creates three arrays x, y and z with such coordinates.
##It takes in the coordinates, and calls a subroutine that computes the centroid of the atoms
##It takes in the coordinates, centroid position, and calls a subroutine that computes the distance of
##every C_alpha atom of the protein from its centroid
##It takes in the distances, and calls a subroutine that prints the histogram of the distances of the
##atoms from the centroid.
##Usage: perl Centroid.pl protein_name.pdb

sub centroid {
my($x_ref,$y_ref,$z_ref)=@_;#passing reference of x,y,z coordinates
my @x=@{$x_ref};
my @y=@{$y_ref};
my @z=@{$z_ref};
my ($sum_x,$sum_y,$sum_z)=0;
my $count=0;
  for (my $i = 0; $i < scalar @x; $i++) { #sum each x/y/z
  $count++;
  $sum_x+=$x[$i];
  $sum_y+=$y[$i];
  $sum_z+=$z[$i];
  }
my $aver_x=$sum_x/$count;
my $aver_y=$sum_y/$count;
my $aver_z=$sum_z/$count;
print "\nCoordinates of the centroid of the protein is","\t","x=",$aver_x,"\t","y=",$aver_y,"\t","z=",$aver_z,"\n\n";
return (\$aver_x,\$aver_y,\$aver_z); #return references of average numbers of x,y,z coordinates
}


sub distance {#passing references of each x,y,z coordinates and the centroid x,y,z coordinates
  my($x_ref,$y_ref,$z_ref,$centerx_ref,$centery_ref,$centerz_ref)=@_;
  my @x=@{$x_ref};
  my @y=@{$y_ref};
  my @z=@{$z_ref};
  my $cx=${$centerx_ref};
  my $cy=${$centery_ref};
  my $cz=${$centerz_ref};
  my $dist=();
  my @dist=();
  for(my $index=0; $index<scalar @x; $index++){
    $dist=((($x[$index]-$cx)**2)+(($y[$index]-$cy)**2)+(($z[$index]-$cz)**2))**0.5;# count distance
    push(@dist,$dist);
  }
  return @dist;#return a array contain all distance from centroid coordinates
}


sub histogram {
  my @dist=@_;
  my ($max,$min)=();
  for(my $index=0; $index<scalar @dist; $index++){
    if (!defined($min) || $dist[$index] > $max) {
    $max = $dist[$index];
  }
   if (!defined($min) || $dist[$index] < $min) {
    $min = $dist[$index];
  }} #find the max and the min distance


  my %histogram; #make a hash
  my $width=5; #set 5 as $width
  for (my $index=0; $index<scalar @dist; $index++) {
    my $range=floor($dist[$index]/$width);
    if(exists $histogram{$range}){
      $histogram{$range}+= 1;
    }
    else{
      $histogram{$range}= 1;
    }} # range as key of hash; numbers of distance in certian range as vaule of hash


  for (my $i = floor($min/$width); $i <= floor($max/$width); $i++){
    my $frequency = $histogram{$i};
    print $i,"\t",$frequency = "#" x $frequency,"\n";
  }#print the range number and display the frequency of distance in this range with "#" counts

  print "======================================================\n\n";
  print     "Width:",$width,"\n";
  print     "Range:",$min,"-",$max,"\n\n";
}


my $pdbfilename;
for(my $i=0;$i<@ARGV; $i++){
	$pdbfilename=$ARGV[$i];}
open(PDBFILE,$pdbfilename)or die "Unable to open $pdbfilename: $!";
chomp(my @pdbdata = <PDBFILE>);
close PDBFILE;#read pdbfile put contents of pdb file into array

my (@x,@y,@z)=();
my ($centroidx_ref,$centroidy_ref,$centroidz_ref)=();
for (my $line = 0; $line < scalar @pdbdata; $line++) {
  if ($pdbdata[$line] =~ m/^ATOM.+/) {
      my $atom= substr($pdbdata[$line],13,2);
      if ($atom eq "CA") {
        my $vaule_x= substr($pdbdata[$line],31,8);
        my $vaule_y= substr($pdbdata[$line],39,8);
        my $vaule_z= substr($pdbdata[$line],47,8);
        push(@x,$vaule_x);
        push(@y,$vaule_y);
        push(@z,$vaule_z);
        }
    }
          last if ($pdbdata[$line] =~/^TER.+/);}#get the coordinates of each ca atom from the first chain(until first "TER")

my($cx_ref,$cy_ref,$cz_ref)=centroid(\@x,\@y,\@z);#passing references of arrays of coordinates into subroutine centroid

my @dist=distance(\@x,\@y,\@z,$cx_ref,$cy_ref,$cz_ref);#passing references of arrays coordinates and centroid coordinates into subroutine distance

histogram(@dist);#passing array of distances into subroutine histogram to generate histogram

exit;
