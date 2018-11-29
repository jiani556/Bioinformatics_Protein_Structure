#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

sub fraction;
##################################################################################################
#It takes as arguments the C_alpha atoms of the amino acids of each of the two chains and        #
#computes the distances Dist_Calpha of all pairs of C_alpha atoms, with one C_alpha in chain A   #
#and the other in chain B. It outputs the interfaced amino acids                                 #
##################################################################################################

sub interfa;
##################################################################################################
# for each chain computes the fraction of the interface amino acids lying on the secondary       #
#                        structures alpha helices and beta sheets                                #
##################################################################################################

sub distance;
##################################################################################################
# for each interface amino acid of a chain determines the closest interface atom of the same     #
# chain to the right in the primary sequence. Then it determines the difference in position of   #
# the two amino acids in the sequence.                                                           #
##################################################################################################

my ($pdbfilename, $chain1, $chain2, $tod) = @ARGV;
open(PDBFILE,$pdbfilename)or die "Unable to open $pdbfilename: $!";
chomp(my @pdbdata = <PDBFILE>);
close PDBFILE;#read pdbfile put contents of pdb file into array
chomp $tod;# threshold for the distance between two amino acids

my %chain1;
my %chain2;#two hash to get the residue residue number atom name and coordinates of two chain
for (my $line = 0; $line < scalar @pdbdata; $line++) {
			if($pdbdata[$line] =~ m/^ATOM/){
				my $atom_name = substr($pdbdata[$line],13,4);
				$atom_name =~ s/^\s+|\s+$//g;
        my $chain = substr($pdbdata[$line],21,2);
        $chain =~ s/^\s+|\s+$//g;
				if($atom_name eq "CA" && ($chain eq $chain1 || $chain eq $chain2)){
					my $res = substr($pdbdata[$line],17,10);
					$res =~ s/^\s+|\s+$//g;
					my $x = substr($pdbdata[$line],30,8);
					$x =~ s/^\s+|\s+$//g;
					my $y = substr($pdbdata[$line],38,8);
					$y =~ s/^\s+|\s+$//g;
					my $z = substr($pdbdata[$line],46,8);
					$z =~ s/^\s+|\s+$//g;
					my $cord = $x."_".$y."_".$z;
					if($chain eq $chain1){
						if(exists($chain1{$res})) {
							$chain1{$res}{$atom_name} = $cord;
						}
						else{
							$chain1{$res}{$atom_name} = $cord;
					}}
					if($chain eq $chain2){
						if(exists($chain1{$res})) {
							$chain2{$res}{$atom_name} = $cord;
						}
						else{
							$chain2{$res}{$atom_name} = $cord;
						}
					}
				}}}

my($inter_ref1,$inter_ref2)=interfa(\%chain1,\%chain2,\$tod);
fraction(\@pdbdata,$inter_ref1);
fraction(\@pdbdata,$inter_ref2);
my @inter1=@{$inter_ref1};
my @inter2=@{$inter_ref2};
distance(@inter1);
distance(@inter2);

sub interfa{     #take the reference of two hash and the threshold
my ($a,$b,$c)=@_;
my %chain1=%{$a};
my %chain2=%{$b};
my $tod=${$c};
my(@interface1,@interface2);
my $distance;
	foreach my $key1(sort keys %chain1){
		my $x1; my $y1; my $z1;
		foreach my $atom1(keys $chain1{$key1}){
			my @ar1 = split("_",$chain1{$key1}{$atom1});
	  				$x1 = $ar1[0];
						$y1 = $ar1[1];
						$z1 = $ar1[2];
			foreach my $key2(sort keys %chain2){
				my $x2;my $y2;my $z2;
				foreach my $atom2(keys $chain2{$key2}){
					my @ar2 = split("_",$chain2{$key2}{$atom2});
		  				$x2 = $ar2[0];
							$y2 = $ar2[1];
							$z2 = $ar2[2];
							$distance = (($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2)**0.5;
		if ($distance<=$tod){
		push(@interface1,$key1."\n");
		push(@interface2,$key2."\n");
		last;}}}
		if ($distance<=$tod){
		last;}}}#calculated the distance and push into array if less than threshold
	my %seen = ();
	@interface1 = grep { ! $seen{ $_ }++ } @interface1;
	@interface2 = grep { ! $seen{ $_ }++ } @interface2;
	print "\nInterfaced amino acid for the first chain:"."\n"."aa"." "."chain"." "."No."." "."\n";
	print @interface1;
	print "\nInterfaced amino acid for the second chain:"."\n"."aa"." "."chain"." "."No."." "."\n";
	print @interface2;
	return (\@interface1,\@interface2);
}     #return the reference of two chain's interface aa



sub fraction(){    #take the array of two chain's interface aa
my($ref_pdb,$ref_inter)=@_;
my @pdbdata=@{$ref_pdb};
my @interface=@{$ref_inter};
my $hcount=0;
my $scount=0;
my $chain;
for(@interface){
$chain=substr($_,4,1);
my $id=substr($_,5,4);
chomp $id;
for (my $line = 0; $line < scalar @pdbdata; $line++) {
			if($pdbdata[$line] =~ m/^HELIX/){
				my $hchain = substr($pdbdata[$line],19,1);
				my $hstart = substr($pdbdata[$line],21,4);
				$hstart =~ s/^\s+|\s+$//g;
        my $hend= substr($pdbdata[$line],33,4);
        $hend =~ s/^\s+|\s+$//g;
				if($chain eq $hchain && $id >= $hstart && $id <= $hend){
					$hcount+=1;}}
			elsif($pdbdata[$line] =~ m/^SHEET/){
				my $schain = substr($pdbdata[$line],21,1);
				my $sstart = substr($pdbdata[$line],22,4);
				$sstart =~ s/^\s+|\s+$//g;
				my $send= substr($pdbdata[$line],33,4);
				$send =~ s/^\s+|\s+$//g;
				if($chain eq $schain && $id >= $sstart && $id <= $send){
				$scount+=1;}}
				}}
my $hfra=$hcount/scalar(@interface);
my $sfra=$scount/scalar(@interface);
print "\nFor chain $chain \n";
print "the fraction of the interface aas lying on the alpha helices is \n $hfra \n";
print "the fraction of the interface aas lying on the beta sheets is \n $sfra \n\n";
}   #count the fraction of interface aa in helices and sheets

sub distance(){   #take the interface array
	my @interface=@_;
	my $idc;
	my @dist;
	my %hash;
	my $chain;

	foreach(@interface){
		my $aa=substr($_,0,3);
		$chain=substr($_,4,1);
		my $id=substr($_,5,4);
		$hash{$id} = $aa;
		}

	foreach my $id(sort keys %hash){
			if (! defined $idc){
				$idc=0;
				push(@dist,$hash{$id}." ".$idc."\n");
			}
			else{
				my $cc=$id-$idc;
				push(@dist,$hash{$id}," ".$cc."\n");
			}
	$idc=$id;
	}

	print "\nDistance of the closest interface atom of the chain $chain \n";

	foreach(@dist){ print $_;}
}#Count the distance of the closest interface atom
exit;
