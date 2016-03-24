#!/usr/bin/perl -w
use strict;
use warnings;

my $alignpath="./CDS_unaligned_NU";    # unaligned Nu files folder
my $outfilefoder='CDS_unaligned_AA';       # output unaligned AA files folder


opendir(FOLDER,"$alignpath");
my @array = grep(/fa/,readdir(FOLDER));
close FOLDER;

system ("mkdir $outfilefoder"); 


foreach my $filename ( @array ){

my $DNAfilename = $filename;

open (F, "$alignpath/$DNAfilename");
     my $species;
     my %sequences;
     my @order;
	
        while( <F> )
	{
		    my $line=$_;  
                if( $line =~ /^(>.+\n)/ )
                {
                        $species = $1;  push (@order,$1);
                        $sequences{$species} = '';   
                }
                else
                   {
                            $sequences{$species} .= $line;                      
                   }

	 }
	close F;


my $outname=$filename; $outname=~s/\.fas/_AA\.fas/;
open (OUT, ">$outfilefoder/$outname");

foreach my $spe (@order){
my $DNA=$sequences{$spe};
my $protein='';

my $codon;
for(my $i=0;$i<(length($DNA)-2);$i+=3)
{
$codon=substr($DNA,$i,3);
$protein.=&codon2aa($codon);
}

print OUT "$spe$protein\n";

}

}

sub codon2aa{
my($codon)=@_;
$codon=uc $codon;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
if(exists $g{$codon})
{
return $g{$codon};
}
else
{

return 'X';
#print STDERR "Bad codon \"$codon\"!!\n";
#exit;
}
}
                                           
