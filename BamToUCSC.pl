#######################################################################################
###     		     					  BamToUCSC	  		   					    ###
#######################################################################################

use Getopt::Long;
use strict;
use warnings;

###################################################################################
### 	  	 				    	INITIALIZE    	 	 						###
###################################################################################

my	$BamPath		= "";
my	$Sorted			= "";
my	$Index			= "";
my	$Genome			= "";
my	$Organism		= "";
my	$UCSCDir		= "";
my	$URLDir			= "";
my 	$Info			= "";

###################################################################################
### 	  	 				    	GET OPTIONS    	 	 						###
###################################################################################

GetOptions ("Bam=s"			=> 	\$BamPath,
			"Sorted=s"		=>	\$Sorted,
			"Index=s"		=>	\$Index,
			"Genome=s"		=>	\$Genome,
			"Organism=s"	=>	\$Organism,
			"UCSCDir=s"		=>	\$UCSCDir,
			"URLDir=s"		=>	\$URLDir,
			"Info=s"		=>	\$Info);

unless 	(-f $BamPath && $BamPath =~ m/.bam$/)					{die ("$BamPath is not a .bam file!\n");}
unless 	($Sorted eq "yes" || $Sorted eq "no")					{die ("$BamPath should either be sorted or unsorted!\n");}
unless	($Index eq "yes" || $Index eq "no")						{die ("$BamPath should either be Indexed or not indexed!\n");}
unless	($Genome ne "" && $Organism ne "")						{die ("Genome and Organism is a mandatory argument\n");}
unless	($UCSCDir =~ m/^\/home\/www\/UCSC\//)					{die ("UCSCDir should be a /home/www/UCSC/ subdirectory\n");}
unless	(-d $URLDir)											{die ("URLDir should be a directory\n");}

###################################################################################
### 	  	 				    	PROCESS BAM    	 	 						###
###################################################################################

# Users and rights
my 	$user = `who am I`;
	$user =~ s/ .*$//g;
	$user =~ s/\n//g;

### UCSCDir ###

unless (-d $UCSCDir){system ("mkdir -p $UCSCDir");
					 system ("chmod 775 $UCSCDir");
					 system ("chown $user\:researchers $UCSCDir");}

### FileNames ###

my	$BamFileName		= $BamPath;
	$BamFileName		=~ s/^.+\///g;
my 	$TrackFileName		= $BamFileName;
	$TrackFileName 		=~ s/\.bam/\.track/g;
my 	$TrackName			= $UCSCDir;
	$TrackName			=~ s/^\/home\/www\/UCSC\///g;
	$TrackName			=~ s/\/$//g;
my 	$BaiPath			= $BamPath;
	$BaiPath			=~ s/\.bam/\.bai/g;

### Copy ###

					 system ("cp $BamPath /home/www/UCSC/$TrackName\/");	
if (-f $BaiPath)	{system ("cp $BaiPath /home/www/UCSC/$TrackName\/");	}
	
### URL.txt ###

my $URLFile = $URLDir . "URL.txt";

if (-f $URLFile){
	
	open 	OUTPUT, ">>$URLFile" or die ("Can't open $URLFile\n");
}
else {

	open 	OUTPUT, ">>$URLFile" or die ("Can't open $URLFile\n");
	print 	OUTPUT 	"TrackName\tOrganism\tGenome\tInfo\tURL\n";
}



### Sort ###

my 	$SortedBam = "/home/www/UCSC/$TrackName\/$BamFileName";
	
if ($Sorted eq "no"){
	
	$SortedBam =~ s/\.bam/\.sorted/g;
	
	system ("./samtools sort /home/www/UCSC/$TrackName\/$BamFileName -o $SortedBam");
}

### Index ###

if ($Index eq "no"){

	system ("./samtools index $SortedBam");

}

### Sorted BamFile ###

my	$SortedBamFileName		= $SortedBam;
	$SortedBamFileName		=~ s/^.+\///g;

### Track ###
	
open TRACK, ">$UCSCDir$TrackFileName" or die ("DIE $UCSCDir$TrackFileName\n");		
print TRACK "track type=bam name=\"$TrackName $Info $Organism\" visibility=squish db=$Genome\ bigDataUrl=http://medgenpr.ugent.be/UCSC/$TrackName\/$SortedBamFileName\n";
close TRACK;
		
print OUTPUT "$TrackName\t$Organism\t$Genome\t$Info\thttp://genome.ucsc.edu/cgi-bin/hgTracks?db=$Genome\&hgt.customText=http://medgenpr.ugent.be/UCSC/$TrackName\/$TrackFileName\n";

close OUTPUT;

