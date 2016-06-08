###################################################################################
### 	  					 Prepare Genomes - Wouter Steyaert     				###
###################################################################################
###################################################################################
###     		     			   Load Libraries   	  		   				###
###################################################################################

use strict;
use warnings;
use Getopt::Long;

###################################################################################
### 	  	 				    	INITIALIZE    	 	 						###
###################################################################################

my $PREFIX			= "";
my $Genome 			= "";

###################################################################################
### 	  	 		  	 	  		   PREFIX   	 	 						###
###################################################################################

if ($0 =~ m/^\//){
	
	($PREFIX) = $0 =~ m/(.*)\/GetGenome\.pl/;
}
else {

	$PREFIX = `pwd`;
	$PREFIX	=~ s/\n//g;
}

###################################################################################
### 	  	 				    	GET OPTIONS    	 	 						###
###################################################################################

GetOptions ("Genome=s"	=> 	\$Genome);

	$Genome 		=~ s/ //g;

###################################################################################
### 	  	 		  	 	  		PROCES GENOME   	 	 					###
###################################################################################

my	$JAVA_INSTALLATION		= $PREFIX . "/jre1.8.0_91/bin/java";
my 	$SAMTOOLSBINDIRECTORY	= $PREFIX . "/samtools-1.3.1/bin/";
my	$BWADIRECTORY			= $PREFIX . "/bwa-0.7.13/";
my 	$PICARD_INSTALLATION	= $PREFIX . "/picard-tools-2.3.0/";
my 	$GenomeDirectory		= $PREFIX . "/genomes/$Genome\/";
my	$PicardLogFile 			= "$PREFIX\/logs/picard.log";
my	$BWALogFile 			= "$PREFIX\/logs/bwa.log";

unless (-d "$PREFIX\/logs/"){
	mkdir("$PREFIX\/logs/");
}

unless (-d $GenomeDirectory){
	system ("mkdir $GenomeDirectory");
}
system ("mv $PREFIX\/genomes/$Genome\.fa.gz $GenomeDirectory 2>/dev/null");
system ("gunzip --silent $GenomeDirectory\/*");

### CreateSequenceDictionary ###

chdir 	("$PICARD_INSTALLATION");

system 	("$JAVA_INSTALLATION -Xmx4g -jar picard.jar CreateSequenceDictionary R=$PREFIX\/genomes\/$Genome\/$Genome\.fa O=$PREFIX\/genomes\/$Genome\/$Genome\.dict > $PicardLogFile 2>> $PicardLogFile");
	
### Faidx ###

chdir 	("$SAMTOOLSBINDIRECTORY");
	
system	("./samtools faidx $PREFIX\/genomes\/$Genome\/$Genome\.fa");

### Index ###

chdir ($BWADIRECTORY);

system("./bwa index $PREFIX\/genomes\/$Genome\/$Genome\.fa > $BWALogFile 2>> $BWALogFile");

print "Done!\n";