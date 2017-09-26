###################################################################################
### 	 					Install BATCH-GE - Wouter Steyaert     				###
###################################################################################

use strict;
use Getopt::Long;

###################################################################################
### 	  	 		  	 	  		   PREFIX   	 	 						###
###################################################################################

my 	$PREFIX				= "";
my 	$WebAccFolder		= "";

if ($0 =~ m/^\//){
	
	($PREFIX) = $0 =~ m/(.*)BATCH-GE\.pl/;
}
else {

	$PREFIX = `pwd`;
	$PREFIX	=~ s/\n//g;
}

###################################################################################
### 	  	 				    	GET OPTIONS    	 	 						###
###################################################################################

GetOptions ("WebAccFolder=s"	=> 	\$WebAccFolder);

if ($WebAccFolder && ! -d $WebAccFolder){
	die ("$WebAccFolder should be a directory\nYou can also install BATCH-GE without this option\n");
}

###################################################################################
### 	  	 				    	 INSTALL    	 	 						###
###################################################################################

if (-d $WebAccFolder){

	my $WebAccFolderFilePath = $PREFIX . "/.web-acc";
	
	open WEBACC, ">$WebAccFolderFilePath" or die ();
	print WEBACC "$WebAccFolder";	
	close WEBACC;
}

###################################################################################
### 	  	 				    	 EXAMPLE    	 	 						###
###################################################################################

my $EXAMPLE_1_EXPERIMENT_FILE 	= $PREFIX . "/_EXAMPLE_1/" . "ExperimentFile.csv";
my $EXAMPLE_2_EXPERIMENT_FILE 	= $PREFIX . "/_EXAMPLE_2/" . "ExperimentFile.csv";
my $EXAMPLE_3_EXPERIMENT_FILE 	= $PREFIX . "/_EXAMPLE_3/" . "ExperimentFile.csv";
my $EXAMPLE_4_EXPERIMENT_FILE 	= $PREFIX . "/_EXAMPLE_4/" . "ExperimentFile.csv";

ChangePrefix($EXAMPLE_1_EXPERIMENT_FILE);
ChangePrefix($EXAMPLE_2_EXPERIMENT_FILE);
ChangePrefix($EXAMPLE_3_EXPERIMENT_FILE);
ChangePrefix($EXAMPLE_4_EXPERIMENT_FILE);

###################################################################################
### 	  	 				    	INSTALL BWA    	 	 						###
###################################################################################

my	$BWADIRECTORY			= $PREFIX . "/bwa-0.7.13/";

chdir($BWADIRECTORY);
system("make -s 2>/dev/null");

###################################################################################
### 	  	 				    	SUBROUTINES    	 	 						###
###################################################################################

sub ChangePrefix {

	my $EXPERIMENT_FILE_PATH 	= $_[0];
	my @EXP_FILE 				= ();
	my $Line_Nr					= 0;
	
	if (-f $EXPERIMENT_FILE_PATH){
		
		open EXP, "$EXPERIMENT_FILE_PATH" or die ("$!");
		@EXP_FILE = <EXP>;
		close EXP;
		
		system ("rm $EXPERIMENT_FILE_PATH");
		
		open EXP, ">$EXPERIMENT_FILE_PATH" or die ("$!");
		
		foreach my $EXP_LINE (@EXP_FILE){
			
			$EXP_LINE =~ s/\n//g;
			$EXP_LINE =~ s/;$//g;
		
			if ($Line_Nr == 0){
				
				print EXP "$EXP_LINE\;";
			
			}
			else {
				
				my @EXP_FILE_VALUES = split (";", $EXP_LINE);
				
				foreach my $EXP_FILE_VALUE (@EXP_FILE_VALUES){
					
					if 		($EXPERIMENT_FILE_PATH =~ m/_EXAMPLE_1/){
						
						$EXP_FILE_VALUE =~ s/.*\/_EXAMPLE_1/$PREFIX\/_EXAMPLE_1/g;

						print EXP "$EXP_FILE_VALUE\;";
						
					}
					elsif 	($EXPERIMENT_FILE_PATH =~ m/_EXAMPLE_2/){
						
						$EXP_FILE_VALUE =~ s/.*\/_EXAMPLE_2/$PREFIX\/_EXAMPLE_2/g;

						print EXP "$EXP_FILE_VALUE\;";
						
					}
					elsif 	($EXPERIMENT_FILE_PATH =~ m/_EXAMPLE_3/){
						
						$EXP_FILE_VALUE =~ s/.*\/_EXAMPLE_3/$PREFIX\/_EXAMPLE_3/g;

						print EXP "$EXP_FILE_VALUE\;";
						
					}
					elsif 	($EXPERIMENT_FILE_PATH =~ m/_EXAMPLE_4/){
						
						$EXP_FILE_VALUE =~ s/.*\/_EXAMPLE_4/$PREFIX\/_EXAMPLE_4/g;

						print EXP "$EXP_FILE_VALUE\;";
						
					}					
				}
			}
			print EXP "\n";
			$Line_Nr++;
		}
		close EXP;
	}
}
