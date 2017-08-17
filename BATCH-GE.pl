###################################################################################
### 	  Calculate efficiency of CRIPR/CAS mutagenesis - Wouter Steyaert     	###
###################################################################################
###################################################################################
###     		     			   Load Libraries   	  		   				###
###################################################################################

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

###################################################################################
### 	  	 		  	 	  		   PREFIX   	 	 						###
###################################################################################

my 	$PREFIX		= "";

if ($0 =~ m/^\//){
	
	($PREFIX) = $0 =~ m/(.*)BATCH-GE\.pl/;
}
else {

	$PREFIX = `pwd`;
	$PREFIX	=~ s/\n//g;
}

###################################################################################
### 	  	 		  	 	  		DEPENDENCIES    	 	 					###
###################################################################################

my	$LOCATION_BAMTOUCSCSCRIPT 	= $PREFIX . "/BamToUCSC.pl";
my	$LOCATION_BBMAP 			= $PREFIX . "/bbmap/";
my 	$FASTXBINDIRECTORY			= $PREFIX . "/fastx_toolkit-0.0.14/bin/";
my 	$SAMTOOLSBINDIRECTORY		= $PREFIX . "/samtools-1.3.1/bin/";
my	$BWADIRECTORY				= $PREFIX . "/bwa-0.7.13/";
my	$JAVA_INSTALLATION			= $PREFIX . "/jre1.8.0_91/bin/java";
my 	$PICARD_INSTALLATION		= $PREFIX . "/picard-tools-2.3.0/";
my 	$GENOMES_MAINFOLDER			= $PREFIX . "/genomes/";
my 	$WEBACCFOLDERPATH			= $PREFIX . "/.web-acc";

###################################################################################
### 	  	 				    	INITIALIZE    	 	 						###
###################################################################################

my 	$LineNr				= 0;
my	$ExperimentFile		= "/";
my 	$WEBACC				= "";
my 	$Help				= 0;
my 	$Man				= 0;
my 	%ExperimentFile		= ();
my 	%OutputDirs			= ();
my	%FileHeaders		= ();

###################################################################################
### 	  	 				  Web Accessible Folder     	 					###
###################################################################################

if (-f $WEBACCFOLDERPATH){
	
	open WEBACC, "$WEBACCFOLDERPATH" or die ("$!\n");
	while (<WEBACC>){
		$WEBACC = $_;
		$WEBACC =~ s/\n//g;
	}
	close WEBACC;
}

###################################################################################
### 	  	 				    	GET OPTIONS    	 	 						###
###################################################################################

GetOptions ("ExperimentFile=s"	=> 	\$ExperimentFile,
			"help|?"   			=>  \$Help,
			"man"   			=>  \$Man) or
			
pod2usage(2);
pod2usage(1) if $Help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $Man;

unless (-f $ExperimentFile){die ("\nYou should pass an ExperimentFile as an argument to the script\n");}

###################################################################################
### 	  	 	    		READ IN EXPERIMENT FILE    	 	 					###
###################################################################################

print "\nCRISPR-CAS MUTAGENESIS EFFICIENCY CALCULATOR\n";
print "\n\tReading experiment file ...\n";

open EXP, "$ExperimentFile" or die ("Can't open $ExperimentFile!\n");
while (<EXP>){

	my 	$Line 				= $_;
		$Line				=~ s/\n//g;
		$Line				=~ s/\r//g;
	my 	@LineValues			= split(";", $Line);
	my 	$FastqDir			= "/";
	my  $SampleNumbers		= "/";
	my 	$Genome				= "/";
	my 	$CutSite			= "/";
	my 	$OutputDir			= "/";
	my 	$CutSitesFile		= "/";
	my 	$RepairSequence		= "";

	if ($LineNr == 0){
		
		my $ColumnNr = 0;
					
		foreach my $LineValue (@LineValues){
			
			$LineValue 					=~ s/\s/_/g;
			$LineValue 					=~ s/\(|\)//g;
			$FileHeaders{$LineValue} 	= $ColumnNr;
			$ColumnNr++;
		}
	}
	else {
			
		if (exists $FileHeaders{"FastqDir"} 		&& exists $LineValues[$FileHeaders{"FastqDir"}])					
			{$FastqDir			= $LineValues[$FileHeaders{"FastqDir"}];}						else{$FastqDir 	= "/";}
		
		if (exists $FileHeaders{"SampleNumbers"} 	&& exists $LineValues[$FileHeaders{"SampleNumbers"}])					
			{$SampleNumbers		= $LineValues[$FileHeaders{"SampleNumbers"}];}					else{$SampleNumbers 	= "/";}
				
		if (exists $FileHeaders{"Genome"} 			&& exists $LineValues[$FileHeaders{"Genome"}])					
			{$Genome 				= $LineValues[$FileHeaders{"Genome"}];}						else{$Genome 			= "/";}
		
		if (exists $FileHeaders{"CutSite"} 			&& exists $LineValues[$FileHeaders{"CutSite"}])					
			{$CutSite 				= $LineValues[$FileHeaders{"CutSite"}];}					else{$CutSite 			= "/";}
			
		if (exists $FileHeaders{"OutputDir"} 		&& exists $LineValues[$FileHeaders{"OutputDir"}])					
			{$OutputDir 			= $LineValues[$FileHeaders{"OutputDir"}];}					else{$OutputDir 		= "/";}
			
		if (exists $FileHeaders{"CutSitesFile"} 	&& exists $LineValues[$FileHeaders{"CutSitesFile"}])					
			{$CutSitesFile 			= $LineValues[$FileHeaders{"CutSitesFile"}];}				else{$CutSitesFile 		= "/";}
			
		if (exists $FileHeaders{"RepairSequence"} 	&& exists $LineValues[$FileHeaders{"RepairSequence"}])					
			{$RepairSequence 		= $LineValues[$FileHeaders{"RepairSequence"}];}				else{$RepairSequence 	= "";}
			
		### Check genome ###
		
		if (! -f "$PREFIX\/genomes/$Genome\/$Genome\.fa"){
		
			die( "\n\t$Genome is not yet available on your system. Download the appropriate fasta file, run the PrepareGenome.pl and retake the analysis\n");	
		}
		
		my 	$IndividualSampleNumbersRef = DetermineIndividualSampleNumbers	($SampleNumbers);	# Pools are initially in the format 7,11,13,15-20 #
			
		foreach my $SampleNumber (@$IndividualSampleNumbersRef){
		
			if ($FastqDir ne "/" && $SampleNumbers ne "/" && $Genome ne "/" && $CutSite ne "/" && $OutputDir ne "/" && $CutSitesFile ne "/"){
			
				$ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}{$OutputDir}{$CutSitesFile}{$CutSite} = "$RepairSequence";
			
			}
		}
	}
	$LineNr++;
}
close EXP;

###################################################################################
### 	  	 				      GO OVER POOLS    	 	 						###
###################################################################################

foreach my $FastqDir (sort keys %ExperimentFile){
	foreach my $SampleNumber (sort keys %{$ExperimentFile{$FastqDir}}){
		foreach my $Genome 		(sort keys %{$ExperimentFile{$FastqDir}{$SampleNumber}}){		
			foreach my $OutputDir (sort keys %{$ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}}){
							
				# The OutputDir and all of its content will be deleted unless this OutputDir is already used in the same Experiment file. 	#
				# If so, the efficiency file and the variant file are extended with the new analysis, the previous results won't be deleted	#
				
				unless 	(-d $OutputDir)						{system ("mkdir -p $OutputDir");}
				if 		(!exists $OutputDirs{$OutputDir})	{CleanUpLogAndOutput ($OutputDir);}
				
				$OutputDirs{$OutputDir} = undef;
				
				my $EfficiencyFile 				= 	"$OutputDir\/Efficiencies.txt";
				my $VariantFile					= 	"$OutputDir\/Variants.txt";
				my $RepairFile 					= 	"$OutputDir\/RepairReport.txt";
				my $CutSitesFileForSampleFile 	= 	"$OutputDir\/CutSitesPool\.bed";
				my $CutSitesFileForSampleString	= 	"";
									
				open EFF, 				">>$EfficiencyFile" 			or die ("Can't open $EfficiencyFile\n");
				open VAR, 				">>$VariantFile" 				or die ("Can't open $VariantFile\n");										
				open SPECIFIC_CUTSITES, ">$CutSitesFileForSampleFile" 	or die ("Can't open $CutSitesFileForSampleFile\n");
				
				print 					"\n\tSample number $SampleNumber from $FastqDir is being analyzed ...\n";
				print 	EFF				"\n\tSample number $SampleNumber from $FastqDir is being analyzed ...\n";
				print 	VAR				"\n\tSample number $SampleNumber from $FastqDir is being analyzed ...\n";
				print 					"\t[Genome=$Genome\]\n\n";
				print 	EFF				"\t[Genome=$Genome\]\n";
				print 	VAR				"\t[Genome=$Genome\]\n";
				
				# Write the cutsites for this particular sample to a new temporary bed file. The coordinate of the cutsites can be stored in	#
				# several cutsitesfiles. The new bed file will be used later to create a smaller bam file for further processing.				#

				foreach my $CutSitesFile (sort keys %{$ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}{$OutputDir}}){
									
					my $CutSitesAllRef = ReadInCutSitesFile ($CutSitesFile);	# Read in Cutsitesfile #

					foreach my $CutSite (sort keys %{$ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}{$OutputDir}{$CutSitesFile}}){
					
						$CutSitesFileForSampleString .= $CutSite . "~";

						if ($CutSitesAllRef->{$CutSite} && $ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}{$OutputDir}{$CutSitesFile}{$CutSite}){	
						
							(my $CutSiteChr, my $CutSiteStart, my $CutSiteStop) 			= split ("_", $CutSitesAllRef->{$CutSite});
							
							print SPECIFIC_CUTSITES "$CutSiteChr\t$CutSiteStart\t$CutSiteStop\t$CutSite\;RepairSequence=$ExperimentFile{$FastqDir}{$SampleNumber}{$Genome}{$OutputDir}{$CutSitesFile}{$CutSite}\n";
						}
						elsif ($CutSitesAllRef->{$CutSite}){
							
							(my $CutSiteChr, my $CutSiteStart, my $CutSiteStop) 			= split ("_", $CutSitesAllRef->{$CutSite});
							
							print SPECIFIC_CUTSITES "$CutSiteChr\t$CutSiteStart\t$CutSiteStop\t$CutSite\n";
						
						}
						else {
							die ("CutSite $CutSite should be defined in your CutSitesFile: $CutSitesFile\n");
						}
					}
				}
				close SPECIFIC_CUTSITES;
				
				$CutSitesFileForSampleString =~ s/~$//g;
				
				################################################################################################
				# 									Process Sequencing Reads								   #
				################################################################################################

				my $CopiedFilesRef 					= 	CopyFiles						($FastqDir, $SampleNumber, $OutputDir);
														TrimReadsOnQuality				($OutputDir, $CopiedFilesRef);										
														RepairReads						($OutputDir, $CopiedFilesRef, $LOCATION_BBMAP);						
				my $AmBasis							=	MapReadsToReference				($Genome, $OutputDir, $CopiedFilesRef, $GENOMES_MAINFOLDER);
				
				
				my $CutAmBasis						=	CollectCutSiteReadsAndProcess	($FastqDir, 
																						$SampleNumber, 
																						$Genome, 
																						$OutputDir, 
																						$AmBasis, 
																						$CutSitesFileForSampleFile,
																						$LOCATION_BAMTOUCSCSCRIPT, 
																						$CutSitesFileForSampleString,
																						$PICARD_INSTALLATION,
																						$JAVA_INSTALLATION,
																						$GENOMES_MAINFOLDER,
																						$WEBACC);
				(my $ConflictReadPairsRef, 
				 my $FilteredReadsRef, 
				 my $CutSiteReadGroupsRef) 			=	AnalyzeSAM						($FastqDir, 
																						 $SampleNumber, 
																						 $Genome, 
																						 $OutputDir, 
																						 $CutAmBasis,
																						 $CutSitesFileForSampleFile,
																						 $RepairFile);															
														CleanUpFolders					($OutputDir, $CutSitesFileForSampleFile);

				close EFF;
				close VAR;
			}				
		}
	}
}

###################################################################################
### 	  	 				     	 SUBROUTINES    	 	 					###
###################################################################################

sub DetermineIndividualSampleNumbers{

	my 	@IndividualSampleNumbers	= ();
	my 	@SampleNumbersCommaSep		= split (",", $_[0]);
	
	foreach my $SampleNrCommaSep (@SampleNumbersCommaSep){
		if ($SampleNrCommaSep =~ m/^\d+$/){
			push (@IndividualSampleNumbers, $SampleNrCommaSep);
		}
		elsif ($SampleNrCommaSep =~ m/^\d+-\d+$/){
			(my $FirstSampleNr, my $LastSampleNr) = split ("-", $SampleNrCommaSep);
			
			for (my $SampleNr = $FirstSampleNr; $SampleNr <= $LastSampleNr; $SampleNr++){
				push (@IndividualSampleNumbers, $SampleNr);
			}
		}
		else {
			@IndividualSampleNumbers = @SampleNumbersCommaSep;
		}
	}
	
	return \@IndividualSampleNumbers;
}

sub CleanUpLogAndOutput {

	my 	$OutputDir = $_[0];
			
	if (-d "$OutputDir\/1_OriginalReads")				{system("rm -rf $OutputDir\/1_OriginalReads");}
	if (-d "$OutputDir\/2_QualityFilteredReads")		{system("rm -rf $OutputDir\/2_QualityFilteredReads");}
	if (-d "$OutputDir\/3_Mappings")					{system("rm -rf $OutputDir\/3_Mappings");}
	if (-d "$OutputDir\/4_CutSiteReads")				{system("rm -rf $OutputDir\/4_CutSiteReads");}
	if (-f "$OutputDir\/EfficienciesCRISPRCas.txt")		{system("rm -rf $OutputDir\/EfficienciesCRISPRCas.txt");}
	if (-f "$OutputDir\/RepairReport.txt")				{system("rm -rf $OutputDir\/RepairReport.txt");}
	if (-f "$OutputDir\/URL.txt")						{system("rm -rf $OutputDir\/URL.txt");}
	if (-f "$OutputDir\/Variants.txt")					{system("rm -rf $OutputDir\/Variants.txt");}
}

sub ReadInCutSitesFile {
							
	my 	$CutSitesFile 	= $_[0];
	my 	%CutSitesAll	= ();
			
	 open CUTSITES, "$CutSitesFile" or die $!;
	 while (<CUTSITES>){
	 
		my 	$Line = $_;
			$Line =~ s/\n//g;
			$Line =~ s/\r//g;
		(my $Chr, my $Start, my $Stop, my $Cutsite) 		= split ("\t", $Line);
					
		$CutSitesAll{"$Cutsite"} 							= "$Chr\_$Start\_$Stop";
		
	 }
	 close CUTSITES;
	
	return \%CutSitesAll;
}

sub CopyFiles {

	print		 "\t\t\tCopy Files To Working Directory\n";

	(my $FastqDir, my $SampleNumber, my $Outputdir) = @_;

	 my $SubdirToCopyTo 	= "$Outputdir\/1_OriginalReads\/";
	 my $SampleNumberQM		= quotemeta ($SampleNumber);
	 my %CopiedFiles		= ();
	 
	 unless (-d $SubdirToCopyTo) {system("mkdir -p $SubdirToCopyTo");}
	 
	 opendir (ORIGINAL_DATA_DIR, $FastqDir) or die ("Can't open $FastqDir!\n");
	 
	 while (my $FastQGZFile = readdir(ORIGINAL_DATA_DIR)){
	 		
		if ($FastQGZFile =~ m/.*_S$SampleNumber\_.*fastq\.gz$/i){
			
			system("cp $FastqDir$FastQGZFile $SubdirToCopyTo");
			$CopiedFiles{$FastQGZFile} = undef;
		}
	 }
	 closedir(ORIGINAL_DATA_DIR);
	 
	 if (scalar keys %CopiedFiles != 2){
		die ("ERROR $FastqDir, SAMPLE:$SampleNumber\tThere should be 2 sequencing files for each sample!\n");
	 }
	 
	 return \%CopiedFiles;
}

sub TrimReadsOnQuality {

	print 		"\t\t\tTrim Reads On Quality\n";

	(my $Outputdir, my $CopiedFilesRef) = @_;
	(my $in_1, my $in_2) 				= sort keys %$CopiedFilesRef;
	 my $FastXLogFile 					= "$Outputdir\/0_Logs/fastx.log";

	unless (-d "$Outputdir\/0_Logs/") 					{system("mkdir -p $Outputdir\/0_Logs/");}
	unless (-d "$Outputdir\/2_QualityFilteredReads/") 	{system("mkdir -p $Outputdir\/2_QualityFilteredReads/");}
		
	chdir("$FASTXBINDIRECTORY");
	
	system ("gunzip $Outputdir\/1_OriginalReads/$in_1");
	system ("gunzip $Outputdir\/1_OriginalReads/$in_2");
	
	$in_1	=~ s/\.gz$//g;
	$in_2	=~ s/\.gz$//g;
	
	system ("./fastq_quality_trimmer -Q33 -t 30 -z -i $Outputdir\/1_OriginalReads/$in_1 -o $Outputdir\/2_QualityFilteredReads/$in_1\.gz > $FastXLogFile");
	system ("./fastq_quality_trimmer -Q33 -t 30 -z -i $Outputdir\/1_OriginalReads/$in_2 -o $Outputdir\/2_QualityFilteredReads/$in_2\.gz >> $FastXLogFile");
	
}

sub RepairReads {

	print 		"\t\t\tRepair reads\n";

	(my $Outputdir, my $CopiedFilesRef, my $LOCATION_BBMAP) = @_;
	(my $in_1, my $in_2) 				= sort keys %$CopiedFilesRef;
	 my $BBMapLogFile 					= "$Outputdir\/0_Logs/bbmap.log";
	
	unless (-d "$Outputdir\/3_RepairedReads/") 	{system("mkdir -p $Outputdir\/3_RepairedReads/");}
	
	chdir($LOCATION_BBMAP);
	
	system("./repair.sh --overwrite=t in1=$Outputdir\/2_QualityFilteredReads/$in_1 in2=$Outputdir\/2_QualityFilteredReads/$in_2 out1=$Outputdir\/3_RepairedReads/$in_1 out2=$Outputdir\/3_RepairedReads/$in_2 >> $BBMapLogFile 2>&1");

}

sub MapReadsToReference {

	print 		"\t\t\tMap Reads To Reference\n";
	
	(my $Genome, my $Outputdir, my $CopiedFilesRef, my $GENOMES_MAINFOLDER) = @_;
	(my $in_1, my $in_2) 							= sort keys %$CopiedFilesRef;
	(my $AmBasis)									= $in_1 =~ m/^(.*)\.fastq\.gz$/i;
	
	 my $BwaLogFile 		= "$Outputdir\/0_Logs/bwa.log";
	 my $SamLogFile 		= "$Outputdir\/0_Logs/sam.log";
	 my $BamOut				= $AmBasis . ".bam";
	 my	$SamOut				= $AmBasis . ".sam";
	
	unless (-d "$Outputdir\/4_Mappings/"){system("mkdir -p $Outputdir\/4_Mappings/");}
	
	chdir ($BWADIRECTORY);

	system("./bwa mem -t 2 -M $GENOMES_MAINFOLDER\/$Genome/$Genome\.fa $Outputdir\/3_RepairedReads/$in_1 $Outputdir\/3_RepairedReads/$in_2 > $Outputdir\/4_Mappings/$SamOut 2>> $BwaLogFile");
	
	chdir ($SAMTOOLSBINDIRECTORY);
	
	system("./samtools view -Sb  $Outputdir\/4_Mappings/$SamOut > $Outputdir\/4_Mappings/$BamOut 2>> $SamLogFile");
	
	return $AmBasis;
}

sub CollectCutSiteReadsAndProcess {

	print 		"\t\t\tCollect CutSite Reads\n";
		
	(my $FastqDir, my $SampleNumber, my $Genome, my $Outputdir, my $AmBasis, my $CutSitesFileForSampleFile, my $LOCATION_BAMTOUCSCSCRIPT, my $CutSitesFileForSampleString, my $PICARD_INSTALLATION, my $JAVA_INSTALLATION, my $GENOMES_MAINFOLDER, my $WEBACC) = @_;

	 my $CutAmBasis			= $AmBasis . "_" . $CutSitesFileForSampleString;	
	 my $PicardLogFile 		= "$Outputdir\/0_Logs/picard.log";
	 my $GATKLogFile 		= "$Outputdir\/0_Logs/gatk.log";

	unless 	(-d "$Outputdir\/5_CutSiteReads/")		{system("mkdir -p $Outputdir\/5_CutSiteReads/");}
	
	chdir ($SAMTOOLSBINDIRECTORY);
	
	system 	("./samtools view -b -L  $CutSitesFileForSampleFile $Outputdir\/4_Mappings/$AmBasis\.bam > $Outputdir\/5_CutSiteReads/$CutAmBasis\.bam");
	
	### SortSam ###
	
	chdir 	("$PICARD_INSTALLATION");
	
	system	("$JAVA_INSTALLATION -Xmx4g -jar picard.jar SortSam I=$Outputdir\/5_CutSiteReads/$CutAmBasis\.bam O=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.bam SO=coordinate 2>> $PicardLogFile");
	
	### MarkDuplicates ###
	
	system 	("$JAVA_INSTALLATION -Xmx4g -jar picard.jar MarkDuplicates I=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.bam O=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.bam M=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.metrics REMOVE_DUPLICATES=true 2>> $PicardLogFile");
	
	### AddOrReplaceReadGroups ###
	
	system 	("$JAVA_INSTALLATION -Xmx4g -jar picard.jar AddOrReplaceReadGroups I=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.bam O=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.bam SO=coordinate LB=Nextera PL=Illumina PU=$FastqDir\~$SampleNumber SM=$CutAmBasis 2>> $PicardLogFile");
	
	### Index ###
	
	chdir ($SAMTOOLSBINDIRECTORY);
	
	system 	("./samtools index $Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.bam");
		
	### bamToSam ###
		
	system 	("./samtools view -h -o $Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.sam $Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.bam");
	
	### bamToUCSC ###
	
	if ($WEBACC){
		system ("perl $LOCATION_BAMTOUCSCSCRIPT --Bam=$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.bam --Sorted=no --Index=no --Organism=$Genome --Genome=$Genome	--UCSCDir=$WEBACC --URLDir=$Outputdir --Info=Sample$SampleNumber");
	}
		
	return $CutAmBasis;
}

sub AnalyzeSAM {
	
	print 		"\t\t\tAnalyzeSAM\n";
	
	(my $FastqDir, my $SampleNumber, my $Genome, my $Outputdir, my $CutAmBasis, my $CutSitesFileForSampleFile, my $RepairFile) = @_;
	 
	 my $FilteredReadPairsRef				= {};
	 my $NrOfFilteredReadsRef				= {};
	 my $CutSitesPoolRef					= {};
	 my	$RepairSequenceRef					= {};
	 my $CutSiteReadGroupsRef				= {};
	 my	$ConflictReadPairsRef				= {};
	 my $PossibleRepairSequencesRef			= {};
	 my $RepairSequencesInIndelsRef			= {};
	 my $RepairSequencesIndelsConflictsRef	= {};
	 my $InDelsRef							= {};
	 
	#############################################
	### 	  Read in Cut Sites For Pool	  ###
	#############################################
	 
	 open SPECIFIC_CUTSITES, "$CutSitesFileForSampleFile" or die $!;
	 while (<SPECIFIC_CUTSITES>){
		 my $Line = $_;
			$Line =~ s/\n//g;
		(my $ReadChr, my $RegionStart, my $RegionStop, my $Sites) 			= split ("\t", $Line);
		
		if ($Sites =~ m/RepairSequence\=/){
		
			(my $CutSite, my $RepairSequence) = split (";", $Sites);
			
			$RepairSequence =~ s/RepairSequence\=//;
			
			$RepairSequenceRef	->{$CutSite} 								= $RepairSequence;
			$RepairSequenceRef	->{$CutSite} 								=~ s/\(|\)|\[|\]//g;
			$CutSitesPoolRef	->{$ReadChr}{$RegionStart}{$RegionStop} 	= "$CutSite";

			$PossibleRepairSequencesRef = CreatePossibleRepairSequences ($RepairSequence, $PossibleRepairSequencesRef);
		}
		else {
		
			$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop} 		= "$Sites";
		}
	 }
	 close SPECIFIC_CUTSITES;
	
	#############################################
	### 		     Go over SAM			  ###
	#############################################
	
	open SAM, "$Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.sam" or die ("Can't open $Outputdir\/5_CutSiteReads/$CutAmBasis\.sorted.remdup.rg.sam\n");
	while (<SAM>){
	
		next if ($_ =~ m/^\@/);
			
		my 	$Line 				= $_;
			$Line 				=~ s/\n//g;
		my 	@LineValues			= split ("\t", $Line);
		
		my 	$ReadPairName 		= $LineValues[0];
		my 	$ReadChr 			= $LineValues[2];
		my 	$ReadStart 			= $LineValues[3];
		my 	$Cigar 				= $LineValues[5];
		my 	$Sequence 			= $LineValues[9];
		my	$ReadStop			= CalculateReadStopFromReadStartAndCigar ($ReadStart, $Cigar);

		if ($CutSitesPoolRef->{$ReadChr}){
			foreach my $RegionStart (keys %{$CutSitesPoolRef->{$ReadChr}}){
				foreach my $RegionStop (keys %{$CutSitesPoolRef->{$ReadChr}{$RegionStart}}){
					
					# if ((($ReadStart) <= $RegionStart && $ReadStop >= $RegionStart) || 							This condition states that there is an overlap
					#(($ReadStart) >= $RegionStart && ($ReadStart) <= $RegionStop)){								between read and cutsite - we dont use this
					
					if (($ReadStart) <= $RegionStart && $ReadStop >= $RegionStop){											# 	This condition states that the cutsite region is
																															# 	fully part of the read - read goes fully throug it
						if 	(exists $NrOfFilteredReadsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}})		# 	Determine total number of reads - total means that		
							{$NrOfFilteredReadsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}++;}			# 	the possible conflicts not yet excluded
						else 																								# 	Determine total number of fragments i.e. readpairs		
							{$NrOfFilteredReadsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}} = 1;}			# 	with 1 or 2 reads		
							
						
						$FilteredReadPairsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{$ReadPairName} = undef;
						
						if ($Cigar =~ m/[N|P|X]/g){die("$_");}
						if ($Cigar =~ m/[D|I]/g){																# 	Read with INDEL - but not necessary in cutsiteregion 	#
																												# 	Thus, the position of the INDEL needs to be checked		#
							
							my	$IndelInCutSite = "no";										
							my 	$PositionInRead	= $ReadStart;
							my 	$CigarRef	 	= SplitCigar($Cigar);
							my 	@Cigar			= @$CigarRef;

							for (my $I = 0; $I < scalar @Cigar; $I+=2){
								
								if ($Cigar[$I+1] =~ m/[D|I]/g){
									
									if (($PositionInRead < $RegionStart && ($PositionInRead + $Cigar[$I] ) > $RegionStart) ||					# RegionStart is inside I/D 	#
										($PositionInRead >= $RegionStart && $PositionInRead <= $RegionStop)){									# I/D starts after RegionStart 	#
														
										if (!exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"NO_INDEL"}{$ReadPairName} &&
											!exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"REPAIR"}{$ReadPairName}){
										
											$CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL"}{$ReadPairName} = undef;
																						
											# ($RepairSequencesInIndelsRef, $RepairSequencesIndelsConflictsRef) = CheckForRepairInIndelRead ($CutSitesPoolRef, $RepairSequenceRef, $PossibleRepairSequencesRef , $ReadChr, $RegionStart, $RegionStop, $Sequence, $ReadPairName, $RepairSequencesInIndelsRef, $RepairSequencesIndelsConflictsRef);
																						
											$InDelsRef = DetermineAndStoreIndelCoordinate ( $CigarRef, $I, $PositionInRead, $ReadChr, $ReadPairName, 
																							$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}, $InDelsRef);
											$IndelInCutSite = "yes";
										}
										else {
										
											$ConflictReadPairsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{$ReadPairName} = undef;
										}
									}
								}
								
								$PositionInRead += $Cigar[$I];
							}
							
							if ($IndelInCutSite eq "no"){																								# Check Status #
							
								 ($CutSiteReadGroupsRef, $ConflictReadPairsRef) 	= StoreAsNoIndelOrRepair ($CutSiteReadGroupsRef, $RepairSequenceRef, $CutSitesPoolRef, $ConflictReadPairsRef, $ReadChr, $RegionStart, $RegionStop, $ReadPairName, $Sequence, $PossibleRepairSequencesRef);	
							}
						}
						else { 																															# Read without INDEL #
							
							($CutSiteReadGroupsRef, $ConflictReadPairsRef) 			= StoreAsNoIndelOrRepair ($CutSiteReadGroupsRef, $RepairSequenceRef, $CutSitesPoolRef, $ConflictReadPairsRef, $ReadChr, $RegionStart, $RegionStop, $ReadPairName, $Sequence, $PossibleRepairSequencesRef);	
							
						}
					}
				}
			}
		}
	}
	close SAM;
	
	#############################################
	### 		   Remove conflicts 		  ###
	#############################################	
	
	foreach my $CutSite (keys %$ConflictReadPairsRef){
	
		foreach my $ReadPairName (keys %{$ConflictReadPairsRef->{$CutSite}}){
		
			### Remove in CutSiteReadGroups ###
			
			if (exists $CutSiteReadGroupsRef->{$CutSite}{"INDEL"}{$ReadPairName}){
				delete $CutSiteReadGroupsRef->{$CutSite}{"INDEL"}{$ReadPairName};
			}
			elsif (exists $CutSiteReadGroupsRef->{$CutSite}{"NO_INDEL"}{$ReadPairName}){
				delete $CutSiteReadGroupsRef->{$CutSite}{"NO_INDEL"}{$ReadPairName};
			}
			elsif (exists $CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}{$ReadPairName}){
				delete $CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}{$ReadPairName};
			}
			
			### Remove in InDels ###
			
			if ($InDelsRef->{$CutSite}){
				
				foreach my $Chr (keys %{$InDelsRef->{$CutSite}}){
					foreach my $Coorindates (keys %{$InDelsRef->{$CutSite}{$Chr}}){
						foreach my $Type (keys %{$InDelsRef->{$CutSite}{$Chr}{$Coorindates}}){
							foreach my $Length (keys %{$InDelsRef->{$CutSite}{$Chr}{$Coorindates}{$Type}}){
								
								if (exists $InDelsRef->{$CutSite}{$Chr}{$Coorindates}{$Type}{$Length}{$ReadPairName}){
									delete $InDelsRef->{$CutSite}{$Chr}{$Coorindates}{$Type}{$Length}{$ReadPairName};
								}
							}
						}
					}
				}
			}
		}
	}
	
	### Remove in RepairSequencesInIndels ###
	
	foreach my $CutSite (keys %$RepairSequencesIndelsConflictsRef){
		foreach my $ReadPairName (keys %{$RepairSequencesIndelsConflictsRef->{$CutSite}}){
		
			if (exists $RepairSequencesInIndelsRef->{$CutSite}{"INDEL_NO_REPAIR"}{$ReadPairName}){
				delete $RepairSequencesInIndelsRef->{$CutSite}{"INDEL_NO_REPAIR"}{$ReadPairName};
			}
			elsif (exists $RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}{$ReadPairName}){
				delete $RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}{$ReadPairName};
			}
		}
	}
	
	#############################################
	### 	     Write out efficiency		  ###
	#############################################
	
	foreach my $CutSite (sort keys %$CutSiteReadGroupsRef){
	
		print 			"\t[Cutsite=$CutSite\]\n";
		print 	EFF		"\t[Cutsite=$CutSite\]\n";
		print 	VAR		"\t[Cutsite=$CutSite\]\n";
	
		my 	$Nondels 		= 	scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"NO_INDEL"}};
			$Nondels		+=	scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}} if ($CutSiteReadGroupsRef->{$CutSite}{"REPAIR"});
		my 	$Indels			=	scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"INDEL"}};
		my 	$Efficiency		= Round($Indels/($Indels + $Nondels));
		
		print 		"\n\t\t\t\tMutagenesis efficiency for $CutSite is $Efficiency ($Indels readpairs with indel(s) versus $Nondels readpairs without indel(s))\n";		
		print EFF 	"\n\t\t\t\tMutagenesis efficiency for $CutSite is $Efficiency ($Indels readpairs with indel(s) versus $Nondels readpairs without indel(s))\n";
		
		if ($CutSiteReadGroupsRef->{$CutSite}{"REPAIR"} && 
			keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}}){	
			
			my $Repairs				= scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}};
			my $RepairEfficiency	= $Repairs/$Indels;
			my $TotalNrOfReadPairs	= $Indels + $Nondels;
			my $RepairFraction		= Round($Repairs/$TotalNrOfReadPairs);
			
			open 	REP, 	">>$RepairFile" or die ("Can't open $RepairFile\n");		
			print 	REP		"\n\tFastqDir $FastqDir Sample $SampleNumber\n";
			print 	REP	 	"\t[Cutsite=$CutSite\]\n";
			
			print 			"\t\t\t\tRepair efficiency for $CutSite is $RepairFraction ($Repairs readpairs with repair versus $TotalNrOfReadPairs readpairs in total)\n\n";		
			print 	REP 	"\t\t\t\tRepair efficiency for $CutSite is $RepairFraction ($Repairs readpairs with repair versus $TotalNrOfReadPairs readpairs in total)\n\n";
			print 	EFF 	"\t\t\t\tRepair efficiency for $CutSite is $RepairFraction ($Repairs readpairs with repair versus $TotalNrOfReadPairs readpairs in total)\n\n";
		
			close 	REP;
		}		
	}
	
	#############################################
	### 	    Write out repair report	  	  ###
	#############################################
		
	foreach my $CutSite (sort keys %$CutSiteReadGroupsRef){
	
		if ($CutSiteReadGroupsRef->{$CutSite}{"REPAIR"} && 
			keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}}){
		
			my 	$RepairSequence 				= $RepairSequenceRef->{$CutSite};
			my 	$Repairs 						= scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}};
			my 	%RepairSequencesNoIndelCount	= ();
			my 	%RepairSequencesIndelCount		= ();
		
			foreach my $ReadPairName (keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}}){
				
				if (exists $RepairSequencesNoIndelCount{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}{$ReadPairName}}){
				
					$RepairSequencesNoIndelCount{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}{$ReadPairName}}++;
				
				}
				else {
					$RepairSequencesNoIndelCount{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}{$ReadPairName}} = 1;
				
				}			
			}
			
			if ($RepairSequencesNoIndelCount{$RepairSequence}){
			
				my $FullHDRNoInDel = $RepairSequencesNoIndelCount{$RepairSequence};
			
				open 	REP, 	">>$RepairFile" or die ("Can't open $RepairFile\n");
				print 			"\t\t\t\tThere are $FullHDRNoInDel readpairs with FULL HDR (on a total of $Repairs repairs)\n";	
				print 	REP 	"\t\t\t\tThere are $FullHDRNoInDel readpairs with FULL HDR (on a total of $Repairs repairs)\n";	
				close 	REP;
			}				
			
			foreach my $OtherRepairSequence (keys %RepairSequencesNoIndelCount){
			
				unless ($OtherRepairSequence eq $RepairSequence){
					
					open 	REP, 	">>$RepairFile" or die ("Can't open $RepairFile\n");
					print 			"\t\t\t\tThere are $RepairSequencesNoIndelCount{$OtherRepairSequence} readpairs with PARTIAL HDR: $OtherRepairSequence (on a total of $Repairs repairs)\n\n";	
					print 	REP 	"\t\t\t\tThere are $RepairSequencesNoIndelCount{$OtherRepairSequence} readpairs with PARTIAL HDR: $OtherRepairSequence (on a total of $Repairs repairs)\n\n";
					close 	REP;
						
				}
			}
			
			if ($RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}){
			
				foreach my $ReadPairName (keys %{$RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}}){
					
					if (exists $RepairSequencesIndelCount{$RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}{$ReadPairName}}){
					
						$RepairSequencesIndelCount{$RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}{$ReadPairName}}++;
					}
					else {
					
						$RepairSequencesIndelCount{$RepairSequencesInIndelsRef->{$CutSite}{"INDEL_REPAIR"}{$ReadPairName}} = 1;
					}
				}
				
				if ($RepairSequencesIndelCount{$RepairSequence}){
				
					my 	$FullHDRIndel 	= $RepairSequencesIndelCount{$RepairSequence};
					my 	$Indels			=	scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"INDEL"}};
						$Indels			+=	scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}} if ($CutSiteReadGroupsRef->{$CutSite}{"REPAIR"});
			
					open 	REP, 	">>$RepairFile" or die ("Can't open $RepairFile\n");
					print 			"\t\t\t\tThere are $FullHDRIndel readpairs with FULL HDR and INDEL\n";	
					print 	REP 	"\t\t\t\tThere are $FullHDRIndel readpairs with FULL HDR and INDEL\n";

				}
				
				foreach my $OtherRepairSequence (keys %RepairSequencesIndelCount){
			
					unless ($OtherRepairSequence eq $RepairSequence){
					
						
						open 	REP, 	">>$RepairFile" or die ("Can't open $RepairFile\n");
						print 			"\t\t\t\tThere are $RepairSequencesIndelCount{$OtherRepairSequence} readpairs with INDEL and PARTIAL HDR: $OtherRepairSequence\n\n";	
						print 	REP 	"\t\t\t\tThere are $RepairSequencesIndelCount{$OtherRepairSequence} readpairs with INDEL and PARTIAL HDR: $OtherRepairSequence\n\n";	
						close 	REP;
					}
				}
			}
		}
	}
	
	#############################################
	### 	      Write out variants		  ###
	#############################################
	
	print VAR "Chromosome\tPositions\tType\tLength\tFlankingSequence (+ strand)\tAbsoluteFrequency\tRelativeFrequency\tNrOfReads\n";
	
	foreach my $CutSite (sort keys %$InDelsRef){
		foreach my $Chr (sort keys %{$InDelsRef->{$CutSite}}){
			foreach my $Positions (sort keys %{$InDelsRef->{$CutSite}{$Chr}}){
				foreach my $Type (sort keys %{$InDelsRef->{$CutSite}{$Chr}{$Positions}}){
					foreach my $Length (sort keys %{$InDelsRef->{$CutSite}{$Chr}{$Positions}{$Type}}){
					
						 my $AbsoluteFrequencyInNoConflictFragments	= scalar (keys %{$InDelsRef->{$CutSite}{$Chr}{$Positions}{$Type}{$Length}});
						 my $AbsoluteFrequencyInNoConflictReads		= 2*$AbsoluteFrequencyInNoConflictFragments;
						
						 my $NoConflictNrOfFragments					= (scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"NO_INDEL"}}) + 
																	  (scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"INDEL"}}) +
																	  (scalar keys %{$CutSiteReadGroupsRef->{$CutSite}{"REPAIR"}});
						 my $NoConflictNrOfReads						= 2*$NoConflictNrOfFragments;
						
						 my $RelativeFrequencyInNoConflictFragments	= Round ($AbsoluteFrequencyInNoConflictFragments/$NoConflictNrOfFragments);
						
						(my $Start, my $Stop) = split ("_", $Positions);
							$Start 	= $Start-10;
							$Stop 	= $Stop + 10;						
						
						 my $Sequence = qx{wget --quiet --output-document=- http://togows.org/api/ucsc/$Genome\/$Chr\:$Start\-$Stop};
						 
						 if (length($Sequence) > 20){
							
							$Sequence = substr($Sequence, 0, 10) . "[" . substr($Sequence, 10, length($Sequence)-20) . "]"	. substr($Sequence, 10+length($Sequence)-20, 10)
						 
						 }
						 else {
						 
							$Sequence = substr($Sequence, 0, 10) . "[]" . substr($Sequence, 10, 10)
						 }				
						
						 print VAR "$Chr\t$Positions\t$Type\t$Length\t$Sequence\t$AbsoluteFrequencyInNoConflictReads\t$RelativeFrequencyInNoConflictFragments\t$NoConflictNrOfReads\n";
					}
				}
			}
		}
	}	
}

sub CreatePossibleRepairSequences {
				
	(my $RepairSequence, my $PossibleRepairSequencesRef) 	= @_;
	 
	 my $CleanRepairSequence								= $RepairSequence;
		$CleanRepairSequence								=~ s/\(|\)|\[|\]//g;
	 my	$LastOpenBracket									= 0;
	 my	$LastClosedBracket									= 0;
	 my @Nucleotides										= ("A", "C", "G", "T");
	 my @FinalSequences										= ();
	
	$PossibleRepairSequencesRef->{$CleanRepairSequence}{$CleanRepairSequence} = undef;
	
	for (my $I = 0; $I < length ($RepairSequence); $I++){
		
		if (substr ($RepairSequence, $I, 1) eq "("){
				
			my 	$Change 			= "";
			my	@Changes			= ();
				$LastOpenBracket	= $I;
				
			$I++;
			
			do {
											
				$Change .= substr ($RepairSequence, $I, 1);
				$I++;
			
			}
			until (substr ($RepairSequence, $I, 1) eq ")");
			
			$LastClosedBracket = $I;
			
			### Create Hash With All Possible Changes ###
			
			if (length ($Change) == 1){
			
				@Changes = @Nucleotides;
			
			}
			else {
			
				@Changes = @Nucleotides;
				
				for (my $J = 1; $J < length ($Change); $J++){ 
				
					my @NewChanges = ();
					
					foreach my $Change (@Changes){
						foreach my $Nucleotide (@Nucleotides){
						
							push (@NewChanges, "$Change$Nucleotide");
							
						}
					}
					
					@Changes = @NewChanges;
				}
			}
			
			### Create Final Hash ###
			
			if (@FinalSequences){
			
				my @NewFinalSequences = ();
			
				foreach my $FinalChange (@FinalSequences){

					my $InterChangeSequence = substr($RepairSequence, length($FinalChange), $LastOpenBracket - length($FinalChange));

					foreach my $Change (@Changes){
					
						my $NewFinalSequence = $FinalChange . $InterChangeSequence . "(" . $Change . ")";
				
						push (@NewFinalSequences, $NewFinalSequence);
			
					}
				}
				
				@FinalSequences = @NewFinalSequences;
			}
			else {
			
				foreach my $Change (@Changes){
					
					my $FinalSequence = substr($RepairSequence, 0, $LastOpenBracket) . "(" . $Change . ")";
					
					push (@FinalSequences, $FinalSequence);
			
				}
			}
		}
	}
	
	foreach my $FinalSequence (@FinalSequences){
	
		$FinalSequence =~ s/\(|\)|\[|\]//g;

		$FinalSequence .= substr($RepairSequence, $LastClosedBracket+1, length($RepairSequence)-$LastClosedBracket-1);
		
		$FinalSequence =~ s/\(|\)|\[|\]//g;
		
		$PossibleRepairSequencesRef->{$CleanRepairSequence}{$FinalSequence} = undef;
	}
	return $PossibleRepairSequencesRef;
}

sub CalculateReadStopFromReadStartAndCigar {
			
	(my $ReadStart, my $Cigar) = @_;
	
	 my $CigarRef	= SplitCigar($Cigar);
	 my $ReadStop 	= $ReadStart -1;
	 
	 for (my $I = 0; $I < scalar @$CigarRef; $I+=2){ 
		
		if 		($CigarRef->[$I+1] eq "M" || $CigarRef->[$I+1] eq "D"){
		
			$ReadStop += $CigarRef->[$I];
		}
	 }
	 
	 return $ReadStop;
}

sub SplitCigar{
	
	my 	$Cigar 					= $_[0];
	my	@OriginalCigarSplited 	= split (//, $Cigar);
	my	@FinalCigar				= ();
	my 	$NumberOfBases			= "";
	
	foreach my $CigarValue (@OriginalCigarSplited){
		
		if ($CigarValue =~ m/[A-Z]/){
			push (@FinalCigar, $NumberOfBases);
			push (@FinalCigar, $CigarValue);
			$NumberOfBases = "";
		}
		else {
			$NumberOfBases .= $CigarValue;
		}
	}
	return \@FinalCigar;
}

sub CheckForRepairInIndelRead {
												
	(my $CutSitesPoolRef, my $RepairSequenceRef, my $PossibleRepairSequencesRef , my $ReadChr, my $RegionStart, my $RegionStop, my $Sequence, my $ReadPairName, my $RepairSequencesInIndelsRef, my $RepairSequencesIndelsConflictsRef) = @_;
	
	if (exists $RepairSequenceRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}){
		
		my $RepairSequence 	= $RepairSequenceRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}};
		
		foreach my $PossibleRepairSequence (keys %{$PossibleRepairSequencesRef->{$RepairSequence}}){
			
			my 	$PossibleRepairSequenceRC 	= reverse ($PossibleRepairSequence);
				$PossibleRepairSequenceRC 	=~ tr/ACGTacgt/TGCAtgca/;
			my 	$PossibleRepairSequenceQM	= quotemeta ($PossibleRepairSequence);
			my 	$PossibleRepairSequenceRCQM	= quotemeta ($PossibleRepairSequenceRC);
			
			if (! exists $RepairSequencesInIndelsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL_NO_REPAIR"}{$ReadPairName} &&
				($Sequence =~ m/$PossibleRepairSequenceRC/ || $Sequence =~ m/$PossibleRepairSequenceRCQM/)){
				
				$RepairSequencesInIndelsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL_REPAIR"}{$ReadPairName} = $PossibleRepairSequence;
				
			}
			elsif (! exists $RepairSequencesInIndelsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL_REPAIR"}{$ReadPairName} &&
					($Sequence !~ m/$PossibleRepairSequenceRC/ && $Sequence !~ m/$PossibleRepairSequenceRCQM/)){
			
				$RepairSequencesInIndelsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL_NO_REPAIR"}{$ReadPairName} = $PossibleRepairSequence;
			}
			else {
			
				$RepairSequencesIndelsConflictsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{$ReadPairName} = undef;
			}
		}
	}
	return $RepairSequencesInIndelsRef, $RepairSequencesIndelsConflictsRef;
}

sub DetermineAndStoreIndelCoordinate {

	(my $CigarRef, my $I, my $PositionInRead, my $Chr, my $ReadPairName, my $CutSite, my $InDelsRef) = @_;
	
	if ($CigarRef->[$I+1] eq "D"){
		
		my $DeletionStart 	= $PositionInRead;
		my $DeletionEnd		= $PositionInRead+$CigarRef->[$I]-1;
		my $Length			= $CigarRef->[$I];
		
		unless (exists $InDelsRef->{$CutSite}{$Chr}{"$DeletionStart\_$DeletionEnd"}{"DEL"}{$Length}{$ReadPairName}){
			$InDelsRef->{$CutSite}{$Chr}{"$DeletionStart\_$DeletionEnd"}{"DEL"}{$Length}{$ReadPairName} = undef;
		}
	}
	
	elsif($CigarRef->[$I+1] eq "I"){
	
		my $InsertionStart 	= $PositionInRead;
		my $InsertionEnd	= $PositionInRead-1;
		my $Length			= $CigarRef->[$I];
		
		unless (exists $InDelsRef->{$CutSite}{$Chr}{"$InsertionStart\_$InsertionEnd"}{"INS"}{$Length}{$ReadPairName}){
			$InDelsRef->{$CutSite}{$Chr}{"$InsertionStart\_$InsertionEnd"}{"INS"}{$Length}{$ReadPairName} = undef;
		}
	}
	
	return $InDelsRef;
}

sub StoreAsNoIndelOrRepair {
									
	(my $CutSiteReadGroupsRef, my $RepairSequenceRef, my $CutSitesPoolRef, my $ConflictReadPairsRef, my $ReadChr, my $RegionStart, my $RegionStop, my $ReadPairName, my $Sequence, my $PossibleRepairSequencesRef) = @_;
	
	if 		($RepairSequenceRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}){							# This cutsite has repairsite
	
		my $RepairSequence 	= $RepairSequenceRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}};
		my @RepairSequences	= ();
	
		foreach my $PossibleRepairSequence (keys %{$PossibleRepairSequencesRef->{$RepairSequence}}){
				
			my 	$PossibleRepairSequenceRC 	= reverse ($PossibleRepairSequence);
				$PossibleRepairSequenceRC 	=~ tr/ACGTacgt/TGCAtgca/;
			my 	$PossibleRepairSequenceQM	= quotemeta ($PossibleRepairSequence);
			my 	$PossibleRepairSequenceRCQM	= quotemeta ($PossibleRepairSequenceRC);
			
			if 	($Sequence =~ m/$PossibleRepairSequenceQM/ || $Sequence =~ m/$PossibleRepairSequenceRCQM/){
					 
				push (@RepairSequences, $PossibleRepairSequence);
			}
		}

		if	(!exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL"}{$ReadPairName} 		&&
			 !exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"REPAIR"}{$ReadPairName}		&&
			 !@RepairSequences){
		
			$CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"NO_INDEL"}{$ReadPairName} = undef;
		}
		elsif 	(!exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL"}{$ReadPairName} 		&&
				 !exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"NO_INDEL"}{$ReadPairName} 	&&
				 scalar @RepairSequences == 1){
				 
			$CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"REPAIR"}{$ReadPairName} = $RepairSequences[0];
		
		}
		else {
		
			$ConflictReadPairsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{$ReadPairName} = undef;
		}
	}
	elsif (!exists $CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"INDEL"}{$ReadPairName}) {
	
		$CutSiteReadGroupsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{"NO_INDEL"}{$ReadPairName} = undef;
	}
	else {
	
		$ConflictReadPairsRef->{$CutSitesPoolRef->{$ReadChr}{$RegionStart}{$RegionStop}}{$ReadPairName} = undef;
	}
	
	return $CutSiteReadGroupsRef, $ConflictReadPairsRef;
}

sub CleanUpFolders {

	print 		"\t\t\tClean up Folders\n";
	
	(my $Outputdir, my $CutSitesFileForSampleFile) = @_;
		
	if (-d "$Outputdir\/1_OriginalReads")		{system("rm -rf $Outputdir\/1_OriginalReads");}
	if (-d "$Outputdir\/2_QualityFilteredReads"){system("rm -rf $Outputdir\/2_QualityFilteredReads");}
	if (-d "$Outputdir\/3_RepairedReads")		{system("rm -rf $Outputdir\/3_RepairedReads");}
	if (-d "$Outputdir\/4_Mappings")			{system("rm -rf $Outputdir\/4_Mappings");}
	if (-d "$Outputdir\/5_CutSiteReads")		{system("rm -rf $Outputdir\/5_CutSiteReads");}
	if (-f "$CutSitesFileForSampleFile")		{system("rm $CutSitesFileForSampleFile");}
}

sub Round {

	my $Number 		= $_[0];
	my $Count		= 0;
	my $AfterComma	= 0;
	my $NewNumber	= "";
	
	for (my $I=0; $I<length($Number); $I++){
	
		my $Digit = substr($Number,$I,1);
		
		if ($Digit =~ m/\./){
		
			$AfterComma	= 1;
		}
		if ($AfterComma	== 0){
		
			$NewNumber	.= $Digit;
		}
		elsif ($AfterComma	== 1){
		
			$Count++;
			
			if ($Count < 6){
				
				$NewNumber	.= $Digit;
			}
			else {
				last;
			}
		}
	}
			
	return $NewNumber;
}

__END__

=head1 NAME

BATCH-GE.pl - A tool to performs batch analysis of NGS data for the assessment of knock-out and knock-in genome editing experiments

=head1 SYNOPSIS

Use:

    perl BATCH-GE.pl [--help] [--man] --Experimentfile=path_to_experimentfile
	
Examples:

    perl BATCH-GE.pl --help
    perl BATCH-GE.pl --ExperimentFile .../BATCH-GE/_EXAMPLE_1/ExperimentFile.csv

=head1 DESCRIPTION

BATCH-GE is a tool to performs batch analysis of NGS data for the assessment of knock-out and knock-in genome editing experiments. It calculates the efficiency of NHEJ and HDR mediated DSB repair.

=head1 ARGUMENTS

=over 4

=item --help (-?)

(Optional.) Displays the usage message.

=item --ExperimentFile=path_to_experimentfile

(Required.) The experiment file should contain all necessary information in order for the script to be able to run.
An example can be found in the example directory.

=back

=head1 AUTHOR

Wouter Steyaert, E<lt>wasteyae.Steyaert@UGent.beE<gt>

=head1 DATE

06-Jun-2016

=cut