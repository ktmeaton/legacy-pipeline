#!/usr/bin/perl

## Poinar Lab sequencing pipeline for ANCIENT samples
## Ana Duggan's parting gift 
## Last modifed Sept 22 2021

#use strict;
use warnings;
use Data::Dumper;

# ------------------------------
# Print Usage Information

if (scalar(@ARGV) < 2) 
{
    print "\n";
    print "Version : v0.1.2\n";
    print "Date    : 2021-08-25\n";
    print "Usage   : legacy.pl input_fastq_directory reference_genome(s) (modern)\n";
    print "Notes   : Redirect output to file called 'Makefile' and launch with command 'make'\n";
    print "          Reference genomes must be indexed using network aware bwa found at /home/ana/scripts/network-aware-bwa/bwa\n";
    print "          If the final parameter is 'modern' this will disable --ancientdna mode for leeHom.\n";
    print "\n";
    exit;
};


# ------------------------------
# Variable Declaration

my $directory = $ARGV[0];
my $ref = $ARGV[1];
my @refgenomes;
for (my $i=2; $i<scalar(@ARGV); $i++) 
{
    push(@refgenomes, $ARGV[$i]);
};


# The last argument will be checked for the keyword 'modern'
my $modern_flag = $ARGV[$#ARGV];
my $ancient_flag = "";
if($modern_flag ne "modern")
{
    $ancient_flag = "--ancientdna";
}
else
{
    # remove the word "modern" from the list of reference genomes
    pop @refgenomes;
}

my $fq2bam = "/home/ana/scripts/BCL2BAM2FASTQ/fastq2bam/fastq2bam";
my $bwa    = "/home/ana/scripts/network-aware-bwa/bwa";
my $samtools = "/home/ana/scripts/samtools-patched/sam";
my $fixpair = "/usr/local/biohazard/bin/bam-fixpair";
my $retrieve = "/home/gabriel/libbam/retrieveMapped_single_and_ProperlyPair";
my $rmdup = "/usr/local/biohazard/bin/bam-rmdup";
my $fld = "/home/gabriel/libbam/insertSize";
my $plot = "/home/ana/scripts/plotFLD.r";
my $dedup = "/home/keaton/scripts/NGSXRemoveDuplicates";
my @listTargets;
my $stringMakeToPrint="";

## Takes directory names from input 
my @sampleArray = readDirSample($directory);

#print Dumper(@sampleArray);

for my $directorySample (@sampleArray)
{
    my @mapped;
    my %extramapped;
    
    # Split directory name into elements on underscores
    my @fieldSample=split("_",$directorySample);

    # print "my directory sample: $directorySample\nmy filed sapmple: Dumper(@fieldSample)\n";
    ## Ensures there is more than one element or dies
    if($#fieldSample <  1)
    {  
	die "Unexpected number of fields for filename $directorySample, found ".$#fieldSample." fields ... exiting\n";
    }

    # Assigns final element of directory name as sample identifier
    my $sampleName= $fieldSample[$#fieldSample];

    # Reads through file names of directory, pushes them into an array
    my @gzArray = readDirFile($directorySample);

    # --------------------------------------------------
    # Create command to align each file to the reference
    
    for my $gzFile (@gzArray)
    {
        # Split file names on underscores
	my @fieldGZ=split("_",$gzFile);

	# Dies if there aren't the expected 5 fields (Sample name, index pair, Lane, Read, file #)
	if($#fieldGZ <  4)
        {
	    die "Unexpected number of fields for file $gzFile ... exiting\n";
	}	
        # Check if R1 file exists
	if($fieldGZ[$#fieldGZ-1] eq "R1")
        {
            my @fieldsToJoin = @fieldGZ;   

            # Swap R1 for R2 in name only
	    $fieldsToJoin[ $#fieldGZ-1 ] = "R2";

            # Join fields together
	    my $filer2 =  join("_",@fieldsToJoin);

            # If the R2 file doesn't exist (ie. single-end) just make the file name blank
            # Make a note in the Makefile that this sample is single-end
            if ( ! -e $filer2 )
            {
		print "# Single end sample: $gzFile\n"; 
             	$filer2 = "";
            }

	    my $fullRefPath =  $ref;
	    my @fieldREF = split("/", $fullRefPath);
            my @fullRefName = split("\\.", $fieldREF[-1]);
	    my $refName = $fullRefName[0];
	    my $tempOutFile = $gzFile.$refName.".sort.bam";

            # Convert fastq files to bam format, align to reference genome, sort
	    my $cmd =  "$fq2bam -r  $sampleName -o /dev/stdout ".$gzFile."\t".$filer2."  | leeHom $ancient_flag --log ".$gzFile."leehom.txt"." -o /dev/stdout  /dev/stdin | $bwa bam2bam -n 0.01 -o 2 -l 16500 -g $fullRefPath /dev/stdin | $samtools sort -o /dev/stdin ".$gzFile.".sort"." >$tempOutFile";
	    push(@listTargets,  $tempOutFile);
	    push(@mapped,  $tempOutFile);

	    # Add command to variable which will be printed in Makefile
	    $stringMakeToPrint = $stringMakeToPrint.$tempOutFile.":\n\t$cmd\n\n";

            # Repeat the command creation process for (optional) additional references
	    for (my $i=0; $i<scalar(@refgenomes); $i++) 
            {
	        my $extraRefPath = $refgenomes[$i];
	        my @extraFieldRef = split("/", $extraRefPath);
	        my @extraFullRefName = split("\\.", $extraFieldRef[-1]);
	        my $extraRefName = $extraFullRefName[0];
	        my $extraTempOutFile = $gzFile.$extraRefName.".sort.bam";
	        my $extraCmd = "$bwa bam2bam -n 0.01 -o 2 -l 16500 -g $extraRefPath <($fixpair -o /dev/stdout $tempOutFile) | $samtools sort -o /dev/stdin ".$gzFile.".sort"." >$extraTempOutFile";
	        push(@listTargets, $extraTempOutFile);
	        push(@{$extramapped{$extraRefName}}, $extraTempOutFile);
    	        $stringMakeToPrint = $stringMakeToPrint.$extraTempOutFile.": ".$tempOutFile."\n\t$extraCmd\n\n";
	    }
	}
        else
        {
            if($fieldGZ[ $#fieldGZ-1 ] eq "R2")
            {
	        #ignore
            }
            else
            {
	        die "Unexpected 5th field for file $gzFile ... exiting\n";
            }
        }
    }

    # ----------------------------------------------------------
    # Create commands to merge and post-process the aligned files
    
    my $fullRefPath =  $ref;
    my @fieldREF = split("/", $fullRefPath);
    my @fullRefName = split("\\.", $fieldREF[-1]);
    my $refName = $fullRefName[0];
    my $metaReadsFileName = $sampleName."_"."EditDist.min24.fasta";
    my $uniqMetaReadsFileName = $sampleName."_"."AllUniq.EditDist.min24.fasta";
    my $dedupStatsFileName = $sampleName."_"."DeDupStats.txt";
    my $mapFileName = $sampleName."_".$refName.".mapped.bam";
    my $fldFileName = $sampleName."_".$refName.".FLD.txt";
    my $plotFileName = $sampleName."_".$refName.".FLD.pdf";
    my $filterFileName = $sampleName."_".$refName.".min24MQ30.bam";
    my $exhaustionFileName = $sampleName."_".$refName.".Exhaust.txt";
    my $fileList = join(" ", @mapped);
    
    # Output for metagenomics
    #   1. Combine the multiple bam files for a single library, 
    #   2. Apply 24 bp length filter,
    #   3/ Use agrep to filter for remaining similarity to adapter sequences
    my $metaCmd = "$samtools merge $fileList | $samtools view -m 24 /dev/stdin | cut -f 1,10 | agrep -v -1 AGATCGGAA | agrep -v -1 TTCCGATCT | sed 's/^/>/g' | sed 's/\\t/\\n/g'  >$metaReadsFileName";
    push(@listTargets,$metaReadsFileName);
    $stringMakeToPrint = $stringMakeToPrint.$metaReadsFileName.": $fileList\n\t$metaCmd\n\n";

    # Remove duplicates based on string deduplication
    my $uniqCmd = "$dedup --fasta $metaReadsFileName --output $uniqMetaReadsFileName --stats $dedupStatsFileName";
    push(@listTargets,$uniqMetaReadsFileName);
    $stringMakeToPrint = $stringMakeToPrint.$uniqMetaReadsFileName.": $metaReadsFileName\n\t$uniqCmd\n\n";

    # Combine the multiple bam files for a single library 
    # (derived from the multiple fastq files where libraries are split across lanes 
    # and multiple files due to line limits) 
    # into a single bam file for a given library and reference combination
    my $mmfCmd = "$samtools merge $fileList | $retrieve /dev/stdin /dev/stdout | $rmdup -c -o $mapFileName /dev/stdin 1>$exhaustionFileName";  ## Combine the multiple bam files for a single library (derived from teh multiple fastq files where libraries are split across lanes and multiple files due to line limits) into a single bam file for a given library and reference combination
    push(@listTargets,$mapFileName);
    $stringMakeToPrint = $stringMakeToPrint.$mapFileName.": $fileList\n\t$mmfCmd\n\n";

    # Calculate the fragment length distribution
    my $fldCmd = "$fld $mapFileName | sort >$fldFileName ";
    push(@listTargets, $fldFileName);
    $stringMakeToPrint = $stringMakeToPrint.$fldFileName.": $mapFileName\n\t$fldCmd\n\n";

    # Filter for length and map quality
    my $minCmd = "$samtools view -b -m 24 -q 30 -o $filterFileName $mapFileName";
    push(@listTargets, $filterFileName);
    $stringMakeToPrint = $stringMakeToPrint.$filterFileName.": $mapFileName\n\t$minCmd\n\n";
   
    # Repeat the metagenomics prep commands for (optional) additional references 
    for (my $i=0; $i<scalar(@refgenomes); $i++) 
    {
        my $extraRefPath = $refgenomes[$i];
	my @extraFieldRef = split("/", $extraRefPath);
	my @extraFullRefName = split("\\.", $extraFieldRef[-1]);
	my $extraRefName = $extraFullRefName[0];
        my $extraMapFileName = $sampleName."_".$extraRefName.".mapped.bam";
	my $extraFLDFileName = $sampleName."_".$extraRefName.".FLD.txt";
	my $extraPlotFileName = $sampleName."_".$extraRefName.".FLD.pdf";
	my $extraFilterFileName = $sampleName."_".$extraRefName.".min24MQ30.bam";
	my $extraFileList = join(" ", @{$extramapped{$extraRefName}});
    	my $extraExhaustionFileName = $sampleName."_".$extraRefName.".Exhaust.txt";
	
        my $extraMmfCmd = "$samtools merge $extraFileList | $retrieve /dev/stdin /dev/stdout | $rmdup -c -o $extraMapFileName /dev/stdin 1>$extraExhaustionFileName";
	push(@listTargets, $extraMapFileName);
        $stringMakeToPrint = $stringMakeToPrint.$extraMapFileName.": $extraFileList\n\t$extraMmfCmd\n\n";
	my $extraFldCmd = "$fld $extraMapFileName | sort >$extraFLDFileName ";
        push(@listTargets, $extraFLDFileName);
	$stringMakeToPrint = $stringMakeToPrint.$extraFLDFileName.": $extraMapFileName\n\t$extraFldCmd\n\n";
	my $extraMinCmd = "$samtools view -b -m 24 -q 30 -o $extraFilterFileName $extraMapFileName";
	push(@listTargets, $extraFilterFileName);
	$stringMakeToPrint = $stringMakeToPrint.$extraFilterFileName.": $extraMapFileName\n\t$extraMinCmd\n\n";
    }
}

# --------------------------------------
# Write all commands to the Makefile
print "\nSHELL := /bin/bash\n\nDefault:\tall\n\n".$stringMakeToPrint."\n\nall:\t".join(" ",@listTargets)."\n\nclean:\n\trm -fv ".join(" ",@listTargets)."\n\n";



# ---------------------------------
# Subroutines

#Reads in list of directory names from stdin, 
sub readDirSample
{
    my ($dir) =  shift;
    my @arrayToReturn;
    opendir(D, $dir) || die "Can't open directory: $dir\n";
    while (my $f = readdir(D)) {
	#print "\$f = $f\n";
	if($f =~ /^Sample/){
	    push(@arrayToReturn,$dir."/".$f);
	}
    }
    closedir(D);
    
    return @arrayToReturn;
}

sub readDirFile
{
    my ($dir) =  shift;
    my @arrayToReturn;
    opendir(D, $dir) || die "Can't open directory: $dir\n";
    while (my $f = readdir(D)) {
	#print "\$f = $f\n";
	if($f =~ /gz$/){	    
	    push(@arrayToReturn,$dir."/".$f);
	}
    }
    closedir(D);
    
    return @arrayToReturn;
}

