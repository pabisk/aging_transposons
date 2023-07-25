#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
#https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl
#https://www.michaelgerth.net/news--blog/how-to-efficiently-bulk-download-ngs-data-from-sequence-read-databases
# This script automates the download of fastq files from the European Nucleotide Archive. Uses regular NCBI SRA accession numbers as input.
# Download instructions are from http://www.ebi.ac.uk/ena/browse/read-download 
# using a modified aria2c instead of wget https://github.com/aria2/aria2/issues/1039

my $aspera;
my $myid;
my $dl_count;
GetOptions(	'ascp'  => \$aspera,
		'id=s' => \$myid,
		'dl_count=s' => \$dl_count,
		);

print "$dl_count";
	
#Check if id file path was provided by user, if not, use default path from Ubuntu installation
if(!defined $myid){
	#$myid="\$HOME/kpabis/.aspera/connect/etc/asperaweb_id_dsa.openssh";
	$myid="/home/kpabis/.aspera/connect/etc/asperaweb_id_dsa.openssh";
}

# define input files
my $srafile = $ARGV[0];
my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
my $filename = $ymd."_wget.log";
if(!defined $ARGV[0]){
	print "This script will download fastq files from the European Nucleotide Archive database.\n\nUSAGE\tperl sra_dowload.pl [OPTIONS] accessions.txt\nOPTIONS\n  --ascp\tUse ascp (Aspera secure copy) instead of wget to fetch the files\n  --id\t\tPath to Aspera private key file (defaults to ~/.aspera/connect/etc/asperaweb_id_dsa.openssh, only in conjunction with --ascp)\n\nNOTES:\t-Requires EITHER wget version 1.16 & up (https://ftp.gnu.org/gnu/wget/) OR ascp (https://downloads.asperasoft.com/en/downloads/50)\n\t-Input file should contain only a single SRA accession number per line and no white spaces (accession numbers are identical to NCBI's SRA accession numbers)\n\t-This script will access European servers, and thus is probably most efficient when used in Europe.\n\t-The script will sort the file with accession numbers, if an unsorted one was provided.\n";
	}
else{
# Check if input files were provided, else show usage information
	system "sort -o $srafile $srafile";
	open (SRA , '<' , $srafile) or die "\n>FILE MISSING<\n\nplease specify path to SRA file (one accession number per line)!\n\n";
	my $count = 0;
	while( <SRA> ) { $count++; }
	printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
	if(!defined $aspera){print "\tFile apparently contains $count accession numbers. Beginning download with wget. Writing log file to $filename!"};
	if(defined $aspera){print "\tFile apparently contains $count accession numbers. Beginning download with ascp."};
	close(SRA);
	open (SRA , '<' , $srafile);
	my $libcount = 0;
	my $line;
	my $currentlib;
	while ($line = <SRA>) {
	# process each accession number separately
        	chomp $line;
		#count processed lines
                $libcount++;
		# Print date & time (here & in front of any output)
		my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
		printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
		# Print current library name
		print "\tDownloading library $line\.\n";
		# Check if SRA accesion numbers have the right length. If not, print error message and move to next entry. 
		if(length($line)<9 or length($line)>12){
			printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
			print "\tWARNING! $line does not appear to be a valid SRA accesion number. Please check!";
		}
		if(!defined $aspera){
			# When valid SRA accession number is provided, and aspera option is not used, call wget and download library. 
			# Definitions for folders in which fastq files are stored can be found under http://www.ebi.ac.uk/ena/browse/read-download  
			if(length($line)==12){
			print "12";
			    print "execute: data/aria3c -c -j $dl_count -x $dl_count  'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6).substr($line,9,3)."/".$line."/".$line."_1.fastq.gz'";
				system "data/aria3c -c -j $dl_count -x $dl_count  'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6).substr($line,9,3)."/".$line."/".$line."_1.fastq.gz'";
								system "data/aria3c -c -j $dl_count -x $dl_count  'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6).substr($line,9,3)."/".$line."/".$line."_2.fastq.gz'";
				
			}	
			if(length($line)==11){
			print "11";
			    print "execute: data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/0".substr($line,9,2)."/".$line."/".$line."_1.fastq.gz'";
		               	system "data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/0".substr($line,9,2)."/".$line."/".$line."_1.fastq.gz'";
		               			               	system "data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/0".substr($line,9,2)."/".$line."/".$line."_2.fastq.gz'";
			}	
			if(length($line)==10){
			print "10";
			    print "execute: data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/00".substr($line,9,1)."/".$line."/".$line."_1.fastq.gz'";
				system "data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/00".substr($line,9,1)."/".$line."/".$line."_1.fastq.gz'";
				system "data/aria3c -c -j $dl_count -x $dl_count 'https://ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/0".substr($line,9,1)."/".$line."/".$line."_2.fastq.gz'";
			}		
			if(length($line)==9){
			    print system "data/aria3c -c -j $dl_count -x $dl_count 'ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/".$line."/*'";
				system "data/aria3c -c -j $dl_count -x $dl_count 'ftp.sra.ebi.ac.uk/vol1/fastq/".substr($line,0,6)."/".$line."/*'";
			}
			$currentlib = $count-$libcount;
			# Print info on how many libraries were downloaded and how many remain
			my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
			my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
		        printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
		        print "\tRemaining donwloads: $currentlib of $count libraries." ;
			}
		if(defined $aspera){
			# Or, perform download with Aspera client.
			if(length($line)==12){
			system "ascp -QT -l 1000m -P33001 -i ".$myid." era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/".substr($line,0,6).substr($line,9,3)."/".$line."/. .";
			}	
			if(length($line)==11){
                       	system "ascp -QT -l 1000m -P33001 -i ".$myid." era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/".substr($line,0,6)."/0".substr($line,9,2)."/".$line."/. .";
			}	
			if(length($line)==10){
			system "ascp -QT -l 1000m -P33001 -i ".$myid." era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/".substr($line,0,6)."/00".substr($line,9,1)."/".$line."/. .";
			}		
			if(length($line)==9){
			system "ascp -QT -l 1000m -P33001 -i ".$myid." era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/".substr($line,0,6)."/".$line."/. .";
			}
			$currentlib = $count-$libcount;
			my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
			my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
                	printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
                	print "\tRemaining donwloads: $currentlib of $count libraries." ;
                	}
		}
	close(SRA);
	my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
        my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
	printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
	print "\tDownloaded $count libraries. Performing checks.\n";
	# Check if for every accession number in file, at least file/folder with the corresponding name. 
	# If not, print warning. (Could use more sophisticated check?)	
	my $missing= `ls * | cut -f1 -d'.' | cut -f1 -d'_' | sort -u | comm -13 - $srafile | sed "s/ //g"`;	
	if($missing  ne ''){
		my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
                my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
		printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
		print("\tWARNING! The follwing SRA accessions were not downloaded. Please check files!\n\n$missing\n");
		
	}
	else{
		my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
                my $ymd = sprintf("%04d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min);
		printf("\n[%02d:%02d:%02d]", $hour, $min, $sec);
		print("\tAll done!\n");	
	}	
}

