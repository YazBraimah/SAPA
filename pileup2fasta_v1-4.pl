#! /usr/bin/perl

## Program Info:
#
# Name: pileup2fasta
#
# Purpose: Intended to generate a fasta file, with minimum bias to a backbone, from a BAM file of a scaffolded
#          genome assembly.  It also generates the SNPs and INDELS from the reference -> test genome
#
# Author: John Nash
# Copyright (c) Public Health Agency of Canada, 2011,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged 
#   for use, and as long as the author/copyright attributions 
#   are not removed. It remains the property of the copyright holder.
#
# I learnt C in 1986 - I use curly brackets like an old fart. Bite me!
#
# History:
#
# 2011-02-23: v1.0 - Original...
# 2012-12-06: v1.1 - '-s' option fixed
# 2012-12-07: v1.2 - Handling mpileup files from bamfiles with multiple contigs in the reference file.
#                    (GFF output only)
# 2012-12-08: v1.3 - Fixed bug which tried to send data to a GFF file when none existed
# 2012-12-08: v1.4 -  Handling mpileup files from bamfiles with multiple contigs in the reference file.
#                    (FASTA output) 
##

## start of main function

use warnings;
use strict;
use Text::Wrap;
$Text::Wrap::columns = 70;

# Init variables:
my $title = "pileup2fasta";
my $version = "1.4";
my $date = "09 December, 2012";
my $author = 'john.nash@phac-aspc.gc.ca';

# error messages:
my $error_msg = "Type \"$title -h\" for help.";
my $bug_msg = "I've never seen this in a samtools mpileup file, so I would love to see it in order to figure out how to deal with it. Please report to $author.";


# Get and process the command line params:
my ($infile, $outfile, $gff_file, $min_base_coverage, $del_qual_limit, $snp_qual_limit, $verbose) 
  = process_command_line();

print "$title version:$version working...\n" if ($verbose);

## File handling:

## Does the input sequence file exist:  handle errors
if ($infile eq '') {
    die ("Error: (m)pileup file not supplied\n", $error_msg,"\n");
}
else {
  open (INFILE, $infile) 
    or die ("Error: Cannot open $infile\n", $error_msg,"\n");
  print STDERR "mpileup file: $infile\n" if ($verbose);
}

if ($outfile ne '')  {
  open (OUTFILE, ">$outfile") 
    or die ("Error: Cannot open $outfile\n", $error_msg,"\n");
  print STDERR "Output fasta file: $outfile\n" if ($verbose);
}
else {
  print STDERR "No fasta output desired\n" if ($verbose);
}
	
if ($gff_file ne '')  {
  open (GFF_FILE, ">$gff_file") 
    or die ("Error: Cannot open $gff_file\n", $error_msg,"\n");
  print STDERR "Output GFF file: $gff_file\n" if ($verbose);
}
else {
  print STDERR "No GFF output desired\n" if ($verbose);
}

## Run the main pileup parsing routing:
&process();
exit (0);
## end of main function


sub help {
print <<EOHelp;
$title $version (released $date)
Please send bug reports to $author.

Usage: $title -i MPILEUP file (obligatory) -o FASTA file (optional) -g GFF file (optional)
       -b NUMBER (default = 8) -v NUMBER (default = 0.85) -s NUMBER (default = 0.90) -V (optional)
 
PILEUP file is the output from samtools mpileup command run like:
\"samtools mpileup -f [reference sequence] [BAM file] >MPILEUPfile\"

Old-versions of the deprecated tool pileup work too.

FASTA file is the assembled sequence guided by the reference sequence\'s scaffold but including SNPs and indels. Differences are in lower case.

GFF file is a features file containing the SNP, INSERTION and DELETION information

-b (minimum base coverage) is the minimum number of NGS sequence read coverage before a SNP or INDEL call will be made
-v (threshold for insertions) is the fraction of allowed positive calls to enable an INSERTION call 
-s (threshold for SNPs) is the fraction of allowed positive calls to enable a SNP call

-V (verbose output) - lets you know what's happening.

If this scrolls by too fast type: \"$title \| more\".
EOHelp
die ("\n");
} # end of sub help


sub process_command_line {
# Variables:
  my %opts = ();    # command line params, as entered by user
  my @cmd_line;     # returned value
  my @list;         # %opts as an array for handling
  my $cmd_args;	    # return value for getopts()
  my $item;
  
# Holders for command line parameters:
  my $infile = '';
  my $outfile = '';
  my $gff_file = '';
  my $min_base_coverage = 8;
  my $del_qual_limit = 0.85;  
  my $snp_qual_limit = 0.90;
  my $verbose = 0;

# Get the command=line parameters:
  use vars qw($opt_i $opt_o $opt_g $opt_b $opt_v $opt_s);
  use Getopt::Std;
  $cmd_args = getopts('i:o:g:b:v:s:hV', \%opts);
  
# Die on illegal argument list:
  if ($cmd_args == 0) {
    die ("Error: Missing or incorrect command line parameter(s)!\n",
	 $error_msg, "\n");
  }
  
# Make the hashes into an array:
  @list = keys %opts;
  
# Do a quick check for "help" and the compulsory parameters:
  foreach $item (@list)  {

# Input file:
    if ($item eq "i") { $infile = $opts{$item}; }
# Output file:
    elsif ($item eq "o") { $outfile  = $opts{$item}; }
# GFF file:
    elsif ($item eq "g") { $gff_file  = $opts{$item}; }
# Minimum reads for acceptable coverage:
    elsif ($item eq "b") { $min_base_coverage = $opts{$item}; }
# Threshold level for insertions:
    elsif ($item eq "v") { $del_qual_limit = $opts{$item}; }
# Threshold level for SNPs:
    elsif ($item eq "s") { $snp_qual_limit = $opts{$item}; }
# Help:
    elsif ($item eq "h") { help(); }
# Verbose:
    elsif ($item eq "V") { $verbose = 1; }
  }
  
# Put it in an array:
  @cmd_line = ($infile, $outfile, $gff_file, $min_base_coverage, $del_qual_limit, 
	       $snp_qual_limit, $verbose);
  return @cmd_line;
	
} #end of sub process_command_line()


sub process {
# first field in an mpileup entry is the entry name:
  my $entry_name;
  my $old_entryName;

# Are we at the start of the file? 1 = YES 0 = NO
my $firstLine = 1;

# How many lines have we looked at (debug code):
  my $count;
  my $oldCount = 1;
#  my $lineCount = 0;

# final sequence:
  my $final_sequence;

# Where all of the entries are going:
  my @outgoing;

# Read the pileup file, line-by-line:
  while(<INFILE>) {

# Count the lines processed from the mpileup file:
    $count++;

# Fix up screwed up DOS/Windows CR/LF: 
   s/\r\n/\n/g;
    chomp;

# Check which entry is being processed, if the original bam file was referenced to a 
#  multiple fasta file:
    my @line = split /\t/;

# Because of the various ways that different OSes handle EOL,
#  I'm adding my own because I want to split the files downstream based on "\n":
    my $line = "$_\n";

# mpileup files (w/o consensus) have 6 fields per line:
    die ("Not an mpileup consensus file\n") if scalar @line != 6;
     $entry_name = $line[0];

 # Special case for first line:
    if ($firstLine) {

# This is where the first line of where the first entry_name appears:
      print STDERR "\nWorking on: ", $entry_name, "\tLine: ", $count, "\n" if ($verbose);
      $old_entryName = $entry_name;
      $firstLine = 0;
    }

# Work on each entry in a multiple sequence entry, one-by-one:
    if ($old_entryName eq $entry_name) {
      push @outgoing, $line;
    }
    else {
      $final_sequence =  &process_line(@outgoing);
      if ($outfile ne '') {
	print OUTFILE '>', $entry_name, " Output from ", $infile, " by $title\n";
	print OUTFILE wrap '', '', $final_sequence, "\n";
      }

# This is the first element of the new $entry_name. Zero it out.  Start afresh:
      print STDERR "Working on: ", $entry_name, "\tLine: ", $count, "\n" if ($verbose);
      @outgoing = ();
      push @outgoing, $line;
      $oldCount = $count;
    } # end of last else
    
    $old_entryName = $entry_name;
  } # end of while(<INFILE>) 
  
# Deal with last contig...
  $final_sequence =  &process_line(@outgoing);
  if ($outfile ne '') {
    print OUTFILE '>', $entry_name, " Output from ", $infile, " by $title\n";
    print OUTFILE wrap '', '', $final_sequence, "\n";
  }

# Final status report:
  print STDERR "\n", $count, " lines in ", $infile, " were processed.\n" if ($verbose);
}  # end of sub process


sub process_line {

# first field in an mpileup entry is the entry name:
  my $entry_name;

# second field is the position:
  my $pos = 0;
  my $old_pos = 0;

# third field is the base:
  my $ref_sequence;

# fourth field is the number of reads:
  my $num_base_reads;

# fifth field is the encoded bases:
  my $variants;
  my $calc_var_qual = 0.0;

# The calculated sequnce:
  my $final_sequence;

# sixth field is the encoded quality string. I could use it but have had no need to yet.
 
# @incoming is the array of all of the reads from each $entry_name:
  my @incoming = @_;
  foreach (@incoming) {
    chomp;

# This helps me decide how to treat deletions which are subbed by Ns:
    my $dumpbase_flag = 0;
 
# Process each line:
    my @line = split /\t/;

# Assign the various mpileup entries to variables:
    $entry_name = $line[0];

# Positional checks in case the incoming mpileup file has deletions:
    $old_pos = $pos;
    $pos = $line[1];

# To eliminate sequences which do not start at Position 1 being considered as deletions from the start, 
# i.e. for analysis of regions, we have to reset $old_pos for the first run:
    if (($old_pos == 0) && ($pos > 1)) {
      $old_pos = ($pos - 1);
    }

    $ref_sequence = $line[2];
    $num_base_reads = $line[3];
    $variants = $line[4];

## Check position accuracy - part one:
    if ($pos == $old_pos) {
      die ("Error: unknown position error type 1. ", $bug_msg, "\n");
    }

## Check position accuracy - part two:
    if ($pos <= $old_pos) {
      die ("Error: unknown position error type 2. ", $bug_msg, "\n");
    }

# Find deletions and pad with * characters:
    if ($pos != ($old_pos + 1))  {
      for ($old_pos+1 .. $pos-1) {
	$final_sequence .= '*';
      }	
      my $deletion = $pos-1 - $old_pos + 1;
      print GFF_FILE "$entry_name\t$title\tDELETION\t$old_pos\t$pos\t.\t.\t.\tSize=$deletion\n" 
	if ($gff_file ne '');
    } # end of if ($pos != ($old_pos + 1))

    my $sign;
    my $digit;
    my $match;

# If we have the desired base coverage, continue:
    if ($num_base_reads >= $min_base_coverage) {
    
# sub out the characters beginning with ^ as they are irrelevant:
      $variants =~ s/\^.{1}//g;
      $variants =~ s/\$//g;

# do not consider strings where there are no variants, i.e. all . or ,:
      if ($variants =~ /[.,]{$num_base_reads}/) {
	if ($dumpbase_flag == 0) {
	  $final_sequence .= $ref_sequence;
	  $dumpbase_flag = 1;
	}
      }

# Deal with variants:
      if ($variants =~ /([\+\-]*)(\d*)([GATC\*]+)/gi) {
	$sign = $1;
	$digit = $2;
	$match = $3;

# Process INSERTIONS:
	if ($sign) {
	  my $var_count = 0;
	  $var_count++ while $variants =~ /[\+\-]\d[GATC\*]+/gi;
	  $calc_var_qual = $var_count / $num_base_reads;
	  if ($calc_var_qual >= $del_qual_limit) {
	    if ($gff_file ne '') {
	      print GFF_FILE "$entry_name\t$title\tINSERTION\t$pos\t$pos\t";
	      printf GFF_FILE "%0.2f", $calc_var_qual;
	      print GFF_FILE "\t.\t.\tInsertion=", uc($match), "\n";
	    } 
	    $final_sequence .= $match;
	  }
	  else {
	    if ($dumpbase_flag == 0) {
	      $final_sequence .= $ref_sequence;
	      $dumpbase_flag = 1;
	    }
	  }
	} # end of if ($sign)
      
# Process SNPs:
	if (!$sign)  {
	  my $var_count = length $match;
	  $match =~ /([GATC\*])/gi;
	  $match = $1;
	  $calc_var_qual = $var_count / $num_base_reads;
	  if ($calc_var_qual >= $snp_qual_limit) {
	    if ($gff_file ne '') {
	      print GFF_FILE "$entry_name\t$title\tSNP\t$pos\t$pos\t";
	      printf GFF_FILE "%0.2f", $calc_var_qual;
	      print GFF_FILE "\t.\t.\tAlleles=$ref_sequence->", uc $match, "\n";
	    }	      
	    $final_sequence .= $match;
	    $dumpbase_flag = 1;
	  }
	  else {
	    if ($dumpbase_flag == 0) {
	      $final_sequence .= $ref_sequence;
	      $dumpbase_flag = 1;
	    }
	  }
	} # end of if (!$sign)
      } # end of if ($variants =~ /([\+\-]*)(\d*)([GATC\*]+)/gi)
    } # end of if ($num_base_reads >= $min_base_coverage)

# Handle anything processed which failed... call it "reference":
    if ($dumpbase_flag == 0) {
      $final_sequence .= $ref_sequence;
      $dumpbase_flag = 1;
    }
  } # end of foreach (@incoming)

# Return the differential sequence from the contig:
  return $final_sequence;
} # end of sub process_lines
