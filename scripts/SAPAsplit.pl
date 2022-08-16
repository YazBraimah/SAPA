#!/usr/bin/env perl

########################
### Required Modules ###

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;
use File::Basename;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Alignment::Clustalw;
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);
# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;

my $usage = <<__EOUSAGE__;

##########################################################################################################################################################################
#   
#     ######     ###    ########     ###    
#    ##    ##   ## ##   ##     ##   ## ##   
#    ##        ##   ##  ##     ##  ##   ##  
#     ######  ##     ## ########  ##     ## 
#          ## ######### ##        ######### 
#    ##    ## ##     ## ##        ##     ## 
#     ######  ##     ## ##        ##     ## 
#
#           _ _ _   
# ___ _ __ | (_) |_ 
#/ __| '_ \| | | __|
#\__ \ |_) | | | |_ 
#|___/ .__/|_|_|\__|
#    |_|            
#
#
#	SAPAsplit v1.0 (Sequence Alignment and PAML Analysis) aligns coding sequences, translates them into proteins, calculates pairwise Ka/Ks, and performs PAML analyses.
#
#					
#	Required arguments:
#
#	--CDS_file|c <string>				sequence file containing CDS from multiple species.
#
#	--gene_id|i <string>				individual gene_id to be analyzed
#   
#	--output|o					name of directory for outputs 
#
#
#
#   Options:
#
#	--show_samples					(optional) display species IDs in CDS_file and exit (only requires "-c <CDS_file>" argument).
#
# 	--include_samples 				include specified species/lines (specify with --samples_file)
#
#	--exclude_samples 				exclude specified species/lines (specify with --samples_file)
#
#       --samples_file <string>				(use line ID if present, otherwise use species ID).
#
#	--save_alignments				(optional) retain alignment files
#
#	--save_split_alignments				(optional) retain split alignment files
#
#	--split 					window,offset
#		--refseq				reference sequence ID (species_line)
#
#	--run_PAML					Perform PAML codeml analyses
#
#		--tree_file				(required) newick format phylogeny
#		
#		*** Running PAML will execute the default NULL, i.e. model = 0, NSsites = 0. To execute the branch model and/or the branch-site model(s), specify as follows:
#
#		--branch	(optional)
#		--branchSite	(optional)	
#		--restrict_samples			File with species IDs to use in branch-type PAML tests. Required if --branch or --branchSite is specificied. 
#									One ID per line. must match species ID in *.tree file and FASTA file.
#
#	--help						
#                       
##########################################################################################################################################################################


__EOUSAGE__

    ;


#####################################
## Set BioPerl alignment protocol ###

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;


###################################
### Define command line options ###

my $CDS_file;
my $gene_id;

my $output_dir;

my $show_samples_flag = 0;
my $include_samples_flag = 0;
my $exclude_samples_flag = 0;
my $samples_file;

my $save_alignments_flag = 0;
my $save_split_alignments_flag =0;

my $calculate_KaKs_flag = 0;
my $split_param;
my $ref_seq;
my $run_PAML_flag = 0;
my $tree_file;
my $branch_flag = 0;
my $branchSite_flag = 0;
my $restrict_samples;


my $help_flag = 0;


&GetOptions (  

    'CDS_file|c=s' => \$CDS_file,
    'gene_id|i=s' => \$gene_id,
    'output|o=s' => \$output_dir,
    
    'show_samples' => \$show_samples_flag,
    'include_samples' => \$include_samples_flag,
    'exclude_samples' => \$exclude_samples_flag,
    'samples_file=s' => \$samples_file,
    
    'save_alignments' => \$save_alignments_flag,
    'save_split_alignments' => \$save_split_alignments_flag,
    
    'calculate_KaKs' => \$calculate_KaKs_flag,
    'split=s' => \$split_param,
    'refseq=s' => \$ref_seq,
    'run_PAML' => \$run_PAML_flag,
    'tree_file=s' => \$tree_file,
    'branch' => \$branch_flag,
    'branchSite' => \$branchSite_flag,
    'restrict_samples=s' => \$restrict_samples,
    
    'help' => \$help_flag,
    
);


#####################################################
### Check command line arguments and housekeeping ###

# if command line argument are not recognized
if ($help_flag) {
    die $usage;
}

# only show samples in CDS_file
if ($show_samples_flag) {
    system ("grep '>' $CDS_file | sed 's/.*species=/species=/g' | sort -u");
	exit(0);
}

# if Mandatory arguments not present
unless ($gene_id || $output_dir) {
	die "\nNo CDS file or gene ID or gene.\n\nAnd you need to specify output directory!\n $usage";
}

if (@ARGV) {
	die "Error, do not recognize params: @ARGV ";
}

# specify output dir cname
if ($output_dir) {
	mkdir($output_dir);
}

# if calculate KaKs, make tmp dir
if ($calculate_KaKs_flag) {
	system ("mkdir KaKs.tmp");
}

# if run PAML, make tmp dir
if ($run_PAML_flag) {
	unless ($tree_file) {
		die "\nError, Need to provide tree file";
	}
	system ("mkdir PAML.output");
}



######################################
###### Main execution of processes ###

main: {
				
		my @gene_object;
		
		###############################################################
		##### Set the gene object to a single gene or list of genes. ##
		push (@gene_object, $gene_id);
		
		
		###############################################################
		
		my @strain_object;
			
			if ($restrict_samples) {
        	open (FILE1, $restrict_samples) or die ("Could not open file \n");
			while ($restrict_samples = <FILE1>){
				chomp $restrict_samples;
				my @file_bits2 = split(/\s+/,$restrict_samples);
				my $sample_names = $file_bits2[0];
				push (@strain_object, $sample_names);
        		}
			}
		###############################################################
        ### Initiate the loop to analyze individual transcripts #######
    
        foreach my $transcript (@gene_object){
            eval {
            ###########################################################
            ### Add orthologous sequences to the analysis...  #########
            ### Processes FASTA header to include species and #########
            ### line ID, if present.                          #########
            
            
                system ("fastagrep.pl -X $transcript $CDS_file > $transcript.fa");
                system ("sed -i.bak '/>/ s/ line=/_/g' $transcript.fa");
                system ("sed -i.bak '/>/ s/>.*species=/>/' $transcript.fa");
                #######################################################
                ### Restrict to a subset of samples ###################
                
                if ($include_samples_flag) {
                    system ("fastagrep.pl -f $samples_file $transcript.fa > $transcript.trimmed.fa");
                    system ("rm $transcript.fa");
                    system ("mv $transcript.trimmed.fa $transcript.fa");
                }

                if ($exclude_samples_flag) {
                    system ("fastagrep.pl -v -f $samples_file $transcript.fa > $transcript.trimmed.fa");
                    system ("rm $transcript.fa");
                    system ("mv $transcript.trimmed.fa $transcript.fa");
                }
            
            
            
            ##########################################################
            ### BioPerl process for sequence analysis ################
			
			# read in sequence	
			my $seqio = Bio::SeqIO->new(-file => "$transcript.fa",
									   -format => 'fasta');
			my %seqs;
			my @prots;
			
			# process each sequence
			while ( my $seq = $seqio->next_seq ) {
				$seqs{$seq->display_id} = $seq;
				
				# translate them into protein
				my $protein = $seq->translate();
				my $pseq = $protein->seq();
				if( $pseq =~ /\*/ &&
				$pseq !~ /\*$/ ) {
					warn("Oi! $transcript has a stop codon!!!");
				}
				
				# Tcoffee can't handle '*' even if it is trailing
				$pseq =~ s/\*//g;
				$protein->seq($pseq);
				push @prots, $protein;
			}
 			
 			# Warn if only 1 sequence present
			if( @prots < 2 ) {
				warn("$transcript - Need at least 2 CDS sequences to proceed");
			}

			# Align the sequences with clustalw
			my $aa_aln = $aln_factory->align(\@prots);
			
			# project the protein alignment back to CDS coordinates
			my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs); 
			my @each = $dna_aln->each_seq();
			
			# output alignments for downstream analysis
			my $out_dna = Bio::AlignIO->new(-file => ">$transcript.aln_dna.afa" ,
										-format => 'fasta');
			$out_dna -> write_aln($dna_aln);

			my $out_aa = Bio::AlignIO->new(-file => ">$transcript.aln_aa.afa" ,
										-format => 'fasta');
			$out_aa -> write_aln($aa_aln);
			
			# clean-up alignment header files
			system ("sed -i.bak '/>/ s/\\/.*//g' $transcript.*.afa");
			system ("rm *bak");
			
			##
			## Split the cDNA file into chunks
			if ($split_param){
				system ("msa_split $transcript.aln_dna.afa --refseq $ref_seq --in-format FASTA --windows $split_param --out-format FASTA --out-root $transcript.split");
				system ("ls -1 $transcript.split* > splitList.txt");
			}
			

			my @split_object; 

			if ($split_param) {
				my $split_list = 'splitList.txt';
			open (FILE1, $split_list) or die ("Could not open file \n");
			while ($split_list = <FILE1>){
				chomp $split_list;
				my @file_bits = split(/\s+/,$split_list);
				my $file_names = $file_bits[0];
				push (@split_object, $file_names);
			}
		}

		foreach my $splitSeq (@split_object){
            eval {
			############################################################
			### Perform the analyses specified in the command line #####
			
###############################################################################################			
			# Run PAML analyses
			if ($run_PAML_flag) {
				my $codeml_output = "codeml.ctl";
				my $codemlFile .= "seqfile = paml.phy\n";
    			$codemlFile .= "treefile = paml.tree\n";
    			$codemlFile .= "outfile = paml.out\n";
				$codemlFile .= "noisy = 9\n"; 
      			$codemlFile .= "verbose = 1\n"; 
      			$codemlFile .= "runmode = 0\n"; 
      			$codemlFile .= "seqtype = 1\n"; 
    			$codemlFile .= "CodonFreq = 2\n";
    			$codemlFile .= "clock = 0\n";
    			$codemlFile .= "aaDist = 0\n";
				$codemlFile .= "model = 0\n"; 
				$codemlFile .= "NSsites = 0\n"; 
				$codemlFile .= "icode = 0\n"; 
				$codemlFile .= "Mgene = 0\n"; 
				$codemlFile .= "fix_kappa = 0\n"; 
				$codemlFile .= "kappa = 2\n"; 
				$codemlFile .= "fix_omega = 0\n"; 
				$codemlFile .= "omega = 1\n"; 
				$codemlFile .= "fix_alpha = 1\n"; 
				$codemlFile .= "alpha = 0\n"; 
				$codemlFile .= "Malpha = 0\n"; 
				$codemlFile .= "ncatG = 8\n"; 
				$codemlFile .= "getSE = 0 \n"; 
				$codemlFile .= "RateAncestor = 1\n"; 
				$codemlFile .= "Small_Diff = .5e-6\n"; 
				$codemlFile .= "cleandata = 1\n"; 
				$codemlFile .= "method = 0 \n"; 
					open (my $ofh, ">$codeml_output") or die "Error, cannot write to $codeml_output";    
    				print $ofh $codemlFile;
    				close $ofh;

			my $samples_present = `grep \">\" $splitSeq | sed 's/> //g' | tr '\n' ',' | sed 's/,\$//g'`;
			chomp $samples_present;
			
			system ("tree_doctor -P $samples_present -n -N -t $tree_file > tmpPruned.tree");
			system ("fastaSortByName.pl $splitSeq > $splitSeq.dnaP.afa");

			# translate this CDS sequence, 
			system ("transeq -sequence $splitSeq.dnaP.afa -outseq $splitSeq.aaP.afa");
			system ("sed -i.bak '/>/ s/_1//g' $splitSeq.aaP.afa");

			# system ("fastaSortByName.pl $transcript.aln_aa.afa > $transcript.aaP.afa");
			system ("sed -i.bak '/>/ s/\\/.*//g' $splitSeq.*P.afa");
			system ("paml_prep.pl $splitSeq.aaP.afa $splitSeq.dnaP.afa -nogap -output paml > $splitSeq.phy");
			system ("sed -i.bak 's/paml.phy/$splitSeq.phy/g' codeml.ctl");
			system ("sed -i.bak 's/paml.tree/tmpPruned.tree/g' codeml.ctl");
			
			### 1. NULL MODEL
			system ("sed -i.bak 's/paml.out/$splitSeq.out/g' codeml.ctl");
			system ("codeml");
			system ("sed -i.bak 's/$splitSeq.out/paml.out/g' codeml.ctl");
			system ("mv $splitSeq.out PAML.output");
	
			if ($branch_flag || $branchSite_flag){
			
			
			foreach my $strain (@strain_object){
			if ($branch_flag){
			### 2. BRANCH MODELS 
			system ("sed -i.bak 's/model = 0/model = 2/g' codeml.ctl");	
				
				# edit tree
				system ("sed -i.bak 's/$strain/$strain\#1/g' tmpPruned.tree");
				# edit $codemlFile 
				system ("sed -i.bak 's/paml.out/$splitSeq.$strain.br.out/g' codeml.ctl");
				# run
				system ("codeml");
				# back to normal
				system ("sed -i.bak 's/$splitSeq.$strain.br.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/$strain\#1/$strain/g' tmpPruned.tree");
			system ("sed -i.bak 's/model = 2/model = 0/g' codeml.ctl");
			system ("mv $splitSeq.$strain.br.out PAML.output");
			}
			
			if ($branchSite_flag) {
			### 3. BRANCH-SITE MODELS
			system ("sed -i.bak 's/NSsites = 0/NSsites = 2/g' codeml.ctl");
			system ("sed -i.bak 's/model = 0/model = 2/g' $codemlFile");
			system ("sed -i.bak 's/fix_omega = 0/fix_omega = 1/g' codeml.ctl");
				
				# edit tree
				system ("sed -i.bak 's/$strain/$strain\#1/g' tmpPruned.tree");
				# edit $codemlFile 
				system ("sed -i.bak 's/paml.out/$splitSeq.$strain.brSt.H0.out/g' codeml.ctl");
				### run
				system ("codeml");
				### swith to H1
				system ("sed -i.bak 's/$splitSeq.$strain.brSt.H0.out/$splitSeq.$strain.brSt.H1.out/g' codeml.ctl");
				system ("sed -i.bak 's/fix_omega = 1/fix_omega = 0/g' codeml.ctl");
				### run 
				system ("codeml");
				### back to normal
				system ("sed -i.bak 's/$splitSeq.$strain.brSt.H1.out/paml.out/g' codeml.ctl");
				system ("sed -i.bak 's/$strain\#1/$strain/g' tmpPruned.tree");

			system ("sed -i.bak 's/NSsites = 2/NSsites = 0/g' codeml.ctl");
			system ("sed -i.bak 's/model = 2/model = 0/g' codeml.ctl");
			system ("mv $splitSeq.$strain.brSt.H*.out PAML.output");	
			}
			}
			}
			## Clean-up
 			# system ("rm $splitSeq.phy");
			system ("rm tmpPruned.tree");
			
			}
			
			
			# option to save the alignment files
			# if ($save_alignments_flag ) {
			# 	system ("mv $transcript.aln_dna.afa $output_dir/$transcript.DNA_alignment.afa");
			# 	system ("mv $transcript.aln_aa.afa $output_dir/$transcript.Protein_alignment.afa");
			# } else {
			# 	system ("rm $transcript.aln*.afa");
			# }
			
			
			if ($save_split_alignments_flag) {
				system ("mkdir $output_dir/split_slignment");
				system ("mv $splitSeq*afa $output_dir/split_slignment");
				# system ("mv $transcript.aln_aa.afa $output_dir/$transcript.Protein_alignment.afa");
			} else {
				system ("rm $splitSeq*");
			}

			# system("rm $splitSeq")
			
			##########################################################
			############### END OF LOOP ##############################
		}		
		}
		# option to save the alignment files
			if ($save_alignments_flag) {
				system ("mv $transcript.aln_dna.afa $output_dir/$transcript.DNA_alignment.afa");
				system ("mv $transcript.aln_aa.afa $output_dir/$transcript.Protein_alignment.afa");
			} else {
				system ("rm $transcript.aln*.afa");
			}


		# clean-up
		system ("rm $transcript.fa");
	}
}
		# Format final output files
		if ($run_PAML_flag){
				system ("rm *bak");
				system ("mv PAML.output $output_dir/");
				system ("rm r* 2N* 4* lnf");
				system ("rm codeml.ctl");
			}

		if ($split_param){
				system ("rm splitList.txt");
			}

		
	
}

