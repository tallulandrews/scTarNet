#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil);

my $VERSION = "0.0";
my $NAME = "Interactions Pipeline";

##### Read input #####

my ($graph, $write, $help, $version, $sample) = '';
my $expression_metafile = undef();
my $background_file = undef();
my $inputGene_file = undef();
my $prefix = "dCor_Interactions";
my $recurr_thresh = 0;
my $threads = 1;
my $multitest="bon";
my $R = "R";

GetOptions('g'=>\$graph, 'w'=>\$write, 'h'=>\$help,'v'=>\$version,'s'=>\$sample,'x=s'=>\$expression_metafile,'p=s'=>\$prefix,'r=i'=>\$recurr_thresh, 't=i'=>\$threads,'b=s'=>\$background_file,'i=s'=>\$inputGene_file,'graph'=>\$graph, 'write'=>\$write, 'help'=>\$help,'version'=>\$version,'sample'=>\$sample,'expression=s'=>\$expression_metafile,'prefix=s'=>\$prefix,'repeated=i'=>\$recurr_thresh, 'threads=i'=>\$threads,'background=s'=>\$background_file,'input=s'=>\$inputGene_file, 'multitest=s'=>\$multitest, 'R=s'=>\$R);
#print ("h=$help\nv=$version\nx = $expression_metafile\ni=$inputGene_file\np=$prefix\nb=$background_file\nr=$recurr_thresh\nt=$threads\ng=$graph\nw=$write\n\n");
if ($version) {
	print "This is: $NAME\n";
	print "Version: $VERSION\n";
	exit();
}
if ($help) {
	print "$NAME v$VERSION\n";
	print 'Summary: 
	This program uses the dcor statistic to identify regulatory interactions between a set of candidate
	transcription factors (or a random sampling of genes with similar expression levels). It compares results
	across a set of gene expression matrices (eg. different normalizations of the same original data) in parallel
	(if possible) and filters putative interactions based on their robustness. This program does NOT map/modify 
	gene names so please ensure all files use identically formatted gene names.
';
	print 'Arguments (short-form given priority in case of conflicts):
	-h/--help             print this summary to screen
	-v/--version          print the version of this program
	-x/--expression       file containing a list of all expression-matrix files to use
	-i/--input            file containing original list of candidate genes
	-p/--prefix           string prefix for all outputfiles & temporary files. (default="dCor_Interactions")
	-b/--background       file containing list of genes to sample from (optional)
	-r/--repeated         integer, the number of expression files a dependency/interaction must be identified in 
			      to be considered valid (default=1/2*number of expression files, rounded up)
	-s/--sample	      if present will randomly sample genes with similar expression levels as the provided input
			      and run the analysis on the sample only.
	-t/--threads	      integer, the maximum number of threads (parallel processes) to run on. (default=1)
	-g/--graph            if present will plot a graph(s) of the networks
	-w/--write            if present will write the adjacency list(s) of the plotted graph(s)
	--multitest	      multiple testing correction to apply. options: "fdr=VAL", "bon" (default="bon")
	-R		      path to R version to use
';
	exit();
}
if (!defined($expression_metafile)) {die "Must provide a file containing full path to each expression matrix file to consider (-x)\n"}
if (!defined($inputGene_file)) {die "Must provide a file containing genes to consider/to base sampling on (-i)\n"}
if (!($multitest eq "bon" || $multitest =~ /fdr=[\d\.]+/)) {die "Error: $multitest is not a valid multiple testing correction.\n";}
my $Rscript = $R."script";

my @expressionMatrix_files = ();
open(my $ifh, $expression_metafile) or die "Error: Cannot open expression-matrix metafile\n";
while (<$ifh>){
	chomp;
	if ($_ !~ /^\s*$/) { #Ignore blank lines
		push(@expressionMatrix_files,$_);
	}
} close($ifh);
if (scalar(@expressionMatrix_files) == 0) {
	die "Error: Expression-matrix metafile is empty\n"
}
# Default value of recurr_thresh/repeated/r
if ($recurr_thresh == 0) {
	$recurr_thresh = ceil(scalar(@expressionMatrix_files)/2)
}

##### Generate random sample #####
if ($sample) {
	my $sample_output_file = "$prefix\_Sampled_Genes.tmp\n";
	my $str_expression_files = join(" ", @expressionMatrix_files);
	if (defined($background_file)) {
		system("$Rscript Sample_Genes.R $inputGene_file $sample_output_file $background_file $str_expression_files");
	} else {
		system("$Rscript Sample_Genes.R $inputGene_file $sample_output_file \"all\" $str_expression_files");
	}
	if (-e $sample_output_file && -s $sample_output_file) { #If file exists & is non-empty
		$inputGene_file = $sample_output_file;
	} else {
		die "Error: Sampled zero genes.\n";
	}
}
##### Find Interactions in parallel #####
my @children = ();
my @interactions_files = ();
for (my $i = 0; $i < scalar(@expressionMatrix_files); $i++) {
	my $outputfile = "$prefix\_$i\_IntDepend_Output.Rd";
	my $max_children = $threads;
	if (scalar(@children) < $threads) {
		my $pid = fork();
		if (!defined $pid) {die "Error: fork() failed!";}
		if ($pid) {
			#Parent - keep track of children
			push(@children,$pid);
			push(@interactions_files,$outputfile);
		} else {
			#Child
			my $expr_file = $expressionMatrix_files[$i];
			system("$R CMD BATCH --no-save --no-restore '--args $expr_file $inputGene_file $prefix\_$i $multitest' Find_Interactions.R  Find_Interactions_$prefix\_$i.Rout");
			print("rm Find_Interactions_$prefix\_$i.Rout");
			exit();
		}
	} else {
		# wait for one to finish then do this one again.
		waitpid(shift(@children),0);
		$i--;
	}
}
foreach my $pid (@children) {waitpid($pid,0);} #Wait for all children to be finished.
##### Plot/Write Networks #####
my $analysis_args = "$prefix";
if ($write) {
	$analysis_args = $analysis_args." 1";
} else {
	$analysis_args = $analysis_args." 0";
}
if ($graph) {
	$analysis_args = $analysis_args." two-panel";
} else {
	$analysis_args = $analysis_args." none";
}
my $interactions_files_str = join(" ", @interactions_files);
$analysis_args = $analysis_args." $recurr_thresh $interactions_files_str";

system("$R CMD BATCH --no-save --no-restore '--args $analysis_args' Combine_Output.R $prefix\_Combine_Output.Rout");
print("rm $prefix\_Combine_Output.Rout");

##### Clean-Up intermediate files #####
my @garbage = glob("$prefix\_*\_IntDepend_Output.Rd");
print("rm @garbage");
