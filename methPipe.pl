use lib qw(/home/skhan/local_lib/);
use Sort::Naturally;
use Cwd;

my $num_args = $#ARGV + 1;
if ($num_args != 3) {
  print "\nUsage: format_tsv.pl directory prefix file\n";
  # exit;
}
my $dir = $ARGV[0];
my $prefix = $ARGV[1];
my $file = $ARGV[2];

do_all($file,$prefix,$dir);

# chr1 839541 - CHH 0.142857 7
# chr1 839544 + CHH 0 10
# chr1 839548 - CHH 0 7
# chr1 839550 - CHH 0 8
# chr1 839551 + CG 0.75 8
# chr1 839552 - CG 0.75 8
# chr1 839554 + CHH 0 8
# chr1 839557 - CHH 0 7
# chr1 839566 + CG 0.857143 7
# chr1 839567 - CG 0.714286 7
# ./hmr -o Human_ESC.hmr Human_ESC.meth
# ./hmr -p Human_ESC.hmr.params -o Human_ESC.hmr Human_ESC.meth
# ./hmr -P Human_ESC.hmr.params -o Human_NHFF_ESC_params.hmr Human_NHFF.meth
# awk ’{$5=1-$5; print $0}’ Col0.meth > Col0_inverted.meth
# ./hmr -o Col0.hmr Col0_inverted.meth
# This kind of HyperMR analysis produces continuous blocks of hyper-methylated CpGs. However in some regions,
# intragenic regions in particular, such continuous blocks of hyper-methylated CpGs are separated by a few unmethylated
# CpGs, which have distinct sequence preference when compared to those CpGs in the majority of unmethylated
# genome. The blocks of hyper-methylated CpGs and gap CpGs together form composite HyperMRs. The hmr-plant
# program, which implements a three-state HMM, is used to identify such HyperMRs. Suppose the methcounts output
# file is Col0 Meth.bed, to find HyperMRs from this dataset, run
# ./hmr-plant -o Col0.hypermr Col0.meth
# ./hmr -partial -o Human_ESC.pmr Human_ESC.meth
# ./hmr -partial -o Human_ESC.pmr Human_ESC.meth



sub do_all{
my $tsv = shift;
my $chr = (split/\./,(split/_/,$tsv)[-1])[0];
my $prefix = shift;
my $dir = shift;
chdir $dir;
my $count = 0;
my $cg = "CG_".$prefix."_$chr".".txt";
my $chg = "CHG_".$prefix."_$chr".".txt";
my $chh = "CHH_".$prefix."_$chr".".txt";
my $outdir = "/home/skhan/data/stacey_methylation/meth_pipe_analysis/$prefix";
mkdir $outdir;
open my $CG,">$outdir/$cg";
open my $CHG,">$outdir/$chg";
open my $CHH,">$outdir/$chh";

# foreach my $input(@$tsv){
open my $fl ,$tsv;
while(<$fl>){
chomp;
$count++;
if($count==1)
{
# chromosome	position strand context estimated_methylation_level Total_number_of reads
print $CG join("\t",qw(chr position strand context methylation total)),"\n";
print $CHG join("\t",qw(chr position strand context methylation total)),"\n";
print $CHH join("\t",qw(chr position strand context methylation total)),"\n";
}
my @chr = split /\s+/,$_;
my $chr;
if($chr[0]=~/\d{2}/){
$chr = "Gm".$chr[0];
}
else{
$chr = "Gm0".$chr[0];
}
# my $chrbase = $chr.'.'.$chr[1];
my $base = $chr[1];
my $cov =  $chr[5];
my $freqc;
eval {
    $freqc = sprintf("%.2f",($chr[4]/$chr[5]));
    1;
} or do {
    @yNewVal = (0);
};
# my $freqT = sprintf("%.2f", 100*(($chr[5]- $chr[4])/ $chr[5])); 

if($_=~/CG[ATCGN]/){
print $CG "$chr\t$base\t$chr[2]\tCG\t$freqc\t$cov\n";
}
if($_=~/(C[ATC][ATCN])|(CHH)/){
print $CHH "$chr\t$base\t$chr[2]\tCHH\t$freqc\t$cov\n";
}
if($_=~/(C[ATCN]G)|(CHG)/){
print $CHG "$chr\t$base\t$chr[2]\tCHG\t$freqc\t$cov\n";
}
else{next;}
}

# }
}

#current format
# 1       42      +       CG      8       8

#needed format
#chromosome	position strand context estimated_methylation_level Total_number_of reads
# chr1 839537 - CG 0.8 5
