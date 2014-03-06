use lib qw(/home/skhan/local_lib/);
use Sort::Naturally;
use Cwd;

my $num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: format_tsv.pl directory \n";
  # exit;
}

my $dir = $ARGV[0] || getcwd;
# print $dir,"\n";
chdir $dir;
my @all  = grep !/(hyper|hypo|\.pl)/,glob("*h*");
# print join("\n",@all),"\n";
foreach my $d(@all){
chdir $d;
# print getcwd,"\n";
my @tsv  = nsort grep/allc/,(grep !/_(C|L)\.tsv/,glob("*.tsv"));
# print join("\n",@tsv),"\n";
foreach my $ts(sort @tsv){
my $dr = $dir."/".$d."/";
my $bsub = (split/\./,$ts)[0];
print $d,"\n";
# print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $bsub -o $bsub.o%J -e $bsub.e%J -q normal \"perl $dir/methPipe.pl $dr $d $ts\"\n\n";
# system  "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $bsub -o $bsub.o%J -e $bsub.e%J -q normal \"perl $dir/methPipe.pl $dr $d $ts\"";
# system "perl $dir/methPipe.pl $dr $d $ts";
chdir $dir;
# exit;
}
}