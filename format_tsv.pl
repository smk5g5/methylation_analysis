use lib qw(/home/skhan/local_lib/);
use Sort::Naturally;
use Cwd;
use Parallel::ForkManager;


# my $num_args = $#ARGV + 1;
# if ($num_args != 1) {
  # print "\nUsage: format_tsv.pl directory \n";
  # exit;
# }

my $dir = getcwd;
print $dir,"\n";
chdir $dir;
my @all  = grep !/\.(pl|txt)/, glob("*h*");
# print join("\n",@all),"\n";
# exit;
my $pm = Parallel::ForkManager->new(6);

foreach my $d(@all){
$pm->start and next; # do the fork
chdir $d  if (-d $d);
my @tsv  = nsort (grep /allc/, (grep !/(_C|_L|\.pl)\.tsv/,glob("*.tsv")));
print join("\n",@tsv),"\n";
do_all($d,@tsv);
print "\n\n";
chdir $dir;
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;

# @them = nsort(qw(
   # foo12a foo12z foo13a foo 14 9x foo12 fooa foolio Foolio Foo12a
  # ));

## chrBase chr base strand coverage freqC freqT
## 1 chr20.12314 chr20 12314 F 21 95.24 4.76
## 2 chr20.12378 chr20 12378 R 13 76.92 23.08
## 3 chr20.12382 chr20 12382 R 13 38.46 61.54
## 4 chr20.12323 chr20 12323 F 21 85.71 14.29
## 5 chr20.55740 chr20 55740 R 53 67.92 32.08

# 10      71      -       CG      1       1
# 10      78      -       CG      1       1
# 10      85      -       CG      1       1
# 10      92      -       CG      1       1
# 10      99      -       CG      1       1
# 10      106     -       CG      1       1
# 10      393     -       CHG     2       2
# 10      405     +       CHG     1       1
# 10      407     -       CHG     2       2
# 10      414     +       CHH     2       2


sub do_all{
my $prefix = shift;
my @tsv = @_;
chop $prefix;
# my $count = 0;
my $cg = "CG_"."$prefix".".txt";
my $chg = "CHG_"."$prefix".".txt";
my $chh = "CHH_"."$prefix".".txt";
open my $CG,">$cg";
open my $CHG,">$chg";
open my $CHH,">$chh";
print $CG join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
print $CHG join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
print $CHH join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";

foreach my $input(@tsv){
print $input,"\n";
open my $fl ,$input;
while(<$fl>){
chomp;
if($_=~/methylated/){next;}
my @chr = split /\s+/,$_;
my $chr;
if($chr[0]=~/\d{2}/){
$chr = "Gm".$chr[0];
}
else{
$chr = "Gm0".$chr[0];
}
my $chrbase = $chr.'.'.$chr[1];
my $str;
if($chr[2] eq '-'){$str = 'R';} else{$str = 'F';}
my $base = $chr[1];
my $cov =  $chr[5];
my $freqc =  sprintf("%.2f",100*($chr[4]/ $chr[5])); 
my $freqT = sprintf("%.2f", 100*(($chr[5]- $chr[4])/ $chr[5])); 

if($_=~/CG[ATCGN]/){
print $CG "$chrbase\t$chr\t$base\t$str\t$cov\t$freqc\t$freqT\n";
}
if($_=~/(C[ATCN][ATCN])|(CHH)/){
print $CHH "$chrbase\t$chr\t$base\t$str\t$cov\t$freqc\t$freqT\n";
}
if($_=~/(C[ATCN]G)|(CHG)/){
print $CHG "$chrbase\t$chr\t$base\t$str\t$cov\t$freqc\t$freqT\n";
}
}
}
}