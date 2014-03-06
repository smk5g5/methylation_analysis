use lib qw(/home/skhan/local_lib/);
use File::Basename;
use File::Find;
use Data::Dumper;

my @all_files = (); 
my $dir = '.';
find(\&print_name_if_dir, $dir);

my @chr = grep/\.txt/,@all_files;
my @kys = qw(RH-NF-12h  RH-mock-12h  RH_HS_24h_25oc  RH_HS_24h_40oC  SR-NF-12h  SR-mock-12h STR_HS_24h_25oC STR_HS_24h_40oC);
my %hash;


foreach my $ky(@kys) {
$hash{$ky}{'CHH'} = two_in_one('CHH_',$ky);
$hash{$ky}{'CHG'} = two_in_one('CHG_',$ky);
$hash{$ky}{'CG'} = two_in_one('CG_',$ky);
}

# print Dumper(\%hash),"\n";

Bsub('RH-mock-12h','RH-NF-12h','RH-mock_NF');
Bsub('RH-mock-12h','SR-mock-12h','RH-SR_mock');
Bsub('SR-mock-12h','SR-NF-12h','SR-mock_NF');
Bsub('RH-NF-12h','SR-NF-12h','SR-RH_NF');
Bsub('RH_HS_24h_25oc','RH_HS_24h_40oC','RH_HS_25oc_40oc');
Bsub('STR_HS_24h_25oC','STR_HS_24h_40oC','STR_HS_25oc_40oC');
Bsub('RH_HS_24h_40oC','STR_HS_24h_40oC','STR-RH_HS_40oC');
Bsub('RH_HS_24h_25oc','STR_HS_24h_25oC','STR-RH_HS_25oc');

# RH_HS_25oc_40oc

# Bsub('RH-mock-12h','RH-NF-12h','RH-mock_NF');
# Rscript two_sample_comp.r ./RH-NF-12h/CG_RH-NF-12h.txt ./SR-NF-12h/CG_SR-NF-12h.txt RH-NF-12h SR-NF-12h tiled_CG_SR-RH_NF

sub Bsub {
my $sample_1 = shift;
my $sample_2 = shift;
my $out_prefix = shift;
print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CHH_$out_prefix -o CHH_$out_prefix.o%J -e CHH_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CHH'} $hash{$sample_2}{'CHH'} $sample_1 $sample_2 CHH_$out_prefix CHH\"\n\n";
system "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CHH_$out_prefix -o CHH_$out_prefix.o%J -e CHH_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CHH'} $hash{$sample_2}{'CHH'} $sample_1 $sample_2 CHH_$out_prefix CHH\"";

print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CG_$out_prefix -o CG_$out_prefix.o%J -e CG_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CG'} $hash{$sample_2}{'CG'} $sample_1 $sample_2 CG_$out_prefix CG\"\n\n";
system "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CG_$out_prefix -o CG_$out_prefix.o%J -e CG_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CG'} $hash{$sample_2}{'CG'} $sample_1 $sample_2 CG_$out_prefix CG\"";


print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CHG_$out_prefix -o CHG_$out_prefix.o%J -e CHG_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CHG'} $hash{$sample_2}{'CHG'} $sample_1 $sample_2 CHG_$out_prefix CHG\"\n\n";
system "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J CHG_$out_prefix -o CHG_$out_prefix.o%J -e CHG_$out_prefix.e%J -q normal -n 6 \" Rscript two_sample_comp.r $hash{$sample_1}{'CHG'} $hash{$sample_2}{'CHG'} $sample_1 $sample_2 CHG_$out_prefix CHG\"";

# system "bsub -R \"rusage[mem=9000] span[hosts=1]\" -J 
#CHH_RH-mock_NF -o CHH_RH-mock_NF.o%J -e CHH_RH-mock_NF.e%J 
#-q normal -n 6 \"
#Rscript two_sample_comp.r $hash{'RH-mock-12h'}{'CHH'} 
#$hash{'RH-NF-12h'}{'CHH'} RH-mock-12h RH-NF-12h CHH_RH-mock_NF CHH\"";
}


############################

sub print_name_if_dir
{
    push(@all_files,$File::Find::name)if -f && $File::Find::name;
}

sub two_in_one{
my $pat_1  = shift;
my $pat_2  = shift;
my($tmp) = grep /$pat_2/,(grep/$pat_1/,@chr);
# print $tmp,"\n";
return($tmp);
}