my @data = glob("*Data");
# print join "\n",@data,"\n";
foreach my $dta(@data){
my $dat_name = (split/\./,$dta)[0];
print "bsub -R \"span[hosts=1]\" -J circos_$dat_name -o circos_$dat_name.o%J -e circos_$dat_name.e%J Rscript ideoDMC.r $dta $dat_name\n";
system "bsub -R \"span[hosts=1]\" -J circos_$dat_name -o circos_$dat_name.o%J -e circos_$dat_name.e%J Rscript ideoDMC.r $dta $dat_name";
}