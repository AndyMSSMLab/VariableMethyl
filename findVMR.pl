###To find VMRs using sliding window analysis on std dev.
use warnings;
use strict;

my $file = $ARGV[0];    #Input file
my $chr = $ARGV[1];     #Chromosome column (0 based)
my $start = $ARGV[2];   #Start column
my $end = $ARGV[3];     #End column
my $window = $ARGV[4];  #Size of window
my $sd_col = $ARGV[5];  #column containing std dev for each probe
my $sd_thres = $ARGV[6];        #Threshold for std dev to be considered significant
my $minp = $ARGV[7];    #Minimum number of singificant probes in the window
my $min_p_ratio = $ARGV[8];     #Minimum ratio of total probes to be significant in the window
my $prefix = $ARGV[9]; if(!defined($prefix)){$prefix = "";}

### Reading input file
open(IN, $file) or die "***Error*** Can not open the file $file";
my $i = 0; my $j = 0;
#my $header = <IN>;
my @data = ();
my $header = <IN>; chomp $header;
while(<IN>){
    chomp $_;
    my @s = split(/\t/,$_);
    for($j = 0; $j<@s; $j++){
        $data[$i][$j] = $s[$j];
    }
    $i++;
}
close IN;

my $nrow = $i;
my $ncol = $j;
undef $i;
undef $j;

my @sign_col = (0) x $nrow;

###Slinding window analysis
for(my $i = 0; $i<$nrow; $i++){
    my $sign_count=0;
    my $total_probes=0;
    my $j;
    for($j = $i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$i][$chr] eq $data[$j][$chr]; $j++){
        if($data[$j][$sd_col]>=$sd_thres){$sign_count++;}
        $total_probes++;
    }

###A sliding window is considered significant if it contians minimum "sign_count" number of significant probes
###and they make up minimum "min_p_ratio" of total probes in window
    if(($sign_count>=$minp) && ($sign_count/$total_probes>=$min_p_ratio)){
        for(my $k=$i; $k<$j; $k++){
            $sign_col[$k]=1;
        }
    }
}


###Printing output with additional column annotating the significant regions or VMR as 1 otherwise 0
print $header,"\t",$prefix,"_Sign_sd\n";

for(my $i = 0; $i <$nrow; $i++){
    for(my $j = 0; $j<$ncol; $j++){
        print $data[$i][$j],"\t";
    }
    print $sign_col[$i],"\n";
}
