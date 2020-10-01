# /usr/bin/perl

use strict;
use warnings;

#phosphorylation sites -  regex!
#protein kinase c [ST].[RK] --> S|T[A-Z]R|K
#casein kinase ii [ST]..[DE] --> S|T[A-Z]{2}D|E
#cAMP and cGMP [RK][RK].[ST] --> R|K R|K[A-Z]S|T
my $acc_num = shift;

if (!open(TREMBL, "sequencesTrembl.txt")) {
	 die "Check input file\n";
}
my $string = join("", <TREMBL>);
my @array = split(">", $string);



#how many proteins contain at least one occurence of each of the sites?
sub site_find{
  my $search = shift;
  my $length = length($search);
  my $PKC = 0;
  my $CKII = 0;
  my $cNMP = 0;
  if (defined $search){
    while ($search =~ /((S|T)[A-Z](R|K))/g){
    $PKC++;
    if ($PKC > $length){last;}#fail fast
	}
    while ($search =~ /((S|T)([A-Z]{2})(D|E))/g){
  	$CKII++;
  	if ($CKII > $length){last;}
	}
  	while ($search =~ /((R|K)(R|K)[A-Z](S|T))/g){
		$cNMP++;
		if ($cNMP > $length){last;}
}
	}
return ($PKC, $CKII, $cNMP);
}

#using subroutine
my @PKCs;
my @CKIIs;
my @cNMPs;

foreach (@array){
my @x = site_find($_);
push (@PKCs, $x[0]);
push (@CKIIs, $x[1]);
push (@cNMPs, $x[2]);
}

my $counter=0;
foreach (@PKCs){
  if ($_ != 0){
    $counter++;
  }
}
print "Number of proteins containing at least one Protein kinase C site: $counter\n";
$counter=0; #reset
foreach (@CKIIs){
  if ($_ != 0){
    $counter++;
  }
}
print "Number of proteins containing at least one Casein kinase II site: $counter\n";
$counter=0; #reset
foreach (@cNMPs){
  if ($_ != 0){
    $counter++;
  }
}
print "Number of proteins containing at least one cAMP or cGMP site: $counter\n";

#which protein contains the largest number of occurences of each of the sites?
my %pkc;
my %ckii;
my %cnmp;
my $y=0;
foreach (@array){
	if ($_ =~ /(\w+)\s/){
		$pkc{$1} = $PKCs[$y];
		$ckii{$1} = $CKIIs[$y];
		$cnmp{$1} = $cNMPs[$y];
	}
	$y++;
}
my @k_ckii;
my @n_ckii;
foreach my $key (sort {$ckii{$a} <=> $ckii{$b}} keys %ckii) {
     push(@k_ckii, $pkc{$key});
		 push(@n_ckii, $key);
}
my $high_ckii = $k_ckii[-1];
my $ki_ckii = $n_ckii[-1];
print "Highest occurence of ckii: $ki_ckii \t $high_ckii sites\n";
my @k_cnmp;
my @n_cnmp;
foreach my $key (sort {$cnmp{$a} <=> $cnmp{$b}} keys %cnmp) {
     push(@k_cnmp, $cnmp{$key});
		 push(@n_cnmp, $key);
}
my $high_cnmp = $k_cnmp[-1];
my $ki_cnmp = $n_cnmp[-1];
print "Highest occurence of cNMP: $ki_cnmp \t $high_cnmp sites\n";
my @k;
my @n;
foreach my $key (sort {$pkc{$a} <=> $pkc{$b}} keys %pkc) {
     push(@k, $pkc{$key});
		 push(@n, $key);
}
my $high = $k[-1];
my $ki = $n[-1];
print "Highest occurence of pkc: $ki \t $high sites\n";

#take acc num as argument, look for that entry, then find all prot kin C sites
my $num=0;
my $protkinc=0;
my $seq="";

foreach (@array){
if ($_ =~ /($acc_num)\s/){
$protkinc=$PKCs[$num];
$seq =$_;
$num++;
}}
$seq =~ s/($acc_num)\W+//;
my @fen = split(/(S|T)\w(R|K)/m, $seq);
my $len=0;
foreach (@fen){
	 $len = $len + length($_);
	my $end = $len+3;
	print "PKC Start position: $len and end position $end for protein $acc_num \n";
$protkinc++;
}

print "PKCs for entry $acc_num: $protkinc protein kinase c phosphorylation sites";
#and print out start and end positions of each sites
