#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;

die "usage: find_targets.pl gbk predicted_ncRNA.seq output.picky match_len \n" if @ARGV!=4;

print STDERR "Find sRNAs which it has targets located in mRNA and not matched with polyU \n";

# to reverse and complement DNA fragment
sub reverse_complement {
  my $s = reverse($_[0]);
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
}

my @genes;
my $chromosome;
my $complement;
my ($g_name, $left, $right, $type) = ('', 0, 0, '');
#check with mRNA target should be in mRNA
sub check_mRNA {
    my ($start, $stop, $dir) = @_;
    my $rlt = 0;
    for (@genes){
        ($complement, $left, $right) = @{$_}[0..2];
        next if $complement != $dir;
        next if $right < $start;
        last if $left > $stop;
        if (($left <= $start && $right >= $start) ||
            ($left <= $stop && $right >= $stop) ||
            ($left <= $start && $right >= $stop) ||
            ($start <= $left && $stop >= $right)) {
            $rlt = 1; #match with mRNA should save
            last;
        }
    }
    return $rlt;
}

print STDERR "Read in GBK file\n";
die "can't open gbk file $ARGV[0]" unless open GBK, $ARGV[0];

while (<GBK>) {
  if(/^ACCESSION\s+(\.+)$/){
      $chromosome = $1;
      print STDERR "chrom: $chromosome\n";
  }  
  if (/\s+\bgene\b\s+(complement\()?(\d+)\.\.(\d+)/) {
    $complement = defined($1); $left = $2; $right = $3;
    $complement = 0 unless defined($complement);
  }elsif (/\s+\bgene\b\s+(join\()?(\d+)\.\.(\d+)/) {
     $complement = defined($1); $left = $2; $right = $3;
     $complement = 0 unless defined($complement);
  }elsif (/\s+\bgene\b\s+complement\(join\(/) {
    my @pos = split/[\s,\(,\.,\)]/, $_;
    $complement = 1; $left = $pos[19]; $right = $pos[24];
  } elsif (/^ORIGIN/) {
    last;
  } elsif ($left) {
    if (m:/gene=\"(\w+)\":) {
      $g_name = $1;
    } elsif (/\s+\bCDS\b\s+(complement\()?(\d+)\.\.(\d+)/) {
      die "mismatch complement $complement" if $complement != defined($1);
      warn "mismatch left $left != $2" if $left != $2;
      warn "mismatch right $right != $3" if $right != $3;
      $type = "mRNA";
      push @genes, [$complement, $left, $right, $g_name, $type];
      #my @temp = ($left,$right);
      #$gbk{$name} = 1;
      $left = 0;  
    } elsif (/\s+\bCDS\b\s+(join\()?(\d+)\.\.(\d+)/) {
      $type = "mRNA";
      push @genes, [$complement, $left, $right, $g_name, $type];
      $left=0;
    } elsif (/\s+\bCDS\b\s+complement\(join\(/) {
      $type = "mRNA";
      push @genes, [$complement, $left, $right, $g_name, $type];
      $left=0;
    } elsif (/\s+\b\w+\b\s+(\w+\()?(\d+)\.\.(\d+)/) {
#      print STDERR $_;
       $left = 0;
    }
  }
}
close GBK;
@genes = sort {$a->[1] <=> $b->[1]} @genes;
print STDERR scalar @genes, " mRNAs are read in!!\n";

print STDERR "Read in ncRNA sequence\n";
my $count=0;
my %ncRNAs;
my $row;
my %pos;
die "can't open predict_seq file $ARGV[1]" unless open RNA, $ARGV[1];
while (<RNA>) {
  chomp;
  next if /^\s*$/;
  
  if (/^>/){
	$row = $_;
  }else{
	#print STDERR "row: $row\n";
	my @fields = split(' ', $row);
	#print STDERR "4: $fields[4] 5:$fields[5] 6:$fields[6]\n";
 	$fields[0] =~ /^>(\w+)$/;
	$fields[0] = $1;
        my $sRNA = $fields[4]."_".$fields[5]."_".$fields[6];
	$pos{$sRNA} = [ ] unless exists $pos{$sRNA};
	$ncRNAs{$_} = [ ] unless exists $ncRNAs{$_};
	push @{$ncRNAs{$_}}, \@fields;
	$count++;
  }
  
}
close RNA;

print STDERR scalar keys %pos ," predicted ncRNAs positions ared read in \n";
print STDERR "$count ncRNA sequences read in, ",  $count-keys(%ncRNAs), " of these are not unique:\n";

print STDERR "finally read in Picky output\n";
my $line = $count = 0;
my %report;
my $seq;
my %with_targets;
my $no_match = 0;
my $rna_len = 0;
my $rna_temp=0;
my $match_len = $ARGV[3];
die "can't open PICKY file $ARGV[2]" unless open PICKY, $ARGV[2];
while (<PICKY>) {
  chomp;
  if (/^\S/) { #no space
    if ($line==0) { #sequence
	my @fields = split/\t/;
	unless (exists $ncRNAs{$fields[0]}) {
	    print STDERR "undefined seqence $fields[0]\n";
	    $no_match++;
	    $seq = "";
	}
	$seq = $fields[0];
	$rna_len = $fields[1];
	$rna_temp = $fields[4];
    } elsif (/^(U|S|RS|RU)/) { #perfect match just skipped
	next;    
    } else { #unperfect match
	next if length($seq) == 0;
	my @fields = split/\t/;

	my $name = $fields[6];

	my $direc = $fields[0] eq ">"?0:1;

	my ($start, $stop);
	if ($direc) {     
	    $start = $fields[5];
	    $stop = $fields[4];
	} else {
	    $start = $fields[4];
	    $stop = $fields[5];
	}

	my $m_len = $fields[3]-$fields[2] + 1;
        next if $m_len > $match_len;

	for (@{$ncRNAs{$seq}}) {
	    if ($_->[3] eq $name) { #only chromosome are same
	
  		#check with mRNA
		#my $rlt = should_save2($m_s, $m_e, $start, $stop, $direc); #put match position	
		my $mRNA_rlt = check_mRNA($start, $stop, $direc);
		#print STDERR "chiX : $mRNA_rlt\n" if $seq eq 'AGAGTATTGGCAGGATGGTGAGATTGAGCGACAATCGAGTTACACCGTCGCTTAAAGTGACGGCATAATAATAAAAAAATGAAATTCCTCTTTGACGGGC';
		if ($mRNA_rlt){
	            #print STDERR "seq $seq\n";
	            #print STDERR "cont, $_->[0] $_->[4]\n";
		    $count++;
		    my $sRNA = $_->[4]."_".$_->[5]."_".$_->[6];
		    #print STDERR "predict sRNA: $sRNA \n";
		    #get match position
		    my $m_start = $_->[1]+ $fields[2]-1;
		    my $m_end = $m_start + $fields[3];
		    $report{$sRNA} = [ ] unless exists $report{$sRNA};
		    push @{$report{$sRNA}}, [ $start, $stop, $direc, $m_start, $m_end ];
		}
	    }	    
	}		
    }
    $line++;
  }else{
      $line = 0;
  }
}
close PICKY;
print STDERR scalar keys %report, " has $count targets.\n\n";

print STDERR "now generates text output\n";
print "sRNA_S\tsRNA_E\tsRNA_D\tMatch_S\tMatch_E\tTarget_S\tTarget_E\tTarget_D\n";
open ncRNA_F, '>RNA_only_ml.txt';
print ncRNA_F "sRNA\tsRNA_S\tsRNA_E\tsRNA_Dir\tTarget_NO\tEtc\n";
foreach my $key (keys %report){
    my @sRNA  = split('_', $key);
    my $rna = '';
    for(@sRNA){
	$rna .= $_."\t";
    }
    my $list = $report{$key};
    for (@{$list}){
	print $rna,$_->[3],"\t",$_->[4],"\t",$_->[0],"\t",$_->[1],"\t",$_->[2],"\n";
    }
    my $target_no = scalar @{$list};
    print ncRNA_F "RNA\t$rna\t$target_no\tN/A\n";
}
close ncRNA_F;

print STDERR "Job is completed!!!!!\n";
