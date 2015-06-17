#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;

die "usage: analyze_wo_text.pl genome.seq baits.seq output.picky bait_length match_length \n" if @ARGV!=5;

my $bait_length = $ARGV[3];
my $match_length = $ARGV[4];

# this subroutine formats outlink URL to a sequence database based on
# the input sequence name pattern; you very likely will need to modify
# this subroutine to match the specific database related to your genes

# to complement DNA fragment
sub complement {
  my $s = $_[0];
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
}

# to reverse and complement DNA fragment
sub reverse_complement {
  my $s = reverse($_[0]);
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
}

# to compare if two nested data structures have the same data
sub compare {
  my ($a, $b) = @_;
  if (ref($a) eq "ARRAY") {
    return 0 unless ref($b) eq "ARRAY";
    return 0 unless @{$a}==@{$b};
    for (my $i=0; $i<@{$a}; ++$i) {
      return 0 unless compare($a->[$i], $b->[$i]);
    }
    return 1;
  } else {
    return $a eq $b;
  }
}

sub should_save2 {
  #print STDERR "Save $_ \n";
  my ($s, $e) = @_;
  my $len = $e-$s+1;
  return 0 if $len > $match_length;
 
  return 1;
}

print STDERR "first read in genomic sequences\n";
my ($name, $seq);
my %genome;
open GENOME, $ARGV[0];
while (<GENOME>) {
  chomp;
  if (/^>(\w+)/) {
    $genome{$name} = $seq if $name;
    $name = $1;
    $seq = "";
  } else {
    $seq .= $_;
  }
}
close GENOME;
$genome{$name} = $seq if $name;
print STDERR scalar keys %genome, " genomic sequences read in.\n\n";


print STDERR "next read in BAIT information\n";
my $count=0;
my %baits;
my $row = 0;
my $gene;
open BAIT, $ARGV[1];
while (<BAIT>) {
  chomp;
  my @fields = split;
  $fields[0] =~ /^>(.+)$/;
  $fields[0] = $1;
  my $bait = uc <BAIT>;
  chomp($bait);

  # optional code to double check; can be removed if unhappy
  my $operon_start = $fields[4];
  my $operon_end = $fields[5]; 
 
  my $part = substr($genome{$name}, $operon_start, $bait_length);
  #$part = &reverse_complement($part) if $fields[6] eq "reverse";

  die "mismatch BAIT: $operon_start \n$bait\n$part\n " unless $bait eq $part;
    
  # end of optional code for double checking.

  $baits{$bait} = [ ] unless exists $baits{$bait};
  push @{$baits{$bait}}, \@fields;
  $count++;
  
}
close BAIT;

print STDERR "$count BAIT sequences read in, ",  $count-keys(%baits), " of these are not unique:\n";

print STDERR "finally read in Picky output\n";
print STDERR "unanticipated additional UTR matches:\n";
my $line = $count = 0;
my @data;
my %report;
my $list;
$row = 0;
my $high_temp = 0;
open PICKY, $ARGV[2];
while (<PICKY>) {
  chomp;
  if (/^\s*$/) {
    if ($count) {
      #print STDERR "\tknown matches for this UTR are:\n";
      #for (@{$list}) {
	#print STDERR "\t$_->[3] ", $_->[1]<$_->[2]?">$_->[1]":"<$_->[1]", " $_->[0]\n";
      #}
    }
    if ($line>1) {
      if (!exists($report{$data[0][0]})) {
	  $report{$data[0][0]} = [ $list, @data[0..$line-1] ];       
      } else {
	die "same Baits with different Picky screening results?"
	    unless compare($report{$data[0][0]}, [ $list, @data[0..$line-1] ]);
      }
    }
    $line = $count = 0;
  } elsif ($line==0) {
    my @fields = split/\t/;
    die "undefined BAIT sequence $fields[0]" unless exists $baits{$fields[0]};
    $high_temp = $fields[4];
    $list = $baits{$fields[0]};
    $data[$line++] = \@fields;
  } elsif (/^(U|RS|RU)/) { #perfect match
      next;    
  } else { #unperfect match
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
    my $i = 0;
    for (@{$list}) {
	if ($_->[3] eq $name) {	
	    my $interval = (abs($_->[2]-$_->[1])+1)/2;
	    ##skip overlaps;
	    if ($_->[6] eq "reverse") {
		last if $direc==1 && $_->[2]+$interval> $stop && $_->[2]<=$start;
	    } else {
		last if $direc==0 && $_->[1]-$interval<=$start && $_->[1]>$stop;
	    }
	}
	++$i;
    }
    
    unless ($i<@{$list}) {
      my $rlt = &should_save2($start, $stop);    
       #print STDERR " $rlt :$_ \n" if $fields[4] == 3743025;
      if ($rlt) {
	$data[$line++] = \@fields;
      }
    }
  } 
}
close PICKY;
print STDERR scalar keys %report, " potential sRNAs to report.\n\n";

print STDERR "now generates text output\n";

$count = 1;
print "ncRNA\tsRNA_Loc\tncRNA_dir\tTarget_Loc\tTarget_dir\tTemp_Deff\tMatch_Len\tChromosome\n";
for (sort { $a->[0][0][0] cmp $b->[0][0][0] } values %report) {
  #last if $count>100;
  my @target_hits = @{$_}[2..$#{$_}]; ##from picky file
  my $fragment_len = $_->[1][1];
  my $t_temp = $_->[1][4];
  my $gene_name = '';
  # loop through each target gene in a group sharing the same 100 bp fragment
  for (@{$_->[0]}) { ##from baits

    # 1: target gene list
    $gene_name = $_->[0];

    # 2: target locations
    my $sRNA_direc = $_->[6] eq "reverse"?1:0;
    my $sRNA_start = $_->[1];
    my $sRNA_end = $_->[2];
       
    for (@target_hits) {
      # 3: sRNA location
	my $name = $_->[6];
	my $direc = 0;
	$direc = 1 if $_->[0] eq "RU" || $_->[0] eq "R>";

	my ($start, $stop);
      if ($direc) {
	  $start = $_->[5]+1;
	  $stop = $_->[4]+1;
      } else {
	  $start = $_->[4]+1;
	  $stop = $_->[5]+1;
      }
	die "wrong cordination of predicted ncRNAs. \n" if $start > $stop;
     
      print $gene_name, "\t",$sRNA_start+1,"-",$sRNA_end+1,"\t$sRNA_direc\t$start, $stop\t$direc\t",$t_temp-$_->[1],"\t",$_->[3]-$_->[2]+1,"\t",$name,"\n";
      
    }
    $count++; # increases the sRNA counter
  }
}

# reduce counter to actual total
$count--;
##put info into the file

print "\n", scalar keys %report, "potential sRNAs to report.\n";
print "$count sRNA groups identified.\n";

print STDERR "$count sRNA groups identified.\n";
print STDERR "Job is completed!!!!!\n";
