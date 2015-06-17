#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;
use Getopt::Long;

my $usage = <<"USAGE";
usage: $0 seq_file picky_file single_file extracted_info match_length has_name(0/1);
USAGE

die $usage unless @ARGV==6;


print STDERR "Find sRNAs which it has targets located in mRNA and not matched with polyU \n";

# to reverse and complement DNA fragment
sub reverse_complement {
  my $s = reverse($_[0]);
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
}

print STDERR "Read in sequences\n";
my $count=0;
my %ncRNAs;
my $row;
my %pos;
die "can't open seq_file $ARGV[0]" unless open SEQ, $ARGV[0];
while (<SEQ>) {
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
	$pos{$sRNA} = 1 unless exists $pos{$sRNA};
	$ncRNAs{$_} = [ ] unless exists $ncRNAs{$_};
	push @{$ncRNAs{$_}}, \@fields;
	$count++;
  }
  
}
close SEQ;

die "can't open singleLine_file $ARGV[2]" unless open SINGLE, $ARGV[2];
my %singles;
while (<SINGLE>) {
    next if /^sRNA/;
    chomp;
    my @fields = split /\t/;
    $singles{$fields[0]} = $fields[1];
}
close SINGLE;

print STDERR scalar keys %pos ," predicted ncRNAs positions ared read in \n";
print STDERR "$count ncRNA sequences read in, ",  $count-keys(%ncRNAs), " of these are not unique:\n";

print STDERR "finally read in Picky output\n";
my $line = $count = 0;
my %reports;
my $seq;
my $rna_len = 0;
my $rna_temp=0;
my $match_len = $ARGV[4];
die "can't open PICKY file $ARGV[1]" unless open PICKY, $ARGV[1];
while (<PICKY>) {
  chomp;
  if (/^\S/) { #no space
    if ($line==0) { #sequence
	my @fields = split/\t/;
	unless (exists $ncRNAs{$fields[0]}) {
	    print STDERR "undefined seqence $fields[0]\n";
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
       
	#only keep smaller than match lenth
	my $m_len = $fields[3]-$fields[2] + 1;
	next if $m_len > $match_len;
	for (@{$ncRNAs{$seq}}) {
	    if ($_->[3] eq $name) { #only chromosome are same
		#check with singles
		my $start_pos = $_->[1]-$_->[4]-1;
		my $m_s = $start_pos+$fields[2];
		my $m_e = $start_pos+$fields[3];	       
		my $key = $_->[4]."_".$_->[5]."_".$_->[6];
		next if exists $reports{$key};
		next unless exists $singles{$key};
		my $single = $singles{$key};
		#print STDERR "$key $single \n";
		my @lines = split /,/, $single; 
		my $should_save = 0;
		for(@lines){
		    my @temp = split /_/;
		    next if ($temp[1] - $temp[0] + 1) < 7; #7 is seed length
		    if(($m_s <= $temp[0] && $m_e >= $temp[0]) ||
		       ($m_s <= $temp[1] && $m_e >= $temp[1]) ||
		       ($m_s <= $temp[0] && $m_e >= $temp[1]) ||
		       ($m_s >= $temp[0] && $m_e <= $temp[1])) {
			#print STDERR "$key $m_s $m_e $temp[0] $temp[1] \n";
			$should_save = 1;
			last;
		    }
		}
		$reports{$key} = 1 if $should_save;
	    }	    
	}		
    }
    $line++;
  }else{
      $line = 0;
  }
}
close PICKY;
print STDERR scalar keys %reports, " ncRNAs are left.\n\n";

die "can't open INFO file $ARGV[3]" unless open INFO, $ARGV[3];
my $has_name = $ARGV[5];
while (<INFO>){
   chomp;
   if(/^sRNA/) { #title
      print "$_\n";
   }else {
      my @fields = split /\t/;
      my $rna;
      if($has_name){
	  $rna = $fields[1]."_".$fields[2]."_".$fields[3]; 
      }else{
	  $rna = $fields[0]."_".$fields[1]."_".$fields[2];   
      }
      print "$_\n" if exists $reports{$rna};
   }
}
print STDERR "Complete the job!!\n";
