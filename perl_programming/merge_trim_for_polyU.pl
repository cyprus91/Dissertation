#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;
use Getopt::Long;

my $usage = <<"USAGE";
usage: $0 gbk predicted_file max_length ext_length;
USAGE

die $usage unless @ARGV==4;

my $max_length = $ARGV[2];
my $ext_length = $ARGV[3];

die "can't open gbk file $ARGV[0]" unless open GBK, $ARGV[0];
my @genes;
my $chromosome;
my $complement;
my ($g_name, $left, $right, $type) = ('', 0, 0, '');
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
      next if $g_name eq 'sgrT';
    } elsif (/\s+\bCDS\b\s+(complement\()?(\d+)\.\.(\d+)/) {
      die "mismatch complement $complement" if $complement != defined($1);
      warn "mismatch left $left != $2" if $left != $2;
      warn "mismatch right $right != $3" if $right != $3;
      $type = "mRNA";
      next if $g_name eq 'sgrT';
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
print STDERR scalar @genes, " mRNAs are collected.\n";

die "can't open predicted file $ARGV[0]" unless open PREDICT, $ARGV[1];
my $line = 0;
my @predicted;
my %saved;
while (<PREDICT>) {
  if ($line) {
    last if /^\s*$/;
    chomp;
    my @fields = split/\t/;
    my @temp = split ("-", $fields[1]);
    my @target = split(", ", $fields[3]);
    my $sRNA = $temp[0]."_".$temp[1]."_".$fields[2];
    die "wrong ncRNA coordiate orientation" unless $temp[0]<$temp[1];
    #push @predicted, [ $temp[0], $temp[1], $fields[4],$target[0], $target[1], $fields[2] ];   
    push @predicted, [ $temp[0], $temp[1], $fields[2] ] unless exists $saved{$sRNA};
    $saved{$sRNA} = 1;    
  }
  $line++;
}
close PREDICT;

print STDERR scalar @predicted, " predicted ncRANs are read in \n";

my @sorted_predict = sort {$a->[2] <=> $b->[2] || $a->[0] <=> $b->[0]} @predicted;
#my @sorted_predict = sort {$a->[0] <=> $b->[0]} @predicted;

$line = 0;
my @trimed;
my ($p_left, $p_right, $p_dir) = (0, 0, 0);
#my $should_next = 0;
for (@sorted_predict){
    my ($left, $right, $dir) = @{$_};
    #extend    
    if ($dir) { #complement so extend right only
	 $right = $right + $ext_length;
    }else{ #extend left only
	 $left = $left - $ext_length;
	 $left = $left<0?0:$left;
    }
   
    if ($line == 0) {
	($p_left, $p_right, $p_dir) = ($left, $right, $dir);
    }else{
	#check max length
	my $len = $right-$left;
	if($len >= $max_length){
	    
	    push @trimed, [ $p_left, $p_right, $p_dir ];
	    push @trimed, [ $left, $right, $dir ];
	    ($p_left, $p_right, $p_dir) = ($left, $right, $dir);
	    $line = 0; #start as new
	    next;
	}
	if ($p_dir != $dir) { #start complement ones
	    push @trimed, [ $p_left, $p_right, $p_dir ];
	    ($p_left, $p_right, $p_dir) = ($left, $right, $dir);	  
	}else{
	    #next if $p_left == $left && $p_right == $right;
	    next if $p_left <= $left && $p_right >= $right; #old one_enclude new one
	    if ($p_left <= $left && $p_right < $right){ #overlapped
		#check length
		$len = $right-$p_left;
		if($len >= $max_length){	    
		    push @trimed, [ $p_left, $p_right, $p_dir ];		   
		    ($p_left, $p_right, $p_dir) = ($left, $right, $dir);
		   
		    next;
		}
		$p_right = $right; #combined		
	    }elsif ($left > $p_right){ #far from previous one, keep previous.
		push @trimed, [ $p_left, $p_right, $p_dir ];
		($p_left, $p_right, $p_dir) = ($left, $right, $dir);#start new one	    
	    }
	}
    }

    $line ++;
}

print STDERR scalar @trimed, " ncRNAs are collected after merge. \n";
open T, ">after_trim.txt";
print T "Start End Dir\n";
for(@trimed){
    my ($s, $e, $d) = @{$_};
    print T "$s\t$e\t$d\n";
}
close T;
#filtered trimed-predicted ncRNA if it is overapped with mRNA
#my @filtered;

my $hit = 0;

print "sRNA\tsRNA_Start\tsRNA_End\tsRNA_Dir\tLength\tMark\n";
#sort both of them
@genes = sort {$a->[1] <=> $b->[1]} @genes;
@trimed = sort {$a->[0] <=> $b->[0]} @trimed;
my $cnt = 1;
for (@trimed){
    my $should_remove = 0;
    my ($start, $stop, $dir) = @{$_};
    for (@genes){
	($complement, $left, $right) = @{$_}[0..2];
	my $nm = @{$_}[3];
	#print STDERR "There is sgrT in mRNA - $nm \n" if $nm eq 'sgrT';
	next if $complement != $dir;
	next if $right < $start;
	last if $left > $stop;
	if (($left <= $start && $right >= $start) ||
	    ($left <= $stop && $right >= $stop) ||
	    ($left <= $start && $right >= $stop) ||
	    ($start <= $left && $stop >= $right)) {
	    #push @filtered, [ $start, $stop, $dir ];
	    $should_remove = 1;
	    last;
	}
    }
    if ($should_remove == 0) { #can save but check length
	$hit++;
	my $len = $stop-$start+1;
	my $mark= $len > $max_length?"Big":"Decent";
	print "Predict$hit\t$start\t$stop\t$dir\t$len\t$mark\n";
    }
}

print STDERR "$hit ncRNAs are collected after triming with gbk \n";
