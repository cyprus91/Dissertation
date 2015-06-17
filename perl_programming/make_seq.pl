#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;

die "usage: fasta_file merge_file need_chop make_rna known_pos(0/#)\n" if @ARGV != 5;

# to reverse and complement DNA fragment
sub reverse_complement {
  my $s = reverse($_[0]);
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
}

#convert DNA to RNA
sub convert_to_RNA {
    my $s = $_[0];
    $s =~ tr/ATGC/AUGC/;
    return $s;
}

my $need_chop = $ARGV[2];
my $convert_rna = $ARGV[3];

my $genome = "";
my $chromosome;
open FASTA, $ARGV[0];
while (<FASTA>) {
    chomp;
    if (/^>(\w+)$/){
	$chromosome = $1;
    }else{
	$genome .= $_;
    }
}
close FASTA;

open KNOWN, $ARGV[1];
my @knowns;
while (<KNOWN>){
    chomp;
    next if /^(Name|sRNA)/;    
    my @fields = split/\t/;
    die "wrong ncRNA coordiate orientation" unless $fields[1]<$fields[2];

    push @knowns, \@fields;

}
close KNOWN;

print STDERR "Total ", $#knowns+1, " known ncRNAs are collected!! \n";
my $known_pos = $ARGV[4];
for (@knowns){
    my ($name, $start, $end, $dir) = @{$_};
    my $known =  $known_pos == 0?0:@{$_}[$known_pos];

    #print STDERR "check $o_name $l_start $start $l_end $end \n";
    my $seq = substr($genome, $start, $end-$start+1);
    my $len = $end-$start + 1;
    die "sequence length is not matched!! $len vs. ", length($seq), "\n" if $len != length($seq);
    $seq =  &reverse_complement($seq) if $dir;
    $seq =  &convert_to_RNA($seq) if $convert_rna;

    if($need_chop) { 
	my $len = 100; #picky allow upto 100
	my $inc_len = 80; #increase by 80, that there is 20 overlap
	my $sub_ord = 1;
	for(my $i=0; $i<length($seq); $i+= $inc_len) {
	    my $sub_seq = substr($seq, $i, $len);
	    my $sub_s = $start+$i + 1;
	    my $sub_e = $sub_s + length($sub_seq)-1;
	    my $sub_title = $name."_".$sub_ord;
	    print  ">$sub_title $sub_s $sub_e $chromosome $start $end $dir $known\n";
	    print  "$sub_seq\n";
	    $sub_ord ++;
	}
	#print "\n";
    }else{
	print  ">$name $start $end $chromosome $start $end $dir $known\n";
	print  "$seq\n";
    }
}

print STDERR "Complete the job!!\n";
