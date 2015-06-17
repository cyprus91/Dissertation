#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;

die "usage: genome DNA/RNA poly_length bait_length\n" if @ARGV != 4;

my $genome = '';
my $polyU = '';
my $polyA = '';
my $bait_length = $ARGV[3];
my $chromosome;

my $seed =  $ARGV[1] eq 'DNA'?'T':'U';
my $poly_len = $ARGV[2];
for( my $i = 0; $i < $poly_len; $i++){
    $polyU .= $seed;
    $polyA .= 'A';
}
#print STDERR $polyU,"\n";
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


#find every position of polyU
my $count = () = $genome =~ /$polyU/g;
my $offset = 0;
my $position = index($genome, $polyU, $offset);

my @polyUs;
my $prev_pos = 0;
while($position != -1) {
    $prev_pos = $position;
    $offset = $position + 1;
    $position = index($genome, $polyU, $offset);
    push @polyUs, $prev_pos if $prev_pos < $position-$poly_len;
}
print STDERR "$count polyU in the given genome\n";

print STDERR "Working for polyA \n";

#find every position of polyA
$count = () = $genome =~ /$polyA/g;
$offset = 0;
$position = index($genome, $polyA, $offset);

print STDERR "$count polyA in the given genome\n";

my @polyAs; #for complements
$prev_pos = 0;
while($position != -1) {
    $prev_pos = $position;
    $offset = $position + 1;
    $position = index($genome, $polyA, $offset);
    push @polyAs, $prev_pos if $prev_pos < $position-$poly_len;;
}

#make baits for polyU
my $bait_cnt = 1;
for(@polyUs){
    my ($start, $end);
    $end = $_+$poly_len;
    $start = $end-$bait_length+1;
    $start = $start < 0? 0:$start;
    my $seq = substr $genome, $start, $bait_length;
    print ">BAIT$bait_cnt $start $end $chromosome $start $end forward 1\n";
    print "$seq\n";
    $bait_cnt++;
}
for(@polyAs) {
    my ($start, $end);
    $start = $_-1;
    $end = $start+$bait_length-1;
    $end = $end > length($genome)?length($genome):$end;
    my $seq = substr $genome, $start, $bait_length;
    print ">BAIT$bait_cnt $start $end $chromosome $start $end reverse 1\n";
    print "$seq\n";
    $bait_cnt++;
}
print STDERR $bait_cnt, " baits are created!! \n";
print STDERR "Complete the job!!\n";
