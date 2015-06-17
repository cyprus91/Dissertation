#!/usr/bin/perl
#
# Copyright: (c) 2014-2015 Hyejin Cho <hyejin.cho91@gmail.com>
# License: NO WARRANTY; licensed under LGPL <https://www.gnu.org/copyleft/lesser.html>

use strict;
use warnings;
use Getopt::Long;

my $usage = <<"USAGE";
usage: $0 structure_file extracted_info poly_length has_name(0/1);
USAGE

die $usage unless @ARGV==4;

my $poly_len = $ARGV[2];
my $has_name = $ARGV[3];
my $seed_len = 7;
my $binding_len = 8;
my $interval = 5;
my $single_len = 2; #at least 2 single RNA regions

my $seq;
my $strc;
my $rna;
my %structures;
die "can't open STRUCTURE file $ARGV[0]" unless open STRUCTURE, $ARGV[0];
while (<STRUCTURE>){
   chomp;
   if(/^\s*$/){ #save info
       die " dupliated $rna \n " if exists $structures{$rna};
       $structures{$rna} = [ $seq, $strc ] unless exists $structures{$rna};
   }elsif (/^>/){
     my @fields = split /\s/;
     $rna = $fields[1]."_".$fields[2]."_".$fields[6];
   }elsif(/[AUTGC]+/){
       $seq = $_;
   }else{ #sturcture info
       $strc = $_;
   }
}
close STRUCTURE;


print STDERR  scalar keys %structures, " known ncRNA-structures are collected!! \n";

#find loops and last polyU tail
my %loops;
my %singles;
for my $key (keys %structures) {
    $strc = @{$structures{$key}}[1];

    my @struc_arr = split(//, $strc);
   
    my $last_open_bk = 0;
    my $first_cl_bk = 0;
    my $first_dot = 0;
    my $last_dot = -1;
    #print STDERR "$key ", join('', @{$seq}), "\n";
   for(my $i = 0; $i < @struc_arr; $i++){     
       #print STDERR "$i $struc_arr[$i] : " if $key eq '931556_931680_1';  
       if ( $struc_arr[$i] eq '(' ){ #open bracket	   
	   $last_open_bk = $i;
	   $first_cl_bk = 0 if $first_cl_bk > 0;
       }elsif( $struc_arr[$i] eq ')' ){#close bracket
	   if( $first_cl_bk == 0 ){
	       my $loop = $last_open_bk."_".$i;
	       $loops{$key} = [ ] unless exists $loops{$key};
	       push @{$loops{$key}}, $loop;
	   }		
	   $first_cl_bk = $i;	    
	}
       #find singles
       if ( $struc_arr[$i] eq '.' ){
	   $first_dot = $i if $first_dot == 0;
	   $last_dot = $i;
       }else{
	   if (($last_dot-$first_dot) == 0 || $last_dot == -1) {#only sigle dot
              $first_dot = 0;
              $last_dot = 0;
  	      next; 
	   } 
	    
	   my $single = $first_dot."_".$last_dot;
	   $singles{$key} = [ ] unless exists $singles{$key};
	   push @{$singles{$key}}, $single;
	   $first_dot = 0;
           $last_dot = 0;
       }
   }
   #print STDERR "\n" if $key eq '931556_931680_1';
   #last one
    if (($last_dot-$first_dot) > 0){
	my $single = $first_dot."_".$last_dot;
	$singles{$key} = [ ] unless exists $singles{$key};
	push @{$singles{$key}}, $single;
    }

}

#check structres
my %filtered;
foreach my $key (keys %singles) {
    my @list = @{$singles{$key}};
    my $should_save = 0;
    my $has_binding = 0;
    my $big_length = $seed_len + $binding_len + $interval;
      
    for my $i (reverse 0 .. $#list ) {
	my @temp = split /_/, $list[$i];
	die "Wrong coordinates of single line: $list[$i] at $key \n" unless $temp[0] < $temp[1];     
	my $len = $temp[1]-$temp[0]+1;
	#check polyU tail length
	if ($i == $#list) { #check the last one whether it is polyU or not
	   #print STDERR "$len vs. $poly_len gcvB \n" if $key eq '2940798_2940922_0';
	   last if $len < ($poly_len-1);
	}else{
	    if($len >= $big_length) {
		$should_save = 1;
		last;
	    }
	    last if ($i == 0 && $has_binding == 0); #no binding
	    if ($has_binding) { #already have a binding sites
		if($seed_len >= $len){
		    $should_save = 1;
		    last;

		}
	    }
	    $has_binding = 1 if $len >= $binding_len;	
	}
    }
    $filtered{$key} = 1 if $should_save;
}

print STDERR scalar keys %filtered, " ncRNAs are collected after filtering.\n";
#print out saved ones
die "can't open INFO file $ARGV[1]" unless open INFO, $ARGV[1];
open SINGLE, ">single_lines.txt";
print SINGLE "sRNA\tLines\n";
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
      if (exists $filtered{$rna}) {
	  print "$_\n";
	  #get single lines and write them into the file.
	  my @list = @{$singles{$rna}};
	  print SINGLE "$rna\t";
	  for (@list){
	      print SINGLE "$_,";
	  }
	  print SINGLE "\n";
      }
   }
}
close SINGLE;
print STDERR "Complete the job!!\n";
