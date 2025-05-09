#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

#input: tree with both branch lengths and support values
#)sv:br

#output: see below

sub badInput {
  my $message = "Usage: perl $0
	-i=<tree>
	-o=<output tree>";
  print STDERR $message;
  die "\n";
}

GetOptions(
	"i=s"=>\my $tree,
	"o=s"=>\my $outtree,
);

badInput() if not defined $tree;
badInput() if not defined $outtree;

my $tree_contents = `cat $tree`;
#$tree_contents =~ s/E-\d+//g; #) (no sv)                         
#$tree_contents =~ s/\)(\d+(\.\d+)?):(\d+(\.\d+)?)/):$4/g; #):br  # this is for phylonet
#$tree_contents =~ s/\)(\d+(\.\d+)?):(\d+(\.\d+)?)/)/g; #)
$tree_contents =~ s/:(-)?(\d+(\.\d+)?)//g; #) (no sv)                               # eta thaklei hoi
$tree_contents =~ s/\)(\d+(\.\d+)?)/)/g; #) (no sv)       )1.344 ei type jinishgula strip korar jonno                   

#this is ugly:
#some are: taxa:br, )sv:br
#some taxa are just numbers -> ambiguity
#$tree_contents =~ s/\((.*?):(\d+(\.\d+)?)/\($1/g;
#$tree_contents =~ s/\,(.*?):(\d+(\.\d+)?)/\,$1/g;
#$tree_contents =~ s/\)(\d+(\.\d+)?)?:(\d+(\.\d+)?)/)/g;


open(OUT, ">", $outtree) or die "can't open $outtree: $!";
 print OUT $tree_contents;
close(OUT);

print "output at $outtree\n";
print "done.\n";
