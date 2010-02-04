#!/usr/bin/perl

##############################################################################
#   DRAFT_PRODIGAL (Script to process multiple sequences in one file)
#   Copyright (C) 2008-2009 University of Tennessee / UT-Battelle
#
#   Code Author:  Doug Hyatt
#
#   Version 1.20:  September, 2009
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

$| = 1;

$SIG{INT} = 'die_clean';
$SIG{TERM} = 'die_clean';
$SIG{QUIT} = 'die_clean';
$SIG{PIPE} = 'die_clean';
$SIG{KILL} = 'die_clean';

$prodigal = "prodigal";
$tmp = ".";
$trn = "$tmp/dprod.trn.$$";
$trn_inp = "$tmp/dprod.trinp.$$";
$in = "$tmp/dprod.seq.$$";
$out = "$tmp/dprod.out.$$";
$trans = "$tmp/dprod.trans.$$";
$start = "$tmp/dprod.start.$$";

$args = join ' ', @ARGV;
if($args =~ /\-[tT]\s/) { die "error: cannot run $0 with -t argument (temp training file is automatically created)\n"; }
if($args =~ /\-[cC]/) { warn "\n\nDRAFT_PRODIGAL:  (WARNING!) -c option is NOT RECOMMENDED for draft sequences!\n\n"; }
if($args =~ /\-[aA]\s(\S+)/) { 
  $final_trans = $1;
  $args =~ s/\-[aA]\s(\S+)//g;
  $args .= " -a $trans ";
  open WTH, ">$final_trans" or die_clean("couldn't open $final_trans for writing\n");
}
if($args =~ /\-[sS]\s(\S+)/) { 
  $final_start = $1;
  $args =~ s/\-[sS]\s(\S+)//g;
  $args .= " -s $start ";
  open WSH, ">$final_start" or die_clean("couldn't open $final_start for writing\n");
}

$ctr = 0; $state = 0;
while($line = <STDIN>) {
  if($state == 0 && $line =~ /^DEFINITION/) {
    $hdr[$ctr] = $line;
    $hdr[$ctr] =~ s/^DEFINITION\s+//g;
    $hdr[$ctr] = ">".$hdr[$ctr];
  }
  if($line =~ /^ORIGIN/) {
    $seq[$ctr] = "";
    $state = 1;
    $ctr++;
    next;
  }
  if($line =~ /^>/) {
    $hdr[$ctr] = $line;
    $seq[$ctr] = "";
    $state = 1;
    $ctr++;
    next;
  }
  if($state == 1 && $line =~ /^\/\//) {
    $state = 0;
    $hdr[$ctr] = "Prodigal Sequence $ctr";
    next;
  }
  if($state == 1) { 
    $line =~ tr/A-Z/a-z/;
    $line =~ s/[^a-z\n]//g;
    $seq[$ctr-1] .= $line; 
  }
}


print STDERR "\n\nDRAFT_PRODIGAL: Training on $ctr sequences\n\n";
unlink($trn);
open FH, ">$trn_inp" or die_clean("failed to write $trn_inp file\n");
for($i = 0; $i < $ctr; $i++) {
  print FH "$hdr[$i]$seq[$i]";
}
close FH;

$rc = system("$prodigal -t $trn < $trn_inp");
unlink($trn_inp);
if($rc != 0) { die_clean("failed to train\n"); }

for($i = 0; $i < $ctr; $i++) {
  open FH, ">$in" or die_clean("couldn't write tmp input seq\n");
  print FH "$hdr[$i]$seq[$i]";
  close FH;
  print "DEFINITION  $hdr[$i]";
  $cmd = "$prodigal -t $trn $args < $in > $out";
  print STDERR "\n\nDRAFT_PRODIGAL:  Executing $cmd\n\n";
  system($cmd) == 0 or die_clean("failed to run prodigal on sequence ".($i+1)." of $ctr");
  open FH, $out or die_clean("couldn't open tmp output file\n");
  while($line = <FH>) { print $line; }
  close FH;
  print "//\n";

  if(defined($final_trans)) {
    print WTH "DEFINITION  $hdr[$i]";
    if(-e $trans) {
      open FH, "$trans" or die_clean("failed to open $trans file\n");
      while($line = <FH>) { print WTH $line; }
      close FH;
    }
    print WTH "//\n";
  }

  if(defined($final_start)) {
    print WSH "DEFINITION  $hdr[$i]";
    if(-e $start) {
      open FH, "$start" or die_clean("failed to open $start file\n");
      while($line = <FH>) { print WSH $line; }
      close FH;
    }
    print WSH "//\n";
  }

  unlink($out); unlink($trans); unlink($in); unlink($start);
}

close WSH;
close WTH;
unlink($trn);

sub die_clean() {
  $msg = shift @_;
  if($msq eq "") { $msg = "caught a signal\n"; }
  unlink($trn);
  unlink($trn_inp);
  unlink($out); 
  unlink($trans); 
  unlink($in); 
  unlink($start);
  die "$msg\n";
}
