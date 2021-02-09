#!/usr/bin/perl
use strict;
my $usage = q/Usage:
 gtf2cdb.pl stringtie_gtf > gtf.cfa
 
 Converts a stringtie GTF into a cdbfasta indexable file 
 with exon coverage info preserved etc.
/;
my $t; #current transcript ID 
my @td; # (chr, strand, start, end, cov, fpkm, tpm) for current transcript
        #   0      1       2     3   4     5     6
my @ex; #exon data for current exon: list of [start, end, coverage]
die("$usage\n") if @ARGV==0 || $ARGV[0]=~m/^[\-]+h/;
while(<>) {
  next if m/^#/;
  next unless m/^chr/;
  my ($c, $z, $f, $s, $e, $y, $strand, $w, $a)=split(/\t/);
  my ($tid)=($a=~m/transcript_id "([^"]+)/);
  next unless $tid;
  my ($cov)=($a=~m/cov "([^"]+)/);
  if ($tid ne $t) { #new transcript, flush the previous one
    #this MUST be a transcript line 
    tflush() if $t; #writes record
    @ex=();
    $t=$tid;
    my ($fm)=($a=~m/FPKM "([^"]+)/);
    my ($tm)=($a=~m/TPM "([^"]+)/);
    @td=($c,$strand, $s, $e, sprintf('%.2f',$cov), 
         sprintf('%.2f',$fm),sprintf('%.2f',$tm));
  } else {
    #must be an exon
    die("Error: exon expected at line: $_") unless $f eq 'exon';
    push(@ex, [$s, $e, sprintf('%.2f',$cov)]);
  }
}

tflush() if $t;

#-----------------
sub tflush {
  #uses globals $t, @td, @ex
  print ">$t ".join(' ',@td)."\n";
  print join(',', (map { $$_[0].'-'.$$_[1].':'.$$_[2] } @ex))."\n";
}
