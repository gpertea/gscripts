#!/usr/bin/perl
#basic email sender utility with msmtp
#proper configuration expected in ~/.msmtprc
use strict;
use Getopt::Std;
my $usage=q{Usage:
 email.pl [-s '<subject>'] [-b '<body_text>'] \
   [-f <file_to_include_in_body>] <to@email.com>
};
my %v=();
getopts('s:b:f:', \%v) || die "$usage\n";
my $to=shift(@ARGV) || 
   die("${usage}Error: desination email required!\n");
my $subj=$v{s} || '[no subject]';
my $body=$v{b};
my $fname=$v{f};
if ($fname && $fname ne '-') {
  die("Error: file $fname not found!\n")
       unless -f $fname;
}
my $etxt="To: <$to>\n";
$etxt.="Subject: $subj\n\n";
$etxt.=$body if $body;
open(MSMTP, "| msmtp --tls-certcheck=off -t '$to'");
print MSMTP "$etxt\n";
if ($fname) {
 @ARGV=();
 print MSMTP "\n" if $body;
 if ($fname ne '-') {
     @ARGV=($fname);
     #print MSMTP "-------- content of $fname follows: --------\n";
 }
 print MSMTP $_ while (<>);
}
close(MSMTP);
