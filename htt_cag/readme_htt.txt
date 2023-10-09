I'll generate all the HTT variants with CAG counts from 17 to 34 CAGs 
(17 synthetic transcripts), then rerun the whole mapping vs this 
"extended transcriptome".  This will allow us to specifically count the 
unambiguous reads that fully cross the CAG region -- as any partial 
mappings are ambiguous.

I pulled the original HTT-001 transcript from Gencode 25 transcriptome file 
into HTT-001.tfa and renamed that long FASTA ID to just "HTT-001" (this is 
the equivalent of HTT-019 in this naming scheme here.

Reminders:
To get the coordinates of a CAG expansion in a sequence in a FASTA file:
  seqmanip -XL HTT-001.tfa | perl -lne 'print( ($-[0]+1)," $+[0]") if m/((CAG){3,})/' 

To get the CAG count in the expansion:
  seqmanip -XL HTT-001.tfa | perl -nle 'print (length($1)/3) if m/((CAG){13,})/' 

Generating the synthetic transcripts:

for v in {17..34}; do 
  seqmanip -L HTT-001.tfa | perl -ple 'BEGIN { $r="CAG"x'$v' } s/((CAG){13,})/$r/;s/HTT-0\d\d/HTT-0'$v'/'  > HTT-0$v.tfa
done 

This code generates files HTT-017.tfa to HTT-034.tfa with increasing CAG expansions. 
HTT-019 will have the same sequence as HTT-001 so we'll delete that file: 
 rm HTT-019.tfa

we'll rename the entry in gencode25pri.sel.tfa as simply HTT-001 and create the new gencode25_xHTT.tfa file with all the 
extra HTT synthetic transcripts:

 perl -pe '$_=">HTT-001\n" if m/\|HTT-001\|/' gencode25pri.sel.tfa > gencode25_xHTT.tfa
 cat HTT-0[123]?.tfa >> gencode25_xHTT.tfa
 cdbfasta -D gencode25_xHTT.tfa
