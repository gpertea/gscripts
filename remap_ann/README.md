These shell scripts were used for the procedure of remapping the Oryza sativa reference annotation, as found in RAP-DB in May 2018, 
to the new Carolina Gold genome assembly.

* **prep_chrs.sh** : prepares the per-chromosome data (indexed FASTA 
  files with DNA and peptide sequences as found on the reference rice genome)
* **proc_chrs.sh** : processes the per-chromosome data, transferring the annotation to the new assembled chromosome sequences; 
  uses other scripts from the __gscripts__ parent directory and repository:
  * paftools.js
  * gff_liftover.pl 
  * gffcount
  * cds2rna_blat.pl
  
  These scripts also require other programs like **gffread**, **cbdfasta/cdbyank**, **samtools**, **gmap**, **minimap2** 
  and **blat**.
  
