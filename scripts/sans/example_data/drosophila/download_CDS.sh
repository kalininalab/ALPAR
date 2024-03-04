mkdir -p CDS;
cd CDS;
for f in `cat ../download_CDS.list`; do wget $f; done;
for f in *.fasta.gz; do nf=$(echo $f | sed -E 's/d(...).*\.fasta.gz/\1.fasta.gz/g'); mv $f $nf; done;
cd ..
ls CDS/*.fasta.gz > CDS_list.txt;
