mkdir -p WG;
cd WG;
for f in `cat ../download_WG.list`; do wget $f; done;
for f in *.fasta.gz; do nf=$(echo $f | sed -E 's/d(...).*\.fasta.gz/\1.fasta.gz/g'); mv $f $nf; done;
cd ..
ls WG/*.fasta.gz > WG_list.txt;
