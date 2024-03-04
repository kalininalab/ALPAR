alias sans='../../SANS'

# run ALF

alfsim 100_5-params.drw


# extract genomes

cd 100_5/results/100_5/DB
mkdir ../../../fa/
for f in SE*_dna.fa; do ../../../../alf.genomes2dna.py $f > ../../../fa/$f; done
cd ../../../fa/
ls *.fa > list.txt


# get reference tree

../../../../scripts/newick2sans.py ../results/100_5/RealTree.nwk > ../RealTree.sans
sed -r -i "s/(SE...)/\1_dna\.fa/g" ../RealTree.sans


# run SANS-serif

sans -v -i list.txt -f strict -o ../noNs.sans


# add N's

mkdir ../iupac
for f in `cat list.txt`; do ../../IUPACize.py $f 0.001 -n > ../iupac/$f; done
cd ../iupac
ls *.fa > list.txt


# run SANS-serif

sans -v -i list.txt -f strict -o ../skipNs.sans
sans -v -i list.txt -f strict -x 16 -o ../x16.sans


# compare to reference

cd ..
for f in *.sans; do ../../scripts/comp.py $f RealTree.sans fa/list.txt 2 > $f.comp; done
