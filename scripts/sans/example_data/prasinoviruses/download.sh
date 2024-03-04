SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p "$SCRIPT_DIR/fa";

rm -f "$SCRIPT_DIR/list.txt";
cat "$SCRIPT_DIR/IDs.txt" | while read l;
do
	NC=$(echo $l | cut -f1 -d" "); 
	TAX=$(echo $l | cut -f2 -d " ");
	curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$NC&rettype=fasta" > "$SCRIPT_DIR/fa/$TAX.fasta"
        echo "fa/$TAX.fasta" >> "$SCRIPT_DIR/list.txt"
done;
