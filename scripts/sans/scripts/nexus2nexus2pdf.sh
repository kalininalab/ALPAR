# modify the SplitsTree command or make it available via your PATH variable
SPLITSTREE="SplitsTree"
echo "OPEN FILE=$1" > splitstreecommands.tmp
echo "UPDATE" >> splitstreecommands.tmp
echo "UPDATE" >> splitstreecommands.tmp
echo "UPDATE" >> splitstreecommands.tmp
echo "EXPORTGRAPHICS format=PDF TEXTASSHAPES=YES file=$1.pdf REPLACE=yes" >> splitstreecommands.tmp
echo "QUIT" >> splitstreecommands.tmp
eval "$SPLITSTREE -g -c splitstreecommands.tmp"
rm splitstreecommands.tmp
