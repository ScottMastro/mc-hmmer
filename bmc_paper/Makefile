all: e.tsv a.tsv

ss.txt:
	wget http://www.rcsb.org/pdb/files/$@
	sed -i 's/ /N/g' $@
	./process_ss.py > tmp.txt
	mv tmp.txt ss.txt

e.tsv: e_paper.tsv ss.txt
	./process_e.py > $@

a.tsv: a_paper.tsv
	./process_a.py > $@
