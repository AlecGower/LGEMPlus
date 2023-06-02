python3 experiments/create_theory_files.py "yeastGEM"
python3 experiments/create_theory_files.py "iMM904"
python3 experiments/create_theory_files.py "iND750"

for model in $(ls experiments/theories/); do
	python3 experiments/single_gene_deletions.py \
        experiments/theories/$model/$(ls experiments/theories/$model | tail -n 1)
done
