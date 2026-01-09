from Bio import Phylo
import sys

# Redirect output to log
sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

tree = Phylo.read(snakemake.input.tree, "newick")
ref_terminals = [t for t in tree.get_terminals() if t.name == snakemake.params.reference]

if ref_terminals:
    tree.prune(ref_terminals[0])
    print(f"Removed reference from tree: {snakemake.params.reference}")
else:
    print(f"Warning: Reference {snakemake.params.reference} not found in tree")

Phylo.write(tree, snakemake.output.tree, "newick")
print("Done!")