# python rf.py true_37.tre out.tre

#!/usr/bin/env python3
import sys

try:
    import dendropy
    from dendropy.calculate import treecompare
except ImportError:
    sys.exit("ERROR: this script requires DendroPy. Install with:\n    pip install dendropy")

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <tree1.tre> <tree2.tre>\n")
        sys.exit(1)

    t1_path, t2_path = sys.argv[1], sys.argv[2]

    # load both trees into the same taxon namespace
    taxa = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(path=t1_path, schema="newick",
                              taxon_namespace=taxa,
                              preserve_underscores=True)
    tree2 = dendropy.Tree.get(path=t2_path, schema="newick",
                              taxon_namespace=taxa,
                              preserve_underscores=True)

    # compute unrooted RF distance using the new API
    rf = treecompare.symmetric_difference(tree1, tree2)

    # maximum possible RF for an unrooted binary tree with n leaves is 2*(n-3)
    n_leaves = len(tree1.leaf_nodes())
    max_rf = 2 * (n_leaves - 3) if n_leaves >= 3 else 0

    # normalize
    rf_rate = rf / max_rf if max_rf > 0 else 0.0

    print(f"RF distance: {rf}")
    print(f"Max possible RF: {max_rf}")
    print(f"RF rate (normalized): {rf_rate:.4f}")

if __name__ == "__main__":
    main()
