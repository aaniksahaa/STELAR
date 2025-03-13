import matplotlib.pyplot as plt
from collections import Counter
import re

def split_top_level(s):
    """
    Split a string s (without outer parentheses) at commas that are at top level.
    """
    parts = []
    level = 0
    current = []
    for char in s:
        if char == '(':
            level += 1
        elif char == ')':
            level -= 1
        if char == ',' and level == 0:
            parts.append(''.join(current).strip())
            current = []
        else:
            current.append(char)
    if current:
        parts.append(''.join(current).strip())
    return parts

def parse_newick(s):
    """
    Recursively parse a Newick-style string representing a binary tree.
    Returns either an integer (or string) for a leaf, or a tuple ("internal", left, right).
    If more than two elements are found at one level, they are combined left-associatively.
    """
    s = s.strip()
    # Remove trailing semicolon if present
    if s.endswith(";"):
        s = s[:-1]
    # If the string is wrapped in parentheses, remove them
    if s.startswith("(") and s.endswith(")"):
        # Check for matching outer parentheses by counting
        # (This avoids stripping if the parentheses are not the true outer ones)
        level = 0
        proper = True
        for i, char in enumerate(s):
            if char == '(':
                level += 1
            elif char == ')':
                level -= 1
            if level == 0 and i < len(s)-1:
                proper = False
                break
        if proper:
            s = s[1:-1].strip()
            parts = split_top_level(s)
            # If only one part is inside, parse it directly.
            if len(parts) == 1:
                return parse_newick(parts[0])
            elif len(parts) == 2:
                return ("internal", parse_newick(parts[0]), parse_newick(parts[1]))
            else:
                # For more than two parts, combine them left-associatively.
                tree = parse_newick(parts[0])
                for part in parts[1:]:
                    tree = ("internal", tree, parse_newick(part))
                return tree
    # If no surrounding parentheses, then it is a leaf.
    # Try converting to integer, otherwise leave as string.
    try:
        return int(s)
    except ValueError:
        return s

def get_leaves(tree):
    """
    Recursively collect all leaves (as a set) in the tree.
    """
    if isinstance(tree, tuple) and tree[0] == "internal":
        return get_leaves(tree[1]) | get_leaves(tree[2])
    else:
        return {tree}

def traverse_and_collect(tree, biparts, subtree_sizes):
    """
    Recursively traverse the tree.
    For each internal node, compute the bipartition (as a pair of frozensets of leaves)
    and record the size (number of leaves) of each child subtree.
    Append these to biparts (list) and subtree_sizes (list).
    """
    if isinstance(tree, tuple) and tree[0] == "internal":
        left_tree, right_tree = tree[1], tree[2]
        left_leaves = get_leaves(left_tree)
        right_leaves = get_leaves(right_tree)
        # Canonicalize bipartition (order by sorted list of labels)
        if sorted(left_leaves) <= sorted(right_leaves):
            bip = (frozenset(left_leaves), frozenset(right_leaves))
        else:
            bip = (frozenset(right_leaves), frozenset(left_leaves))
        biparts.append(bip)
        subtree_sizes.append(len(left_leaves))
        subtree_sizes.append(len(right_leaves))
        traverse_and_collect(left_tree, biparts, subtree_sizes)
        traverse_and_collect(right_tree, biparts, subtree_sizes)

def process_trees(filename):
    """
    Process a file of trees (one tree per line) and return:
    - A Counter of bipartition frequencies.
    - A list of subtree sizes.
    """
    bipartitions = []
    subtree_sizes = []
    with open(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            # Parse the tree
            tree = parse_newick(line)
            # Collect bipartitions and subtree sizes from the tree.
            traverse_and_collect(tree, bipartitions, subtree_sizes)
    # Count frequencies of each canonical bipartition.
    bip_counter = Counter(bipartitions)
    return bip_counter, subtree_sizes

def plot_bipartition_frequencies(bip_counter):
    """
    Plot a bar plot showing the distribution of bipartition frequencies.
    The x-axis represents the number of times a bipartition occurs,
    and the y-axis represents the count of bipartitions with that frequency.
    """
    freq_counts = Counter(bip_counter.values())
    x = sorted(freq_counts.keys())
    y = [freq_counts[f] for f in x]

    plt.figure(figsize=(8, 6))
    plt.bar(x, y, width=0.6, edgecolor='black')
    plt.xlabel("Bipartition occurrence count")
    plt.ylabel("Number of unique bipartitions")
    plt.title("Distribution of Bipartition Frequencies")
    plt.xticks(x)
    plt.tight_layout()
    plt.show()

def plot_subtree_size_distribution(subtree_sizes):
    """
    Plot a histogram of the distribution of subtree sizes.
    """
    plt.figure(figsize=(8, 6))
    plt.hist(subtree_sizes, bins=range(min(subtree_sizes), max(subtree_sizes)+2), edgecolor='black', align='left')
    plt.xlabel("Subtree size (number of leaves)")
    plt.ylabel("Frequency")
    plt.title("Distribution of Subtree Sizes")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Replace 'trees.txt' with your actual file name.
    filename = "./Data-set/test-dataset/test-100-gt/gene_trees.txt"
    filename = "./Data-set/37-taxon/noscale.100g.500b/R20/R20.genetrees"
    filename = "./in.tre"
    filename = "./in2.tre"
    
    bip_counter, subtree_sizes = process_trees(filename)
    
    print("Number of unique bipartitions found:", len(bip_counter))
    # Uncomment the next line to see the raw frequency counts:
    # for bip, count in bip_counter.items():
    #     print(bip, ":", count)
    
    # Plot bipartition frequency distribution
    plot_bipartition_frequencies(bip_counter)
    
    # Plot subtree size distribution
    plot_subtree_size_distribution(subtree_sizes)


    # aa
    # bb
