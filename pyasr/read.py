import re 
import dendropy
from phylopandas import DataFrame

def read_codeml_output(filename):
    """Read codeml output and get ancestors as DataFrame. 
    Returns DataFrame and tree as a string.
    """
    # Read paml output.
    with open(filename, 'r') as f:
        data = f.read()
    
    # Rip all trees out of the codeml output. 
    regex = re.compile('\([()\w\:. ,]+;')
    trees = regex.findall(data)

    # First tree in codeml file is the original input tree
    main_tree = dendropy.Tree.get(data=trees[0], schema='newick')

    # Third tree in codeml fule is ancestor tree.
    anc_tree = dendropy.Tree.get(data=trees[2], schema='newick')

    # Map ancestors onto main tree object
    ancestors = anc_tree.internal_nodes()
    for i, node in enumerate(main_tree.internal_nodes()):
        node.label = ancestors[i].label

    # Compile a regular expression to find blocks of data for internal nodes
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+""")
    # Strip the node number from this block of data.
    node_num_regex = re.compile("[0-9]+")

    # Iterate through each block of internal node data
    index, sequences, posteriors = [], [], []
    for node in node_regex.findall(data):
        # Initialize a dictionary for site data
        site_data = {}

        # Compile regex for matching site data
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")

        site_num = 0
        # Iterate through each match for site data.
        seq, post = [], []
        for site in site_regex.findall(node):
            # Iterate through residues
            scores = [float(site[i+2:i+7]) for i in range(0,len(site), 9)]
            j, p = max(enumerate(scores), key=lambda item: item[1])
            seq.append(site[j*9])
            post.append(p)

        # Add site data to ancestor_data
        index.append(node_num_regex.search(node).group(0)) # Get node_number defined by PAML
        sequences.append("".join(seq))
        posteriors.append(sum(post)/len(post))
        
    df = DataFrame({'id':index, "sequence":sequences, "posterior":posteriors}, index=index)
    return main_tree, df

    
