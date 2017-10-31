import re 
import dendropy
from phylopandas import DataFrame

def read_codeml_output(filename):
    """Parse the 'rst' file returned by codeml. Returns Maximum Likelihood (ML) and ML+1  
    residues for all ancestors, and their posterior probability.

    Returns
    -------
    ancestors : dict of DataFrames
        The key is the node label in the tree. The value is a DataFrame with 
        4 columns: ml_residue, ml_posterior, alt_residue, and alt_posterior 
    tree : dendropy.Tree
        Tree object with ancestors and tips labelled. Keeps branch length as well.
    """
    # Read paml output.
    with open(filename, 'r') as f:
        data = f.read()
    
    # Rip all trees out of the codeml output. 
    regex = re.compile('\([()\w\:. ,]+;')
    trees = regex.findall(data)

    # First tree in codeml file is the original input tree
    tip_tree = dendropy.Tree.get(data=trees[0], schema='newick')

    # Third tree in codeml fule is ancestor tree.
    anc_tree = dendropy.Tree.get(data=trees[2], schema='newick')

    # Main tree to return
    tree = tip_tree

    # Map ancestors onto main tree object
    ancestors = anc_tree.internal_nodes()
    for i, node in enumerate(tree.internal_nodes()):
        node.label = ancestors[i].label

    # Compile a regular expression to find blocks of data for internal nodes
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+\n\n""")
    
    # Strip the node number from this block of data.
    node_num_regex = re.compile("[0-9]+")

    # Get dataframes for all ancestors.
    ancestors = {}
    for node in node_regex.findall(data):
        data = {'ml_residue': [], 'ml_posterior':[], 'alt_residue':[], 'alt_posterior':[]}

        # Compile regex for matching site data
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")
        
        # Iterate through each match for site data.
        for site in site_regex.findall(node):
            # Iterate through residues
            scores = [float(site[i+2:i+7]) for i in range(0,len(site), 9)]            
            residues = [site[i] for i in range(0, len(site), 9)]
            
            # Get the indices of sorted scores
            sorted_score_index = [i[0] for i in sorted(enumerate(scores), key=lambda x:x[1], reverse=True)]
            ml_index = sorted_score_index[0]
            alt_index = sorted_score_index[1]
            
            # Store ML sequence and Alternative sites.
            data['ml_residue'].append(residues[ml_index])
            data['ml_posterior'].append(scores[ml_index])
            data['alt_residue'].append(residues[alt_index])
            data['alt_posterior'].append(scores[alt_index])

        # Convert to dataframe. 
        df = DataFrame(data)
                
        # Node number
        node_label = node_num_regex.search(node).group(0)
        ancestors[node_label] = df
        
    return ancestors, tree
