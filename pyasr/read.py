import re
import dendropy
from phylopandas import DataFrame


def read_codeml_output(
    filename,
    df,
    altall_cutoff=0.2,
    ):
    """Read codeml file.
    """
    # Read paml output.
    with open(filename, 'r') as f:
        data = f.read()

    # Rip all trees out of the codeml output.
    regex = re.compile('\([()\w\:. ,]+;')
    trees = regex.findall(data)
    anc_tree = trees[2]

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

    # Map nodes onto dataframe.
    df['reconstruct_label'] = None
    for node in tree.postorder_node_iter():

        # Ignore parent node
        if node.parent_node is None:
            pass

        elif node.is_leaf():
            node_label = node.taxon.label
            parent_label = node.parent_node.label
            # Set node label.
            df.loc[df.uid == node_label, 'reconstruct_label'] = node_label

            # Set parent label.
            parent_id = df.loc[df.uid == node_label, 'parent'].values[0]
            df.loc[df.id == parent_id, 'reconstruct_label'] = node.parent_node.label

        elif node.is_internal():
            label = node.label
            parent_id = df.loc[df.reconstruct_label == label, 'parent'].values[0]
            df.loc[df.id == parent_id, 'reconstruct_label'] = node.parent_node.label


    # Compile a regular expression to find blocks of data for internal nodes
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+\n""")

    # Strip the node number from this block of data.
    node_num_regex = re.compile("[0-9]+")

    # Get dataframes for all ancestors.
    df['ml_sequence'] = None
    df['ml_posterior'] = None
    df['alt_sequence'] = None
    df['alt_posterior'] = None

    for node in node_regex.findall(data):
        # Get node label
        node_label = node_num_regex.search(node).group(0)

        # Compile regex for matching site data
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")

        # Iterate through each match for site data.
        ml_sequence, ml_posterior, alt_sequence, alt_posterior = [], [], [], []

        for site in site_regex.findall(node):

            # Iterate through residues
            scores = [float(site[i+2:i+7]) for i in range(0,len(site), 9)]
            residues = [site[i] for i in range(0, len(site), 9)]

            # Get the indices of sorted scores
            sorted_score_index = [i[0] for i in sorted(
                enumerate(scores),
                key=lambda x:x[1],
                reverse=True)]

            ml_idx = sorted_score_index[0]
            alt_idx = sorted_score_index[1]

            # Should we keep alterative site.
            ml_sequence.append(residues[ml_idx])
            ml_posterior.append(scores[ml_idx])

            if scores[alt_idx] < altall_cutoff:
                alt_idx = ml_idx

            alt_sequence.append(residues[alt_idx])
            alt_posterior.append(scores[alt_idx])

        keys = [
            "ml_sequence",
            "ml_posterior",
            "alt_sequence",
            "alt_posterior"
        ]

        vals = [
            "".join(ml_sequence),
            sum(ml_posterior) / len(ml_posterior),
            "".join(alt_sequence),
            sum(alt_posterior) / len(alt_posterior),
        ]

        df.loc[df.reconstruct_label == node_label, keys] = vals

    return df
