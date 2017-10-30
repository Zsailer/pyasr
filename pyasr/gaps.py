import dendropy

def add_gaps_to_ancestors(tree, df_seq, df_anc, id_col='id'):
    """Estimate gaps for ancestral nodes.
    Handling Gaps
    -------------
    Uses Fitch Parsimony to determine sites in ancestral sequences that
    are likely gaps. The posterior probability of these sites are ignored
    when calculating the average posterior probability of the ancestor.
    """
    taxa = tree.taxon_namespace

    # Get alignment as fasta
    alignment = df_seq.to_fasta(id_col=id_col, id_only=True)

    # Build a Sequence data matrix from Dendropy
    data = dendropy.ProteinCharacterMatrix.get(
        data=alignment,
        schema="fasta",
        taxon_namespace=taxa)
    
    # Get the alphabet of Dendropy's ProteinCharacterMatrix
    alphabet = data.state_alphabets[0].symbols
    gap = alphabet.index('-')

    # Construct a map object between sequence data and tree data.
    taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=False)

    # Fitch algorithm to determine placement of gaps
    dendropy.model.parsimony.fitch_down_pass(tree.postorder_node_iter(),
            taxon_state_sets_map=taxon_state_sets_map)
    dendropy.model.parsimony.fitch_up_pass(tree.preorder_node_iter())

    # Interate through ancestors and insert gaps.
    for node in tree.internal_nodes():
        # get ancestor sequence
        row = df_anc.loc[df_anc['id']==node.label]

        anc_seq = list(row['sequence'][0])
        
        for i, site in enumerate(node.state_sets):
            # Get a list of possible residues at each site.
            sites = list(node.state_sets[i])
            
            # Check if site should be a gap.
            if len(sites) == 1 and sites[0] == gap:
                anc_seq[i] = "-"
                
        # Set value in ancestor dataframe
        df_anc.loc[df_anc['id']==node.label, 'sequence'] = ''.join(anc_seq)
        #df_anc.set_value(node.label, 'sequence', ''.join(anc_seq))
        
    return tree, df_seq, df_anc
