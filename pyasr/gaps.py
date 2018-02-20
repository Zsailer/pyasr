import dendropy

def infer_gaps_in_tree(df_seq, tree, id_col='id', sequence_col='sequence'):
    """Adds a character matrix to DendroPy tree and infers gaps using
    Fitch's algorithm.

    Infer gaps in sequences at ancestral nodes.
    """
    taxa = tree.taxon_namespace

    # Get alignment as fasta
    alignment = df_seq.phylo.to_fasta(id_col=id_col, id_only=True,
                                      sequence_col=sequence_col)

    # Build a Sequence data matrix from Dendropy
    data = dendropy.ProteinCharacterMatrix.get(
        data=alignment,
        schema="fasta",
        taxon_namespace=taxa)

    # Construct a map object between sequence data and tree data.
    taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=False)

    # Fitch algorithm to determine placement of gaps
    dendropy.model.parsimony.fitch_down_pass(tree.postorder_node_iter(),
            taxon_state_sets_map=taxon_state_sets_map)
    dendropy.model.parsimony.fitch_up_pass(tree.preorder_node_iter())
    return tree
