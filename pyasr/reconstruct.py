import os
import shutil
import subprocess
import pkg_resources

from Bio.Phylo.PAML import codeml
from phylopandas import DataFrame

from .read import read_codeml_output
from .gaps import infer_gaps_in_tree

def reconstruct(df_seq, tree, id_col='id', sequence_col='sequence', working_dir='', save_ancestors=False, altall_cutoff=0.2, infer_gaps=True, aaRatefile='lg', **kwargs):
    """Use PAML to contruct ancestral sequences by Maximum Likelihood.
    
    Parameters
    ----------
    df_seq : phylopandas.DataFrame
        dataframe containing information (including sequences) about tips of tree.
    tree : dendropy.Tree
        Tree object that includes leaf nodes for each sequence in the dataframe.
    id_col: str (default = 'id')
        column in df_seq that maps the sequences in df_seq to the tree tips.
    sequence_col : str (default = 'sequence')
        column of aligned sequences to use for reconstruction.
    working_dir : str (default = '')
        directory to spew PAML output.
    save_ancestors : bool (default=False)
        save ancestors dataframes to file.
    altall_cutoff : float (default=0.2)
        probability cutoff for selecting alternative sites to flip in AltAll sequence. Flips any
        sites whose posterior is greater than or equal to altall_cutoff.
    infer_gaps : bool (default='true')
        If tree, uses Fitch's algorithm to infer gaps in the ancestral sequences.
    aaRatefile : str (default='lg')
        Evolutionary model to use for reconstruction. Read more about these models
        in the PAML documentation.
        
    Returns
    -------
    df_seq : phylopandas.DataFrame
        dataframe containing information (including sequences) about tips of the tree.
    df_ancs : phylopandas.DataFrame
        dataframe containing information (including sequences) about ancestors of the tree.
    tree_ancs : dendropy.Tree
        Tree updated with ancestral nodes labelled according to df_ancs
    """
    # Construct default arguments
    default_options = dict(verbose=9, CodonFreq=None, cleandata=0,
        fix_blength=2, NSsites=None, fix_omega=None, clock=None,
        ncatG=8, runmode=0, fix_kappa=None, fix_alpha=1, Small_Diff=1.0e-6,
        method=0, Malpha=None, aaDist=None, RateAncestor=2, icode=None,
        alpha=None, seqtype=2, omega=None, getSE=None, noisy=3, Mgene=None,
        kappa=None, model=3, ndata=None)
        
    # Update default arguments in place.
    default_options.update(**kwargs)

    # Write out alignment for PAML
    n, m = len(df_seq), len(df_seq['sequence'][0])
    alignment_file = 'alignment.phy'
    alignment_path = os.path.join(working_dir, alignment_file)
    alignment_str = df_seq.to_fasta(sequence_col=sequence_col, id_col=id_col, id_only=True)
    
    # FORMATTING ALIGNMENT ANNOYINGLY. This is a hack for now.
    lines = alignment_str.strip().split('\n')
    header = "{}\n".format(lines[0])
    seq = ''
    fasta_str = '' 
    for line in lines[1:]:
        if line[0] == '>':
            fasta_str += "{}\n".format(header)
            header = "{}\n".format(line)
        else:
            header += line
    fasta_str += "{}\n".format(header)
    alignment_str = "{} {}\n".format(n,m) + fasta_str
    
    # Write to file
    with open(alignment_path, 'w') as f:
        f.write(alignment_str)
    
    # Write tree to file
    tree_file = 'tree-to-reconstruct.newick'
    tree_path = os.path.join(working_dir, tree_file)
    tree.write(path=tree_path, schema='newick', 
        suppress_internal_taxon_labels=True,
        suppress_internal_node_labels=True)
    
    output_file = 'results.txt'
    output_path = os.path.join(working_dir, output_file)
    
    # copy model from package to project directory. 
    path_to_model = pkg_resources.resource_filename('pyasr', os.path.join('dat', '{}.dat'.format(aaRatefile)))
    model_file = '{}.dat'.format(aaRatefile)
    model_path = os.path.join(working_dir, model_file)
    shutil.copyfile(path_to_model, model_path)

    # Build control file.
    cml = codeml.Codeml(alignment=alignment_path, 
        tree=tree_path,
        out_file=output_path,
        working_dir=working_dir)
    cml.set_options(aaRatefile=model_file, **default_options)

    # Write out control file.
    cml.ctl_file = os.path.join(working_dir, 'codeml_options.ctl')
    cml.write_ctl_file()     

    # Run PAML: 1. change directory, 2. run codeml, 3. parse results and 4. switch back to current dir.
    current_dir = os.getcwd()
    project_dir = os.path.join(current_dir, working_dir)
    os.chdir(project_dir)
    output = subprocess.run(['codeml', 'codeml_options.ctl'])
    os.chdir(current_dir)
    
    # Parse output.
    rst_file = os.path.join(working_dir, 'rst')
    ancestors, anc_tree = read_codeml_output(rst_file)

    # Map ancestors (and support) onto main tree object.
    anc_nodes = anc_tree.internal_nodes()
    for i, node in enumerate(tree.internal_nodes()):
        node.support = node.label
        node.label = anc_nodes[i].label

    # Save ancestors to 'ancestors' directory in working_dir
    if save_ancestors:
        os.makedirs(os.path.join(working_dir, 'ancestors'))
        for anc, df in ancestors.items():
            path = os.path.join(working_dir, 'ancestors', '{}'.format(anc))
            df.to_csv(path)

    # Summarize data in a single dataframe
    data = {'id':[], 'ml_sequence':[], 'ml_posterior':[], 'alt_sequence':[], 'alt_posterior':[], 'support':[]}
    if infer_gaps:
        # Infer gaps in tree object
        tree = infer_gaps_in_tree(df_seq, tree, id_col=id_col, sequence_col=sequence_col)
        gap = 21 # Gap index in state_sets
        
        for anc, df in ancestors.items():
            # Get ML_sequence
            ml_seq = list(df['ml_residue'])
            alt_seq = list(df['alt_residue'])
            ml_p, alt_p = [], [] 
            
            # Get node in tree object.
            node = tree.find_node_with_label(anc)
            
            # Build sequences and probability statistics
            for i, site in enumerate(ml_seq):
                # Get a list of possible residues at each site.
                sites = list(node.state_sets[i])
                
                # Check if site should be a gap.
                if len(sites) == 1 and sites[0] == gap:
                    ml_seq[i] = "-"
                    alt_seq[i] = "-"
                else:
                    # If no gap, append ML site to ML sequence
                    ml_p.append(df['ml_posterior'][i]) 
                    
                    # Should we flip this site for the altall seq? 
                    if df['alt_posterior'][i] < altall_cutoff:
                        # Don't flip, keep the ML value.
                        alt_seq[i] = df['ml_residue'][i]
                        alt_p.append(df['ml_posterior'][i])
                    else:
                        # Flip and get posterior value.
                        alt_p.append(df['alt_posterior'][i])
            
            # Start everything in data.
            ml_sequence = ''.join(ml_seq)
            alt_sequence = ''.join(alt_seq)
            ml_posterior = sum(ml_p) / len(ml_p)
            alt_posterior = sum(alt_p) / len(alt_p)
            data['id'].append(anc)
            data['ml_sequence'].append(ml_sequence)
            data['ml_posterior'].append(ml_posterior)
            data['alt_sequence'].append(alt_sequence)
            data['alt_posterior'].append(alt_posterior)
            data['support'].append(node.support)
    else:

        for anc, df in ancestors.items():
            # Get ML_sequence
            ml_seq = list(df['ml_residue'])
            alt_seq = list(df['alt_residue'])
            ml_p, alt_p = [], [] 
            
            # Get node in tree object.
            node = tree.find_node_with_label(anc)
            
            # Build sequences and probability statistics
            for i, site in enumerate(ml_seq):
                # Get a list of possible residues at each site.
                sites = list(node.state_sets[i])
                
                # If no gap, append ML site to ML sequence
                ml_p.append(df['ml_posterior'][i]) 
                
                # Should we flip this site for the altall seq? 
                if df['alt_posterior'][i] < altall_cutoff:
                    # Don't flip, keep the ML value.
                    alt_seq[i] = df['ml_residue'][i]
                    alt_p.append(df['ml_posterior'][i])
                else:
                    # Flip and get posterior value.
                    alt_p.append(df['alt_posterior'][i])
            
            # Start everything in data.
            ml_sequence = ''.join(ml_seq)
            alt_sequence = ''.join(alt_seq)
            ml_posterior = sum(ml_p) / len(ml_p)
            alt_posterior = sum(alt_p) / len(alt_p)
            data['id'].append(anc)
            data['ml_sequence'].append(ml_sequence)
            data['ml_posterior'].append(ml_posterior)
            data['alt_sequence'].append(alt_sequence)
            data['alt_posterior'].append(alt_posterior)
            data['support'].append(node.support)
            
    df_ancs = DataFrame(data)
    return df_seq, df_ancs, tree
