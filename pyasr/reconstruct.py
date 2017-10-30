import os
import shutil
import subprocess
import pkg_resources

from Bio.Phylo.PAML import codeml

from .read import read_codeml_output
from .gaps import add_gaps_to_ancestors

def reconstruct(df, tree, id_col='id', sequence_col='sequence', working_dir='', infer_gaps=True, aaRatefile='lg', **kwargs):
    """Use PAML to build a phylogenetic tree."""
    # Construct default arguments
    default_options = dict(verbose=9, CodonFreq=None, cleandata=0,
        fix_blength=2, NSsites=None, fix_omega=None, clock=None,
        ncatG=8, runmode=0, fix_kappa=None, fix_alpha=1, Small_Diff=1.0e-6,
        method=0, Malpha=None, aaDist=None, RateAncestor=2, icode=None,
        alpha=None, seqtype=2, omega=None, getSE=None, noisy=3, Mgene=None,
        kappa=None, model=3, ndata=None)
        
    # Update default arguments in place.
    default_options.update(**kwargs)

    # Write file to disk
    ### NEED TO FIX SEQUENCE COL!!
    n, m = len(df), len(df['sequence'][0])
    alignment_file = 'alignment.phy'
    alignment_path = os.path.join(working_dir, alignment_file)
    alignment_str = df.to_fasta(sequence_col=sequence_col, id_col=id_col, id_only=True)
    alignment_str = "{} {}\n".format(n,m) + alignment_str
    
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
        # verbose = verbose,
        # CodonFreq = CodonFreq,
        # cleandata = cleandata,
        # fix_blength = fix_blength,
        # NSsites = NSsites,
        # fix_omega = fix_omega,
        # clock = clock,
        # ncatG = ncatG,
        # runmode = runmode,
        # fix_kappa = fix_kappa,
        # fix_alpha = fix_alpha,
        # Small_Diff = Small_Diff,
        # method = method,
        # Malpha = Malpha,
        # aaDist = aaDist,
        # RateAncestor = RateAncestor,
        # aaRatefile = model_file,
        # icode = icode,
        # alpha = alpha,
        # seqtype = seqtype,
        # omega = omega,
        # getSE = getSE,
        # noisy = noisy,
        # Mgene = Mgene,
        # kappa = kappa,
        # model = model,
        # ndata = ndata)
        
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
    tree_ancs, df_ancs = read_codeml_output(rst_file)

    # Infer gaps
    if infer_gaps == True:
        tree_ancs, df_seq, df_ancs = add_gaps_to_ancestors(tree_ancs, df, df_ancs, id_col=id_col)

    return tree_ancs, df_seq, df_ancs
