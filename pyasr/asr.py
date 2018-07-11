import os
import shutil
import subprocess
import pkg_resources

from Bio.Phylo.PAML import codeml

from .read import read_codeml_output

# Run paml
def reconstruct(
    df,
    id_col='uid',
    sequence_col='sequence',
    working_dir='',
    save_ancestors=False,
    altall_cutoff=0.2,
    infer_gaps=True,
    aaRatefile='lg',
    **kwargs
    ):

    df = df.copy()

    # Construct default arguments
    default_options = dict(verbose=9, CodonFreq=None, cleandata=0,
        fix_blength=2, NSsites=None, fix_omega=None, clock=None,
        ncatG=8, runmode=0, fix_kappa=None, fix_alpha=1, Small_Diff=1.0e-6,
        method=0, Malpha=None, aaDist=None, RateAncestor=2, icode=None,
        alpha=None, seqtype=2, omega=None, getSE=None, noisy=3, Mgene=None,
        kappa=None, model=3, ndata=None)

    # Update default arguments in place.
    default_options.update(**kwargs)

    # ---------------- Prepare model ----------------
    # copy model from package to project directory.
    path_to_model = pkg_resources.resource_filename(
        'pyasr', os.path.join('dat', '{}.dat'.format(aaRatefile)))

    model_file = '{}.dat'.format(aaRatefile)
    model_path = os.path.join(working_dir, model_file)
    shutil.copyfile(path_to_model, model_path)

    # ----------------------

    curr_path = os.getcwd()
    proj_path = os.path.join(curr_path, working_dir)
    ali_path = os.path.join(working_dir, 'ali-to-reconstruct.phy')
    tree_path = os.path.join(working_dir, 'tree-to-reconstruct.phy')
    out_path = os.path.join(working_dir, 'results.txt')
    ctl_path = os.path.join(working_dir, 'codeml_options.ctl')
    rst_path = os.path.join(working_dir, 'rst')

    df.phylo.to_fasta(
        filename=ali_path,
        id_col=id_col,
        sequence_col=sequence_col,
    )

    df.phylo.to_newick(
        filename=tree_path,
        taxon_col=id_col,
        node_col=id_col,
        suppress_internal_node_labels=True,
    )

    df.phylo.to_newick(
        taxon_col=id_col,
        node_col=id_col,
        suppress_internal_node_labels=True,
    )

    # Build and write out control file.
    cml = codeml.Codeml(alignment=ali_path,
        tree=tree_path,
        out_file=out_path,
        working_dir=working_dir)
    cml.set_options(aaRatefile=model_file, **default_options)
    cml.ctl_file = ctl_path
    cml.write_ctl_file()

    # ----------------------

    os.chdir(proj_path)
    output = subprocess.run(['codeml', 'codeml_options.ctl'])
    os.chdir(curr_path)

    # ----------------------

    return read_codeml_output(rst_path, df)
