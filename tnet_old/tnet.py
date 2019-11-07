#!/usr/bin/python3
import sys, os
import subprocess as sp
from phylo import Phylogeny, TransmissionNetwork

def rename(node):
    if node.isleaf:
        node.name = 'outbreak{}'.format(node.name)
    if node.brlen is None or node.brlen < 0.01:
        node.brlen = 0.01

def rename_back(node):
    if not node.name == 'root':
        node.name = node.name.split('outbreak')[1]

def netgen(file_path, tnet_outfile, verbose=False, tnet_path='TNet/lib/tnet.linux'):
    # Prepare input tree file
    tree = Phylogeny(file_path)
    # if verbose:
        # print('Renaming source tree...')
    tree.apply_function(rename)
    tree.write_newick('{}.temp'.format(file_path))

    # Run tnet
    # if verbose:
        # print('Running tnet...')
    cmd = '{} {}.temp {}'.format(tnet_path, file_path, tnet_outfile)
    child = sp.Popen(cmd.split(), stdout=sp.PIPE)
    out, err = child.communicate()
    if err:
        if verbose:
            print('Failure')
            print(err.decode('utf-8'))
        sys.exit(1)
    if verbose:
        print(out.decode('utf-8'))

    # Parse tnet outfile into FAVITES transmission network
    tnetwork = TransmissionNetwork()
    tnetwork.read_tnet_outfile(tnet_outfile)
    tnetwork.apply_function(rename_back)
    return tnetwork

def cli():
    # Parse Args
    if not len(sys.argv) == 3:
        raise IndexError("Usage: python3 tnet_gen.py [input phylogeny file] [desired output file]")
    TREE_FILE = os.path.abspath(sys.argv[1])
    OUT_FILE = os.path.abspath(sys.argv[2])
    TNET_OUTFILE = '{}.tnet.log'.format(TREE_FILE)

    net = netgen(TREE_FILE, TNET_OUTFILE,verbose=True)
    net.write_favites_net(OUT_FILE)

    # Clean up
    os.remove('stage2.txt')

def api(tree_file, out_file, log_file,  tnet_path):
    
    netgen(tree_file, log_file, tnet_path=tnet_path).write_favites_net(out_file)

# main
if __name__ == '__main__':
    cli()
