import sys
import random
import decimal
import os

# Classes
class TreeNode:
    def __init__(self, name=None, parent=None):
        # Structure
        self.name = name
        self.children = []
        self.parent = parent
        self.brlen = None
        if type(self.parent) == TreeNode:
            self.parent.add_child(self)

        # Default Attributes
        self.internal_node = False
        self.isroot = False
        self.isleaf = False
        self.rchild = False
        self.attributes = {}

    # Getters and Setters
    def add_child(self, node):
        if node not in self.children:
            self.children.append(node)

    def rm_child(self, node):
        self.children.remove(node)
    
    def get_children(self):
        return self.children

    def set_parent(self, node):
        self.parent = node

    def get_parent(self):
        return self.parent

    def set_attribute(self, key, value):
        self.attributes[key] = value

    def rm_attribute(self, key):
        del self.attributes[key]

    def get_attribute(self, key):
        return self.attributes[key]

    def has_attribute(self, key):
        try:
            self.attributes[key]
            return True
        except KeyError:
            return False

    def __repr__(self):
        return '<TreeNode {}>'.format(self.name)

class PhyloNode(TreeNode):
    '''
    To Do:
    Probably can remove some calls to parent.add_child now that the Node __init__ does that
    
    '''
    def __init__(self, name=None, parent=None, supressWarnings=False):
        super().__init__(name, parent)
        if type(self.parent) == PhyloNode:
            self.parent.add_child(self)
        self.supressWarnings = supressWarnings

    # I/O
    def read_newick(self, newick_tree):
        # Set Up
        current_node = self
        buff = []
        reading_name = 1
        next_internal = False
        try:
            f = open(newick_tree)
            newick_tree = f.read()
        except (FileNotFoundError, OSError):
            if not self.supressWarnings:
                print('Note: Using string {} as raw input. If this is not what was intended, there may be an error in the file name'.format(newick_tree))

        # Read String
        for c in newick_tree:
            c = c.strip()
            if c == ';':
                break
            elif c == '(':
                new_node = PhyloNode(parent=current_node)
                current_node = new_node
                reading_name = 1
            elif c == ')':
                if reading_name:
                    current_node.name = ''.join(buff)
                else:
                    current_node.brlen = float(''.join(buff))
                current_node.rchild = True
                if next_internal:
                    current_node.internal_node = True
                    next_internal = False
                else:
                    current_node.internal_node = False
                if current_node.parent:
                    current_node.parent.add_child(current_node)
                    current_node = current_node.parent
                buff=[]
                reading_name=1
                next_internal = True
            elif c == ':':
                current_node.name = ''.join(buff)
                reading_name = 0
                buff=[]
            elif c == ',':
                if reading_name:
                    current_node.name = ''.join(buff)
                else:
                    current_node.brlen = float(''.join(buff))
                current_node.parent.add_child(current_node)
                if next_internal:
                    current_node.internal_node = True
                    next_internal = False
                else:
                    current_node.internal_node = False
                new_node = PhyloNode(parent=current_node.parent)
                current_node = new_node
                reading_name = 1
                buff=[]
            elif c == ' ':
                pass
            else:
                buff.append(c)
        self.empty = False
        for n in self.iter_nodes():
            if len(n.children) == 0:
                n.isleaf = True

    def write_newick(self):
        buff = []
        if self.internal_node:
            buff.append('(')
        for child in self.children:
            buff.append(child.write_newick())
            if not child.rchild:
                buff.append(',')
        if not self.isroot:
            buff.append(self.name)
            if self.brlen:
                buff.append(':')
                buff.append('{:.20f}'.format(self.brlen))
            if self.rchild:
                buff.append(')')
        else:
            buff.append(';\n')
        return ''.join(buff)

    def get_ascii(self, level, prestring):
        buff = []
        if self.isroot:
            buff.append('\u255F\u2500\u2500\u2500\u2510\n')
        elif self.internal_node:
            if self.name and not self.isroot:
                print_name = '\u2510{}'.format(self.name)
            else:
                print_name = '\u2510'
            if self.rchild:
                buff.append('{}\u2514\u2500\u2500\u2500{}\n'.format(prestring, print_name))
                buff.append('{}    \u2502\n'.format(prestring))
            else:
                buff.append('{}\u251C\u2500\u2500\u2500{}\n'.format(prestring, print_name))
                buff.append('{}\u2502   \u2502\n'.format(prestring))
        else:
            if self.rchild:
                buff.append('{}\u2514\u2500\u2500\u2500{}\n'.format(prestring, self.name))
            else:
                buff.append('{}\u251C\u2500\u2500\u2500{}\n'.format(prestring, self.name))
        for child in self.children:
            if self.rchild:
                buff.append(child.get_ascii(level+1, '{}    '.format(prestring)))
            else:
                buff.append(child.get_ascii(level+1, '{}\u2502   '.format(prestring)))
        if not self.internal_node:
            if self.rchild:
                buff.append('{}   \n'.format(prestring))
            else:
                buff.append('{}\u2502   \n'.format(prestring))
        return ''.join(buff)

    # Functions and Traversals
    def iter_nodes(self, mode='preorder'):
        '''
        Iterates nodes in tree. Note that inorder traversal is only well defined for binary trees.
        '''
        MODES = [   'preorder',
                    'postorder',
                    'inorder'   ]
        assert mode in MODES, 'Invalid mode'
        if mode == 'preorder':
            yield self
            for child in self.children: yield from child.iter_nodes(mode)
        elif mode == 'postorder':
            for child in self.children: yield from child.iter_nodes(mode)
            yield self
        elif mode == 'inorder':
            if len(self.children) == 0:
                yield self
            else:
                yield from self.children[0].iter_nodes(mode)
                yield self
                yield from self.children[1].iter_nodes(mode)

    def iter_leaves(self):
        for node in self.iter_nodes():
            if node.isleaf: yield node
    
    def apply_function(self, f, mode='preorder'):
        '''
        Takes function f(node) and applies it at all nodes of the tree
        '''
        for n in self.iter_nodes(mode=mode):
            f(n)

    def __repr__(self):
        return '<PhyloNode {}>'.format(self.name)

class Phylogeny(PhyloNode):
    def __init__(self, t_file=None, name='root', supressWarnings=False):
        super().__init__(name=name, supressWarnings=supressWarnings)
        self.default_brlen = None
        self.empty = True
        self.rchild = True
        self.internal_node = True
        self.isroot = True
        if t_file:
            self.read_newick(t_file)
        self.isleaf = False

    def read_newick(self, newick_tree):
        super().read_newick(newick_tree)

    def write_newick(self, stream=None):
        string = super().write_newick()
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(string)
        f.flush()
        if f != sys.stdout:
            f.close()

    def get_newick(self):
        return super().write_newick()

    def write_ascii(self, stream=None):
        string = super().get_ascii(0, '')
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(string)
        f.flush()
        if f != sys.stdout:
            f.flush()
    
    def __repr__(self):
        if self.empty:
            return '<Empty Phylogeny Object>'
        else:
            return '<Phylogeny: {}>'.format(self.get_newick().strip())

    def __str__(self):
        if self.empty:
            return self.__repr__()
        else:
            return super().get_ascii(0, '')

class PhyloGroup():
    def __init__(self):
        self.trees = []

    def add_tree(self, tree):
        assert isinstance(tree, Phylogeny), 'Tree is not a Phylogeny object'
        self.trees.append(tree)

    def rm_tree(self, tree):
        self.trees.remove(tree)

    def read_nexus(self, stream):
        try:
            f = open(stream, 'r')
        except:
            print('File {} was unable to be opened'.format(stream))
            sys.exit()
        tree_block = 0
        for line in f:
            if line.lower().startswith('begin trees;'):
                tree_block = 1
            elif tree_block and line.strip().lower().startswith('tree'):
                line_list = line.split()
                name = line_list[1]
                newick = line_list[3]
                new_tree = Phylogeny(name=name, supressWarnings=True)
                new_tree.read_newick(newick)
                self.add_tree(new_tree)
            elif line.lower().startswith('end;'):
                tree_block = 0
            else: pass
        f.close()

    def get_nexus(self):
        buffer = ["Begin trees;\n"]
        for t in self.trees:
            buffer.append("\ttree {} = {}\n".format(t.name, t.get_newick().strip()))
        buffer.append("End;")
        return ''.join(buffer)

    def write_nexus(self, stream=None):
        string = self.get_nexus()
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write('#NEXUS\n')
        f.write(string)
        if f != sys.stdout:
            f.close()

    def __repr__(self):
        return '<PhyloGroup of size {}>'.format(len(self.trees))

class TransmissionNetwork(PhyloNode):
    def __init__(self, t_file=None, name='root'):
        super().__init__(name=name)
        self.empty = True
        if t_file:
            self.read_favites_net(t_file)
        self.isroot = True
        self.internal_node = True
        self.rchild = True
        self.isleaf = False

    def hosts(self):
        return [node.name for node in self.iter_nodes()]

    # I/O Methods
    def read_tnet_outfile(self, tnet_file):
        try:
            f = open(tnet_file, 'r')
        except:
            print('File {} was unable to be opened'.format(tnet_file))
            sys.exit()
        
        transmission_lines = []
        hosts = {'root': self}

        for line in f:
            if '->' in line:
                transmission_lines.append(line.split())

        for line in reversed(transmission_lines):
            dest2 = None
            if len(line) > 3:
                dest2 = line[3]
            source = line[0]
            dest = line[2]
            time = 0

            if source in hosts.keys():
                new_node = PhyloNode(name=dest, parent=hosts[source])
                new_node.set_attribute('infection_time', time)
                hosts[source].add_child(new_node)
                hosts[dest] = new_node
            else:
                new_source = PhyloNode(name=source, parent=self)
                self.add_child(new_source)
                hosts[source] = new_source
                new_dest = PhyloNode(name=dest, parent=new_source)
                new_source.add_child(new_dest)
                hosts[dest] = new_dest

            if dest2:
                new_node = PhyloNode(name=dest2, parent=hosts[source])
                new_node.set_attribute('infection_time', time)
                hosts[source].add_child(new_node)
                hosts[dest2] = new_node

        self.empty = False
        self.apply_function(_finish_tnet)
        f.close()
        
    def read_favites_net(self, favites_net):
        try:
            f = open(favites_net, 'r')
        except:
            print('File {} was unable to be opened'.format(favites_net))
            sys.exit()
        hosts = {'root': self}
        for line in f:
            line = line.strip().split('\t')
            if len(line) != 3:
                line = ''.join(line).split(' ')
            assert len(line) == 3, "Malformed FAVITES network file"
            source = line[0]
            dest = line[1]
            time = line[2]
            if source == dest:
                continue
            elif source in hosts.keys():
                new_node = PhyloNode(name=dest, parent=hosts[source])
                new_node.set_attribute('infection_time', time)
                hosts[source].add_child(new_node)
                hosts[dest] = new_node
            else:
                new_node = PhyloNode(name=dest, parent=self)
                new_node.set_attribute('infection_time', time)
                hosts['root'].add_child(new_node)
                hosts[dest] = new_node
        self.empty = False
        self.apply_function(_finish_tnet)
        f.close()

    def read_phyloscanner_net(self, favites_net):
        try:
            f = open(favites_net, 'r')
        except:
            print('File {} was unable to be opened'.format(favites_net))
            sys.exit()
        hosts = {'root': self}
        i = 0
        for line in f:
            line = line.strip().split('\t')
            assert len(line) == 3, "Malformed FAVITES network file"
            source = line[0]
            dest = line[1]
            time = line[2]
            if source in hosts.keys() and dest in hosts.keys() and hosts[source] in hosts[dest].get_children():
                pass
            elif source in hosts.keys():
                if dest in hosts.keys():
                    hosts[dest].get_parent().rm_child(hosts[dest])
                    hosts[dest].set_parent(hosts[source])
                    hosts[source].add_child(hosts[dest])
                else:
                    new_node = PhyloNode(name=dest, parent=hosts[source])
                    hosts[source].add_child(new_node)
                    hosts[dest] = new_node
                new_node.set_attribute('infection_time', time)
            elif dest in hosts.keys() and source not in hosts.keys():
                new_node = PhyloNode(name=source, parent=hosts[dest])
                new_node.set_attribute('infection_time', time)
                hosts[dest].add_child(new_node)
                hosts[source] = new_node
            else:
                new_source = PhyloNode(name=source, parent=self)
                self.add_child(new_source)
                new_node = PhyloNode(name=dest, parent=new_source)
                new_node.set_attribute('infection_time', time)
                hosts[source] = new_source
                hosts[source].add_child(new_node)
                hosts[dest] = new_node
            i += 1
        # Final pass to stitch any duplicate nodes
        self.empty = False
        self.apply_function(_finish_tnet)
        f.close()

    def read_beastlier_net(self, beastlier_net):
        try:
            f = open(beastlier_net, 'r')
        except:
            print('File {} was unable to be opened'.format(beastlier_net))
            sys.exit()
        hosts = {'root': self}
        for line in f:
            if 'Child' in line:
                continue
            line = line.strip().split(',')
            assert len(line) == 3, "Malformed Beastlier network file"
            dest = line[0]
            source = line[1]
            time = line[2]
            if source in hosts.keys():
                new_node = PhyloNode(name=dest, parent=hosts[source])
                new_node.set_attribute('infection_time', time)
                hosts[source].add_child(new_node)
                hosts[dest] = new_node
            else:
                new_node = PhyloNode(name=dest, parent=self)
                new_node.set_attribute('infection_time', time)
                hosts['root'].add_child(new_node)
                hosts[dest] = new_node
        self.empty = False
        self.apply_function(_finish_tnet)
        f.close()

    def write_favites_net(self, stream=sys.stdout):
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(self.get_favites_net())
        f.close()

    def get_favites_net(self):
        buff = []
        nodes = list(self.iter_nodes(mode='preorder'))
        nodes.remove(self)
        if nodes[0].has_attribute('infection_time'):
            def getTime(node):
                return float(node.get_attribute('infection_time'))
            nodes.sort(key=getTime)
            for node in nodes:
                if node.get_parent().isroot:
                    parent_name = 'None'
                else:
                    parent_name = node.get_parent().name
                buff.append('{}\t{}\t{}\n'.format(parent_name, node.name, getTime(node)))
        else:
            for node in nodes:
                if node.get_parent().isroot:
                        parent_name = 'None'
                else:
                    parent_name = node.get_parent().name
                buff.append('{}\t{}\t{}\n'.format(parent_name, node.name, float(0)))
        return ''.join(buff)

    def write_newick(self, stream=sys.stdout):
        string = super().write_newick()
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(string)
        f.flush()
        if f != sys.stdout:
            f.close()

    def write_ascii(self, stream=sys.stdout):
        string = super().get_ascii(0, '')
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(string)
        f.flush()
        if f!= sys.stdout:
            f.close()

    def __repr__(self):
        if self.empty:
            return '<Empty TransmissionNetwork Object>'
        else:
            return '<TransmissionNetwork: {}>'.format(super().write_newick().strip())

    def __str__(self):
        if self.empty:
            return self.__repr__()
        else:
            return super().get_ascii(0, '')

class Sequence():
    def __init__(self, seq, taxa=None, num=None, datatype='dna'):
        self.seq = seq
        t_list = list(taxa)
        for i in range(len(t_list)-1):
            if t_list[i] == ' ':
                t_list[i] == '_'
        taxa = ''.join(t_list)
        self.taxa = taxa
        self.len = len(seq)
        if not num:
            self.num = self.taxa
        self.datatype = datatype

    def __repr__(self):
        return 'Sequence of type {} with name {}'.format(self.datatype, self.num)

    def __str__(self):
        return '{}\t{}'.format(self.taxa, self.seq)

class SeqGroup():
    def __init__(self, datatype='dna'):
        self.seqs = []
        self.datatype = datatype
        self.length = None

    def add_seq(self, sequence):
        if not self.length:
            self.length = sequence.len
        else:
            assert sequence.len == self.length, 'Sequence length ({})does not match group length ({})'.format(sequence.len, self.length)
        assert sequence.datatype == self.datatype, 'Sequence datatype ({})does not match group datatype ({})'.format(sequence.datatype, self.datatype)
        self.seqs.append(sequence)

    def rm_seq(self, sequence):
        self.seqs.remove(sequence)

    def read_fasta(self, stream=None):
        try:
            f = open(stream, 'r')
        except:
            print('File {} was unable to be opened'.format(stream))
            sys.exit()
        
        firstseq = 1
        buffer = []
        seqname = ''
        for line in f:
            if line.startswith('>'):
                seqname = line.strip().split('>')[1]
                if not firstseq:
                    self.seqs.append(Sequence(''.join(buffer), datatype=self.datatype, taxa=seqname))
                    buffer = []
                else:
                    firstseq = 0
            else:
                buffer.append(line.strip())
        self.add_seq(Sequence(''.join(buffer), datatype=self.datatype, taxa=seqname))

    def get_fasta(self):
        buffer = []
        for s in self.seqs:
            buffer.extend(['>', s.taxa, '\n', s.seq, '\n'])
        return ''.join(buffer)

    def write_fasta(self, stream=None):
        string = self.get_fasta()
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write(string)

    def read_nexus(self, stream):
        try:
            f = open(stream, 'r')
        except:
            print('File {} was unable to be opened'.format(stream))
            sys.exit()
        data_block = 0
        matrix = 0
        for line in f:
            if line.lower().startswith('begin data;'):
                data_block = 1
            elif line.lower().strip().startswith('matrix'):
                matrix = 1
            elif line.lower().strip().startswith(';'):
                data_block = 0
                matrix = 0
            elif data_block and matrix:
                line_list = line.split()
                name = line_list[0]
                seq = line_list[1]
                self.add_seq(Sequence(seq, taxa=name))
            else: pass
        f.close()

    def get_nexus(self):
        buffer = ['Begin data;\n']
        buffer.append('Dimensions ntax={} nchar={};\n'.format(len(self.seqs), self.length))
        buffer.append('Format datatype={} missing=? gap=-;\n'.format(self.datatype))
        buffer.append('Matrix\n')
        for s in self.seqs:
            buffer.extend(['\t', s.taxa, ' ', s.seq, '\n'])
        buffer.append(';\nEnd;')
        return ''.join(buffer)

    def write_nexus(self, stream=None):
        string = self.get_nexus()
        try:
            f = open(stream, 'w+')
        except:
            f = sys.stdout
        f.write('#NEXUS\n')
        f.write(string)
        f.close()

    def __repr__(self):
        return 'SeqGroup with {} seqs of type {}'.format(len(self.seqs), self.datatype)

# Transmission Network Inference Functions
def label_internal_nodes(tree):
    if isinstance(tree, Phylogeny):
        t = tree
    else:
        t = Phylogeny(tree)
    t.apply_function(_prep_node)
    _fitch_recurse_up(t)
    _fitch_recurse_down(t)
    return t

def construct_transmission_network(labeled_tree):
    hosts = {}
    hosts['root'] = TransmissionNetwork()
    _tnet_down(labeled_tree, hosts)
    tnet = hosts['root']
    tnet.empty = False
    tnet.rchild = True
    tnet.internal_node = True
    tnet.isroot = True
    tnet.apply_function(_finish_tnet)
    return tnet

def fitch_score(tree):
    if isinstance(tree, Phylogeny):
        t = tree
    else:
        t = Phylogeny(tree)
    t.apply_function(_prep_node)
    return _fitch_recurse_up(t)

def compare_networks(refnet, testnet):
    if not isinstance(refnet, TransmissionNetwork):
        refnet = TransmissionNetwork(refnet)
    if not isinstance(testnet, TransmissionNetwork):
        testnet = TransmissionNetwork(testnet)
    
    ref_netfile = refnet.get_favites_net()
    ref_edgelist = [set([j[0], j[1]]) for j in [i.split('\t') for i in ref_netfile.strip().split('\n')]]
    try:
        test_netfile = testnet.get_favites_net()
        test_edgelist = [set([j[0], j[1]]) for j in [i.split('\t') for i in test_netfile.strip().split('\n')]]
    except IndexError:
        return(0,len(ref_edgelist),0)
    false_negative = 0
    false_positive = 0
    common_edges = 0
    for i in ref_edgelist:
        if i not in test_edgelist:
            false_negative += 1
        else:
            common_edges += 1
    for j in test_edgelist:
        if j not in ref_edgelist:
            false_positive += 1
    return(false_positive, false_negative, common_edges)
    
# Helper Functions for Transmission Network Inference
def _prep_node(node):
    if node.isleaf:
        node.set_attribute('host', node.name.split('_')[0])
        node.set_attribute('strain', node.name.split('_')[1])
        host_set = set()
        host_set.add(node.get_attribute('host'))
        node.set_attribute('host_set', host_set)
    else:
        node.set_attribute('host_set', set())

def _fitch_recurse_up(node):
    if node.isleaf: return 0
    assert len(node.children) == 2, "Tree must be binary"
    lchild = node.children[0]
    rchild = node.children[1]
    s = _fitch_recurse_up(lchild) + _fitch_recurse_up(rchild)
    host_set = lchild.get_attribute('host_set') & rchild.get_attribute('host_set')
    if not len(host_set):
        host_set = lchild.get_attribute('host_set') | rchild.get_attribute('host_set')
        s+=1
    node.set_attribute('host_set', host_set)
    return s

def _fitch_recurse_down(node):
    if node.isleaf: return
    assert len(node.get_children()) == 2, "Tree must be binary"
    if node.isroot:
        node.set_attribute('host', _choose_host(node))
        node.name = node.get_attribute('host')
    elif node.get_parent().get_attribute('host') in node.get_attribute('host_set'):
        node.set_attribute('host', node.get_parent().get_attribute('host'))
    else:
        node.set_attribute('host', _choose_host(node))
    _fitch_recurse_down(node.get_children()[0])
    _fitch_recurse_down(node.get_children()[1])
    if not node.name:
        node.name = node.get_attribute('host')

def _choose_host(node):
    return random.choice(list(node.get_attribute('host_set')))

def _tnet_down(node, hosts):
    node_host = node.get_attribute('host')
    if node.isroot:
        hosts[node_host] = PhyloNode(name=node_host, parent=hosts['root'])
        hosts['root'].add_child(hosts[node_host])
        node.internal_node = True
    elif node_host in hosts.keys():
        pass
    else:
        parent_host = hosts[node.get_parent().get_attribute('host')]
        new_host = PhyloNode(name=node_host, parent=parent_host)
        parent_host.add_child(new_host)
        new_host.internal_node = False
        hosts[node_host] = new_host
    for child in node.get_children():
        _tnet_down(child, hosts)
      
def _finish_tnet(node):
    if len(node.get_children()) == 0:
        node.isleaf = True
    else:
        node.internal_node = True
        for child in node.get_children():
            child.rchild = False
        node.get_children()[-1].rchild = True
