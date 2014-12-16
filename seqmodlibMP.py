#!/usr/bin/python

import sys
import time
import math
import random
import itertools
import cPickle
import operator
import multiprocessing
import datetime
import randsiggen
import gc

def load_proteome(filename, silent=True):
    """
    Unpickles proteome and returns it as a dictionary.

    Args:
        filename: String indicating proteome pickle file location.
        silent: Boolean indicating whether progress should be printed to
            standard output.

    Returns:
        Dictionary containing proteome.
        {'PROTEIN NAME 1': 'AMINO ACID SEQUENCE 1',
         'PROTEIN NAME 2': 'AMINO ACID SEQUENCE 2',
         ...}
    """
    load_proteome_time = time.clock()
    if not silent:
        print("")
        print("unpickling proteome from " + str(filename))
    #the expected proteome dictionary format in the pickled file is
    #{'PROTEIN NAME': 'AMINO ACID SEQUENCE'}
    try:
        proteome = cPickle.load(open(filename))
    except:
        print("Proteome unpickling error.")
        raise
    if not silent:
        #using time.clock() - starting_time is prone to underflow after ~30 min
        print("proteome loaded in " + str(time.clock() - load_proteome_time) + 
              " sec")
        print("...")
    return proteome

def homogenize(peptides, substitute_acid, target_acids):
    """
    Replace each instance of every acid in target_acids with a substitute acid.
    Its original intended use is to simplify generating multichroic signals
    where one color labels more than one amino acid.

    Args:
        peptides: Dictionary of peptides to be homogenized. It is in the same
            format as the dictionary returned by load_proteome.
            {'PROTEIN 1': 'AMINO ACID SEQUENCE 1',
             'PROTEIN 2': 'AMINO ACID SEQUENCE 2',
             ...}
        substitute_acid: One-letter string indicating which amino acid will
            substitute for the target acids.
        target_acids: Tuple or list of one-letter strings indicating which
            amino acids will be substituted with substitute_acid. Example:
            Aspartic acid 'D' will be labeled with the same color as glutamic
            acid 'E'. First call this function with arguments
            substitute_acid='E', target_acids=['D'].

    Returns:
        Dictionary in identical to the peptides dictionary given, except with
        all target_acids replaced with the substitute_acid. This dictionary is
        generated fresh, and does not side-affect the original argument.
    """
    return_peptides = {}
    for protein in peptides:
        sequence = peptides[protein]
        for acid in target_acids:
            homogenized_sequence = sequence.replace(acid, substitute_acid)
        return_peptides.setdefault(protein, homogenized_sequence)
    return return_peptides

def cleave(peptides, cleave_acid, silent=True):
    """
    Cleave each peptide passed to the function at the peptide bond following a
    given amino acid.

    Args:
        peptides: Dictionary of peptides to be cleaved. It is in the same
            format as the dictionary returned by load_proteome.
            {'PROTEIN 1': 'AMINO ACID SEQUENCE 1',
             'PROTEIN 2': 'AMINO ACID SEQUENCE 2',
             ...}
        cleave_acid: One-letter string indicating after which acid to cleave
            the bond, e.g. passing a 'K' cleaves after every lyseine.
        silent: Boolean indicating whether progress should be printed to
            standard output.

    Returns:
        Returns a dictionary mapping the given proteins to the peptides
        resulting from cleaving.
            {'PROTEIN 1': ('CLEAVED SEQUENCE 1', 'CLEAVED SEQUENCE 2', ...),
             'PROTEIN 2': ('CLEAVED SEQUENCE 1', 'CLEAVED SEQUENCE 2', ...),
             ...}
        Empty sequences, i.e. '', are not added. If a protein has no non-empty
        subsequences, its entry is not added. If cleaving a peptide results in
        more than one identical subsequence, they are all added anyways. The
        returned dictionary is generated de novo; the original dictionary
        argument is not affected.
    """
    #progress tracking
    cleave_time = time.clock()
    cleave_progress = 0
    #the new dictionary that will be returned
    return_peptides = {}
    for protein in peptides:
        #disregard empty sequences
        if not peptides[protein]:
            continue
        subsequences = peptides[protein].split(cleave_acid)
        #split(cleave_acid) omits the cleave_acid itself from every resulting
        #substring and adds an empty substring to the end if cleave_acid is at
        #the end of a given sequence, e.g.
        #'ABCABCABCCCC'.split('C') returns ['AB', 'AB', 'AB', '', '', '', '']
        for index, string in enumerate(subsequences[:-1]):
            #add omitted cleave_acid to all but the last gap; if the last acid
            #in the sequence is cleave_acid itself, it leaves an empty
            #subsequence last which needs to be removed anyways; if the last
            #acid is not cleave_acid, then cleave_acid does not need to be
            #added to it
            subsequences[index] += cleave_acid
        if subsequences[-1] == '':
            #if last subsequence resulting from split is empty, this means that
            #the last item in the original sequence was cleave_acid, which is 
            #readded by the loop directly above to the penultimate member
            subsequences.pop()
        #eliminate all empty subsequences
        subsequences = [subsequence for subsequence in subsequences
                        if subsequence]
        #if subsequences is empty list, do not add
        if subsequences:
            return_peptides.setdefault(protein, tuple(subsequences))
        if not silent:
            cleave_progress += 1
            sys.stdout.write("%d of %d peptides cleaved\r" %
                             (cleave_progress, len(peptides)))
    if not silent:
        #using time.clock() - starting_time is prone to underflow after ~30 min
        print("")
        print("proteome cleaved in " + str(time.clock() - cleave_time) + "sec")
        print("...")
    return return_peptides

def attach(peptides, attach_acid, silent=True):
    """
    Attach peptides to substrate via an acid.

    Args:
        peptides: Dictionary of peptides to attach. It is in the same format as
            the return value of cleave()
            {'PROTEIN 1': ('CLEAVED SEQUENCE 1', 'CLEAVED SEQUENCE 2'),
             'PROTEIN 2': ('CLEAVED SEQUENCE 1', 'CLEAVED SEQUENCE 2'), ...}
        attach_acid: One-letter string indicating which acid is used to attach
            the peptides, e.g. 'C' will attach all peptides using a cysteine.
            If the string is 'cterm', then all peptides are attached using
            their carboxyl terminus.
        silent: Boolean indicating whether progress should be printed to
            standard output.

    Returns:
        Returns a dictionary mapping proteins to all attached peptides
        associated with them. Those peptides that do not contain the attaching
        acid are omitted. Proteins with no peptides that can attach are
        omitted. Each attached peptide is partitioned in two: the peptide head,
        which represents the portion of the peptide preceding the first
        attaching acid and hence accessible to Edman chemistry, and peptide
        tail, which represents all amino acids after the first attachment,
        where they are inaccessible to Edman chemistry. Returned dictionary
        format is as follows:
            {'PROTEIN 1': (('PEPTIDE HEAD 1', 'PEPTIDE TAIL 1'),
                           ('PEPTIDE HEAD 2', 'PEPTIDE TAIL 2'),
                            ...),
             'PROTEIN 2': (('PEPTIDE HEAD 1', 'PEPTIDE TAIL 1'),
                           ('PEPTIDE HEAD 2', 'PEPTIDE TAIL 2'),
                            ...),
             ...}
        If a protein has more than one identical attachment pair, they are all
        included. This function does not affect the original peptide argument.
    """
    #progress tracking
    attach_time = time.clock()
    attach_progress = 0
    if not silent:
        print("attaching peptides using " + str(attach_acid))
    #the dictionary that will be returned
    return_peptides = {}
    #special case for cterm
    if attach_acid == 'cterm':
        for protein in peptides:
            for index, sequence in enumerate(peptides[protein]):
                return_peptides.setdefault(protein, []).append((sequence,''))
            return_peptides[protein] = tuple(return_peptides[protein])
        return return_peptides
    for protein in peptides:
        for index, sequence in enumerate(peptides[protein]):
            #if sequence does not contain attaching acid, it is omitted
            if attach_acid in sequence:
                attach_point = sequence.find(attach_acid)
                return_peptides.setdefault(protein, [])
                return_peptides[protein].append((sequence[:attach_point],
                                                 sequence[attach_point:]))
        if protein in return_peptides:
            #if protein has successful attachments, tuplefy result
            return_peptides[protein] = tuple(return_peptides[protein])
        if not silent:
            attach_progress += 1
            sys.stdout.write("peptides for %d of %d proteins attached\r" %
                             (attach_progress, len(peptides)))
    if not silent:
        #using time.clock() - starting_time is prone to underflow after ~30 min
        print("")
        print("peptides attached in " + str(time.clock() - attach_time) +
              " sec")
        print("...")
    return return_peptides

def homogenize_attached(peptides, substitute_acid, target_acids):
    """
    Same as homogenize, but operates on peptides in the format yielded by
    attached.
    """
    return_peptides = {}
    for protein, sequences in peptides.iteritems():
        for head, tail in sequences:
            for acid in target_acids:
                head = head.replace(acid, substitute_acid)
                tail = tail.replace(acid, substitute_acid)
            return_peptides.setdefault(protein, []).append((head, tail))
    for protein, sequences in return_peptides.iteritems():
        return_peptides[protein] = tuple(sequences)
    return return_peptides

class SignalTrie:
    """
    An implementation of a trie (prefix tree) structure to store large numbers
    of signals, and for each of these signals to track the number of times they
    were generated by particular souce proteins. The trie's root node is null,
    representing an empty signal. All descendant nodes identify themselves by a
    letter representing an amino acid and an integer for the its position in
    the sequence. A node represents the signal composed by concatenating nodes'
    gaps and amino acids transversed to reach it, itself included. Each node
    contains a dictionary whose keys are (position, amino acid) pairs pointing
    to further sequence members; each node is likewise pointed to by its
    ancestor. To count the number of times a protein generated the signal
    represented by the node, each node contains a dictionary mapping the source
    protein of its signal to an integer representing the number its signals
    from that protein.
    """
    def __init__(self, (pg, aa)):
        """
        Initializes this node to be made of amino acid aa with position pg.
        This node's descendant dictionary is empty, and it has no signal counts
        at initialization.
        """
        #(position, amino acid)
        self.signal_block = (pg, aa)
        #{next_signal_block: SignalTrie}
        self.descendants = {}
        #{protein: count}
        self.signal_count = {}
    def add_descendant(self, subsignal, source_protein):
        """
        Increments the count of the signal represented in the tree by this node
        followed by subsignal by one, with source_protein being the source of
        the signal. If the signal is not yet present in the tree, recursively
        adds subsignal as a descendant of this node.

        Args:
            subsignal: The remainder following this node of the signal to be
                incremented. It is a tuple as for all signals
                ((position, aa), (position, aa), ...).
            source_protein: String, name of protein that generated the signal.

        Returns:
            Self.
        """
        if len(subsignal) == 0:
            return
        elif self.signal_block == (None, None):
            self.descendants.setdefault(subsignal[0], SignalTrie(subsignal[0]))
            self.descendants[subsignal[0]].add_descendant(subsignal,
                                                          source_protein)
        elif len(subsignal) == 1:
            self.signal_count.setdefault(source_protein, 0)
            self.signal_count[source_protein] += 1
        else:
            self.descendants.setdefault(subsignal[1], SignalTrie(subsignal[1]))
            self.descendants[subsignal[1]].add_descendant(subsignal[1:],
                                                          source_protein)
        return self
    def set_descendant(self, subsignal, count):
        """
        Sets the signal_count in the leaf representing subsignal to a copy of
        count.

        Args:
            subsignal: The remainder following this node of the signal whose
                signal count is to be changed to count's copy. It is a tuple as
                for all signals
                ((position, aa), (position, aa), ...).
            count: Set subsignal's signal count to a copy of this dictionary.

        Returns:
            Self.
        """
        if len(subsignal) == 0:
            return
        elif self.signal_block == (None, None):
            self.descendants.setdefault(subsignal[0], SignalTrie(subsignal[0]))
            self.descendants[subsignal[0]].set_descendant(subsignal, count)
        elif len(subsignal) == 1:
            self.signal_count = count.copy()
        else:
            self.descendants.setdefault(subsignal[1], SignalTrie(subsignal[1]))
            self.descendants[subsignal[1]].set_descendant(subsignal[1:], count)
        return self
    def get_descendant(self, subsignal):
        """
        Return pointer to node represented by subsignal.

        Args:
            subsignal: The remainder following this node of the subsignal whose
                pointer to return. It is a tuple as for all signals
                ((position, aa), (position, aa), ...).

        Returns:
            Pointer to node represented by subsignal if it exists, None
                otherwise.
        """
        if len(subsignal) == 0:
            return
        elif self.signal_block == (None, None):
            if subsignal[0] in self.descendants:
                return self.descendants[subsignal[0]].get_descendant(subsignal)
            else:
                return None
        elif len(subsignal) == 1:
            return self
        else:
            if subsignal[1] in self.descendants:
                return\
                   self.descendants[subsignal[1]].get_descendant(subsignal[1:])
            else:
                return None
    def node_iterator(self):
        """
        An iterator over ALL descendant nodes and this node itself, whether
        with empty signal counts or not. See leaf_iterator to iterate over only
        nonempty descendant nodes. Iteration yields a tuple for each node:
        (signal represented by the node, signal counter, pointer to node)

        #??#This cannot be used while adding or removing nodes from the tree as
        #??#dictionaries cannot change size while iterating.
        """
        for d_trie in self.descendants.itervalues():
            for node in d_trie.node_iterator():
                if self.signal_block == (None, None):
                    yield node
                else:
                    #the decomposition is necessary to prepend
                    #self.signal_block
                    yield ((self.signal_block,) + node[0], node[1], node[2])
        #the following line is the critical difference between node_iterator
        #and leaf_iterator, as leaf_iterator is conditional on the node being
        #nonempty
        yield ((self.signal_block,), self.signal_count, self)
    def pop_node(self, prefix_signal=()):
        """
        Pops one terminal node and returns it and its signal. Cannot pop self.

        Args:
            prefix_signal: Tuple representing this node.

        Returns:
            (signal, node) where signal represents the node, and node has been
                popped.
        """
        d_gap, d_trie = self.descendants.items()[0]
        if len(d_trie.descendants) == 0:
            del self.descendants[d_gap]
            return prefix_signal + (d_gap,), d_trie
        else:
            return d_trie.pop_node(prefix_signal + (d_gap,))
    def leaf_iterator(self):
        """
        An iterator over all NONEMPTY descendant nodes and this node itself.
        See node_iterator to iterate over all descendant nodes whether empty or
        not. Iteration yields a tuple for each node:
        (signal represented by the node, signal counter, pointer to node)

        #??#This cannot be used while adding or removing nodes from the tree as
        #??#dictionaries cannot change size while iterating.
        """
        for d_trie in self.descendants.itervalues():
            for leaf in d_trie.leaf_iterator():
                if self.signal_block == (None, None):
                    yield leaf
                else:
                    #the decomposition is necessary to prepend
                    #self.signal_block
                    yield ((self.signal_block,) + leaf[0], leaf[1], leaf[2])
        #the following line is the critical difference between node_iterator
        #and leaf_iterator, as leaf_iterator is conditional on the node being
        #nonempty
        if len(self.signal_count) > 0:
            yield ((self.signal_block,), self.signal_count, self)
    def find_uniques(self, worst_ratio, absolute_min, maximum_secondary=None):
        """
        Returns all signals and their most responsible proteins where these
        protein are

        A. either unique (to force this requirement , set worst_ratio to None)
        OR ratio of the most responsible protein to the second most responsible
        protein is at least worst_ratio
        
        --AND--
        
        B. there are at least absolute_min signals from the most responsible
        protein

        Args:
            worst_ratio: Floating point ratio of protein that most frequently
                generated this signal to second most frequent protein source.
                If None, only those signals that are absolutely unique to a
                protein, i.e. there is only one protein that produced them,
                are returned. Otherwise, only those signals where the ratio
                exceeds the worst_ratio arguments are returned.
            absolute_min: Integer indicating for a signal to be returned it
                must have been generated at least an absolute_min number of
                times by the most responsible protein.
            maximum_secondary: If not None, integer indicating the maximum
                quantity of the second_worst protein.

        Returns:
            Dictionary of
                {signal:
           ((protein that most frequently generated this signal, its quantity),
             (tuples of proteins that are second most frequent source and their
              quantities), (total number of tertiary proteins)}
        """
        #this is an old implementation not using the iterators, but seems to
        #work fine
        uniques = {}
        #first see if self is a unique
        if len(self.signal_count) > 0:
            best = (None, 0)
            second = (None, 0)
            #find the best and second-best proteins
            for protein, count in self.signal_count.iteritems():
                if count > best[1]:
                    best = (protein, count)
                elif count > second[1]:
                    second = (protein, count)
            if (
                (best[1] >= absolute_min)

                and

                        (
                                   (worst_ratio is None and second[0] is None)
                                   or
                                   (worst_ratio is not None and second[1] == 0)
                                   or
                                   (worst_ratio is not None and
                                    float(best[1]) / second[1] >= worst_ratio)
                        )

                and

                        (
                                   maximum_secondary is None or
                                   second[0] is None or
                                   second[1] <= maximum_secondary
                        )
               ):
                #best protein, array of second-best protein,
                #total count of tertiary
                uniques.setdefault((self.signal_block,), [best, [second], 0])
                for protein, count in self.signal_count.iteritems():
                    if count == second[1] and protein != second[0]:
                        uniques[(self.signal_block,)][1].append((protein, count))
                    elif count < second[1]:
                        uniques[(self.signal_block,)][2] += count
        #now recursively check all descendants for uniques
        for block, descendant in self.descendants.iteritems():
            d_u = descendant.find_uniques(worst_ratio,
                                          absolute_min,
                                          maximum_secondary)
            for signal, (best, secondary, tertiary) in d_u.iteritems():
                if self.signal_block != (None, None):
                    uniques.setdefault((self.signal_block,) + signal,
                                       (best, secondary, tertiary))
                else:
                    uniques.setdefault(signal, (best, secondary, tertiary))
        return uniques
    def find_uniques_absolute(self, minimum_best, maximum_secondary):
        """
        Returns all signals and their most responsible and second most
        responsible proteins, as long as there are at least minimum_best
        primary proteins and at most maximum_secondary secondary proteins.
        """
        #this is an old implementation not using the iterators, but seems to
        #work fine
        uniques = {}
        #first see if self is a unique
        if len(self.signal_count) > 0:
            best = (None, 0)
            second = (None, 0)
            #find the best and second-best proteins
            for protein, count in self.signal_count.iteritems():
                if count > best[1]:
                    best = (protein, count)
                elif count > second[1]:
                    second = (protein, count)
            if best[1] >= minimum_best and second[1] <= maximum_secondary:
                #best protein, array of second-best protein,
                #total count of tertiary
                uniques.setdefault((self.signal_block,), [best, [second], 0])
                for protein, count in self.signal_count.iteritems():
                    if count == second[1] and protein != second[0]:
                        uniques[(self.signal_block,)][1].append((protein, count))
                    elif count < second[1]:
                        uniques[(self.signal_block,)][2] += count
        #now recursively check all descendants for uniques
        for block, descendant in self.descendants.iteritems():
            d_u = descendant.find_uniques_absolute(minimum_best,
                                                   maximum_secondary)
            #as above, best = (best protein, its quantity) and
            #second = (second best protein if present, its quantity)
            while len(d_u) > 0:
                signal, (best, second, tertiary) = d_u.popitem()
                if self.signal_block != (None, None):
                    uniques.setdefault((self.signal_block,) + signal,
                                       (best, second, tertiary))
                else:
                    uniques.setdefault(signal, (best, second, tertiary))
        return uniques
    def count_nodes(self):
        """
        Count total number of nodes in the trie rooted at this node.

        Returns:
            (Number of empty nodes, number of non-empty nodes)
        """
        empty, used = 0, 0
        for leaf in self.node_iterator():
            assert len(leaf[1]) >= 0
            if len(leaf[1]) == 0:
                empty += 1
            else:
                used += 1
        return empty, used
    def prune(self, signal):
        """
        Return the signal back along with its signal counts. If the node
        representing the signal has no descendants, remove it from the trie. If
        it has descendants, reset its signal count dictionary to empty.

        Args:
            signal: Signal to remove from the trie. If it is not present, an
                AssertionError. If signal is a prefix of a longer signal, i.e.
                the prefix trie has non-empty descendant nodes from the
                terminal node of the signal given, then the node remains but
                its signal count dictionary is emptied.

        Returns:
            (signal, signal_count dictionary for the signal)
        """
        #all signals to be pruned must be non-empty
        assert len(signal) > 0
        #the recursion uses the parent node to prune its descendant, hence the
        #root (None, None) node must be the only one to prune signals length 1
        if len(signal) == 1:
            assert self.signal_block == (None, None)
        #if the signal is longer than 1 block long, and this block is the root
        #(None, None) node, then the first signal block must be a descendant of
        #the root node
        elif self.signal_block == (None, None):
            assert signal[0] in self.descendants
        #if the signal is longer than 1 block long, and this block is not the
        #root node, then the first signal block must refer to this node and its
        #second signal block must be a descendant
        else:
            assert signal[0] == self.signal_block,\
                   ('self.signal_block ' + str(self.signal_block) +
                    '; signal ' + str(signal))
            assert signal[1] in self.descendants
        #if the signal is length one, then it is just a length one signal that
        #is a direct descendant of the root (None, None) node; recursion would
        #not result in having a non-root node receiving a length one signal
        if len(signal) == 1:
            if len(self.descendants[signal[0]].descendants) == 0:
                #remove descendant node if it has no children of its own
                return (signal, self.descendants.pop(signal[0]).signal_count)
            else:
                #if desendant node has children of its own, then empty signal
                #count dictionary but leave the node itself
                s_c = self.descendants[signal[0]].signal_count
                self.descendants[signal[0]].signal_count = {}
                return (signal, s_c)
        #if the signal is longer than 1 block, and this node is the root (None,
        #None) node, then pass the recursion to the next node
        elif self.signal_block == (None, None):
            return self.descendants[signal[0]].prune(signal)
        #if the signal is longer than 1 block and this is not a root node, then
        #it must have had its function called from its parent. if the signal is
        #length 2, prune the child and return (base case); otherwise, recurse
        else:
            if len(signal) == 2:
                if len(self.descendants[signal[1]].descendants) == 0:
                    #remove descendant node if it has no children of its own
                    return (signal,
                            self.descendants.pop(signal[1]).signal_count)
                else:
                    #if desendant node has children of its own, then empty
                    #signal count dictionary but leave the node itself
                    s_c = self.descendants[signal[1]].signal_count
                    self.descendants[signal[1]].signal_count = {}
                    return (signal, s_c)
            else:
                r = self.descendants[signal[1]].prune(signal[1:])
                return ((self.signal_block,) + r[0], r[1])
    def graft(self, signal, signal_count):
        """
        Add a signal to this trie with given protein signal counts. If the
        signal is already present in the trie, add the protein signal counts.

        Args:
            signal: The signal sequence to add. This function is recursive, so
                the full signal to be added is the concatenation of the signal
                represented by the node called and subsignal passed.
            signal_count: If the signal is new, the signal_count dictionary
                passed will be used as its protein source count. If the signal
                exists, the signal counts passed will be added to it.                

        Returns:
            Self.
        """
        #the signal passed cannot be empty
        assert len(signal) > 0
        #the signal's first block must either match this node (otherwise we
        #should not be here), or this is the root (None, None) node and the
        #first block is a descendant of the root
        assert signal[0] == self.signal_block or\
               self.signal_block == (None, None),\
               ('signal: ' + str(signal) +
                '; self.signal_block: ' + str(self.signal_block))
        #the signal_count must be non-zero, otherwise we will be adding empty
        #leaf nodes
        assert len(signal_count) > 0
        #if this is the root node, pass to descendant
        if self.signal_block == (None, None):
            self.descendants.setdefault(signal[0], SignalTrie(signal[0]))
            self.descendants[signal[0]].graft(signal, signal_count)
        #if this is not the root node, and the signal length is one, then we
        #must be referring to this node
        elif len(signal) == 1:
            for protein in signal_count:
                self.signal_count.setdefault(protein, 0)
                self.signal_count[protein] += signal_count[protein]
        #this is not the root node, and there is still signal to recurse
        #through
        else:
            self.descendants.setdefault(signal[1], SignalTrie(signal[1]))
            self.descendants[signal[1]].graft(signal[1:], signal_count)
        return self
    def merge(self, trie, cycles=None):
        """
        Iterates through all leafs in trie and grafts them onto this trie. This
        function can be used to copy an entire SignalTrie by initializing and
        emtpy trie via copy = SignalTrie((None, None)) and then
        copy.merge(original).

        Args:
            trie: The trie to merge with this one.
            cycles: If not None, merge only those leafs that are within cycles.

        Returns:
            Self.
        """
        #make sure this is called only from the root node
        assert self.signal_block == (None, None),\
               'merge can only be called on the root node'
        for leaf in trie.leaf_iterator():
            if cycles is None:
                self.graft(leaf[0], leaf[1])
            elif leaf[0][-1][0] <= cycles:
                self.graft(leaf[0], leaf[1])
        return self
    def truncating_projection(self, cycles):
        """
        This function projects the signals and their counts to signals that
        would be observed if the number of Edman cycles was truncated to a
        given number.

        Args:
            cycles: Number of Edman cycles to truncate to; integer.

        Returns:
            Self.
        """
        #first iterate through all leafs; for those leafs who require more than
        #'cycles' Edman cycles to observe, project them onto the shorter cycle
        #space and graft
        for leaf in self.leaf_iterator():
            #leaf[0][-1][0] is the number of Edman cycle at which last signal
            #drop was observed
            if leaf[0][-1][0] > cycles:
                #generate the projected signal; 's_b' is short for signal block
                projected_signal = tuple([s_b for s_b in leaf[0]
                                          if s_b[0] <= cycles])
                #this conditional necessary to prevent grafting empty
                #projected_signals. in this context, empty projected_signals
                #are those that are undetectable within 'cycles'
                if projected_signal:
                    self.graft(projected_signal, leaf[1])
        #now iterate and find all nodes and find those that fit within cycles
        #but point to descendants who do not.
        #note that the trie cannot be iterated and have nodes removed at the
        #same time, so it is necessary to break iterating and removing into two
        #steps
        terminal_node_pointers = [(node[2], descendant) #(term node, desc key)
         for node in self.node_iterator() for descendant in node[2].descendants
                        if node[0][-1][0] <= cycles and descendant[0] > cycles]
        #now iterate eliminate all the terminal node pointers
        for terminal_node, descendant_pointer in terminal_node_pointers:
            del terminal_node.descendants[descendant_pointer]
        #finally, we need to eliminate branches of nodes that do not lead to
        #leaves by iterating over all leaves and (special case) root node and
        #removing those descendants that have no leaves whatsoever to iterate
        #over
        terminal_leaf_pointers = []
        for leaf in self.leaf_iterator():
            for descendant, d_pointer in leaf[2].descendants.iteritems():
                has_subleaf = False
                for subleaf in d_pointer.leaf_iterator():
                    has_subleaf = True
                    break
                if not has_subleaf:
                    terminal_leaf_pointers.append((leaf[2], descendant))
        #special case of root node is because it is never a leaf and will never
        #be iterated over by leaf_iterator
        for descendant, d_pointer in self.descendants.iteritems():
            has_subleaf = False
            for subleaf in d_pointer.leaf_iterator():
                has_subleaf = True
                break
            if not has_subleaf:
                terminal_leaf_pointers.append((self, descendant))
        for leaf_pointer, descendant in terminal_leaf_pointers:
            del leaf_pointer.descendants[descendant]
        return self

def monte_carlo_trie(peptides, p, b, u, windows, sample_size=100,
                      random_seed=random.random(), silent=True):
    """
    Generates sample_size number of random signals based on parameters given
    and returns them represented by SignalTrie.

    Args:
        peptides: Dictionary of peptides to perform simulation on. Format as
            retruned by seqmodlibMP.attach()
        p: Probability of an Edman reaction succeeding. All Edman reactions are
            modelled as independent Bernoulli variables. p must be a number in
            [0, 1] inclusive.
        b: Photobleaching survival constant. Photobleaching is modeled as an
            exponential survival function s(k) = e^-bk, modeling the
            probability of a fluor surviving k laser exposures.
        u: Probability of a fluor being photobleached or unattached or
            otherwise nonfunctional and hence unobservable to begin with.
            Probabilities of each fluor being unobservable are independent.
        windows: Laser excitations patterns for each color stored as a
            dictionary. Keys are single-letter strings representing the labeled
            acids, and values are lists or tuples representing the sequence of
            windows during which peptide luminosities are observed. Each member
            of the sequence refers to the number of the Edman cycle for which a
            difference in luminosity is searched for. For example, (3, 4, 7)
            would represent the search for luminosity drops due to the third,
            fourth, and seventh Edman cycle. This means there will be an
            exposure between the second and third Edman, the third and fourth
            Edman, between the fourth and fifth Edman, followed by exposures
            directly before and after the seventh Edman. If position 1 is
            indicated, there is an initial exposure before any Edmans. Thus,
            the dictionary format is e.g. {'E': (3, 4, 7), 'K': (7, 8, 9)}.
            Each amino acid is assumed to be labeled with a distinct color. If
            more than one amino acid is labeled with the same color,
            homogenize() needs to be applied. Having windows corresponding to
            large numbers of Edman reactions means reactions will keep being
            applied that long; the tails need to survive all of these.
        sample_size: Number of samples to take of each protein. If a sample
            generates an empty signal, it is still considered as taken, however
            it is not added to the resulting trie.
        random_seed: Used to initialize random number generator. Strongly
            recommended for this to be randomly generated itself for each
            instance of this function call to ensure non-duplicated results.
        silent: Boolean indicating whether progress should be printed to
            standard output.

    Returns:
        SignalTrie that represents all the signals generated by the peptides.
    """
    if not silent:
        print('monte_carlo_trie starting at ' + str(datetime.datetime.now()))
        sys.stdout.flush()
    #for tracking computational progress; total signals is total number of
    #signals that will be generated
    signal_progress = 0
    last_print = 0
    total_signals = (sum([len(attached)
                          for attached in peptides.itervalues()]) *
                     sample_size)
    update_interval = total_signals / 10
    #initialize the signal tree to be generated
    return_trie = SignalTrie((None, None))
    #initialize random number generator for this call
    random.seed(random_seed)
    #convert windows to randsiggen
    rsg_windows = {}
    for acid, positions in windows.iteritems():
        rsg_windows.setdefault(acid, max(positions))
    for protein in peptides:
        #weigh number of samples taken by how many peptides are yielded by
        #this protein; hence total number of peptide samples for this
        #protein is sample_size * len(peptides[protein])
        #generate array indicating how many of each peptide to sample
        for i, peptide in enumerate(peptides[protein]):
            sample_counter = sample_size
            while sample_counter > 0:
                batch_size = min(10**3, sample_counter)
                randsiggen.random_signal(peptide, protein, p, b, u,
                                         rsg_windows, batch_size,
                                         random.randint(0, 10**8), return_trie)
                sample_counter -= batch_size
                signal_progress += batch_size
            if (not silent and
                (signal_progress - last_print >= update_interval)):
                print(str(signal_progress) + ' of ' + str(total_signals) +
                      ' signals generated by ' + str(datetime.datetime.now()))
                last_print = signal_progress
                sys.stdout.flush()
    return return_trie
