import pandas as pd

from drugex.data.corpus.vocabulary import VocSmiles, VocGraph
from drugex.logs import logger
from drugex.data.interfaces import DataSplitter, FragmentPairEncoder
from drugex.molecules.converters.interfaces import ConversionException
from drugex.molecules.interfaces import MolSupplier
from drugex.parallel.collectors import ListExtend
from drugex.parallel.evaluator import ParallelSupplierEvaluator
from drugex.parallel.interfaces import ParallelProcessor


class SequenceFragmentEncoder(FragmentPairEncoder):
    """
    Encode fragment-molecule pairs for the sequence-based models.

    """

    def __init__(self, vocabulary=VocSmiles(True), update_voc=True, throw = False):
        self.vocabulary = vocabulary
        self.updateVoc = update_voc
        self.throw = throw

    def encodeMol(self, sequence):
        """
        Encode a molecule sequence.

        Args:
            sequence: sequential representation of the molecule (i.e. SMILES)

        Returns:
            a `tuple` containing the obtained tokens from the sequence (if any) and the corresponding sequence of codes
        """

        tokens = None
        if self.updateVoc:
            tokens = self.vocabulary.addWordsFromSeq(sequence)
        elif self.throw:
            tokens = self.vocabulary.removeIfNew(sequence)

        if tokens:
            # Encode all but end tokens
            output = self.vocabulary.encode([tokens[: -1]])
            code = output[0].reshape(-1).tolist()
            return tokens, code
        return tokens, None

    def encodeFrag(self, mol, mol_tokens, frag):
        """Encode a fragment.

        Is called by `FragmentPairsEncodedSupplier` with the `mol`
        argument being the output of the above `encodeMol` method.
        """
        tokens = None
        if self.updateVoc:
            tokens = self.vocabulary.addWordsFromSeq(frag, ignoreConstraints=True)
        elif self.throw:
            tokens = self.vocabulary.removeIfNew(frag, ignoreConstraints=True)

        if tokens:
            # Encode all but end tokens
            output = self.vocabulary.encode([tokens[: -1]])
            code = output[0].reshape(-1).tolist()
            return code

    def getVoc(self):
        return self.vocabulary

class GraphFragmentEncoder(FragmentPairEncoder):
    """
    Encode molecules and fragments for the graph-based transformer (`GraphModel`).
    """

    def __init__(self, vocabulary=VocGraph()):
        """
        Initialize this instance with the vocabulary to use.

        Args:
            vocabulary: used to perform the encoding
        """

        self.vocabulary = vocabulary

    def encodeMol(self, smiles):
        """
        Molecules are encoded together with fragments -> we just pass the smiles back as both tokens and result of encoding.

        Args:
            smiles:

        Returns:
            The input smiles as both the tokens and as the encoded result.
        """

        return smiles, smiles

    def encodeFrag(self, mol, mol_tokens, frag):
        """
        Encode molecules and fragments at once.

        Args:
            mol: parent molecule SMILES (from `encodeMol`)
            mol_tokens: molecule SMILES (from `encodeMol`)
            frag: SMILES of the fragment in the parent molecule

        Returns:
            One line of the graph-encoded data.
        """

        if mol == frag:
            return None
        try:
            output = self.vocabulary.encode([mol], [frag])
            f, s = self.vocabulary.decode(output)

            assert mol == s[0]
            #assert f == frag[0]
            code = output[0].reshape(-1).tolist()
            return code
        except Exception as exp:
            logger.warn(f'The following exception occured while encoding fragment {frag} for molecule {mol}: {exp}')
            return None

    def getVoc(self):
        return self.vocabulary

class FragmentPairsEncodedSupplier(MolSupplier):
    """
    Transforms fragment-molecule pairs to the encoded representation used by the fragment-based DrugEx models.

    """

    class FragmentEncodingException(ConversionException):
        """
        Raise this when a fragment failed to encode.
        """

        pass

    class MoleculeEncodingException(ConversionException):
        """
        Raise this when the parent molecule of the fragment failed to be encoded.
        """

        pass

    def __init__(self, pairs, encoder):
        """
        Initialize from a `DataFrame` containing the fragment-molecule pairs.

        Args:
            pairs (list): list of (fragment, molecule) `tuple`s that each denotes one fragment-molecule pair
            encoder: a `FragmentPairEncoder` handling encoding of molecules and fragments
        """
        self.encoder = encoder
        self.pairs = iter(pairs)

    def next(self):
        """
        Get the next pair and encode it with the encoder.

        Returns:
            `tuple`: (str, str) encoded form of fragment-molecule pair
        """

        pair = next(self.pairs) # (fragment, molecule)

        # encode molecule
        tokens, encoded_mol = self.encoder.encodeMol(pair[1])
        if not tokens:
            raise self.MoleculeEncodingException(f'Failed to encode molecule: {pair[1]}')

        # encode fragment
        encoded_frag = self.encoder.encodeFrag(pair[1], tokens, pair[0])
        if not encoded_frag:
            raise self.FragmentEncodingException(f'Failed to encode fragment {pair[0]} from molecule: {pair[1]}')

        return encoded_frag, encoded_mol


class FragmentPairsSupplier(MolSupplier):
    """
    Produces fragment-molecule pairs from input molecules.

    """

    def __init__(self, molecules, fragmenter, max_bonds=None):
        """

        Args:
            molecules: the input molecules as a `list`-like object or an iterator
            fragmenter: an instance of `Fragmenter
        """
        self.molecules = molecules if hasattr(molecules, "__next__") else iter(molecules)
        self.fragmenter = fragmenter
        self.currentBatch = None
        self.maxBonds = max_bonds

    def next(self):
        """
        Generate the next fragment-molecule pair.

        Returns:
            a (fragment, molecule) `tuple`
        """

        if not self.currentBatch:
            batch = None
            while not batch:
                # the fragmenter generates multiple pairs at once from one molecule, we use batching to return them one by one
                batch = self.fragmenter(next(self.molecules))
            self.currentBatch = iter(batch)
        try:
            frags = next(self.currentBatch)
        except StopIteration:
            self.currentBatch = None
            return None
        return frags

class FragmentCorpusEncoder(ParallelProcessor):
    """
    Fragments and encodes fragment-molecule pairs in parallel. Each encoded pair is used as input to the fragment-based DrugEx models.

    """

    class FragmentPairsCollector(ListExtend):
        """
        A simple `ResultCollector` that extends an internal `list`. It can also wrap another instance of itself.
        """

        def __init__(self, other=None):
            """

            Args:
                other: another instance of `FragmentPairsCollector` to call after extending
            """
            super().__init__()
            self.other = other

        def __call__(self, result):
            self.items.extend(result[0])
            if self.other:
                self.other(result)

    def __init__(self, fragmenter, encoder, pairs_splitter=None, n_proc=None, chunk_size=None):
        """

        Args:
            fragmenter (MolConverter): a `MolConverter` that returns a `list` of (fragment, molecule) `tuple`s for a given molecule supplied as its SMILES string. See the reference implementation in `Fragmenter`.
            encoder:  a `FragmentPairEncoder` that handles how molecules and fragments are encoded in the final result
            pairs_splitter: a `ChunkSplitter` that divides the generated molecule-fragment pairs from the "fragmenter" to splits (i.e. test and train)
            n_proc: number of processes to use for parallel operations
            chunk_size: maximum size of data chunks processed by a single process (can save memory)
        """

        super().__init__(n_proc, chunk_size)
        self.fragmenter = fragmenter
        self.encoder = encoder
        self.pairsSplitter = pairs_splitter

    def getFragmentPairs(self, mols, collector):
        """
        Apply the given "fragmenter" in parallel.

        Args:
            mols: Molecules represented as SMILES strings.
            collector: The `ResultCollector` to apply to fetch the result per process.

        Returns:
            `None`
        """

        evaluator = ParallelSupplierEvaluator(
            FragmentPairsSupplier,
            kwargs={
                "fragmenter" : self.fragmenter
            },
            chunk_size=self.chunkSize,
            chunks=self.chunks,
            n_proc=self.nProc
        )
        evaluator.apply(mols, collector, desc_string="Creating fragment-molecule pairs")

    def splitFragmentPairs(self, pairs):
        """
        Use the "pairs_splitter" to get splits of the calculated molecule-fragment pairs from `FragmentCorpusEncoder.getFragmentPairs()`

        Args:
            pairs: pairs generated by the "fragmenter"

        Returns:
            splits from the specified "splitter"

        """

        return self.pairsSplitter(pairs) if self.pairsSplitter else [pairs]

    def encodeFragments(self, pairs, collector):
        """
        Encodes fragment-pairs obtained from `FragmentCorpusEncoder.getFragmentPairs()` with the specified `FragmentPairEncoder` initialized in "encoder".

        Args:
            pairs: `list` of (fragment, molecule) `tuple`s to encode
            collector: The `ResultCollector` to apply to fetch encoding data from each process.

        Returns:
            `None`
        """

        evaluator = ParallelSupplierEvaluator(
            FragmentPairsEncodedSupplier,
            kwargs={
                'encoder': self.encoder,
            },
            chunk_size=self.chunkSize,
            chunks=self.chunks,
            n_proc=self.nProc
        )
        evaluator.apply(pairs, collector, desc_string="Encoding fragment-molecule pairs.")

    def apply(self, mols, fragmentPairsCollector=None, encodingCollectors=None):
        """
        Apply fragmentation and encoding to the given molecules represented as SMILES strings. Collectors can be used to fetch fragment-molecule pairs and the final encoding with vocabulary.

        Args:
            mols: `list` of molecules as SMILES strings
            fragmentPairsCollector: an instance of `ResultCollector` to collect results of the fragmentation (the generated fragment-molecule `tuple`s from the given "fragmenter").
            encodingCollectors: a `list` of `ResultCollector` instances matching in length the number of splits given by the "pairs_splitter". Each `ResultCollector` receives a (data, `FragmentPairsEncodedSupplier`) `tuple` of the currently finished process.

        Returns:
            `None`
        """

        pairs_collector = self.FragmentPairsCollector(fragmentPairsCollector)
        self.getFragmentPairs(mols, pairs_collector)
        splits = self.splitFragmentPairs(pairs_collector.getList())
        if encodingCollectors and len(encodingCollectors) != len(splits):
            raise RuntimeError(f'The number of encoding collectors must match the number of splits: {len(encodingCollectors)} != {len(splits)}')
        for split_idx in range(len(splits)):
            self.encodeFragments(splits[split_idx], encodingCollectors[split_idx] if encodingCollectors else None)

class FragmentPairsSplitter(DataSplitter):
    """
    A `DataSplitter` to be used to split molecule-fragment pairs into training and test data.
    """

    def __init__(self, ratio=0.2, max_test_samples=1e4, train_collector=None, test_collector=None, unique_collector=None, make_unique=False, seed=None):
        """
        Set settings for the splitter.

        Args:
            ratio: Ratio of fragment-molecule pairs to move to the test set.
            max_test_samples: Maximum number of test samples (to speed up calculations).
            train_collector: a `ResultCollector` to collect the training set
            test_collector: a `ResultCollector` to collect the test set
            unique_collector: a `ResultCollector` to collect the 'unique' data set (only one example per unique fragment)
            make_unique: make the training set with only unique fragments in addition
            seed: fix the random seed to always get the same split
        """

        self.ratio = ratio
        self.maxTestSamples = max_test_samples
        self.uniqueCollect = unique_collector
        self.trainCollect = train_collector
        self.testCollect = test_collector
        self.makeUnique = make_unique
        self.seed = seed

    def __call__(self, pairs):
        """
        Split the input `list` to the desired data sets. The split is done on the fragments so that no single fragment is contained in both the training and test split.

        Args:
            pairs: `list` of (fragment, molecule) tuples

        Returns:
            a `tuple` with three pandas `DataFrame` instances corresponding to the test, train and uniqe sets, respectively if "unique_only" is `False`.  If "unique_only" is `True`, only the unique data set is created.
        """

        df = pd.DataFrame(pairs, columns=["Frags", "Smiles"])
        frags = set(df.Frags)
        test_len = int(len(frags) * self.ratio)
        if self.seed:
            test_in = df.Frags.drop_duplicates().sort_values()
        else:
            test_in = df.Frags.drop_duplicates()
        if test_len > int(self.maxTestSamples):
            logger.warning(f'To speed up the training, the test set size was automatically capped at {self.maxTestSamples} fragments instead of the default 10% of original data, which would have been: {test_len}.')
            test_in = test_in.sample(int(self.maxTestSamples), random_state=self.seed)
        else:
            test_in = test_in.sample(test_len, random_state=self.seed)
        test = df[df.Frags.isin(test_in)].values.tolist()
        train = df[~df.Frags.isin(test_in)].values.tolist()
        unique = None
        if self.makeUnique:
            unique = df.drop_duplicates(subset="Frags").values.tolist()

        if self.trainCollect:
            self.trainCollect(train)
        if self.testCollect:
            self.testCollect(test)
        if self.uniqueCollect:
            self.uniqueCollect(unique)

        if unique:
            return test, train, unique
        else:
            return test, train