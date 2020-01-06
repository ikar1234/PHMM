import numpy as np
from typing import List
from warnings import warn

from PHMM.src.Records import Record


class PHMM:
    record: Record
    # transition probability matrix
    P: np.ndarray
    # emission probability matrix
    Q: List
    # initial probabilities
    columns: List[int]
    # alphabet
    alph: List[str]
    # length of the model
    l: int

    def __init__(self, record=None, alph="dna"):
        """
        Construct the transition and emission probability matrices,
        which are the main parameters of the model.
        :param record: record with the position-weight matrix
        """
        if record is not None:
            self.record = record
            self._annotate()
            # infer from the record
            if alph == "dna" or self.record.meta.get('alength', 20) == 4:
                self.alph = "A C G T".split(' ')
            else:
                self.alph = "A C D E F G H I K L M N P Q R S T V W Y".split(' ')
            self.P = self._mle_P()
            self.Q = self._mle_Q()
        else:
            self.record = None
            self.P = None
            self.Q = None
            self.columns = None
            self.alph = "A C G T".split(' ')

    def _annotate(self):
        """
        Infer the matching and inserting states
        from the alignment  by counting the number of gaps.
        :return: list of states (M=0,I=1,D=2)
        """
        # if there are exactly 0 gaps, the state is matching
        # if the gaps are <50%, the state is matching (=0)
        # else the state is inserting (=1)

        gaps_func = lambda x: 0 if x <= 0.5 else 1
        self.columns = [gaps_func(x) for x in self.record.matrix[0, :]]

        # the length of the PHMM is the number of matching states
        self.l = sum([x == 0 for x in self.columns])

    def _mle_P(self) -> np.ndarray:
        """
        Maximum likelihood estimator for the transition matrix.
        :return:
        """
        # states
        # TODO: matrices become too sparse
        st = self.columns

        n_match = sum([x == 0 for x in st])
        n_insert = len(st) - n_match
        # sequences from alignment
        seqs = self.record.seqs

        seq_paths = list(PHMM.make_paths(self.columns, seqs))


        match_trans = np.zeros(())
        insert_trans = np.zeros(())
        delete_trans = np.zeros(())
        P_hat = match_trans[match_trans, insert_trans, delete_trans]

        # a long sequence of cases

        P_hat[0] = self._normalize(P_hat[0], pseudo=True)
        P_hat[1] = self._normalize(P_hat[1], pseudo=True)
        P_hat[2] = self._normalize(P_hat[2], pseudo=True)

        return P_hat

    def _mle_Q(self):
        """
        Maximum likelihood estimator for the emission matrix.
        :return:
        """
        # count the number of emissions

        # list of list of dictionary objects
        st = self.columns
        n_match = sum([x == 0 for x in st])
        n_insert = len(st) - n_match
        # emission probabilities for all match and insert states
        match_emission = np.zeros((n_match, len(self.alph)))
        insert_emission = np.zeros((n_insert, len(self.alph)))
        c_m = c_i = 0

        for i in range(len(st)):
            # track the number of matching and inserting states
            if st[i] == 0:
                # dont count the gaps and normalize
                match_emission[c_m, :] = self.record.matrix[1:, i]
                c_m += 1
            else:
                insert_emission[c_i, :] = self.record.matrix[1:, i]
                c_i += 1
        assert c_m == n_match and c_i == n_insert

        Q_hat = [match_emission, insert_emission]

        return Q_hat

    def evaluate(self, seq: str, log: bool = False) -> float:
        """
        Given an observation sequence and a model, find the likelihood of the sequence with respect to the model.
        This value can be estimated from the profile itself.
        :param seq: sequence
        :param log: take the log likelihood. This could also result in the likelihoood being -inf.
        :return: (log-)likelihood of the sequence with respect to the model
        """
        m = self.record.matrix
        # protein/dna alphabet, characters in the same order as the PWM
        if log:
            with np.errstate(divide='ignore'):
                lhd = np.sum([np.log(m[self.alph.index(seq[i]), i]) for i in range(len(seq))])
        else:
            lhd = np.product([m[self.alph.index(seq[i]), i] for i in range(len(seq))])
        return lhd

    def viterbi_decoding(self, seq: str):
        """
        Compute the most probable path of a sequence in
        the HMM using dynamic programming.
        The run-time of algorithm is O(|Q|^2*n),  where |Q| is the number of
        emission-states.
        @:param seq: input sequence
        @:param P: transition probability matrix
        @:param Q: emission probability matrix
        :return: most probable path and the probability of it
        """
        l_seq = len(seq)
        # we have 3 states - match,insert,delete
        trellis = np.zeros((3, l_seq))
        # in order to obtain the sequence itself
        backpointers = np.zeros((3, l_seq), dtype=np.int)
        # initialize
        trellis[:, 0] = np.array([1, 0, 0])

        # since we have only 3 states, the run-time is O(n)
        for n in range(1, l_seq):
            for t in range(3):
                last_column = trellis[:, n - 1] * self.P[t, :]
                ind = self.alph.index(seq[n])
                print(self.Q[t][n])
                trellis[t, n] = self.Q[t][n].get(seq[n]) * np.max(last_column)
                backpointers[t, n] = np.argmax(last_column)
        last = np.argmax(trellis[:, l_seq - 1])

        yield last
        for i in range(l_seq - 1, 0, -1):
            yield backpointers[last, i]
            last = backpointers[last, i]

    def viterbi_training(self, seqs: List[str]):
        """
        Train each sequence using the models parameters
        and concatenate sequences to a multiple alignment.
        :param seqs:
        :return:
        """
        # paths of states (insert=0, match=1, delete=2) for every sequence
        paths = [list(self.viterbi_decoding(seq)) for seq in seqs]
        print(paths)

        for s in range(len(seqs)):
            new_seq = ''
            for i in range(len(seqs[s])):
                # if in matching or inserting state
                if paths[s][i] < 2:
                    new_seq += seqs[s][i]
                # if any other state is inserting
                elif any([paths[other_seq][i] == 2 for other_seq in range(len(seqs))]):
                    new_seq += '-'
                # deleting state
                else:
                    new_seq += '-'
            yield new_seq

    def _forward(self, seq: str) -> np.ndarray:
        """
        Forward algorithm.
        :param seq:
        :return: Probability that the prefix
        sequence of symbols is generated and the system is in
        state k at time i.
        """
        if self.P is None or self.Q is None:
            raise ValueError("This algorithm is not applicable if both matrices are unknown.")

        # number of states
        n = self.P.shape[0]
        l = len(seq)

        fwd = np.zeros((n, l), dtype=np.float)
        fwd[0, 0] = 1

        for i in range(1, l):
            for st in range(n):
                letter = self.alph.index(seq[i])
                fwd[st, i] = self.Q[st][i, letter] * np.dot(fwd[:, i - 1], self.P[:, st])

        return fwd

    def _backward(self, seq: str) -> np.ndarray:
        """
        Backward algorithm.
        :param seq:
        :return: Probability that the prefix sequence of symbols
        starts in state k at time i and then generates the
        sequence of symbols x_{i+1}..x_L .
        """
        if self.P is None or self.Q is None:
            raise ValueError("This algorithm is not applicable if both matrices are unknown.")
        # number of states
        n = self.P.shape[0]
        l = len(seq)
        bwd = np.zeros((n, l))

        # initialize
        bwd[:, l - 1] = self.P[:, 0]

        for i in range(l - 2, 1, -1):
            for k in range(n):
                # What if l > len(seq)?
                letter = self.alph.index(seq[i + 1])
                bwd[k, i] = np.sum(
                    self.P[k, :] * np.array([self.Q[i][st, letter] for st in range(n)]) * bwd[:, i + 1])

        return bwd

    def obs_prob(self, seq):
        """
        Probability that the sequence is generated from the PHMM.
        We use the forward algorithm to compute that.
        :param seq: sequence
        """
        return np.sum(self._forward(seq)[:, -1])

    def _forwards(self, seqs: List[str]):
        for s in seqs:
            yield self._forward(s)

    def _backwards(self, seqs: List[str]):
        for s in seqs:
            yield self._backward(s)

    def baum_welch(self, seqs: List[str], n_iter: int = 400):
        """
        Baum-Welch algorithm to estimate the parameters of the HMM.
        :param seqs: list of sequences
        :param n_iter: number of iterations
        :return:
        """

        # mean length of the sequences
        l = int(np.mean([len(x) for x in seqs]))

        # number of training sequences
        n_seq = len(seqs)

        self.P = np.random.rand(3, 3)
        # TODO
        self.Q = [np.random.rand(3, l), np.random.rand(3, 3)]

        P_hat = self.P
        Q_hat = self.Q

        for i in range(n_iter):
            # forward probabilities for each sequence
            fwd = list(self._forwards(seqs))
            # backward probabilities for each sequence
            bwd = list(self._backwards(seqs))

            for k in range(3):
                for l in range(3):
                    P_hat[k, l] = \
                        np.sum([1 / self.obs_prob(seqs[j]) *
                                sum([fwd[j][k, i] * self.P[k, l] * self.Q[
                                    l, self.alph.index(seqs[j][i])] * bwd[j][l, i + 1] for i in
                                     range(len(seqs[j]) - 1)]) for j in
                                range(n_seq)])
                for b in range(len(self.alph)):
                    Q_hat[k, b] = sum(
                        [1 / self.obs_prob(seqs[j]) * sum(
                            [fwd[j][k, i] * bwd[j][l, i] for i in range(len(seqs[j])) if seqs[j][i] == self.alph[b]])
                         for j in range(n_seq)])

            self.P = self._normalize(P_hat)
            self.Q = self._normalize(Q_hat, cols=True)

    def train(self, seqs: List[str], method="baum_welch") -> str:
        """
        Produces a multiple sequence alignment from a set
        of sequences.
        :param method: Training method. Most common is the
        Baum-Welch algorithm.
        :param seqs: List of sequences
        :return: multiple sequence alignment
        """
        # Viterbi-Training - uses the Viterbi algorithm for every input sequence and
        # aligns the sequence using the generated observed states.
        if method == 'viterbi':
            if self.P is None and self.Q is None:
                raise ValueError("Viterbi training cannot be used if the state path is unknown.\n"
                                 "Use 'baum-welch' instead.")
        # Baum-Welch is used when the state path is unknown.
        elif method == 'baum_welch':
            if self.P is not None and self.Q is not None:
                warn("Baum-Welch assumes that the state path is unknown, but a profile was given.\n"
                     "The profile will be ignored. Use 'viterbi' if you want to use the profile.")
            self.baum_welch(seqs)
        else:
            raise ValueError("Parameter 'method' must be either 'viterbi' or 'baum_welch'.")

        new_seqs = list(self.viterbi_training(seqs))

        return '\n'.join(new_seqs)

    @classmethod
    def compare(cls, obj1, obj2, seq: str, verbose: bool = False):
        """
        Compare two PHMM using the forward algorithm, which computes
        the probability P(seq | M).
        :param obj1: first PHMM
        :param obj2: second PHMM
        :param seq: sequence
        :param verbose: output explanation
        :return: ratio of probabilities
        """
        if not isinstance(obj2, PHMM):
            raise ValueError("The second argument must be a PHMM.")
        ratio = cls.obs_prob(obj1, seq) / cls.obs_prob(obj2, seq)
        if ratio > 1:
            if verbose:
                print(f'The first PHMM is better than the second '
                      f'with ratio {ratio}')
        elif ratio == 1:
            if verbose:
                print(f'Both PHMM are equally well.')
        else:
            if verbose:
                print(f'The second PHMM is better than the first '
                      f'with ratio {1 / ratio}')
        return ratio

    @staticmethod
    def _normalize(m: np.ndarray, cols=False, pseudo=True):
        # TODO: Laplace rule
        """
        Normalize a matrix such that each row sums up to 1.
        @:param cols: normalize the columns instead
        @:param pseudo: add pseudocounts
        """
        if cols:
            if pseudo:
                return (m + 1) / (np.sum(m, axis=0) + m.shape[0])
            else:
                return m / np.sum(m, axis=0)
        else:
            # use broadcasting
            row_sums = m.sum(axis=1)
            if pseudo:
                new_matrix = (m + 1) / (row_sums[:, None] + m.shape[1])
            else:
                new_matrix = m / row_sums[:, None]
            return new_matrix

    @staticmethod
    def make_paths(states, seqs):
        """
        Infer the states of each sequence in the multiple alignment.
        :param states: states of the PHMM
        :param seqs: sequences with gaps
        :return:
        """
        for s in seqs:
            path = []
            for i in range(len(states)):
                if s[i] == '-':
                    # gap in a matching state
                    if states[i] == 0:
                        path.append(2)
                    else:
                        pass
                # symbol
                else:
                    # match
                    if states[i] == 0:
                        path.append(0)
                    # insert
                    else:
                        path.append(1)

            yield path
