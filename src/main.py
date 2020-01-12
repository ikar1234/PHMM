import numpy as np
from numpy import log, argmax
from typing import List
from warnings import warn
from pandas import DataFrame
from PHMM.src.Records import Record
from itertools import zip_longest


class PHMM:
    record: Record
    # transition probability matrix
    P: DataFrame
    # emission probability matrix
    Q: DataFrame
    # alphabet
    alph: List[str]
    # length of the model
    l: int
    # number of matching states
    k: int
    # states
    states: List[str]

    def __init__(self, record=None, alph="dna", threshold: float = 0.25, pscount: int = 0.01):
        """
        Create either a PHMM from multiple sequence alignment or an empty PHMM object.
        :param record: record with the the sequences
        :param alph: alphabet
        :param threshold: threshold for inferring which state is matching
        and which inserting; usually taken to be 0.5 or 0.25
        :param pscount: pseudocount for normalization
        """
        if record is not None:
            self.record = record
            # infer from the record
            if alph == "dna" or self.record.alph == "dna":
                self.alph = "A C G T".split(' ')
            else:
                self.alph = "A C D E F G H I K L M N P Q R S T V W Y".split(' ')
            # initialize parameters
            self.initialize(threshold, pscount)
        else:
            self.record = None
            self.P = None
            self.Q = None
            self.columns = None
            self.alph = "A C G T".split(' ')

    def initialize(self, threshold, pscount):
        """
        Derive the MLE for the transition and emission matrices using
        the multiple sequence alignment and the parameters.
        :param threshold: threshold for inferring which state is matching
        and which inserting; usually taken to be 0.5 or 0.25
        :param pscount: pseudocount for normalization
        :return:
        """

        # multiple sequence alignment
        aln = self.record.seqs

        aln_cols = list(zip(*(aln)))
        m, n = len(aln), len(aln_cols)  # m sequences, n columns

        # indices of columns where '-' count is below threshold
        match_cols = [i for i in range(n) if aln_cols[i].count('-') / m < threshold]

        # state names
        self.k = len(match_cols)  # there k M-states
        states_ = [(f'M{i}', f'D{i}', f'I{i}') for i in range(1, self.k + 1)]
        self.states = ['S', 'I0'] + [i for j in states_ for i in j] + ['E']

        # initialize matrices
        self.P = DataFrame(data=0.0, index=self.states, columns=self.states)
        self.Q = DataFrame(data=0.0, index=self.states, columns=self.alph)

        for seq in aln:
            state_ix = 0
            last_state = 'S'
            for i in range(n):
                if i in match_cols:
                    state_ix += 1
                    if seq[i] != '-':
                        current_state = 'M' + str(state_ix)
                        self.Q.loc[current_state, seq[i]] += 1
                    else:
                        current_state = 'D' + str(state_ix)

                    self.P.loc[last_state, current_state] += 1
                    last_state = current_state

                elif seq[i] != '-':
                    current_state = 'I' + str(state_ix)
                    self.P.loc[last_state, current_state] += 1
                    self.Q.loc[current_state, seq[i]] += 1
                    last_state = current_state

            self.P.loc[last_state, 'E'] += 1

        # normalize
        self.P = self._normalize(self.P, pseudo=False).round(3)
        self.Q = self._normalize(self.Q, pseudo=False).round(3)

        # add pseudocounts
        self.P.iloc[:2, 1:4] += pscount
        self.P.iloc[-4:-1, -2:] += pscount
        for i in range(self.k):
            self.P.iloc[i * 3 - 1:i * 3 + 2, i * 3 + 1:i * 3 + 4] += pscount
            self.Q.iloc[i * 3 + 1:i * 3 + 3, :] += pscount
        self.Q.iloc[-2, :] += pscount

        # normalize again
        self.P = PHMM._normalize(self.P)
        self.Q = PHMM._normalize(self.Q)

    def viterbi_training(self, seqs: List[str]):
        """
        Train each sequence using the models parameters
        and concatenate sequences to a multiple alignment.
        :param seqs:
        :return: multiple sequence alignment
        """
        # paths of states (insert=0, match=1, delete=2) for every sequence
        paths = [self.viterbi_decoding(seq) for seq in seqs]

        cols = [[x[0] for x in l if x is not None] for l in zip_longest(*paths)]

        # construct the multiple alignment
        msa = ""
        for p in range(len(paths)):
            seq = ""
            n_chars = 0
            for char in range(len(paths[p])):
                # delete state
                if paths[p][char][0] == 'D':
                    seq += '-'
                # insert state
                elif paths[p][char][0] == 'I':
                    seq += seqs[p][n_chars]
                    n_chars += 1
                # some other sequence in insert state -> place a gap
                elif 'I' in cols[char]:
                    seq += '-'
                # match state
                else:
                    seq += seqs[p][n_chars]
                    n_chars += 1
            msa += seq
            msa += '\n'

        return msa

    def viterbi_decoding(self, seq):
        """
        Compute the most probable path of a sequence in
        the HMM using dynamic programming.
        The run-time of algorithm is O(|Q|^2*n),  where |Q| is the number of
        emission-states.
        @:param seq: input sequence
        :return: most probable path and the probability of it
        """
        n = len(seq)

        # initialize Viterbi graph with each node contains [score, predecessor]
        V = {(i, j): [-np.inf, ''] for i in self.states for j in range(n + 1)}

        # column 0, all Ds
        V['D1', 0] = [np.log(self.P.loc['S', 'D1']), ('S', 0)]

        for j in range(2, self.k + 1):
            V['D' + str(j), 0] = [V['D' + str(j - 1), 0][0] +
                                  log(self.P.loc['D' + str(j - 1), 'D' + str(j)]), ('D' + str(j - 1), 0)]

        # top 2 rows, I0, M1 of the remaining columns
        V['I0', 1] = [log(self.P.loc['S', 'I0']) + log(self.Q.loc['I0', seq[0]]),
                      ('S', 0)]
        V['M1', 1] = [log(self.P.loc['S', 'M1']) + log(self.Q.loc['M1', seq[0]]),
                      ('S', 0)]

        for i in range(2, n + 1):
            V['I0', i] = [V['I0', i - 1][0] + log(self.P.loc['I0', 'I0']) +
                          log(self.Q.loc['I0', seq[i - 1]]), ('I0', i - 1)]

            V['M1', i] = [V['I0', i - 1][0] + log(self.P.loc['I0', 'M1']) +
                          log(self.Q.loc['M1', seq[i - 1]]), ('I0', i - 1)]

            V['D1', i] = [V['I0', i][0] + log(self.P.loc['I0', 'D1']) +
                          log(self.Q.loc['D1', seq[i - 1]]), ('I0', i)]

        # main recurrence
        for i in range(1, n + 1):
            for j in range(2, self.k + 1):
                prev_states = 'M{} I{} D{}'.format(j - 1, j - 1, j - 1).split()

                edges = [V[s, i - 1][0] + log(self.P.loc[s, 'I' + str(j - 1)]) for s in prev_states]
                V['I' + str(j - 1), i] = [max(edges) + log(self.Q.loc['I' + str(j - 1), seq[i - 1]]),
                                          (prev_states[argmax(edges)], i - 1)]

                edges = [V[s, i - 1][0] + log(self.P.loc[s, 'M' + str(j)]) for s in prev_states]
                V['M' + str(j), i] = [max(edges) + log(self.Q.loc['M' + str(j), seq[i - 1]]),
                                      (prev_states[argmax(edges)], i - 1)]

                edges = [V[s, i][0] + log(self.P.loc[s, 'D' + str(j)]) for s in prev_states]
                V['D' + str(j), i] = [max(edges), (prev_states[argmax(edges)], i)]

        # the last row (I{k_matches})
        for i in range(1, n + 1):
            prev_states = f'M{self.k} I{self.k} D{self.k}'.split()

            edges = [V[s, i - 1][0] + log(self.P.loc[s, 'I' + str(j)]) for s in prev_states]
            V['I' + str(j), i] = [max(edges) + log(self.Q.loc['I' + str(j), seq[i - 1]]),
                                  (prev_states[argmax(edges)], i - 1)]

        # End node
        node = (prev_states[argmax([V[s, n][0] + log(self.P.loc[s, 'E']) for s in prev_states])], n)

        # path reconstruction by backtracing
        path = []
        while node[0] != 'S':
            path.append(node[0])
            node = V[node][1]
        return path[::-1]

    def obs_prob(self, seq):
        """
        Probability that the sequence is generated from the PHMM.
        We use the forward algorithm to compute that.
        :param seq: sequence
        # TODO: test
        """
        return np.sum(self._forward(seq)[:, -1])

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
            for st in range(len(self.states)):
                fwd[st, i] = self.Q[seq[i]][st] * np.dot(fwd[:, i - 1], self.P[self.states[st]])

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
        n = len(self.states)
        l = len(seq)
        bwd = np.zeros((n, l))

        # initialize
        bwd[:, l - 1] = self.P['S']
        for i in range(l - 2, 1, -1):
            for k in range(n):
                bwd[k, i] = np.sum(self.P[self.states[k]] * self.Q[seq[i + 1]] * bwd[:, i + 1])

        return bwd

    def _forwards(self, seqs: List[str]):
        for s in seqs:
            yield self._forward(s)

    def _backwards(self, seqs: List[str]):
        for s in seqs:
            yield self._backward(s)

    def baum_welch(self, seqs: List[str], n_iter: int = 1):
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

        # number of matching states
        self.k = l
        states_ = [(f'M{i}', f'D{i}', f'I{i}') for i in range(1, self.k + 1)]
        self.states = ['S', 'I0'] + [i for j in states_ for i in j] + ['E']

        n_states = len(self.states)

        # initialize both matrices randomly
        self.P = DataFrame(data=0.0, index=self.states, columns=self.states)
        self.Q = DataFrame(data=0.0, index=self.states, columns=self.alph)

        for p in self.P.keys():
            self.P.loc[:, p] += np.random.rand(1, n_states)[0]

        for q in self.Q.keys():
            self.Q.loc[:, q] += np.random.rand(1, n_states)[0]

        P_hat = self._normalize(self.P)
        Q_hat = self._normalize(self.Q)

        self.P = P_hat
        self.Q = Q_hat

        for i in range(n_iter):
            # forward probabilities for each sequence
            fwd = list(self._forwards(seqs))
            # backward probabilities for each sequence
            bwd = list(self._backwards(seqs))
            for k in range(n_states):
                state_k = self.states[k]
                for j in range(n_seq):
                    for l in range(n_states):
                        obs_prob = self.obs_prob(seqs[j])
                        # prevent division by 0
                        if obs_prob == 0:
                            obs_prob = 1
                        P_hat[state_k][l] += sum(
                            [fwd[j][k, i] * self.P[state_k][l] * self.Q[seqs[j][i]][l] * bwd[j][
                                l, i + 1] for i in range(len(seqs[j]) - 1)]) / obs_prob
                    for b in self.alph:
                        obs_prob = self.obs_prob(seqs[j])
                        # prevent division by 0
                        if obs_prob == 0:
                            obs_prob = 1
                        Q_hat[b][k] = sum([fwd[j][k, i] * bwd[j][l, i] for i in range(len(seqs[j])) if
                                           seqs[j][i] == b]) / obs_prob

        self.P = self._normalize(P_hat)
        self.Q = self._normalize(Q_hat)

    def train(self, seqs: List[str], method="baum-welch") -> str:
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
                raise ValueError("Viterbi training cannot be used if the state path is unknown. "
                                 "Use 'baum-welch' instead.")
        # Baum-Welch is used when the state path is unknown.
        elif method == 'baum-welch' or method == "baum_welch":
            if self.P is not None and self.Q is not None:
                warn("Baum-Welch assumes that the state path is unknown, but a profile was given.\n"
                     "The profile will be ignored. Use 'viterbi' if you want to use the profile.")
            self.baum_welch(seqs)
        else:
            raise ValueError("Parameter 'method' must be either 'viterbi' or 'baum-welch'.")

        new_seqs: str = self.viterbi_training(seqs)

        return new_seqs

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
    def _normalize(m, cols=False, pseudo=True):
        """
        Normalize a matrix such that each row sums up to 1.
        @:param cols: normalize the columns instead
        @:param pseudo: add pseudocounts
        """
        if isinstance(m, DataFrame):
            if pseudo:
                return m.div(m.sum(1) + 1e-10, axis=0) + 1e-100
            else:
                return m.div(m.sum(1) + 1e-10, axis=0)
        else:
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
