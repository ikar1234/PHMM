import numpy as np
from typing import List
from warnings import warn

from PHMM.src.Records import Record


class PHMM():
    record: Record
    # transition probability matrix
    P: np.ndarray
    # emission probability matrix
    Q: np.ndarray

    def __init__(self, record=None):
        """
        Construct the transition and emission probability matrices,
        which are the main parameters of the model.
        :param record: record with the position-weight matrix
        """
        self.record = record
        self.P = record.matrix  # TODO:?
        self.Q = record.matrix

    def evaluate(self, seq: str) -> float:
        """
        Given an observation sequence and a model, find the likelihood of the sequence with respect to the model.
        This value can be estimated from the profile itself.
        :param seq: sequence
        :return: log-likelihood of the sequence with respect to the model
        """
        m = self.record.matrix
        # protein alphabet, characters in the same order as the PWM
        al = "A C D E F G H I K L M N P Q R S T V W Y".split(' ')
        return np.sum([np.log(m[al.index(seq[i]), i]) for i in range(len(seq))])[0]

    def forward(self, q: str, i: int) -> float:
        """

        :param q: state (one of I,M,D)
        :param i:
        :return:
        """
        if self.P is None:
            raise ValueError("State path must be known.")

    def backward(self, q: str, i: int) -> float:
        """

        :param q: state (one of I,M,D)
        :param i:
        :return:
        """

        if self.P is None:
            raise ValueError("State path must be known.")

    def viterbi_decoding(self, seq: str):
        # TODO: use log?
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
        # we have 3 states - insert,match,delete
        trellis = np.zeros((3, len(seq)))
        # in order to obtain the sequence itself
        backpointers = np.zeros((3, len(seq)))
        # initialize
        for t in range(3):
            trellis[0][t] = ...
        # since we have only 3 states, the run-time os O(n)
        for n in range(1, len(seq)):
            for t in range(3):
                # TODO: log prob
                last_column = [trellis[n - 1][i] * self.P[i][t] for i in range(3)]
                trellis[n][t] = max(last_column)
                # TODO: check
                backpointers[n][t] = np.argmax(last_column)
                # multiply by emission probabilities
                trellis[n][t] *= self.Q[n][t]

        return backpointers

    def viterbi_training(self, seqs: List[str]):
        """
        Viterbi training.
        :param seqs:
        :return:
        """
        # paths of states (match, insert, delete) for every sequence
        paths = [self.viterbi(seq) for seq in seqs]
        for p in range(len(seqs)):
            yield

    def train(self, seqs: List[str], method="baum_welch", n_iter=1000) -> str:
        """
        Produces a multiple sequence alignment from a set
        of sequences.
        :param method: Training method. Most common is the
        Baum-Welch algorithm.
        :param n_iter: number of iterations for the Baum-Welch algorithm
        :return: multiple sequence alignment
        """
        # Viterbi-Training - uses the Viterbi algorithm for every input sequence and
        # aligns the sequence using the generated observed states.
        if method == 'viterbi':
            if self.P is None and self.Q is None:
                warn("Viterbi training cannot be used if the state path is unknown.\n"
                     "Use 'baum-welch' instead.")
                return ""
            return self.viterbi_training(seqs)
            # multiple sequence alignment
            msa = '\n'.join(path_seqs)

            return msa
        # Baum-Welch is used when the state path is unknown.
        elif method == 'baum_welch':
            if self.P is not None and self.Q is not None:
                warn("Baum-Welch assumes that the state path is unknown, but a profile was given.\n"
                     "The profile will be ignored. Use 'viterbi' if you want to use the profile.")

            # we initialize P and Q to be random Markov matrices
            self.P = np.random.rand(self.P.shape[0], self.P.shape[1])
            self.Q = np.random.rand(self.Q.shape[0], self.P.shape[1])
            # normalize (divide by the column sums and transpose)
            self.P = (self.P / np.sum(self.P, axis=0)).T
            self.Q = (self.Q / np.sum(self.Q, axis=0)).T

            for n in range(n_iter):
                alpha = ...  # self.forward(V, a, b)
                beta = ...  # self.backward(V, a, b)

                xi = np.zeros((M, M, T - 1))
                for t in range(T - 1):
                    denominator = np.dot(np.dot(alpha[t, :].T, a) * b[:, V[t + 1]].T, beta[t + 1, :])
                    for i in range(M):
                        numerator = alpha[t, i] * a[i, :] * b[:, V[t + 1]].T * beta[t + 1, :].T
                        xi[i, :, t] = numerator / denominator

                gamma = np.sum(xi, axis=1)
                a = np.sum(xi, 2) / np.sum(gamma, axis=1).reshape((-1, 1))

                # Add additional T'th element in gamma
                gamma = np.hstack((gamma, np.sum(xi[:, :, T - 2], axis=0).reshape((-1, 1))))

                K = b.shape[1]
                denominator = np.sum(gamma, axis=1)
                for l in range(K):
                    b[:, l] = np.sum(gamma[:, V == l], axis=1)

                b = np.divide(b, denominator.reshape((-1, 1)))

            ...
        else:
            raise ValueError("Parameter 'method' must be either 'viterbi' or 'baum_welch'.")
