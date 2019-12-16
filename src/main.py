# TODO: 1. Evaluation
#   - using the forward and the backward algorithm to estimate the probability of a sequence given a model


# TODO: 2. Decoding Viterbi algorithm - Dynamic programming O(|Q|*l) space and  O(|Q|^2*l) time
#  - backtracing to find most probable path
#    forward-backward = forward 1..i + backward n..i computation of probabilities for a given state i
#    don't forget to use logs for numerical stability/prevent underflow/
#    utils?

# TODO: 3. Training
#  - Given a model structure and a set of sequences, find the model that best fits the data.
#  - 3 possibilities:
#   - MLE (maximum likelihood estimation)
#   - Viterbi training(DO NOT confuse with Viterbi decoding)
#   - Baum Welch =EM; uses forward-backward algorithm


# TODO: parsing data

"""
Wenn wir fur eine Menge von Sequenzen ¨ ~y = (y(1), . . . , y(m)) jetzt allgemein ein mehrfaches Sequenzen-Alignment
berechnen wollen, mussen wir zuerst den Parameter ℓ des Hidden Markov Models festlegen,das ist die Anzahl der
eigentlichen (nichtstillen) Match-Zustande. Diese kann beispielsweise als die mittele Sequenzl¨ange von
~y gew¨ahlt werden. Dann trainieren wir mit den Sequenzen ~y mithilfe des Baum-Welch-Algorithmus die
Zustandsubergangs- und Emissionswahrscheinlichkeiten des Hidden Markov Models. Letztendlich berechnen wir dann wider
den wahrscheinlichsten Pfad jeder Sequenz mithilfe des Viterbi-Algorithmus
(oder des vorhin angepassten Scoring-Algorithmus) durch dieses Hidden Markov Model. Durch die Zustandsfolgen wissen
wir nun fur jedes Zeichen der Sequenz, ob Sie durch ein Match-, Insertion- oder Deletion-Zustand gelaufen sind und
fugen die Sequenzen im Wesentlichen gem¨aß der Match-Zustande zu einem mehrfachen Sequenz-Alignment zusammen.
Zeichen, die zu einem Insert-Zustand geh¨oren, werden in neuen Spalten aufgenommen, so dass die anderen Sequenzen
mit Gaps gefullt werden (außer sie enthalten an dieser Stelle auch Insertionen). Deletionszust¨ande fuhren zu Gaps
in der Sequenz beim Alignieren mit den anderen Sequenzen an den entsprechenden Match-Positionen.
"""