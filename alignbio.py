from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Definir las secuencias
seq2 = "MISIPPLISTRLPVGVARLSKITGPLALPTTGRAHYDVYSCLPITLLHLPAHFQKFSQPAEISHIRYRELLGYSHQRPRLQKGTHSSRQVAALPLVPRSSTLDKYVAFFTAVFFILLVGSFRFLDVAAGTKIPLHLVKSLLLSKIRKPLEVRSSTLFQTFLSANKIIKKGDWKLPYFVFLLLGRIIKGEHPPLMGLRAAFLAWHFH"
seq1 = "MDTGHTVSVDHPKVTSHKNHPKQQLLLHDIHPTTYFY"

# Realizar el alineamiento global
alignments = pairwise2.align.globalxx(seq1, seq2)

# Tomar el primer alineamiento (puedes considerar otros alineamientos si es necesario)
best_alignment = alignments[0]

# Calcular el porcentaje de identidad
identical_count = sum(1 for a, b in zip(best_alignment.seqA, best_alignment.seqB) if a == b)
percentage_identity = (identical_count / len(seq1)) * 100

print(f"Porcentaje de similitud: {percentage_identity:.2f}%")
