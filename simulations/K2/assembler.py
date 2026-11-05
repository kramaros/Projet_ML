
def lire_fasta(fichier):
    sequences = {}
    nom = None

    with open(fichier, "r") as f:
        for ligne in f:
            ligne = ligne.strip()

            if ligne.startswith(">"):
                nom = ligne[1:].split()[0]
                sequences[nom] = ""
            elif ligne and nom is not None:
                sequences[nom] += ligne

    return sequences


segment_A = lire_fasta("segment_A.fasta")
segment_B = lire_fasta("segment_B.fasta")

# Vérifier que les deux segments ont les mêmes espèces
assert set(segment_A) == set(segment_B)
assert len(segment_A) == 4
assert all(len(seq) == 600 for seq in segment_A.values())
assert all(len(seq) == 400 for seq in segment_B.values())

# Assembler les deux segments
alignement = {}

for espece in segment_A:
    alignement[espece] = segment_A[espece] + segment_B[espece]

# Enregistrer au format PHYLIP
with open("simulation_K2.phy", "w") as f:
    f.write("4 1000\n")

    for espece, sequence in alignement.items():
        f.write(f"{espece:<10}{sequence}\n")

print("Assemblage terminé !")
print("Nombre d'espèces :", len(alignement))
print("Longueur des séquences :", len(next(iter(alignement.values()))))
print("K réel : 2")
print("Point de changement : entre 600 et 601")
