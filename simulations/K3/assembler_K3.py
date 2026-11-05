
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


# Lire les trois segments
A = lire_fasta("segment_A.fasta")
B = lire_fasta("segment_B.fasta")
C = lire_fasta("segment_C.fasta")

# Vérifier les espèces et les longueurs
assert len(A) == 4
assert set(A) == set(B) == set(C)
assert all(len(seq) == 350 for seq in A.values())
assert all(len(seq) == 350 for seq in B.values())
assert all(len(seq) == 300 for seq in C.values())

# Assembler les séquences des mêmes espèces
alignement = {}

for espece in A:
    alignement[espece] = A[espece] + B[espece] + C[espece]

# Enregistrer au format PHYLIP
with open("simulation_K3.phy", "w") as f:
    f.write("4 1000\n")

    for espece, sequence in alignement.items():
        f.write(f"{espece:<10}{sequence}\n")

print("Assemblage terminé !")
print("Nombre d'espèces :", len(alignement))
print("Longueur :", len(next(iter(alignement.values()))))
print("K réel : 3")
print("Points de changement : 350/351 et 700/701")

