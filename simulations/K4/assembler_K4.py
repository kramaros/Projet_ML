
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


# Lire les quatre segments
A = lire_fasta("segment_A.fasta")
B = lire_fasta("segment_B.fasta")
C = lire_fasta("segment_C.fasta")
D = lire_fasta("segment_D.fasta")

# Vérifier les espèces et les longueurs
assert len(A) == 4
assert set(A) == set(B) == set(C) == set(D)

for segment in [A, B, C, D]:
    assert all(len(seq) == 250 for seq in segment.values())

# Assembler les quatre segments
alignement = {}

for espece in A:
    alignement[espece] = (
        A[espece] + B[espece] + C[espece] + D[espece]
    )

# Enregistrer au format PHYLIP
with open("simulation_K4.phy", "w") as f:
    f.write("4 1000\n")

    for espece, sequence in alignement.items():
        f.write(f"{espece:<10}{sequence}\n")

print("Assemblage terminé !")
print("Nombre d'espèces :", len(alignement))
print("Longueur :", len(next(iter(alignement.values()))))
print("K réel : 4")
print("Points de changement : 250/251, 500/501 et 750/751")
