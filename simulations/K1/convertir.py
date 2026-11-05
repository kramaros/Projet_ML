
sequences = {}
nom = None

with open("simulation_K1.fasta", "r") as f:
    for ligne in f:
        ligne = ligne.strip()

        if ligne.startswith(">"):
            nom = ligne[1:].split()[0]
            sequences[nom] = ""
        elif ligne:
            sequences[nom] += ligne

# Vérification des séquences
assert len(sequences) == 4
assert all(len(seq) == 1000 for seq in sequences.values())

# Conversion au format PHYLIP
with open("simulation_K1.phy", "w") as f:
    f.write("4 1000\n")

    for nom, sequence in sequences.items():
        f.write(f"{nom:<10}{sequence}\n")

print("Conversion terminée !")
print("Nombre d'espèces :", len(sequences))
print("Longueur :", len(next(iter(sequences.values()))))
print("K réel : 1")
