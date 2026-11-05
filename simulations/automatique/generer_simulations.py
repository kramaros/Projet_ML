
import subprocess
from pathlib import Path
import csv

# Dossier de travail
dossier = Path(__file__).resolve().parent
sortie = dossier / "alignements"
sortie.mkdir(exist_ok=True)

# Nos trois topologies
arbres = {
    "A": "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
    "B": "((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);",
    "C": "((A:0.1,D:0.1):0.1,(B:0.1,C:0.1):0.1);"
}

# Longueurs et arbres des segments
configurations = {
    1: [(1000, "A")],
    2: [(600, "A"), (400, "B")],
    3: [(350, "A"), (350, "B"), (300, "C")],
    4: [(250, "A"), (250, "B"),
        (250, "C"), (250, "A")]
}

nombre_repetitions = 2


def lire_fasta(chemin):
    sequences = {}
    nom = None

    with open(chemin) as f:
        for ligne in f:
            ligne = ligne.strip()

            if ligne.startswith(">"):
                nom = ligne[1:].split()[0]
                sequences[nom] = ""
            elif ligne:
                sequences[nom] += ligne

    return sequences


# Éviter d'écraser des simulations existantes
for k in configurations:
    for rep in range(1, nombre_repetitions + 1):
        identifiant = f"K{k}_rep{rep:03d}"
        if (sortie / identifiant).exists():
            raise FileExistsError(
                f"Simulation déjà présente : {identifiant}"
            )

resultats = []

for k, segments in configurations.items():

    for rep in range(1, nombre_repetitions + 1):

        identifiant = f"K{k}_rep{rep:03d}"
        dossier_sim = sortie / identifiant
        dossier_sim.mkdir()

        sequences_finales = None
        position = 0
        frontieres = []

        for i, (longueur, arbre) in enumerate(segments, 1):

            # Écrire l'arbre
            (dossier_sim / f"arbre_{i}.dnd").write_text(
                arbres[arbre] + "\n"
            )

            # Paramètres de bpp-seqgen 2.4
            parametres = f"""alphabet = DNA
input.tree.file = arbre_{i}.dnd
input.tree.format = Newick
model = HKY85(kappa=2, theta=0.5, theta1=0.5, theta2=0.5)
rate_distribution = Gamma(n=4, alpha=0.5)
number_of_sites = {longueur}
output.sequence.file = segment_{i}.fasta
output.sequence.format = Fasta
"""

            fichier_param = f"segment_{i}.bpp"
            (dossier_sim / fichier_param).write_text(parametres)

            graine = 10000 + k * 1000 + rep * 100 + i

            execution = subprocess.run(
                [
                    "bppseqgen",
                    f"param={fichier_param}",
                    f"--seed={graine}"
                ],
                cwd=dossier_sim,
                capture_output=True,
                text=True,
                check=True
            )

            (dossier_sim / f"segment_{i}.log").write_text(
                execution.stdout + execution.stderr
            )

            sequences = lire_fasta(
                dossier_sim / f"segment_{i}.fasta"
            )

            assert set(sequences) == {"A", "B", "C", "D"}
            assert all(
                len(seq) == longueur
                for seq in sequences.values()
            )

            if sequences_finales is None:
                sequences_finales = sequences.copy()
            else:
                for espece in sequences_finales:
                    sequences_finales[espece] += sequences[espece]

            position += longueur

            if i < len(segments):
                frontieres.append(position)

        # Écrire l'alignement final PHYLIP
        assert position == 1000

        chemin_phy = dossier_sim / f"{identifiant}.phy"

        with open(chemin_phy, "w") as f:
            f.write("4 1000\n")

            for espece in sorted(sequences_finales):
                sequence = sequences_finales[espece]
                assert len(sequence) == 1000
                f.write(f"{espece:<10}{sequence}\n")

        resultats.append({
            "alignement": identifiant,
            "K_reel": k,
            "longueur": 1000,
            "frontieres": ";".join(map(str, frontieres)),
            "fichier": str(chemin_phy.relative_to(dossier))
        })

        print(f"{identifiant} : OK")

# Enregistrer les vraies étiquettes
with open(sortie / "labels.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "alignement", "K_reel", "longueur",
            "frontieres", "fichier"
        ]
    )

    writer.writeheader()
    writer.writerows(resultats)

print("\nSimulations terminées :", len(resultats))
