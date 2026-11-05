
import csv
import random
import subprocess
from pathlib import Path

# Dossiers
dossier = Path(__file__).resolve().parent
sortie = dossier / "donnees_v2"

# Ne pas écraser une génération précédente
if sortie.exists():
    raise FileExistsError(
        "Le dossier donnees_v1 existe déjà."
    )

# Nombre de simulations par classe
repetitions = {
    "train": 10,
    "test": 3
}

# Paramètres variables
longueurs_possibles = [800, 1000, 1200, 1500]
longueur_min_segment = 150

# Topologies pour 4 espèces
topologies = {
    "A": (("A", "B"), ("C", "D")),
    "B": (("A", "C"), ("B", "D")),
    "C": (("A", "D"), ("B", "C"))
}


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


def generer_longueurs(total, k, rng):
    """Générer K segments dont la somme vaut total."""

    minimum = longueur_min_segment
    reste = total - k * minimum

    if reste < 0:
        raise ValueError("Alignement trop court pour K")

    # Répartir aléatoirement la longueur restante
    poids = [rng.random() + 0.01 for _ in range(k)]
    somme_poids = sum(poids)

    longueurs = [
        minimum + int(reste * p / somme_poids)
        for p in poids
    ]

    # Corriger les arrondis
    longueurs[-1] += total - sum(longueurs)

    rng.shuffle(longueurs)

    return longueurs


def generer_topologies(k, rng):
    """Éviter deux topologies identiques consécutives."""

    choix = []

    for i in range(k):
        disponibles = [
            t for t in topologies
            if not choix or t != choix[-1]
        ]

        choix.append(rng.choice(disponibles))

    return choix


def creer_newick(topologie, branche):
    groupes = topologies[topologie]

    (a, b), (c, d) = groupes

    return (
        f"(({a}:{branche:.5f},{b}:{branche:.5f})"
        f":{branche:.5f},"
        f"({c}:{branche:.5f},{d}:{branche:.5f})"
        f":{branche:.5f});"
    )


resultats = []

# Graines distinctes pour entraînement et test
graines_groupes = {
    "train": 20261103,
    "test": 20271103
}

for groupe, nb_repetitions in repetitions.items():

    rng = random.Random(graines_groupes[groupe])

    for k in [1, 2, 3, 4]:

        for rep in range(1, nb_repetitions + 1):

            identifiant = f"{groupe}_K{k}_rep{rep:03d}"
            dossier_sim = sortie / groupe / identifiant
            dossier_sim.mkdir(parents=True)

            # Tirage des paramètres
            longueur_totale = rng.choice(longueurs_possibles)
            longueurs = generer_longueurs(
                longueur_totale, k, rng
            )
            arbres_choisis = generer_topologies(k, rng)

            branche = rng.uniform(0.05, 0.20)
            kappa = rng.uniform(1.5, 4.0)
            alpha = rng.uniform(0.5, 2.0)

            sequences_finales = None
            position = 0
            frontieres = []

            for i, (longueur, arbre) in enumerate(
                zip(longueurs, arbres_choisis), start=1
            ):

                # Arbre au format Newick
                nom_arbre = f"arbre_{i}.dnd"

                (dossier_sim / nom_arbre).write_text(
                    creer_newick(arbre, branche) + "\n"
                )

                # Paramètres bpp-seqgen
                parametres = f"""alphabet = DNA
input.tree.file = {nom_arbre}
input.tree.format = Newick
model = HKY85(kappa={kappa:.4f}, theta=0.5, theta1=0.5, theta2=0.5)
rate_distribution = Gamma(n=4, alpha={alpha:.4f})
number_of_sites = {longueur}
output.sequence.file = segment_{i}.fasta
output.sequence.format = Fasta
"""

                fichier_param = f"segment_{i}.bpp"

                (dossier_sim / fichier_param).write_text(
                    parametres
                )

                graine = rng.randrange(1, 2**31)

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

                if i < k:
                    frontieres.append(position)

            assert position == longueur_totale

            # Créer l'alignement PHYLIP
            chemin_phy = dossier_sim / f"{identifiant}.phy"

            with open(chemin_phy, "w") as f:
                f.write(f"4 {longueur_totale}\n")

                for espece in sorted(sequences_finales):
                    sequence = sequences_finales[espece]
                    assert len(sequence) == longueur_totale

                    f.write(f"{espece:<10}{sequence}\n")

            resultats.append({
                "alignement": identifiant,
                "groupe": groupe,
                "K_reel": k,
                "longueur": longueur_totale,
                "frontieres": ";".join(map(str, frontieres)),
                "topologies": ";".join(arbres_choisis),
                "branche": round(branche, 5),
                "kappa": round(kappa, 4),
                "alpha": round(alpha, 4),
                "fichier": str(chemin_phy.relative_to(sortie))
            })

            print(
                f"{identifiant} : OK"
                f" | longueur = {longueur_totale}"
                f" | frontieres = {frontieres}"
            )

# Enregistrer les métadonnées
with open(sortie / "labels.csv", "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "alignement", "groupe", "K_reel",
            "longueur", "frontieres",
            "topologies", "branche",
            "kappa", "alpha", "fichier"
        ]
    )

    writer.writeheader()
    writer.writerows(resultats)

print("\nSimulations terminées :", len(resultats))
print("Entraînement :", sum(r["groupe"] == "train" for r in resultats))
print("Test :", sum(r["groupe"] == "test" for r in resultats))
print("Dossier :", sortie)
