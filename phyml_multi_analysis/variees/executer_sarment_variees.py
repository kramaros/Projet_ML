
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

racine = Path(__file__).resolve().parents[2]

dossier = racine / "phyml_multi_analysis/variees/donnees_v1"
fichier_phyml = dossier / "resultats_phyml.csv"

sarment = Path.home() / "Projet_10_ML/SARMENT"
script_hmm = Path.home() / "Projet_10_ML/phyml_multi/PartitioningHMM.py"

# Lire les résultats PhyML-Multi
with open(fichier_phyml, newline="") as f:
    analyses = list(csv.DictReader(f))

# Premier test sur un alignement, puis tous les alignements
selection = sys.argv[1] if len(sys.argv) > 1 else "train_K2_rep001"

if selection != "all":
    analyses = [
        ligne for ligne in analyses
        if ligne["alignement"] == selection
    ]

    if not analyses:
        raise ValueError("Alignement introuvable : " + selection)

resultats = []
echecs = []

for analyse in analyses:

    identifiant = analyse["alignement"]
    groupe = analyse["groupe"]
    candidat = int(analyse["arbres_candidats"])

    if candidat == 1:
        continue

    chemin = dossier / groupe / identifiant / f"candidat_{candidat}"

    fichier_lk = chemin / f"{identifiant}.phy_phyml_siteLks.txt"
    journal = chemin / "partition_HMM.log"

    autocorr = analyse["autocorrelation"]

    if not fichier_lk.exists() or not autocorr:
        print("Données manquantes :", identifiant, candidat)
        echecs.append((identifiant, candidat))
        continue

    # Réutiliser les résultats déjà complets
    termine = (
        journal.exists()
        and "Viterbi :" in journal.read_text()
        and "Forward-Backward :" in journal.read_text()
    )

    if not termine:

        environnement = os.environ.copy()
        environnement["PYTHONPATH"] = str(sarment)

        commande = [
            "python3",
            str(script_hmm),
            fichier_lk.name,
            autocorr
        ]

        with open(journal, "w") as f:
            execution = subprocess.run(
                commande,
                cwd=chemin,
                env=environnement,
                stdout=f,
                stderr=subprocess.STDOUT
            )

        if execution.returncode != 0:
            print("Erreur SARMENT :", identifiant, candidat)
            echecs.append((identifiant, candidat))
            continue

    texte = journal.read_text()

    ligne_resultat = {
        "alignement": identifiant,
        "groupe": groupe,
        "arbres_candidats": candidat
    }

    try:
        for methode, prefixe in [
            ("Viterbi", "viterbi"),
            ("Forward-Backward", "fb")
        ]:

            bloc = re.search(
                re.escape(methode) + r"\s*:\s*(.*?)--->",
                texte,
                re.DOTALL
            )

            if bloc is None:
                raise ValueError("Partition introuvable")

            segments = re.findall(
                r"<(\d+)-(\d+)>(\d+)",
                bloc.group(1)
            )

            if not segments:
                raise ValueError("Aucun segment trouvé")

            longueurs = [
                int(fin) - int(debut) + 1
                for debut, fin, etat in segments
            ]

            frontieres = [
                int(fin) + 1
                for debut, fin, etat in segments[:-1]
            ]

            ligne_resultat[f"nb_segments_{prefixe}"] = len(segments)

            ligne_resultat[f"longueur_moy_{prefixe}"] = (
                sum(longueurs) / len(longueurs)
            )

            ligne_resultat[f"longueur_min_{prefixe}"] = min(longueurs)

            ligne_resultat[f"frontieres_{prefixe}"] = ";".join(
                map(str, frontieres)
            )

        resultats.append(ligne_resultat)

        print(
            identifiant,
            "| candidat :", candidat,
            "| Viterbi :", ligne_resultat["nb_segments_viterbi"],
            "| Forward-Backward :", ligne_resultat["nb_segments_fb"],
            flush=True
        )

    except ValueError as erreur:
        print("Erreur lecture :", identifiant, candidat, erreur)
        echecs.append((identifiant, candidat))

# Enregistrer les résultats
fichier_sortie = dossier / "resultats_partitionnement.csv"

colonnes = [
    "alignement",
    "groupe",
    "arbres_candidats"
]

for prefixe in ["viterbi", "fb"]:
    colonnes += [
        f"nb_segments_{prefixe}",
        f"longueur_moy_{prefixe}",
        f"longueur_min_{prefixe}",
        f"frontieres_{prefixe}"
    ]

with open(fichier_sortie, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=colonnes)
    writer.writeheader()
    writer.writerows(resultats)

print("\nPartitionnements réussis :", len(resultats))
print("Échecs :", len(echecs))
print("Fichier :", fichier_sortie)

if echecs:
    print("Analyses à vérifier :", echecs)
