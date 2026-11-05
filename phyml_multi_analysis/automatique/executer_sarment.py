
import csv
import re
import subprocess
import sys
from pathlib import Path

dossier = Path(__file__).resolve().parent
sarment = Path.home() / "Projet_10_ML/SARMENT"
script = Path.home() / "Projet_10_ML/phyml_multi/PartitioningHMM.py"

resultats = []

with open(dossier / "resultats_phyml.csv") as f:
    analyses = list(csv.DictReader(f))

for analyse in analyses:
    identifiant = analyse["alignement"]
    candidat = int(analyse["arbres_candidats"])

    # Le cas à un arbre est traité séparément
    if candidat == 1:
        continue

    autocorrelation = analyse["autocorrelation"]

    if not autocorrelation:
        print("Autocorrélation manquante :", identifiant, candidat)
        sys.exit(1)

    chemin = dossier / identifiant / f"candidat_{candidat}"
    fichier = chemin / f"{identifiant}.phy_phyml_siteLks.txt"
    journal = chemin / "partition_HMM.log"

    if not fichier.exists():
        print("Fichier manquant :", fichier)
        sys.exit(1)

    # Exécuter SARMENT si le résultat n'existe pas
    if not journal.exists():
        with open(journal, "w") as f:
            execution = subprocess.run(
                [
                    "python3",
                    str(script),
                    fichier.name,
                    autocorrelation
                ],
                cwd=chemin,
                env={
                    **__import__("os").environ,
                    "PYTHONPATH": str(sarment)
                },
                stdout=f,
                stderr=subprocess.STDOUT
            )

        if execution.returncode != 0:
            print("Erreur SARMENT :", journal)
            sys.exit(1)

    texte = journal.read_text()

    ligne_resultat = {
        "alignement": identifiant,
        "arbres_candidats": candidat
    }

    for methode, prefixe in [
        ("Viterbi", "viterbi"),
        ("Forward-Backward", "fb")
    ]:
        correspondance = re.search(
            re.escape(methode) + r"\s*:\s*(.*?)--->",
            texte,
            re.DOTALL
        )

        if not correspondance:
            print("Résultat manquant :", methode, journal)
            sys.exit(1)

        segments = re.findall(
            r"<(\d+)-(\d+)>(\d+)",
            correspondance.group(1)
        )

        if not segments:
            print("Aucun segment :", journal)
            sys.exit(1)

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
        "| Forward-Backward :", ligne_resultat["nb_segments_fb"]
    )

with open(dossier / "resultats_partitionnement.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=resultats[0].keys())
    writer.writeheader()
    writer.writerows(resultats)

print("\nAnalyses terminées :", len(resultats))
