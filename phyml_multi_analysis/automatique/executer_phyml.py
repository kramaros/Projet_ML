
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Chemins du projet
racine = Path(__file__).resolve().parents[2]
dossier_sortie = Path(__file__).resolve().parent

dossier_simulations = (
    racine / "simulations" / "automatique"
)

logiciel = Path.home() / "Projet_10_ML/phyml_multi/phyml_multi"

# Lire les alignements simulés
with open(dossier_simulations / "alignements/labels.csv") as f:
    alignements = list(csv.DictReader(f))

# Premier test : 2 alignements
# Pour tout analyser : python3 executer_phyml.py all
if len(sys.argv) > 1 and sys.argv[1] == "all":
    alignements_a_traiter = alignements
else:
    alignements_a_traiter = alignements[:2]

resultats = []

for ligne in alignements_a_traiter:

    identifiant = ligne["alignement"]
    fichier_source = dossier_simulations / ligne["fichier"]

    for candidat in [1, 2, 3]:

        dossier = (
            dossier_sortie / identifiant / f"candidat_{candidat}"
        )
        dossier.mkdir(parents=True, exist_ok=True)

        fichier_phy = dossier / fichier_source.name
        journal = dossier / "phyml_run.log"

        # Conserver une analyse déjà terminée
        deja_termine = (
            journal.exists()
            and "Final Log-likelihood" in journal.read_text()
        )

        if not deja_termine:
            shutil.copy2(fichier_source, fichier_phy)

            commande = [
                str(logiciel),
                fichier_phy.name,
                "0", "s", "1", "0",
                "HKY", "e", "e", "4", "e",
                "BIONJ", "y", "y", "y",
                str(candidat)
            ]

            with open(journal, "w") as f:
                execution = subprocess.run(
                    commande,
                    cwd=dossier,
                    stdout=f,
                    stderr=subprocess.STDOUT
                )

            if execution.returncode != 0:
                print("ERREUR :", identifiant, candidat)
                print("Consulter :", journal)
                sys.exit(1)

        texte = journal.read_text()

        valeurs_lk = re.findall(
            r"Final Log-likelihood of the multiTrees model\s*:\s*([-\d.]+)",
            texte
        )

        if not valeurs_lk:
            print("Log-vraisemblance manquante :", journal)
            sys.exit(1)

        autocorrelations = re.findall(
            r"Autocorrelation parameter before optimization\s*:\s*"
            r"[\d.]+\s+and then after\s*:\s*([\d.]+)",
            texte
        )

        log_lk = float(valeurs_lk[-1])
        autocorr = (
            float(autocorrelations[-1])
            if autocorrelations else ""
        )

        resultats.append({
            "alignement": identifiant,
            "arbres_candidats": candidat,
            "log_likelihood": log_lk,
            "autocorrelation": autocorr
        })

        print(
            identifiant,
            "| candidat :", candidat,
            "| log-lk :", log_lk
        )

# Créer le tableau des résultats
fichier_csv = dossier_sortie / "resultats_phyml.csv"

with open(fichier_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "alignement",
            "arbres_candidats",
            "log_likelihood",
            "autocorrelation"
        ]
    )

    writer.writeheader()
    writer.writerows(resultats)

print("\nAnalyses terminées :", len(resultats))
print("Résultats enregistrés dans :", fichier_csv)
