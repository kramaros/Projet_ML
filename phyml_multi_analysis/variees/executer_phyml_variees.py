
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

racine = Path(__file__).resolve().parents[2]
simulations = racine / "simulations/variees/donnees_v1"
sortie = racine / "phyml_multi_analysis/variees/donnees_v1"
logiciel = Path.home() / "Projet_10_ML/phyml_multi/phyml_multi"

with open(simulations / "labels.csv", newline="") as f:
    alignements = list(csv.DictReader(f))

selection = sys.argv[1] if len(sys.argv) > 1 else "train_K2_rep001"

if selection == "all":
    a_traiter = alignements
else:
    a_traiter = [
        ligne for ligne in alignements
        if ligne["alignement"] == selection
    ]

if not a_traiter:
    raise ValueError("Alignement introuvable : " + selection)

sortie.mkdir(parents=True, exist_ok=True)
fichier_csv = sortie / "resultats_phyml.csv"

colonnes = [
    "alignement",
    "groupe",
    "arbres_candidats",
    "log_likelihood",
    "autocorrelation"
]

resultats = []
echecs = []


def analyse_complete(dossier, identifiant):
    journal = dossier / "phyml_run.log"
    fichier_lk = dossier / f"{identifiant}.phy_phyml_siteLks.txt"

    return (
        journal.exists()
        and fichier_lk.exists()
        and fichier_lk.stat().st_size > 0
        and "Final Log-likelihood" in journal.read_text()
    )


def sauvegarder():
    temporaire = fichier_csv.with_suffix(".csv.tmp")

    with open(temporaire, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=colonnes)
        writer.writeheader()
        writer.writerows(resultats)

    temporaire.replace(fichier_csv)


for ligne in a_traiter:

    identifiant = ligne["alignement"]
    groupe = ligne["groupe"]
    fichier_source = simulations / ligne["fichier"]

    for candidat in [1, 2, 3]:

        dossier = sortie / groupe / identifiant / f"candidat_{candidat}"

        if not analyse_complete(dossier, identifiant):

            # Deux tentatives maximum
            for tentative in [1, 2]:

                # Archiver une analyse incomplète
                if dossier.exists():
                    numero = 1
                    while True:
                        archive = dossier.parent / (
                            f"candidat_{candidat}_echec_{numero}"
                        )
                        if not archive.exists():
                            break
                        numero += 1

                    dossier.rename(archive)

                dossier.mkdir(parents=True, exist_ok=True)

                shutil.copy2(
                    fichier_source,
                    dossier / fichier_source.name
                )

                commande = [
                    str(logiciel),
                    fichier_source.name,
                    "0", "s", "1", "0",
                    "HKY", "e", "e", "4", "e",
                    "BIONJ", "y", "y", "y",
                    str(candidat)
                ]

                print(
                    f"{identifiant} | candidat {candidat}"
                    f" | tentative {tentative}",
                    flush=True
                )

                journal = dossier / "phyml_run.log"

                with open(journal, "w") as f:
                    execution = subprocess.run(
                        commande,
                        cwd=dossier,
                        stdout=f,
                        stderr=subprocess.STDOUT
                    )

                if (
                    execution.returncode == 0
                    and analyse_complete(dossier, identifiant)
                ):
                    break

                print(
                    "Échec de cette tentative :",
                    journal,
                    flush=True
                )

        if not analyse_complete(dossier, identifiant):
            echecs.append((identifiant, candidat))
            print("ANALYSE NON TERMINÉE :", identifiant, candidat)
            continue

        texte = (dossier / "phyml_run.log").read_text()

        valeurs = re.findall(
            r"Final Log-likelihood of the multiTrees model\s*:"
            r"\s*([-\d.]+)",
            texte
        )

        autocorrelations = re.findall(
            r"Autocorrelation parameter before optimization\s*:"
            r"\s*[\d.]+\s+and then after\s*:\s*([\d.]+)",
            texte
        )

        if not valeurs:
            echecs.append((identifiant, candidat))
            continue

        loglk = float(valeurs[-1])

        autocorr = (
            float(autocorrelations[-1])
            if autocorrelations else ""
        )

        resultats.append({
            "alignement": identifiant,
            "groupe": groupe,
            "arbres_candidats": candidat,
            "log_likelihood": loglk,
            "autocorrelation": autocorr
        })

        # Sauvegarde après chaque analyse réussie
        sauvegarder()

        print(
            identifiant,
            "| candidat :", candidat,
            "| log-lk :", loglk,
            flush=True
        )

# Rapport des échecs éventuels
fichier_echecs = sortie / "echecs_phyml.csv"

with open(fichier_echecs, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["alignement", "arbres_candidats"])
    writer.writerows(echecs)

print("\nAnalyses réussies :", len(resultats))
print("Analyses échouées :", len(echecs))
print("Fichier :", fichier_csv)

if echecs:
    sys.exit(1)

