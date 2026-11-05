# Premier test - PhyML-Multi

## Objectif
Tester PhyML-Multi sur un alignement simulé contenant
un événement de recombinaison connu.

## Données
- Alignement : simSeqs.phy
- Point de recombinaison attendu : position 599
- Modèle de substitution : HKY
- Nombre de composantes phylogénétiques : 2
- Méthode : phylo-HMM

## Commande utilisée
../phyml_multi simSeqs.phy 0 i 1 0 HKY e e 4 e BIONJ y y y 2

## Résultats
- Log-vraisemblance finale : -19811.596047
- Paramètre d'autocorrélation : 0.997973
- Temps d'exécution : 32 secondes
- Optimisation convergente
- Exécution terminée sans erreur

## Partitionnement
Le point de recombinaison attendu est à la position 599.
Sa détection avec PartitioningHMM.py reste à vérifier.

## Logiciels
- PhyML-Multi : https://github.com/lgueguen/phyml_multi
- SARMENT : https://github.com/lgueguen/SARMENT

## Remarque
Les résultats ont été obtenus avec le programme PhyML-Multi
compilé sous Ubuntu WSL avec GCC (norme gnu17).
