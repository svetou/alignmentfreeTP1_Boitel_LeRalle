# Alignment free - TP 1
**Svetlana Le Ralle - Léo Boitel**

## Matrice des distances de Jaccard

Pour $k=21$ :

|                   |  GCA_000069965.1  |  GCA_000005845.2  |  GCA_000013265.1  |  GCA_030271835.1  |  GCA_000008865.2  |  
|  ---------------  |  ---------------  |  ---------------  |  ---------------  |  ---------------  |  ---------------  |  
|  GCA_000069965.1  |  1.0000000000000  |  0.0025674665312  |  0.0024370151490  |  0.0311339446756  |  0.0023137873667  | 
|  GCA_000005845.2  |  0.0025674665312  |  1.0000000000000  |  0.3410085892814  |  0.0025766983532  |  0.4364844461428  | 
|  GCA_000013265.1  |  0.0024370151490  |  0.3410085892814  |  1.0000000000000  |  0.0024341020865  |  0.3070500187181  | 
|  GCA_030271835.1  |  0.0311339446756  |  0.0025766983532  |  0.0024341020865  |  1.0000000000000  |  0.0023181149305  | 
|  GCA_000008865.2  |  0.0023137873667  |  0.4364844461428  |  0.3070500187181  |  0.0023181149305  |  1.0000000000000  | 

## Commentaire sur la matrice des distances

On voit dans la matrice que les séquences GCA_000005845.2, GCA_000013265.1 et GCA_000008865.2 sont similaires deux à deux (30 à 40% de k-mers en commun). En revanche, les séquences GCA_000069965.1 et GCA_030271835.1 sont chacunes très isolées (moins de 1% de k-mers en commun avec les autres, 3% entre elles).

Pour tenter d'expliquer ces résultats, nous avons cherché ces séquences dans l'European Nucleotide Archive. Les trois séquences proches proviennent toutes de souches d'E.Coli, ce qui explique leur forte similarité. Les séquences plus isolées viennent respectivement de Proteus appendicitidis et Proteus mirabilis - elles n'ont donc rien à voir avec E. Coli - et sont toutes deux dans le taxon Proteus, d'où leur légère similarité de 3% de k-mers.

## Description des méthodes
### Objectif
On cherche à calculer une mesure de similarité entre séquences. Pour ce faire, on compare leur distribution de k-mers (sous-mots de longueur k). Cela a l'avantage d'être une opération locale, contrairement à une comparaison par alignement qui nécessiterait de travailler avec toute la séquence d'un coup - c'est donc plus rapide (linéaire en la longueur de la séquence avec une petite constante), et potentiellement parallélisable.

### Encodage des k-mers
Pour manipuler les k-mers rapidement, on les encode comme des entiers de 64 bits (puisqu'il y a 4 bases, un nucléotide est encodé sur 2 bits, on travaille donc avec $k<=32$). Cela nous permet par exemple de les comparer en un cycle de CPU au lieu de $k$.

Quand on encode les k-mers d'une longue séquence, on lit la séquence dans l'ordre et on peut générer l'encodage d'un k-mer à partir du k-mer précédent en temps constant : on réutilise les $k-1$ nucléotides qu'ils ont en commun qui sont déjà encodés, et on se contente de les décaler de 2 bits à gauche pour ajouter l'encodage du nouveau nucléotide qu'on lit. On utilise un masque pour éviter de "déborder" de notre espace d'encodage de $2k$ bits.

### K-mers canoniques
On cherche à mesurer la similarité biologique entre nos deux séquences - nous voulons donc traiter de la même manière un k-mer et son reverse complément. Nous ne pouvons pas déterminer le sens de chaque brin d'ADN une fois pour toute, car il y a pu avoir des inversions au cours de l'évolution, et nous pouvons donc avoir deux séquences très proches dont le sens relatif varie.

On encode donc en parallèle chaque k-mer avec son reverse complément, et on utilise le minimum des deux entiers résultants pour représenter le k-mer. Il sera donc encodé de la même manière que son reverse complément le serait. 

On peut encoder rapidement le reverse complément des k-mers d'une longue séquence : on passe d'un k-mer à l'autre en temps constant en décalant les $k-1$ bits de l'encodage du dernier k-mer, à droite cette fois, et on complète les deux bits à l'index $2k$ et $2k-1$ avec le nucléotide qu'on vient de lire.

### Fonctions sur les k-mers
Les fonctions qui manipulent les kmers sont implémentées dans `kmers.py`.

`encode_kmer_forward_sense` encode un k-mer seul dans le sens positif selon l'encodage décrit ci-dessus. 

`encode_kmer_reverse_complement` encode son reverse complément.

`encode_kmer` encode un k-mer vers sa représentation canonique (en prenant le min des deux fonctions précédentes).

`stream_kmers` streame les k-mers encodés d'une séquence, en réutilisant l'encodage du k-mer précédent selon l'optimisation décrite ci-dessus pour générer chaque encodage en temps constant.

### Distance de Jaccard
La fonction `jaccard` du fichier principal du module calcule la distance de jaccard entre deux multi-ensembles de séquences, A et B. Elle utilise nos fonctions sur les k-mers pour générer des streams efficaces de k-mers encodés canoniquement.

Pour calculer la distance de Jaccard, on commence par compter le nombre d'occurrences de chaque k-mer dans A, en streamant tous les k-mers de A et en mettant à jour un dictionnaire de comptes au fur et à mesure.

On stream ensuite tous les k-mers de B ; à chaque fois que nous en trouvons un dans le dictionnaire, nous réduisons d'un son compte dans le dictionnaire et nous incrémentons une variable qui compte les k-mers en commun (intersection des multi-ensembles).

Finalement, on divise le nombre de kmers dans l'intersection qu'on a comptés par l'union de A et B (la somme de leurs tailles respectives moins leur intersection), et on obtient la distance de Jaccard entre eux.

### Tests
Des tests pour chaque fonction implémentée sont disponibles disponibles (`python tests_TP.py`). Ils testent des résultats calculés à la main et la cohérence des sorties d'une fonction à l'autre.

## Seconde partie du TP
Dans cette seconde partie, nous voulons réduire le temps de calcul ainsi que l'espace de stockage utilisé par notre code implémenté la dernière fois, dans le but d'une utilisation sur des séquences beaucoup plus longues. Nous utilisons donc des itérateurs pour manipuler les séquences d'origine (pour ne jamais les avoir en entier en mémoire d'un coup), et nous en dérivons des objets plus petits - des sketches.

### Création de sketches
Nous crééons donc des sketches de nos kmers, c'est à dire des listes d'un certain nombre de plus petits kmers.
Pour cela, nous avons implémenté la fonction 'filter_smallests' qui filtre les kmers et ne garde que les s plus petits kmers parmis ceux générés par la fonction stream_kmers.
Dans la fonction 'filter_smallests', nous utilisons en premier lieu la fonction 'hash_kmers' qui réarrange l'ordre des kmers aléatoirement pour éviter les biais biologiques des minimums, par exemple le biai des queues poly A. 

Après le hashing, nous utilisons la la library heapq, une library standard basée sur les arbres binaires et qui réduit considérablement la complexité de l'algorithme. 
Elle nous permet de passer d'une complixté de x à log(x) en utilisant une heapq au lieu d'une liste.

Nous créons tout d'abord une heap vide, puis, on boucle sur les kmers. Si la longueur de notre heap est inférieure à notre taille de sketches souhaitée, alors on ajoute notre élément à notre heap.
Si le heap est déjà de la longueur souhaité de notre sketch, alors obtient le plus grand élément actuel dans le heap (en prenant le négatif du premier élément du heap, car les éléments sont stockés en valeurs négatives). 
Si le kmer actuel est plus petit que le plus grand élément du heap, alors on le remplace.
On remplace le kmer maximal actuel par le nouveau kmer maximal (encore une fois en utilisant la valeur négative pour maintenir la structure max-heap).
Enfin, on retourne la liste des éléments du heap, inversant les signes pour revenir aux valeurs originales des k-mers, puis les trie pour renvoyer une liste ordonnée des s plus petits kmers.

### Calcul de la distance de jaccard, nouvelle version 
Dans cette nouvelle version, nous utilisons deux pointeurs aux débuts des deux listes ordonnées, nous prenons le plus petit kmer.
Deux autres compteurs, count_union et count_intersection, sont initialisés à zéro pour suivre la taille de l'union et de l'intersection.

Dans la boucle principale,  les deux sketches en parallèle sont explorés.
Si l'élément actuel de sketch1 est le meme que de sketch2, alors on incrémente à la fois l'intersection et l'union, puis on passe au kmer suivant dans les deux sketches.
Si l'élément actuel de sketch1 est plus petit que celui de sketch2, alors cet élément n'existe que dans sketch1, et il est ajouté à l'union (mais pas à l'intersection). Le compteur i est incrémenté pour passer à l'élément suivant de sketch1.
Si l'élément actuel de sketch2 est plus petit, la même logique s'applique, mais pour sketch2.

A la fin de la boucle, la fonction ajoute les éléments restants de l'autre sketches à l'union, car ils n'ont pas été comparés et sont donc uniques à cette esquisse.
Enfin, la fonction retourne le ratio de l'intersection par rapport à l'union, qui correspond à la similarité de jaccard. 
