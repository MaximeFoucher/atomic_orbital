# Simulateur d'Orbitales Atomiques — Atome d'Hydrogène

> Ce document accompagne le simulateur Python interactif des orbitales de l'atome d'hydrogène.
> Il reprend de zéro toute la théorie nécessaire pour comprendre ce que le code calcule et visualise.

---

## Table des matières

1. [Contexte physique](#1-contexte-physique)
2. [Les nombres quantiques n, l, m](#2-les-nombres-quantiques-n-l-m)
3. [La fonction d'onde ψ](#3-la-fonction-donde-ψ)
4. [La partie radiale R(r)](#4-la-partie-radiale-rnlr)
5. [Les polynômes de Laguerre associés](#5-les-polynômes-de-laguerre-associés)
6. [La partie angulaire Y(θ, φ) — harmoniques sphériques](#6-la-partie-angulaire-yθ-φ--harmoniques-sphériques)
7. [La densité de probabilité |ψ|²](#7-la-densité-de-probabilité-ψ)
8. [Les niveaux d'énergie et les raies spectrales](#8-les-niveaux-dénergie-et-les-raies-spectrales)
9. [Catalogue des orbitales principales](#9-catalogue-des-orbitales-principales)
10. [Fun fact — L'analogie avec la vibration d'une membrane](#10-fun-fact--lanalogie-avec-la-vibration-dune-membrane)
11. [Utilisation du simulateur](#11-utilisation-du-simulateur)

---

## 1. Contexte physique

L'atome d'hydrogène est le système quantique le plus simple : un proton au centre, un électron en orbite. Il est remarquable parce que c'est **le seul atome dont l'équation de Schrödinger admet une solution analytique exacte** — c'est-à-dire que l'on peut écrire la réponse sous forme de formules fermées, sans approximation.

L'équation de Schrödinger indépendante du temps s'écrit :

```
Ĥ ψ = E ψ
```

où `Ĥ` est l'opérateur hamiltonien (énergie cinétique + énergie potentielle coulombienne) :

```
Ĥ = -ℏ²/(2mₑ) ∇² - e²/(4πε₀ r)
```

En passant en **coordonnées sphériques** (r, θ, φ) — coordonnées naturelles pour un potentiel à symétrie centrale — l'équation se sépare en deux parties indépendantes : une radiale et une angulaire. C'est cette séparation qui donne naissance aux nombres quantiques.

---

## 2. Les nombres quantiques n, l, m

Les nombres quantiques émergent naturellement des conditions aux limites imposées lors de la résolution de l'équation de Schrödinger. Ce ne sont pas des hypothèses arbitraires : **ils apparaissent parce que la fonction d'onde doit être continue, bornée et normalisable**.

### n — Nombre quantique principal

| Propriété | Valeur |
|-----------|--------|
| Valeurs autorisées | n = 1, 2, 3, 4, … |
| Rôle physique | Détermine le **niveau d'énergie** et la **taille** de l'orbitale |
| Nombre de noeuds radiaux | n − l − 1 |

Plus n est grand, plus l'électron est en moyenne loin du noyau, et plus son énergie est élevée (moins négative). Toutes les orbitales de même n forment une **couche électronique** (K, L, M, N, …).

### l — Nombre quantique azimutal (ou orbital)

| Propriété | Valeur |
|-----------|--------|
| Valeurs autorisées | l = 0, 1, 2, …, n−1 &nbsp;&nbsp; (**dépend de n**) |
| Rôle physique | Détermine le **moment cinétique orbital** et la **forme** de l'orbitale |
| Notation spectroscopique | l = 0 → s, l = 1 → p, l = 2 → d, l = 3 → f |

Le moment cinétique orbital vaut `L = ℏ √(l(l+1))`. Pour n = 3, l peut valoir 0, 1 ou 2 — ce qui correspond aux sous-couches 3s, 3p et 3d.

### m — Nombre quantique magnétique

| Propriété | Valeur |
|-----------|--------|
| Valeurs autorisées | m = −l, −l+1, …, 0, …, l−1, l &nbsp;&nbsp; (**dépend de l**) |
| Rôle physique | Détermine l'**orientation** dans l'espace et la projection du moment cinétique sur l'axe z |
| Nombre de valeurs possibles | 2l + 1 |

La projection sur l'axe z vaut `Lz = m ℏ`. En l'absence de champ magnétique, tous les états de même (n, l) ont la même énergie : c'est la **dégénérescence**.

### Résumé des dépendances

```
n = 1, 2, 3, ...
      │
      └─► l = 0, 1, ..., n−1
                │
                └─► m = −l, ..., 0, ..., +l
```

**Exemple pour n = 3 :**

```
n=3
├── l=0 (3s)  →  m = 0                         → 1 état
├── l=1 (3p)  →  m = −1, 0, +1                 → 3 états
└── l=2 (3d)  →  m = −2, −1, 0, +1, +2         → 5 états
                                          Total : 9 états (dégénérés en énergie)
```

---

## 3. La fonction d'onde ψ

La séparation de variables en coordonnées sphériques conduit à une factorisation exacte :

```
ψ_{n,l,m}(r, θ, φ) = R_{n,l}(r) · Y_l^m(θ, φ)
```

où :
- `R_{n,l}(r)` est la **partie radiale** — elle ne dépend que de la distance au noyau
- `Y_l^m(θ, φ)` sont les **harmoniques sphériques** — elles ne dépendent que de la direction

Cette factorisation est fondamentale : elle signifie que la forme de l'orbitale (angulaire) et son extension radiale sont des propriétés séparées et indépendantes.

### Normalisation

La fonction d'onde doit vérifier :

```
∫∫∫ |ψ|² r² sin(θ) dr dθ dφ = 1
```

Ce qui impose des coefficients de normalisation précis dans R et Y.

---

## 4. La partie radiale R_{n,l}(r)

```
         ┌─────────────────────────────────────────────┐  -(r/na₀)    ( 2r  )^l
R_{n,l} = │  (2/na₀)³ · (n−l−1)!                       │ · e          · ─────    · L^{2l+1}_{n−l−1}(2r/na₀)
         │  ─────────────────────────────              │               (na₀)
         └  2n · [(n+l)!]³                             ┘
```

En pratique, on pose `ρ = 2r / (na₀)` pour simplifier :

```
         ┌──────────────────────────────────┐
R_{n,l} = │  (2/na₀)³ · (n−l−1)!           │^(1/2)  · e^{-ρ/2} · ρ^l · L^{2l+1}_{n−l−1}(ρ)
         └  2n · (n+l)!                     ┘
```

Avec :
- `a₀ = 0.529 Å` : rayon de Bohr (en unités atomiques, a₀ = 1)
- `L^k_n(ρ)` : polynôme de Laguerre associé (voir section suivante)

**Interprétation des trois facteurs :**

| Facteur | Rôle |
|---------|------|
| `e^{-ρ/2}` | Décroissance exponentielle — l'électron est confiné près du noyau |
| `ρ^l` | Comportement en r=0 : nul pour l > 0, non nul pour l = 0 (orbitales s) |
| `L^{2l+1}_{n−l−1}(ρ)` | Introduit des **noeuds radiaux** — sphères où la probabilité est nulle |

Le nombre de **noeuds radiaux** (sphères nodales) est exactement `n − l − 1`.

---

## 5. Les polynômes de Laguerre associés

### Origine

Les polynômes de Laguerre émergent de la résolution de l'équation différentielle radiale de Schrödinger. Après substitution et simplification, on obtient l'**équation de Laguerre associée** :

```
x · y''(x) + (k+1−x) · y'(x) + n · y(x) = 0
```

Les seules solutions polynomiales (bornées, normalisables) de cette équation sont les **polynômes de Laguerre associés** `L^k_n(x)`, qui existent uniquement pour n entier non négatif. C'est précisément cette contrainte qui **quantifie** le nombre quantique principal n.

### Formule de récurrence

```
L^k_0(x) = 1
L^k_1(x) = 1 + k − x
L^k_n(x) = [(2n−1+k−x) · L^k_{n-1}(x) − (n−1+k) · L^k_{n-2}(x)] / n
```

### Premiers polynômes explicites

```
L^k_0(x) = 1

L^k_1(x) = −x + (k+1)

L^k_2(x) = x²/2 − (k+2)x + (k+1)(k+2)/2

L^k_3(x) = −x³/6 + (k+3)x²/2 − (k+2)(k+3)x/2 + (k+1)(k+2)(k+3)/6
```

### Rôle physique

Le polynôme `L^{2l+1}_{n−l−1}(ρ)` est de degré `n−l−1`. Il possède exactement `n−l−1` racines réelles positives, ce qui correspond aux **n−l−1 noeuds radiaux** de la fonction d'onde.

> **Note pour les illustrations :** Les graphes de R_{n,l}(r) montrent clairement ces noeuds.  
> Par exemple, R_{3,0} (orbitale 3s) a 2 noeuds radiaux, R_{3,1} (3p) en a 1, et R_{3,2} (3d) n'en a aucun.
> *(Emplacement prévu pour les images des premières orbitales radiales)*

---

## 6. La partie angulaire Y_l^m(θ, φ) — harmoniques sphériques

### Harmoniques sphériques réelles

Pour la visualisation, on préfère des **harmoniques sphériques réelles** (obtenues par combinaison linéaire de Y^{+m} et Y^{-m}) :

```
Y_l^0    (θ,φ) = N_{l,0} · P_l^0(cos θ)
Y_l^{m+} (θ,φ) = N_{l,m} · P_l^m(cos θ) · cos(mφ)    pour m > 0
Y_l^{m-} (θ,φ) = N_{l,m} · P_l^m(cos θ) · sin(|m|φ)  pour m < 0
```

### Facteur de normalisation

```
         ┌──────────────────────────────────────┐
N_{l,m} = │  (2l+1)   (l−|m|)!                 │^(1/2)     ⎧ 1       si m = 0
         │  ──────  · ────────                 │         × ⎨
         └  4π       (l+|m|)!                  ┘           ⎩ √2      si m ≠ 0
```

### Polynômes de Legendre associés P_l^m

Les `P_l^m(cos θ)` sont les polynômes de Legendre associés, qui décrivent la dépendance en θ :

```
P_0^0(x) = 1
P_1^0(x) = x
P_1^1(x) = −√(1−x²)
P_2^0(x) = (3x²−1)/2
P_2^1(x) = −3x√(1−x²)
P_2^2(x) = 3(1−x²)
```

### Interprétation géométrique

Les harmoniques sphériques décrivent la **forme** de l'orbitale :
- `l = 0` (s) : sphère — aucune dépendance angulaire
- `l = 1` (p) : deux lobes — selon z (m=0), ou orientés en x/y (m=±1)
- `l = 2` (d) : quatre lobes ou tore selon m
- La valeur de m sélectionne l'**orientation** dans le plan.

---

## 7. La densité de probabilité |ψ|²

Ce que l'on visualise dans le simulateur n'est pas ψ lui-même (qui peut être négatif ou complexe), mais **|ψ|²** :

```
|ψ_{n,l,m}(r,θ,φ)|² = |R_{n,l}(r)|² · |Y_l^m(θ,φ)|²
```

**Interprétation :** `|ψ|² dV` est la probabilité de trouver l'électron dans le volume élémentaire `dV` situé en (r, θ, φ).

### Conversion vers les coordonnées cartésiennes

Pour tracer dans le plan XZ (φ = 0) :

```
r     = √(x² + z²)
θ     = arccos(z / r)   ← angle polaire depuis l'axe z
φ     = 0               ← dans le plan XZ
```

Pour une visualisation 3D en nuage de points, on tire des coordonnées (x, y, z) aléatoires et on calcule :

```
r     = √(x² + y² + z²)
θ     = arccos(z / r)
φ     = arctan2(y, x)
```

### Densité radiale de probabilité

La probabilité de trouver l'électron **entre r et r+dr** (toutes directions confondues) est :

```
P(r) dr = r² |R_{n,l}(r)|² dr
```

Le facteur r² provient de l'élément de volume en sphérique. La distance la plus probable (maximum de P(r)) vaut `n² a₀` pour les orbitales s, ce qui redonne le rayon de Bohr pour n=1.

---

## 8. Les niveaux d'énergie et les raies spectrales

### Niveaux d'énergie

La résolution de l'équation radiale impose une quantification de l'énergie :

```
Eₙ = − 13.6 eV / n²
```

Les valeurs numériques pour les premiers niveaux :

| n | Énergie (eV) | Couche |
|---|-------------|--------|
| 1 | −13.6 | K |
| 2 | −3.4  | L |
| 3 | −1.51 | M |
| 4 | −0.85 | N |
| ∞ | 0     | ionisation |

**Points importants :**
- L'énergie ne dépend que de n, pas de l ni de m — c'est la **dégénérescence** de l'hydrogène.
- Le signe négatif signifie que l'électron est **lié** : il faut fournir 13.6 eV pour ioniser l'hydrogène depuis l'état fondamental.

### Transitions et raies spectrales

Quand un électron passe d'un niveau n_i à un niveau n_f < n_i, il émet un **photon** d'énergie :

```
ΔE = Eₙᵢ − Eₙf = 13.6 eV · (1/nf² − 1/nᵢ²)
```

La longueur d'onde du photon émis est donnée par la **formule de Rydberg** :

```
1/λ = R_H · (1/nf² − 1/nᵢ²)
```

avec `R_H = 1.097 × 10⁷ m⁻¹` (constante de Rydberg).

### Séries spectrales

| Série | nf | nᵢ | Domaine | Exemple |
|-------|----|----|---------|---------|
| Lyman | 1 | 2, 3, 4, … | Ultraviolet | 121 nm (Lyman-α) |
| Balmer | 2 | 3, 4, 5, … | Visible | 656 nm (rouge, Hα) |
| Paschen | 3 | 4, 5, 6, … | Infrarouge proche | 1875 nm |
| Brackett | 4 | 5, 6, 7, … | Infrarouge moyen | 4051 nm |

La série de Balmer est historiquement importante : ses raies sont visibles à l'œil nu et ont été mesurées avant même que la théorie quantique n'existe. Balmer avait trouvé empiriquement la formule `λ = B · n²/(n²−4)` en 1885 — la mécanique quantique l'a expliquée 40 ans plus tard.

### Lien avec les orbitales

Les raies spectrales sont une **empreinte directe** de la structure des orbitales. Chaque raie correspond à une transition entre deux niveaux d'énergie, eux-mêmes déterminés par les fonctions d'onde ψ_{n,l,m}. Mesurer le spectre d'un atome, c'est lire indirectement la géométrie de ses orbitales.

---

## 9. Catalogue des orbitales principales

### Orbitales s (l = 0, m = 0) — sphériques

| Orbitale | n | l | m | Noeuds radiaux | Forme |
|----------|---|---|---|---------------|-------|
| 1s | 1 | 0 | 0 | 0 | Sphère concentrée |
| 2s | 2 | 0 | 0 | 1 | Sphère avec un noeud sphérique |
| 3s | 3 | 0 | 0 | 2 | Sphère avec deux noeuds sphériques |

Les orbitales s sont les seules à avoir une densité de probabilité **non nulle au noyau** (r = 0).

### Orbitales p (l = 1) — deux lobes

| Orbitale | m | Orientation |
|----------|---|-------------|
| 2p_z | 0 | Lobes le long de l'axe z |
| 2p_x | +1 | Lobes le long de l'axe x |
| 2p_y | −1 | Lobes le long de l'axe y |

Les trois orbitales p sont équivalentes en énergie et forment un ensemble orthogonal.

### Orbitales d (l = 2) — quatre lobes ou tore

| Orbitale | m | Forme |
|----------|---|-------|
| 3d_{z²} | 0 | Deux lobes sur z + tore équatorial |
| 3d_{xz} | +1 | Quatre lobes dans le plan xz |
| 3d_{yz} | −1 | Quatre lobes dans le plan yz |
| 3d_{x²-y²} | +2 | Quatre lobes sur les axes x et y |
| 3d_{xy} | −2 | Quatre lobes à 45° des axes x et y |

### Noeuds angulaires

En plus des noeuds radiaux, les orbitales possèdent `l` **noeuds angulaires** (plans ou cônes où ψ = 0). Le nombre total de noeuds est toujours `n − 1` :

```
noeuds totaux = noeuds radiaux + noeuds angulaires = (n−l−1) + l = n−1
```

---

## 10. Fun fact — L'analogie avec la vibration d'une membrane

Les orbitales atomiques et les **modes de vibration d'une membrane circulaire** (comme un tambour) sont régis par des équations mathématiquement très similaires — la même famille d'équations aux valeurs propres.

### La correspondance

| Membrane circulaire | Atome d'hydrogène |
|--------------------|------------------|
| Équation des ondes ∇²u = λu | Équation de Schrödinger Ĥψ = Eψ |
| Fréquences de résonance discrètes | Niveaux d'énergie discrets |
| Noeuds nodaux (lignes immobiles) | Noeuds de la fonction d'onde |
| Nombre de noeuds circulaires ~ n | Nombre de noeuds radiaux ~ n−l−1 |
| Nombre de diamètres nodaux ~ m | Orientation ~ m |
| Modes fondamental + harmoniques | Orbitale 1s + orbitales excitées |

### Ce que ça signifie physiquement

Frapper un tambour produit un son composé de multiples **modes propres de vibration** — chaque mode a une forme géométrique précise et une fréquence bien définie. De même, l'électron dans l'atome n'est pas une particule qui "tourne" autour du noyau, mais une **onde stationnaire tridimensionnelle** dont les différentes formes sont les orbitales.

Les noeuds des orbitales (surfaces où ψ = 0) sont l'analogue exact des **lignes nodales** visibles sur les figures de Chladni : des zones où la membrane ne vibre jamais. On peut les rendre visibles en saupoudrant du sable sur une plaque métallique que l'on fait vibrer — le sable s'accumule précisément sur les noeuds.

La différence principale : la membrane vibre dans un espace 2D à temps réel, l'électron existe dans un espace 3D sous forme d'une distribution de probabilité — mais la mathématique sous-jacente est la même.

---

## 11. Utilisation du simulateur

### Prérequis

```bash
pip install numpy scipy plotly dash
```

### Lancement

```bash
python orbitales.py
```

Ouvrir ensuite **http://127.0.0.1:8050** dans un navigateur.

> **Important :** ne pas nommer le fichier `code.py` — conflit avec le module Python standard.

### Contrôles

| Slider | Plage | Contrainte |
|--------|-------|-----------|
| n | 1 → 5 | Aucune |
| l | 0 → n−1 | Automatiquement limité par n |
| m | −l → +l | Automatiquement limité par l |

### Paramètres du nuage de points

Dans la fonction `generate_cloud`, deux paramètres permettent d'ajuster qualité et performance :

```python
N = 30000      # nombre de points tirés aléatoirement — plus = plus dense mais plus lent
mask = prob > 0.001   # seuil minimal — plus bas = plus de volume apparent
```

Dans le scatter, `opacity` contrôle la transparence :
```python
opacity=0.06   # entre 0.01 (très transparent) et 0.1 (plus dense)
```

### Unités

Toutes les distances sont en **rayons de Bohr** (a₀ ≈ 0.529 Å). En posant `a_0 = 1` dans le code, on travaille en **unités atomiques** — c'est la convention standard en physique atomique.

---

## Références

- Griffiths, D. J. — *Introduction to Quantum Mechanics*, Cambridge University Press
- Cohen-Tannoudji, C. — *Mécanique Quantique*, Hermann
- NIST — Atomic Spectra Database : https://physics.nist.gov/asd
- Scipy documentation — `scipy.special.genlaguerre`, `scipy.special.lpmv`
