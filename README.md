# Analiza profilului genetic în R Shiny

Aplicație R Shiny pentru încărcarea, procesarea și interpretarea fișierelor ADN personale. Aplicația extrage rsID-uri din fișiere genetice, le compară cu GWAS Catalog și afișează statistici, interpretări genetice, profil individual, istoric de încărcări și comparații între fișiere.

Proiectul a fost realizat pentru disertație și include atât o variantă completă pentru rulare locală, cât și o variantă demo pentru publicare pe `shinyapps.io`.

## Funcționalități

- încărcare fișiere ADN în formate `23andMe`, `Ancestry`, `FTDNA`, `VCF`, `TXT` și `CSV`;
- extragere și normalizare `rsID`, cromozom și genotip;
- identificarea sexului biologic estimativ pe baza markerilor X/Y;
- integrare cu GWAS Catalog prin fișier local și API `gwasrapidd`;
- cache SQLite pentru evitarea interogărilor API repetate;
- autentificare cu utilizatori și roluri;
- cereri pentru rol admin, aprobare, respingere și revocare rol admin;
- salvarea fișierelor încărcate și a rezultatelor în SQLite;
- istoric upload-uri;
- comparare 2-3 fișiere prin diagramă Venn;
- raport PDF pentru profilul individual;
- interfață adaptată pentru desktop și demo pe mobil.

## Flux GWAS

Aplicația folosește două moduri de lucru:

```text
Local:
cache SQLite -> GWAS Catalog local TSV -> API gwasrapidd pentru lipsuri

shinyapps.io demo:
cache SQLite -> API gwasrapidd limitat
```

În rularea locală, fișierul `gwas_catalog_associations.tsv` permite căutarea rapidă în GWAS Catalog. În demo-ul online, acest fișier este exclus deoarece este foarte mare, iar aplicația folosește API-ul pentru un număr limitat de rsID-uri.

## Fișiere excluse de pe GitHub

Din motive de dimensiune și confidențialitate, următoarele fișiere nu se urcă pe GitHub:

- `gwas_catalog_associations.tsv`;
- baze de date SQLite, precum `genetic_app.sqlite`;
- fișiere ADN reale;
- fișiere R locale, precum `.RData` și `.Rhistory`;
- cache-uri generate local.

Aceste reguli sunt definite în `.gitignore`.

## Fișiere pentru demo

Repository-ul poate conține fișiere demo mici, fără date reale:

```text
test_sample.csv
test_sample.txt
test_sample.vcf
test_sample_female.csv
test_sample_XYlipsa.txt
```

Acestea sunt utile pentru testarea aplicației fără utilizarea unor fișiere ADN reale.

## Rulare locală

1. Instalează pachetele necesare în R:

```r
install.packages(c(
  "shiny", "tidyverse", "DT", "shinythemes", "plotly",
  "shinycssloaders", "shinyjs", "DBI", "RSQLite",
  "gwasrapidd", "httr", "sodium"
))
```

2. Descarcă fișierul GWAS Catalog de pe pagina oficială:

[GWAS Catalog - File Downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads)

Varianta folosită:

```text
All associations v1.0.2 - with added ontology annotations,
GWAS Catalog study accession numbers and genotyping technology
```

3. Pune fișierul `gwas_catalog_associations.tsv` în folderul aplicației, lângă `app.R`.

4. Rulează aplicația:

```r
shiny::runApp()
```

La prima rulare, baza SQLite `genetic_app.sqlite` este creată automat în folderul aplicației.

## Publicare pe shinyapps.io

Pentru demo public, aplicația poate fi publicată pe `shinyapps.io`.

Fișierul `gwas_catalog_associations.tsv` este exclus din deploy prin `.rscignore`, deoarece este foarte mare și poate depăși resursele disponibile pe planul gratuit.

Comandă exemplu:

```r
rsconnect::deployApp(
  appDir = "D:/analiza-profil-genetic-shiny",
  appName = "analiza-profil-genetic-shiny"
)
```


## Roluri utilizatori

Aplicația are două tipuri de admin:

- `admin principal`: contul implicit `admin`, care poate aproba/respinge cereri de admin și revoca roluri admin;
- `admin promovat`: utilizator aprobat de adminul principal, care poate vizualiza și gestiona upload-uri, dar nu poate aproba alți admini.

Utilizatorii obișnuiți pot trimite o cerere pentru rol admin din interfață.

## Confidențialitate

Fișierele ADN conțin date personale sensibile. Nu încărcați pe GitHub fișiere ADN reale, baze de date SQLite cu date reale sau rezultate asociate unor persoane reale.

Pentru date reale sau fișiere mari, se recomandă rularea locală a aplicației.
