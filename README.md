# Analiza profilului genetic în R Shiny

Aplicație R Shiny pentru încărcarea și analiza fișierelor ADN personale în formate precum 23andMe, Ancestry, FTDNA, VCF și CSV. Aplicația extrage rsID-urile, le compară cu GWAS Catalog și afișează statistici, interpretări genetice, profil individual, istoric și comparații între fișiere.

## Funcționalități principale

- autentificare utilizator și rol admin;
- încărcare fișiere ADN în formate diferite;
- extragere rsID, cromozom și genotip;
- căutare GWAS inteligentă: cache SQLite, fișier local GWAS, apoi API prin `gwasrapidd`;
- salvare upload-uri și rezultate în SQLite;
- istoric upload-uri;
- comparare 2-3 fișiere prin diagrama Venn;
- raport PDF pentru profil individual.

## Fișiere care NU sunt incluse pe GitHub

Din motive de dimensiune și confidențialitate, următoarele fișiere nu trebuie urcate în repository:

- `gwas_catalog_associations.tsv`, deoarece este foarte mare;
- fișiere ADN reale, de exemplu 23andMe, Ancestry, FTDNA sau VCF;
- baza de date SQLite locală;
- fișiere generate automat de R, precum `.RData` sau `.Rhistory`.

Acestea sunt excluse prin `.gitignore`.

Important: faptul că `gwas_catalog_associations.tsv` nu este urcat pe GitHub nu înseamnă că aplicația nu îl folosește. Pentru rulare locală sau pentru publicare pe `shinyapps.io`, fișierul trebuie pus manual în folderul aplicației, lângă `app.R`.

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

Varianta folosită în aplicație este:

`All associations v1.0.2 - with added ontology annotations, GWAS Catalog study accession numbers and genotyping technology`

După descărcare, pune fișierul `gwas_catalog_associations.tsv` în folderul aplicației, lângă `app.R`.

Structura locală recomandată este:

```text
app.R
README.md
.gitignore
.rscignore
gwas_catalog_associations.tsv
test_sample.csv
test_sample.txt
test_sample.vcf
```

Fișierul `gwas_catalog_associations.tsv` rămâne local și nu este trimis pe GitHub, dar poate fi inclus în deploy-ul către `shinyapps.io`.

3. Rulează aplicația:

```r
shiny::runApp()
```

## Publicare

Pentru demo public se poate folosi `shinyapps.io`. Linkul generat de `shinyapps.io` poate fi transformat ulterior în cod QR pentru prezentarea disertației.

Atenție: fișierul GWAS Catalog local este mare. Pentru un repository GitHub public se recomandă să nu fie încărcat direct, ci gestionat separat.

Pentru publicare pe `shinyapps.io`, fișierul `gwas_catalog_associations.tsv` trebuie să existe în folderul aplicației înainte de rularea comenzii `rsconnect::deployApp()`. Fișierul este exclus doar de la GitHub prin `.gitignore`, nu și de la deploy-ul `shinyapps.io`.

Baza SQLite nu trebuie creată manual. La prima rulare, aplicația creează automat fișierul `genetic_app.sqlite` în folderul aplicației.

## Confidențialitate

Fișierele ADN conțin date personale sensibile. Nu încărcați pe GitHub fișiere ADN reale, baze de date SQLite cu upload-uri reale sau rezultate asociate unor persoane reale.
