# Analiza profilului genetic în R Shiny

Aplicație R Shiny pentru încărcarea, procesarea și interpretarea fișierelor ADN personale. Aplicația extrage variante genetice de tip `rsID`, le compară cu GWAS Catalog și afișează statistici, interpretări genetice, profil individual, istoric de încărcări și comparații între fișiere.

Proiectul a fost realizat pentru disertație și include atât o variantă completă pentru rulare locală, cât și o variantă demo pentru publicare pe `shinyapps.io`.

## Funcționalități

- încărcare fișiere ADN în formate `23andMe`, `Ancestry`, `FTDNA`, `VCF`, `TXT` și `CSV`;
- extragere și normalizare `rsID`, cromozom și genotip;
- identificarea estimativă a sexului biologic pe baza markerilor X/Y;
- căutare în GWAS Catalog local, cache SQLite și API `gwasrapidd`;
- filtrarea rezultatelor GWAS după alela de risc prezentă în genotip;
- cache SQLite pentru evitarea interogărilor API repetate;
- autentificare cu utilizatori și roluri;
- cereri pentru rol admin, aprobare, respingere și revocare rol admin;
- salvarea fișierelor încărcate și a rezultatelor în SQLite;
- istoric upload-uri;
- comparare 2-3 fișiere prin diagramă Venn;
- raport PDF pentru profilul individual.

## Ce caută aplicația

Aplicația nu caută doar numele SNP-ului. Fluxul de analiză este:

1. Din fișierul ADN se extrag `rsID`, cromozomul și genotipul utilizatorului.
2. `rsID`-urile sunt normalizate, de exemplu `RS4877963` devine `rs4877963`.
3. Aplicația caută aceste `rsID`-uri în cache SQLite, apoi în fișierul local GWAS Catalog, iar pentru lipsuri poate folosi API-ul `gwasrapidd`.
4. Pentru fiecare asociere GWAS se citește coloana `STRONGEST SNP-RISK ALLELE`, de exemplu `rs4877963-T`.
5. Aplicația extrage alela de risc pentru rsID-ul curent, de exemplu `T`.
6. Rezultatul este păstrat doar dacă genotipul utilizatorului conține acea alelă de risc.
7. Dacă același rsID apare în mai multe studii GWAS, aplicația păstrează asocierea reprezentativă cu `p_value` cel mai mic.

Exemplu: dacă GWAS raportează `rs4877963-T`, dar genotipul utilizatorului este `AA`, rândul nu este păstrat. Dacă genotipul este `AT` sau `TT`, asocierea poate fi păstrată.

## Flux GWAS

Aplicația folosește două moduri de lucru:

```text
Local:
cache SQLite -> GWAS Catalog local TSV -> API gwasrapidd pentru lipsuri

shinyapps.io demo:
cache SQLite -> API gwasrapidd limitat
```

În rularea locală, fișierul `gwas_catalog_associations.tsv` permite căutarea rapidă în GWAS Catalog. În demo-ul online, acest fișier este exclus deoarece este foarte mare, iar aplicația folosește API-ul pentru un număr limitat de rsID-uri.

În varianta publicată pe `shinyapps.io`, interogarea API este limitată la maximum `20` de rsID-uri noi per analiză, pentru a evita blocarea aplicației și depășirea resurselor disponibile. Pentru analiza completă se recomandă rularea locală, cu fișierul GWAS Catalog TSV descărcat lângă `app.R`.

## Fișiere pentru demo

Repository-ul conține fișiere demo mici, fără date reale:

```text
demo_23andme_femeie_44.txt
demo_ancestry_femeie_47.txt
demo_csv_femeie_43.csv
demo_ftdna_barbat_45.txt
demo_ftdna_barbat_46.csv
demo_vcf_fara_markeri_42.vcf
demo_vcf_txt_fara_markeri_48.txt
```

Aceste fișiere sunt utile pentru testarea aplicației fără utilizarea unor fișiere ADN reale. Numerele din numele fișierelor indică numărul de variante incluse.

## Fișiere excluse de pe GitHub

Din motive de dimensiune și confidențialitate, următoarele fișiere nu se urcă pe GitHub:

- `gwas_catalog_associations.tsv`;
- baze de date SQLite, precum `genetic_app.sqlite`;
- fișiere ADN reale;
- fișiere R locale, precum `.RData` și `.Rhistory`;
- cache-uri generate local.

Aceste reguli sunt definite în `.gitignore`.

## Rulare locală

Pentru rularea locală se recomandă instalarea `R` și `RStudio Desktop`. Aplicația poate rula și din consola R, dar RStudio este mai ușor de folosit pentru deschiderea proiectului și pornirea aplicației.

1. Instalează R:

[Download R](https://cran.r-project.org/)

2. Instalează RStudio Desktop:

[Download RStudio Desktop](https://posit.co/download/rstudio-desktop/)

3. Deschide folderul aplicației în RStudio.

4. Instalează pachetele necesare în consola R:

```r
install.packages(c(
  "shiny", "tidyverse", "DT", "shinythemes", "plotly",
  "shinycssloaders", "shinyjs", "DBI", "RSQLite",
  "gwasrapidd", "httr", "sodium"
))
```

5. Descarcă fișierul GWAS Catalog de pe pagina oficială:

[GWAS Catalog - File Downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads)

Varianta folosită:

```text
All associations v1.0.2 - with added ontology annotations,
GWAS Catalog study accession numbers and genotyping technology
```

6. Pune fișierul `gwas_catalog_associations.tsv` în folderul aplicației, lângă `app.R`.

7. Rulează aplicația:

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

Linkul generat de `shinyapps.io` poate fi transformat în cod QR pentru prezentarea disertației.

## Roluri utilizatori

Aplicația are două tipuri de admin:

- `admin principal`: contul implicit `admin`, care poate aproba sau respinge cereri de admin și poate revoca roluri admin;
- `admin promovat`: utilizator aprobat de adminul principal, care poate vizualiza și gestiona upload-uri, dar nu poate aproba alți admini.

Utilizatorii obișnuiți pot trimite o cerere pentru rol admin din interfață.

## Confidențialitate

Fișierele ADN conțin date personale sensibile. Nu încărcați pe GitHub fișiere ADN reale, baze de date SQLite cu date reale sau rezultate asociate unor persoane reale.

Pentru date reale sau fișiere mari, se recomandă rularea locală a aplicației.
