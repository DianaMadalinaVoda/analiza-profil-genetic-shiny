# Analiza profilului genetic în R Shiny

Aplicație R Shiny pentru încărcarea, procesarea și interpretarea fișierelor ADN personale. Aplicația extrage variante genetice de tip `rsID`, le compară cu GWAS Catalog și afișează statistici, interpretări genetice, profil individual, istoric de încărcări și comparații între fișiere.

Proiectul a fost realizat pentru disertație și este gândit pentru rulare locală, cu fișierul GWAS Catalog descărcat lângă `app.R` și cu salvarea analizelor într-o bază de date SQLite.

## Funcționalități

- încărcare fișiere ADN în formate `23andMe`, `Ancestry`, `FTDNA`, `VCF`, `TXT` și `CSV`;
- extragere și normalizare `rsID`, cromozom și genotip;
- identificarea estimativă a sexului biologic pe baza markerilor X/Y;
- căutare în GWAS Catalog local, cache SQLite și API `gwasrapidd`;
- filtrarea rezultatelor GWAS după alela de risc prezentă în genotip;
- mesaje de progres pentru pregătirea alelelor de risc, filtrarea genotipurilor și salvarea rezultatelor;
- cache SQLite pentru evitarea interogărilor API repetate;
- autentificare cu utilizatori și roluri;
- cereri pentru rol admin, aprobare, respingere și revocare rol admin;
- schimbarea parolei de către utilizatorul autentificat;
- solicitarea resetării parolei prin cerere adresată administratorului principal, 
  cu generare de parolă temporară și obligativitatea schimbării la prima autentificare;
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

În timpul procesării, interfața afișează separat etapele de pregătire a alelelor de risc și de filtrare a genotipurilor, deoarece aceste operații pot dura mai mult pentru fișiere ADN mari.

## Flux GWAS

Aplicația folosește un flux local cu mecanism de cache și completare prin API:

```text
cache SQLite -> GWAS Catalog local TSV -> API gwasrapidd pentru lipsuri -> filtrare după alela de risc
```

În rularea locală, fișierul `gwas_catalog_associations.tsv` permite căutarea rapidă în GWAS Catalog. Dacă anumite rsID-uri nu sunt găsite local sau în cache, aplicația poate interoga suplimentar API-ul prin `gwasrapidd`, cu limită pentru rsID-uri noi astfel încât analiza să nu blocheze aplicația.

Pentru analiză completă și rapidă se recomandă rularea locală, cu fișierul GWAS Catalog TSV descărcat lângă `app.R`.

## Observație pentru fișiere VCF

În cazul fișierelor VCF, aplicația păstrează doar variantele care au identificator `rsID` valid, deoarece integrarea cu GWAS Catalog se realizează pe baza acestui identificator. Unele fișiere VCF de tip exome conțin poziții genomice fără identificatori `rsID` în coloana `ID`, de exemplu `.`. Aceste poziții sunt valide structural în VCF, dar sunt ignorate în etapa de analiză GWAS deoarece nu pot fi corelate direct cu GWAS Catalog.

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

La compararea a două fișiere, aplicația afișează numărul de rsID-uri comune, numărul total de rsID-uri unice și rsID-urile specifice fiecărui fișier. La compararea a trei fișiere, sunt afișate suplimentar intersecțiile pereche și intersecția comună tuturor celor trei fișiere.

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
GWAS Catalog study accession numbers and genotyping technology (Download (full list))
```

6. Pune fișierul `gwas_catalog_associations.tsv` în folderul aplicației, lângă `app.R`.

7. Rulează aplicația:

```r
shiny::runApp()
```

La prima rulare, baza SQLite `genetic_app.sqlite` este creată automat în folderul aplicației.

Contul implicit de administrator este:

- **Utilizator:** `admin`
- **Parolă inițială:** `genetica123`

La prima autentificare, administratorul este obligat să își schimbe parola din motive de securitate.

> **Notă:** Același mecanism se aplică și utilizatorilor care primesc o parolă 
> temporară prin resetare — la prima autentificare sunt obligați să seteze 
> o parolă nouă.

## Roluri utilizatori

Aplicația are două tipuri de admin:

- `admin principal`: contul implicit `admin`, care poate aproba sau respinge cereri de admin și poate revoca roluri admin;
- `admin promovat`: utilizator aprobat de adminul principal, care poate vizualiza și gestiona upload-uri, dar nu poate aproba alți admini.

Utilizatorii obișnuiți pot:
- trimite o cerere pentru rol admin din interfață;
- solicita resetarea parolei dacă nu o mai cunosc;
- schimba parola din contul propriu după autentificare.

## Resetare parolă

Dacă un utilizator nu își mai cunoaște parola, poate trimite o cerere de resetare 
din ecranul de autentificare, prin butonul **"Am uitat parola"**. Cererea ajunge 
la administratorul principal în tab-ul **Istoric**.

Fluxul este:
1. Utilizatorul introduce numele de utilizator și trimite cererea.
2. Administratorul principal vede cererea și apasă **"Generează parolă temporară"**.
3. Aplicația afișează parola temporară — administratorul o comunică verbal utilizatorului.
4. Utilizatorul se autentifică cu parola temporară.
5. La prima autentificare, sistemul îl obligă să seteze o parolă nouă.

Schimbarea parolei este disponibilă și pentru utilizatorii autentificați, 
din butonul **"Schimbă parola"** din sidebar.

## Confidențialitate

Fișierele ADN conțin date personale sensibile. Nu încărcați pe GitHub fișiere ADN reale, baze de date SQLite cu date reale sau rezultate asociate unor persoane reale.

Pentru date reale sau fișiere mari, se recomandă rularea locală a aplicației.
