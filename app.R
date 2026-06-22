# Pachetele sunt încărcate la pornirea aplicației.
# Mesajele de startup sunt ascunse pentru ca consola să rămână ușor de urmărit.
suppressPackageStartupMessages(
  suppressWarnings({
    library(shiny)
    library(tidyverse)
    library(DT)
    library(shinythemes)
    library(plotly)
    library(shinycssloaders)
    library(shinyjs)
    library(DBI)
    library(RSQLite)
    library(grid)
    library(gwasrapidd)
    library(httr)
    library(sodium)
  })
)

options(timeout = 600)


options(shiny.maxRequestSize = 1024 * 1024 * 1024)

# Setări principale ale aplicației: fișierul GWAS local, baza de date și limitele pentru API.
# Fișierul GWAS Catalog poate fi descărcat de aici:
# https://www.ebi.ac.uk/gwas/docs/file-downloads
gwas_file_path = "gwas_catalog_associations.tsv"
gwas_file_path_candidates = c(
  gwas_file_path,
  file.path(getwd(), "gwas_catalog_associations.tsv")
)

# Baza SQLite este creată automat în folderul aplicației, ca proiectul să fie portabil.
db_file_path = file.path(getwd(), "genetic_app.sqlite")
default_app_username = "admin"
default_app_password = "genetica123"
gwas_source_mode = "hybrid"
gwas_query_chunk_size = 25
gwas_cache_lookup_chunk_size = 400
gwas_cache_fetch_chunk_size = 150
gwas_max_workers = 4
gwas_api_max_new_rsids = 20
sex_detection_large_file_min_snps = 1000
sex_detection_y_min_count = 1000
sex_detection_y_min_ratio = 0.005
sex_detection_y_noise_ratio = 0.0005
sex_detection_x_min_count = 1000
sex_detection_x_min_ratio = 0.003
enable_risk_allele_filter = TRUE

# Culori pastel folosite în grafice și în elementele vizuale ale aplicației.
pastel_colors = c("#FFFFCC", "#CCEBC5", "#DECBE4", "#FBB4AE", "#E5D8BD", "#B3CDE3", "#FDDAEC", "#C9E4DE","#EFD3D7", 
                  "#D6E2E9","#F6DFEB", "#E2F0CB", "#B5EAD7","#FFDAC1", "#C7CEEA")

# Mediu intern folosit ca memorie temporară pentru GWAS local, ca să nu recitim TSV-ul la fiecare căutare.
gwas_cache_env = new.env(parent = emptyenv())



# Funcții ajutătoare
clean_column_names = function(x) {
  x = enc2utf8(x)
  x = gsub("^\ufeff", "", x)
  x = gsub("^#\\s*", "", x)
  x = tolower(trimws(x))
  x = gsub("[^a-z0-9]+", "_", x)
  x = gsub("^_+|_+$", "", x)
  x
}

# Împarte vectorii mari în bucăți mai mici, util pentru SQLite și pentru interogările API.
split_into_chunks = function(values, chunk_size) {
  if (length(values) == 0) {
    return(list())
  }
  split(values, ceiling(seq_along(values) / chunk_size))
}

# Găsește prima cale validă către fișierul local GWAS Catalog.
resolve_gwas_file_path = function() {
  existing_paths = gwas_file_path_candidates[file.exists(gwas_file_path_candidates)]
  if (length(existing_paths) == 0) {
    return(NA_character_)
  }
  
  existing_paths[1]
}

# Alege un număr rezonabil de procese paralele, fără să blocheze complet calculatorul.
# Unele exporturi FTDNA/CSV au fiecare rand de date pus intre ghilimele,
# iar virgulele reale raman in interiorul acelui camp. Normalizam randurile
# inainte de citire, ca RSID, CHROMOSOME, POSITION si RESULT sa fie separate.
read_ftdna_data_table = function(file_path) {
  direct_df = tryCatch(
    read.csv(
      file_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE,
      quote = "\"",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  
  if (!is.null(direct_df) && ncol(direct_df) >= 4) {
    first_col = as.character(direct_df[[1]])
    first_col = first_col[!is.na(first_col) & nzchar(first_col)]
    
    if (length(first_col) > 0 && any(grepl("^rs[0-9]+$", tolower(trimws(first_col))))) {
      return(direct_df)
    }
    
    if (length(first_col) > 0 && any(grepl("^rs[0-9]+,", tolower(trimws(first_col))))) {
      reparsed_df = tryCatch(
        read.csv(
          text = paste(c("rsid,chromosome,position,result", first_col), collapse = "\n"),
          header = TRUE,
          stringsAsFactors = FALSE,
          fill = TRUE,
          check.names = FALSE,
          quote = "\"",
          comment.char = ""
        ),
        error = function(e) NULL
      )
      
      if (!is.null(reparsed_df) && ncol(reparsed_df) >= 4) {
        return(reparsed_df)
      }
    }
  }
  
  ftdna_lines = readLines(file_path, warn = FALSE, encoding = "UTF-8")
  ftdna_lines = gsub("^\ufeff", "", ftdna_lines)
  ftdna_lines = trimws(ftdna_lines)
  ftdna_lines = ftdna_lines[nzchar(ftdna_lines)]
  
  normalize_ftdna_line = function(current_line) {
    trimmed_line = trimws(current_line)
    
    if (grepl('^".*"$', trimmed_line)) {
      trimmed_line = sub('^"', "", trimmed_line)
      trimmed_line = sub('"$', "", trimmed_line)
      trimmed_line = gsub('""', '"', trimmed_line, fixed = TRUE)
    }
    
    trimmed_line
  }
  
  if (length(ftdna_lines) < 2) {
    return(data.frame())
  }
  
  ftdna_lines = vapply(ftdna_lines[-1], normalize_ftdna_line, character(1))
  ftdna_rows = strsplit(ftdna_lines, ",", fixed = TRUE)
  ftdna_rows = ftdna_rows[lengths(ftdna_rows) >= 4]
  
  if (length(ftdna_rows) == 0) {
    return(data.frame())
  }
  
  clean_ftdna_field = function(value) {
    value = trimws(as.character(value))
    value = gsub('""', '"', value, fixed = TRUE)
    gsub('^"|"$', "", value)
  }
  
  data.frame(
    rsid = clean_ftdna_field(vapply(ftdna_rows, function(x) x[1], character(1))),
    chromosome = clean_ftdna_field(vapply(ftdna_rows, function(x) x[2], character(1))),
    position = clean_ftdna_field(vapply(ftdna_rows, function(x) x[3], character(1))),
    result = clean_ftdna_field(vapply(ftdna_rows, function(x) x[4], character(1))),
    stringsAsFactors = FALSE
  )
}

get_parallel_worker_count = function(task_count) {
  available = parallel::detectCores(logical = TRUE)
  if (is.na(available) || available < 1) {
    available = 1
  }
  max(1, min(gwas_max_workers, max(1, available - 1), task_count))
}

# Instrumentare usoara pentru testele de evaluare rulate manual din aplicatie.
get_r_memory_metrics = function(reset_peak = FALSE) {
  gc_stats = suppressWarnings(gc(reset = reset_peak))
  max_used_col = match("max used", colnames(gc_stats))
  
  list(
    used_mb = sum(as.numeric(gc_stats[, 2]), na.rm = TRUE),
    peak_mb = if (!is.na(max_used_col) && max_used_col < ncol(gc_stats)) {
      sum(as.numeric(gc_stats[, max_used_col + 1]), na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

start_evaluation_metrics = function(file_names) {
  memory_start = get_r_memory_metrics(reset_peak = TRUE)
  
  list(
    started_at = Sys.time(),
    proc_time = proc.time(),
    memory_start_mb = memory_start$used_mb,
    file_names = as.character(file_names)
  )
}

finish_evaluation_metrics = function(metrics, success) {
  proc_delta = proc.time() - metrics$proc_time
  elapsed_seconds = as.numeric(difftime(Sys.time(), metrics$started_at, units = "secs"))
  cpu_seconds = unname(proc_delta["user.self"] + proc_delta["sys.self"])
  cpu_percent = if (is.finite(elapsed_seconds) && elapsed_seconds > 0) {
    100 * cpu_seconds / elapsed_seconds
  } else {
    NA_real_
  }
  memory_end = get_r_memory_metrics()
  file_label = paste(short_file_label(metrics$file_names, 54), collapse = " | ")
  

  cat("[EVALUARE] Status analiza:", if (isTRUE(success)) "finalizata" else "eroare", "\n")
  cat("[EVALUARE] Fisiere selectate:", length(metrics$file_names), "\n")
  cat("[EVALUARE] Nume fisiere:", file_label, "\n")
  cat(sprintf("[EVALUARE] Timp total: %.2f secunde\n", elapsed_seconds))
  cat(sprintf("[EVALUARE] CPU R: %.2f secunde | CPU mediu R: %.1f%%\n", cpu_seconds, cpu_percent))
  cat(sprintf(
    "[EVALUARE] RAM R: start %.1f MB | dupa analiza %.1f MB | maxim in analiza %.1f MB\n",
    metrics$memory_start_mb,
    memory_end$used_mb,
    memory_end$peak_mb
  ))
}

# Securizare parole
# Dacă în baza de date există parole vechi în clar, ele sunt acceptate doar la primul login
# și apoi sunt transformate automat în hash sodium.

password_hash_type = function(password_value) {
  if (is.na(password_value) || !nzchar(password_value)) {
    return("empty")
  }
  
  if (grepl("^\\$7\\$", password_value)) {
    return("sodium")
  }
  
  "plain"
}

password_is_hashed = function(password_value) {
  password_hash_type(password_value) == "sodium"
}

# Doar contul admin principal poate aproba cereri de promovare.
is_primary_admin = function(username) {
  identical(username, default_app_username)
}

hash_password = function(password_value) {
  if (!requireNamespace("sodium", quietly = TRUE)) {
    stop("Pentru securizarea parolelor instalează pachetul sodium: install.packages('sodium')")
  }
  
  sodium::password_store(password_value)
}

verify_password = function(stored_password, password_value) {
  if (is.na(stored_password) || !nzchar(stored_password)) {
    return(FALSE)
  }
  
  hash_type = password_hash_type(stored_password)
  
  if (hash_type == "sodium") {
    if (!requireNamespace("sodium", quietly = TRUE)) {
      stop("Parola este hash-uita cu sodium, dar pachetul sodium nu este instalat. Ruleaza install.packages('sodium').")
    }
    
    return(isTRUE(sodium::password_verify(stored_password, password_value)))
  }
  
  # Compatibilitate pentru conturile vechi, dacă în DB existau parole în clar.
  identical(stored_password, password_value)
}

# Convertește genotipul din format VCF (0/1, 1/1 etc.) în baze reale.
convert_vcf_gt_to_genotype = function(gt, ref, alt) {
  gt = as.character(gt)
  ref = as.character(ref)
  alt = as.character(alt)
  
  if (is.na(gt) || !nzchar(gt) || gt %in% c(".", "./.", ".|.")) {
    return("--")
  }
  
  gt_main = strsplit(gt, ":", fixed = TRUE)[[1]][1]
  alleles = unlist(strsplit(gt_main, "[/|]"))
  alt_values = unlist(strsplit(alt, ",", fixed = TRUE))
  
  decode_allele = function(idx) {
    if (is.na(idx) || idx == ".") {
      return("-")
    }
    
    idx_num = suppressWarnings(as.integer(idx))
    if (is.na(idx_num)) {
      return("-")
    }
    
    if (idx_num == 0) {
      return(ref)
    }
    
    if (idx_num >= 1 && idx_num <= length(alt_values)) {
      return(alt_values[idx_num])
    }
    
    "-"
  }
  
  genotype = paste0(vapply(alleles, decode_allele, character(1)), collapse = "")
  genotype = toupper(gsub("[^ACGT]", "-", genotype))
  
  if (!nzchar(genotype)) "--" else genotype
}


# CITIRE FIȘIERE ADN


# Detectează formatul fișierului încărcat și îl transformă într-un format comun.
detect_declared_sex_from_filename = function(file_name) {
  file_name = tolower(basename(as.character(file_name)))
  
  sex_match = regexec("sex[_-]?(xy|xx|unknown)", file_name, ignore.case = TRUE)
  sex_parts = regmatches(file_name, sex_match)[[1]]
  
  if (length(sex_parts) == 0) {
    return(NA_character_)
  }
  
  declared_sex = tolower(sex_parts[2])
  
  if (identical(declared_sex, "xy")) {
    return("Masculin")
  }
  
  if (identical(declared_sex, "xx")) {
    return("Feminin")
  }
  
  if (identical(declared_sex, "unknown")) {
    return("Necunoscut")
  }
  
  NA_character_
}

read_adn_data = function(file_path, original_file_name = NULL) {
  # Citim un bloc mai mare, deoarece fișierele VCF pot avea multe linii de metadate înainte de #CHROM.
  header_lines = readLines(file_path, n = 1000, warn = FALSE, encoding = "UTF-8")
  normalized_header_lines = trimws(gsub("^\ufeff", "", header_lines))
  normalized_header_lines = gsub('^"|"$', "", normalized_header_lines)
  normalized_header_lines = gsub('""', '"', normalized_header_lines, fixed = TRUE)
  header_text = paste(header_lines, collapse = "\n")
  
  # Aplicația acceptă formate diferite, dar toate sunt convertite ulterior la aceleași coloane:
  # rsid, cromozom și genotip.
  is_ftdna = any(grepl("^RSID,CHROMOSOME,POSITION,RESULT", normalized_header_lines, ignore.case = TRUE))
  vcf_header_index = grep("^#CHROM\\s+POS\\s+ID\\s+REF\\s+ALT", header_lines, ignore.case = TRUE)
  is_vcf = any(grepl("^##fileformat=VCF", header_lines, ignore.case = TRUE)) ||
    length(vcf_header_index) > 0
  is_23andme = grepl("23andMe", header_text, ignore.case = TRUE) ||
    any(grepl("^#\\s*rsid\\s+chromosome\\s+position\\s+genotype\\s*$", header_lines, ignore.case = TRUE))
  is_ancestry = grepl("AncestryDNA", header_text, ignore.case = TRUE) ||
    any(grepl("^rsid\\tchromosome\\tposition\\tallele1\\tallele2$", header_lines, ignore.case = TRUE))
  first_data_line_index = which(nzchar(trimws(header_lines)) & !grepl("^#", trimws(header_lines)))[1]
  first_data_line = if (is.na(first_data_line_index)) NA_character_ else header_lines[first_data_line_index]
  is_headerless_genotype = !is.na(first_data_line) &&
    grepl("^rs[0-9]+[,[:space:]]+[0-9XYMTA-Za-z]+[,[:space:]]+[0-9]+[,[:space:]]+[A-Za-z-]+", first_data_line, ignore.case = TRUE)
  
  df = tryCatch({
    if (is_vcf) {
      if (length(vcf_header_index) == 0) {
        return(NULL)
      }
      
      vcf_col_names = strsplit(trimws(header_lines[vcf_header_index[1]]), "\\s+")[[1]]
      vcf_col_names[1] = sub("^#", "", vcf_col_names[1])
      vcf_preview_lines = header_lines[(vcf_header_index[1] + 1):length(header_lines)]
      vcf_preview_lines = vcf_preview_lines[nzchar(trimws(vcf_preview_lines))]
      vcf_field_count = if (length(vcf_preview_lines) == 0) {
        length(vcf_col_names)
      } else {
        max(lengths(strsplit(vcf_preview_lines, "\t", fixed = TRUE)))
      }
      if (vcf_field_count > length(vcf_col_names)) {
        vcf_col_names = c(vcf_col_names, paste0("extra_", seq_len(vcf_field_count - length(vcf_col_names))))
      }
      # Unele fișiere VCF au tab-uri finale sau câmpuri suplimentare accidentale.
      # Coloanele extra sunt ignorate ulterior, dar previn respingerea fișierului la citire.
      vcf_col_names = c(vcf_col_names, paste0("extra_padding_", seq_len(10)))
      
      vcf_df = read.table(
        file_path,
        header = FALSE,
        sep = "\t",
        skip = vcf_header_index[1],
        col.names = vcf_col_names,
        stringsAsFactors = FALSE,
        fill = TRUE,
        check.names = FALSE,
        comment.char = ""
      )
      
      if (ncol(vcf_df) < 8) {
        vcf_df = read.table(
          file_path,
          header = FALSE,
          sep = "",
          skip = vcf_header_index[1],
          col.names = vcf_col_names,
          stringsAsFactors = FALSE,
          fill = TRUE,
          check.names = FALSE,
          comment.char = ""
        )
      }
      
      vcf_df
    } else if (is_ftdna) {
      read_ftdna_data_table(file_path)
    } else if (is_23andme) {
      header_index = grep("^#\\s*rsid\\s+chromosome\\s+position\\s+genotype\\s*$", header_lines, ignore.case = TRUE)
      if (length(header_index) == 0) {
        return(NULL)
      }
      read.table(
        file_path,
        header = FALSE,
        sep = "",
        skip = header_index[1],
        col.names = c("rsid", "chromosome", "position", "genotype"),
        stringsAsFactors = FALSE,
        fill = TRUE,
        check.names = FALSE,
        comment.char = ""
      )
    } else if (is_headerless_genotype) {
      sep_guess = if (grepl(",", first_data_line, fixed = TRUE)) "," else ""
      raw_df = read.table(
        file_path,
        header = FALSE,
        sep = sep_guess,
        stringsAsFactors = FALSE,
        fill = TRUE,
        check.names = FALSE,
        quote = "\"",
        comment.char = "#"
      )
      colnames(raw_df)[seq_len(min(4, ncol(raw_df)))] = c("rsid", "chromosome", "position", "genotype")[seq_len(min(4, ncol(raw_df)))]
      raw_df
    } else if (is_ancestry) {
      read.table(
        file_path,
        header = TRUE,
        sep = "\t",
        comment.char = "#",
        stringsAsFactors = FALSE,
        fill = TRUE,
        check.names = FALSE
      )
    } else {
      sep_guess = if (grepl(",", first_data_line, fixed = TRUE)) "," else "\t"
      read.table(
        file_path,
        header = TRUE,
        sep = sep_guess,
        stringsAsFactors = FALSE,
        fill = TRUE,
        check.names = FALSE,
        quote = "\"",
        comment.char = ""
      )
    }
  }, error = function(e) return(NULL))
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  colnames(df) = clean_column_names(colnames(df))
  
  # Fiecare format are structură proprie, de aceea extragerea genotipului se face separat.
  if (is_vcf) {
    sample_col = setdiff(colnames(df), c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"))[1]
    if (is.na(sample_col) || !nzchar(sample_col)) {
      return(NULL)
    }
    
    gt_values = if ("format" %in% colnames(df)) as.character(df[[sample_col]]) else rep(NA_character_, nrow(df))
    
    res = data.frame(
      rsid = trimws(as.character(df$id)),
      chrom = trimws(as.character(df$chrom)),
      genotype = mapply(
        convert_vcf_gt_to_genotype,
        gt = gt_values,
        ref = df$ref,
        alt = df$alt,
        USE.NAMES = FALSE
      ),
      stringsAsFactors = FALSE
    )
  } else if (is_ftdna) {
    if (!all(c("rsid", "chromosome", "result") %in% colnames(df)) && ncol(df) == 1) {
      ftdna_rows = strsplit(as.character(df[[1]]), ",", fixed = TRUE)
      ftdna_rows = ftdna_rows[lengths(ftdna_rows) >= 4]
      
      if (length(ftdna_rows) == 0) {
        return(NULL)
      }
      
      df = data.frame(
        rsid = vapply(ftdna_rows, function(x) x[1], character(1)),
        chromosome = vapply(ftdna_rows, function(x) x[2], character(1)),
        position = vapply(ftdna_rows, function(x) x[3], character(1)),
        result = vapply(ftdna_rows, function(x) x[4], character(1)),
        stringsAsFactors = FALSE
      )
      df[] = lapply(df, function(x) gsub('^"|"$', "", x))
    }
    
    if (!all(c("rsid", "chromosome", "result") %in% colnames(df))) {
      return(NULL)
    }
    
    res = data.frame(
      rsid = trimws(as.character(df$rsid)),
      chrom = trimws(as.character(df$chromosome)),
      genotype = trimws(as.character(df$result)),
      stringsAsFactors = FALSE
    )
  } else if (is_23andme) {
    res = data.frame(
      rsid = trimws(as.character(df$rsid)),
      chrom = if ("chromosome" %in% colnames(df)) trimws(as.character(df$chromosome)) else "NA",
      genotype = if ("genotype" %in% colnames(df)) toupper(trimws(as.character(df$genotype))) else "--",
      stringsAsFactors = FALSE
    )
  } else {
    n_cols = ncol(df)
    rsid_col = intersect(c("rsid", "rs_id", "snp", "snps", "marker", "id"), colnames(df))[1]
    if (is.na(rsid_col) || !nzchar(rsid_col)) {
      return(NULL)
    }
    
    if (identical(rsid_col, "snps")) {
      snp_list = strsplit(as.character(df[[rsid_col]]), ";", fixed = TRUE)
      res = data.frame(
        rsid = unlist(snp_list, use.names = FALSE),
        chrom = "NA",
        genotype = "--",
        stringsAsFactors = FALSE
      )
    } else {
      chrom_col = intersect(c("chromosome", "chrom", "chr", "chromosome_name"), colnames(df))[1]
      res = data.frame(
        rsid = as.character(df[[rsid_col]]),
        chrom = if (!is.na(chrom_col) && nzchar(chrom_col)) as.character(df[[chrom_col]]) else NA_character_,
        stringsAsFactors = FALSE
      )
      
      if ("allele1" %in% colnames(df) && "allele2" %in% colnames(df)) {
        res$genotype = toupper(paste0(df$allele1, df$allele2))
      } else if (n_cols >= 4) {
        res$genotype = toupper(as.character(df[[4]]))
      } else {
        res$genotype = "--"
      }
    }
  }
  
  res$rsid = tolower(trimws(as.character(res$rsid)))
  # Se păstrează doar identificatorii rsID valizi, pentru a evita trimiterile greșite către GWAS.
  res = res %>% filter(grepl("^rs[0-9]+$", rsid))
  res$genotype = gsub("[^A-Z]", "-", res$genotype)
  res$chrom = normalize_chrom_value(res$chrom)
  
  # Detectarea sexului este orientativă. Semnalele X/Y slabe sau amestecate sunt tratate ca necunoscute.
  total_valid_markers = nrow(res)
  y_marker_count = sum(res$chrom == "Y", na.rm = TRUE)
  x_marker_count = sum(res$chrom == "X", na.rm = TRUE)
  y_marker_ratio = y_marker_count / max(1, total_valid_markers)
  x_marker_ratio = x_marker_count / max(1, total_valid_markers)
  has_y_signal = if (total_valid_markers < sex_detection_large_file_min_snps) {
    y_marker_count > 0
  } else {
    y_marker_count >= sex_detection_y_min_count &&
      y_marker_ratio >= sex_detection_y_min_ratio
  }
  has_x_signal = if (total_valid_markers < sex_detection_large_file_min_snps) {
    x_marker_count > 0
  } else {
    x_marker_count >= sex_detection_x_min_count &&
      x_marker_ratio >= sex_detection_x_min_ratio
  }
  has_y_noise = if (total_valid_markers < sex_detection_large_file_min_snps) {
    FALSE
  } else {
    y_marker_count > 0 && y_marker_ratio > sex_detection_y_noise_ratio
  }
  
  declared_sex_from_filename = detect_declared_sex_from_filename(original_file_name)
  
  if (identical(declared_sex_from_filename, "Masculin") || identical(declared_sex_from_filename, "Feminin")) {
    sex_status = declared_sex_from_filename
  } else if (identical(declared_sex_from_filename, "Necunoscut")) {
    # Pentru fișierele marcate explicit ca sex_unknown păstrăm regula simplă:
    # orice marker Y indică profil masculin, altfel prezența X indică profil feminin.
    if (y_marker_count > 0) {
      sex_status = "Masculin"
    } else if (x_marker_count > 0) {
      sex_status = "Feminin"
    } else {
      sex_status = "Necunoscut (lipsă markeri sexuali)"
    }
  } else if (has_y_signal) {
    sex_status = "Masculin"
  } else if (has_x_signal && !has_y_noise) {
    sex_status = "Feminin"
  } else {
    sex_status = "Necunoscut (lipsă markeri sexuali)"
  }
  
  return(list(
    data = distinct(res, rsid, .keep_all = TRUE),
    sex = sex_status
  ))
}



# CITIRE REFERINȚĂ GWAS


# Citește o copie locală a GWAS Catalog și o ține în memorie pentru viteză.
read_gwas_data_local = function(variant_ids) {
  # Normalizarea previne diferențele de forma între RSID-urile din fișiere și cele din GWAS.
  variant_ids = tolower(trimws(as.character(variant_ids)))
  variant_ids = unique(variant_ids[grepl("^rs[0-9]+$", variant_ids)])
  
  resolved_gwas_file_path = resolve_gwas_file_path()
  if (is.na(resolved_gwas_file_path)) {
    stop("Fișierul local GWAS Catalog nu a fost gășit pentru filtrarea rapidă.")
  }
  
  if (!exists("local_gwas_ref", envir = gwas_cache_env, inherits = FALSE)) {
    # TSV-ul local este citit o singură dată pe sesiune și apoi refolosit din memorie.
    gwas = read.delim(
      resolved_gwas_file_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = ""
    )
    
    if (!"SNPS" %in% colnames(gwas)) {
      stop("Fișierul local GWAS trebuie să conțină coloană SNPS.")
    }
    
    snp_list = strsplit(as.character(gwas$SNPS), ";", fixed = TRUE)
    
    # O linie GWAS poate conține mai multe SNP-uri; aici desfacem lista astfel încât fiecare rsID
    # să poată fi căutat individual.
    ref_local = data.frame(
      rsid = unlist(snp_list, use.names = FALSE),
      categoria = rep(gwas$MAPPED_TRAIT, lengths(snp_list)),
      disease_trait = rep(gwas$`DISEASE/TRAIT`, lengths(snp_list)),
      mapped_trait = rep(gwas$MAPPED_TRAIT, lengths(snp_list)),
      strongest_snp_risk_allele = rep(gwas$`STRONGEST SNP-RISK ALLELE`, lengths(snp_list)),
      p_value = rep(gwas$`P-VALUE`, lengths(snp_list)),
      mapped_gene = rep(gwas$MAPPED_GENE, lengths(snp_list)),
      study = rep(gwas$STUDY, lengths(snp_list)),
      pubmedid = rep(gwas$PUBMEDID, lengths(snp_list)),
      stringsAsFactors = FALSE
    )
    
    ref_local$rsid = tolower(trimws(ref_local$rsid))
    ref_local = ref_local %>%
      filter(grepl("^rs[0-9]+$", rsid)) %>%
      distinct()
    
    assign("local_gwas_ref", ref_local, envir = gwas_cache_env)
  }
  
  ref_local = get("local_gwas_ref", envir = gwas_cache_env)
  ref_local %>% filter(rsid %in% variant_ids)
}

# Folosește TSV-ul local doar ca index rapid; adnotările finale vin din API.
get_gwas_candidate_rsids = function(variant_ids) {
  variant_ids = tolower(trimws(as.character(variant_ids)))
  variant_ids = unique(variant_ids[grepl("^rs[0-9]+$", variant_ids)])
  
  if (length(variant_ids) == 0) {
    return(character())
  }
  
  if (is.na(resolve_gwas_file_path())) {
    stop("Fișierul local GWAS Catalog lipsește, deci nu pot prefiltra înainte de API.")
  }
  
  local_ref = read_gwas_data_local(variant_ids)
  unique(local_ref$rsid)
}

# Interoghează GWAS Catalog online prin gwasrapidd pentru rsid-urile din fișierele ADN.
read_gwas_data_online = function(variant_ids) {
  variant_ids = tolower(trimws(as.character(variant_ids)))
  variant_ids = unique(variant_ids[grepl("^rs[0-9]+$", variant_ids)])
  
  if (length(variant_ids) == 0) {
    return(empty_gwas_ref())
  }
  
  variant_chunks = split_into_chunks(variant_ids, gwas_query_chunk_size)
  
  # Funcție defensivă: încearcă de mai multe ori același batch, deoarece API-ul poate răspunde lent.
  safe_get_associations = function(chunk_ids, max_tries = 3) {
    for (i in seq_len(max_tries)) {
      res = tryCatch(
        gwasrapidd::get_associations(
          variant_id = chunk_ids,
          set_operation = "union",
          interactive = FALSE,
          verbose = FALSE,
          warnings = FALSE
        ),
        error = function(e) e
      )
      
      if (!inherits(res, "error")) {
        return(res)
      }
      
      if (i < max_tries) {
        Sys.sleep(1.5 * i)
      }
    }
    
    NULL
  }
  
  workers = get_parallel_worker_count(length(variant_chunks))
  
  # Dacă pachetele future/furrr există, batch-urile API se pot rula în paralel.
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE) && workers > 1) {
    old_plan = future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    
    association_results = furrr::future_map(
      variant_chunks,
      safe_get_associations,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    association_results = lapply(variant_chunks, safe_get_associations)
  }
  
  association_results = Filter(
    function(x) !is.null(x) && nrow(x@associations) > 0,
    association_results
  )
  
  if (length(association_results) == 0) {
    return(empty_gwas_ref())
  }
  
  associations_obj = association_results[[1]]
  if (length(association_results) > 1) {
    for (i in 2:length(association_results)) {
      associations_obj = gwasrapidd::bind(associations_obj, association_results[[i]])
    }
  }
  
  association_ids = unique(associations_obj@associations$association_id)
  
  if (length(association_ids) == 0) {
    return(empty_gwas_ref())
  }
  
  association_id_chunks = split_into_chunks(association_ids, 200)
  
  # După găsirea asocierilor, se cer studiile și trait-urile ca să putem afișa informații complete.
  safe_get_studies = function(chunk_ids, max_tries = 3) {
    for (i in seq_len(max_tries)) {
      res = tryCatch(
        gwasrapidd::get_studies(
          association_id = chunk_ids,
          set_operation = "union",
          interactive = FALSE,
          verbose = FALSE,
          warnings = FALSE
        ),
        error = function(e) e
      )
      
      if (!inherits(res, "error")) {
        return(res)
      }
      
      if (i < max_tries) {
        Sys.sleep(1.5 * i)
      }
    }
    
    NULL
  }
  
  safe_get_traits = function(chunk_ids, max_tries = 3) {
    for (i in seq_len(max_tries)) {
      res = tryCatch(
        gwasrapidd::get_traits(
          association_id = chunk_ids,
          set_operation = "union",
          verbose = FALSE,
          warnings = FALSE
        ),
        error = function(e) e
      )
      
      if (!inherits(res, "error")) {
        return(res)
      }
      
      if (i < max_tries) {
        Sys.sleep(1.5 * i)
      }
    }
    
    NULL
  }
  
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE) && workers > 1) {
    studies_results = furrr::future_map(
      association_id_chunks,
      safe_get_studies,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    studies_results = lapply(association_id_chunks, safe_get_studies)
  }
  studies_results = Filter(Negate(is.null), studies_results)
  
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE) && workers > 1) {
    traits_results = furrr::future_map(
      association_id_chunks,
      safe_get_traits,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    traits_results = lapply(association_id_chunks, safe_get_traits)
  }
  traits_results = Filter(Negate(is.null), traits_results)
  
  studies_obj = if (length(studies_results) == 0) {
    NULL
  } else {
    studies_results[[1]]
  }
  
  if (length(studies_results) > 1) {
    for (i in 2:length(studies_results)) {
      studies_obj = gwasrapidd::bind(studies_obj, studies_results[[i]])
    }
  }
  
  traits_obj = if (length(traits_results) == 0) {
    NULL
  } else {
    traits_results[[1]]
  }
  
  if (length(traits_results) > 1) {
    for (i in 2:length(traits_results)) {
      traits_obj = gwasrapidd::bind(traits_obj, traits_results[[i]])
    }
  }
  
  association_tbl = associations_obj@associations %>%
    transmute(
      association_id = association_id,
      p_value = as.character(pvalue)
    )
  
  risk_tbl = associations_obj@risk_alleles %>%
    transmute(
      association_id = association_id,
      locus_id = locus_id,
      rsid = tolower(variant_id),
      risk_allele = risk_allele
    )
  
  genes_tbl = associations_obj@genes %>%
    transmute(
      association_id = association_id,
      locus_id = locus_id,
      mapped_gene = gene_name
    )
  
  study_map_tbl = gwasrapidd::association_to_study(
    association_id = association_ids,
    verbose = FALSE,
    warnings = FALSE
  )
  
  trait_map_tbl = gwasrapidd::association_to_trait(
    association_id = association_ids,
    verbose = FALSE,
    warnings = FALSE
  )
  
  studies_tbl = if (is.null(studies_obj)) {
    data.frame(study_id = character(), disease_trait = character(), stringsAsFactors = FALSE)
  } else {
    studies_obj@studies %>%
      transmute(
        study_id = study_id,
        disease_trait = reported_trait
      )
  }
  
  publications_tbl = if (is.null(studies_obj)) {
    data.frame(study_id = character(), pubmedid = character(), study = character(), stringsAsFactors = FALSE)
  } else {
    studies_obj@publications %>%
      transmute(
        study_id = study_id,
        pubmedid = as.character(pubmed_id),
        study = title
      )
  }
  
  traits_tbl = if (is.null(traits_obj)) {
    data.frame(efo_id = character(), mapped_trait = character(), stringsAsFactors = FALSE)
  } else {
    traits_obj@traits %>%
      transmute(
        efo_id = efo_id,
        mapped_trait = trait
      )
  }
  
  ref = risk_tbl %>%
    left_join(
      association_tbl,
      by = "association_id",
      relationship = "many-to-many"
    )%>%
    left_join(
      genes_tbl,
      by = c("association_id", "locus_id"),
      relationship = "many-to-many"
    ) %>%
    left_join(study_map_tbl, by = "association_id") %>%
    left_join(studies_tbl, by = "study_id") %>%
    left_join(publications_tbl, by = "study_id") %>%
    left_join(trait_map_tbl, by = "association_id") %>%
    left_join(traits_tbl, by = "efo_id") %>%
    mutate(
      categoria = mapped_trait,
      strongest_snp_risk_allele = ifelse(
        is.na(risk_allele) | risk_allele == "",
        rsid,
        paste0(rsid, "-", risk_allele)
      )
    ) %>%
    select(
      rsid,
      categoria,
      disease_trait,
      mapped_trait,
      strongest_snp_risk_allele,
      p_value,
      mapped_gene,
      study,
      pubmedid
    ) %>%
    distinct()
  
  # Tabelul final este redus la coloanele folosite în interfață și în raport.
  ref
}

read_gwas_data = function(variant_ids) {
  if (identical(gwas_source_mode, "local")) {
    return(read_gwas_data_local(variant_ids))
  }
  
  if (identical(gwas_source_mode, "online")) {
    return(read_gwas_data_online(variant_ids))
  }
  
  if (identical(gwas_source_mode, "hybrid")) {
    candidate_rsids = get_gwas_candidate_rsids(variant_ids)
    return(read_gwas_data_online(candidate_rsids))
  }
  
  if (identical(gwas_source_mode, "auto") && !is.na(resolve_gwas_file_path())) {
    return(read_gwas_data_local(variant_ids))
  }
  
  read_gwas_data_online(variant_ids)
}

empty_gwas_ref = function() {
  data.frame(
    rsid = character(),
    categoria = character(),
    disease_trait = character(),
    mapped_trait = character(),
    strongest_snp_risk_allele = character(),
    p_value = character(),
    mapped_gene = character(),
    study = character(),
    pubmedid = character(),
    stringsAsFactors = FALSE
  )
}

# Uniformizează coloanele GWAS înainte de bind_rows, indiferent dacă vin din cache, TSV sau API.
normalize_gwas_ref = function(ref_df) {
  template = empty_gwas_ref()
  
  if (is.null(ref_df) || nrow(ref_df) == 0) {
    return(template)
  }
  
  missing_cols = setdiff(names(template), names(ref_df))
  for (col_name in missing_cols) {
    ref_df[[col_name]] = NA_character_
  }
  
  ref_df = ref_df[, names(template), drop = FALSE]
  ref_df[] = lapply(ref_df, function(x) as.character(x))
  ref_df
}

# Extrage alela de risc corespunzatoare rsID-ului curent din campul GWAS.
extract_risk_allele_for_rsid = function(risk_field, rsid) {
  risk_field = as.character(risk_field)
  rsid = tolower(trimws(as.character(rsid)))
  allele_result = rep(NA_character_, length(risk_field))
  
  valid_idx = !is.na(risk_field) & nzchar(risk_field) & !is.na(rsid) & nzchar(rsid)
  direct_idx = valid_idx & startsWith(tolower(trimws(risk_field)), paste0(rsid, "-"))
  
  if (any(direct_idx)) {
    allele_result[direct_idx] = toupper(gsub(
      "[^ACGT]",
      "",
      sub("^[^-]+-", "", sub(";.*$", "", risk_field[direct_idx]))
    ))
    allele_result[!nzchar(allele_result)] = NA_character_
  }
  
  remaining_idx = which(valid_idx & is.na(allele_result) & grepl(";", risk_field, fixed = TRUE))
  if (length(remaining_idx) == 0) {
    return(allele_result)
  }
  
  allele_result[remaining_idx] = mapply(function(field, current_rsid) {
    if (is.na(field) || !nzchar(field) || is.na(current_rsid) || !nzchar(current_rsid)) {
      return(NA_character_)
    }
    
    parts = unlist(strsplit(field, ";", fixed = TRUE), use.names = FALSE)
    parts = trimws(parts)
    current_part = parts[tolower(sub("-.*$", "", parts)) == current_rsid][1]
    
    if (is.na(current_part) || !grepl("-", current_part, fixed = TRUE)) {
      return(NA_character_)
    }
    
    allele = sub("^.*-", "", current_part)
    allele = toupper(gsub("[^ACGT]", "", allele))
    
    if (!nzchar(allele)) {
      return(NA_character_)
    }
    
    allele
  }, risk_field[remaining_idx], rsid[remaining_idx], USE.NAMES = FALSE)
  
  allele_result
}

# Pastreaza doar asocierile GWAS pentru care genotipul contine alela de risc raportata.
filter_matches_by_risk_allele = function(df) {
  if (is.null(df) || nrow(df) == 0 || !"strongest_snp_risk_allele" %in% names(df)) {
    return(df)
  }
  
  rsid_col = if ("rsid_join" %in% names(df)) "rsid_join" else "rsid"
  if (!all(c(rsid_col, "genotype") %in% names(df))) {
    return(df)
  }
  
  df$risk_allele_match = extract_risk_allele_for_rsid(df$strongest_snp_risk_allele, df[[rsid_col]])
  df$genotype_clean = toupper(gsub("[^ACGT]", "", as.character(df$genotype)))
  has_risk_allele = mapply(
    function(allele, genotype_value) {
      if (is.na(allele) || !nzchar(allele)) {
        return(TRUE)
      }
      
      grepl(allele, genotype_value, fixed = TRUE)
    },
    df$risk_allele_match,
    df$genotype_clean,
    USE.NAMES = FALSE
  )
  
  df = df %>%
    filter(has_risk_allele) %>%
    select(-any_of(c("risk_allele_match", "genotype_clean")))
  
  df
}

# Pregateste o singura data alela de risc in tabelul GWAS.
# Astfel evitam recalcularea ei dupa fiecare join si in fiecare tabel/grafic.
prepare_gwas_ref_for_join = function(ref_df) {
  if (is.null(ref_df) || nrow(ref_df) == 0) {
    return(ref_df)
  }
  
  if (!isTRUE(enable_risk_allele_filter)) {
    return(ref_df)
  }
  
  if ("risk_allele_match" %in% names(ref_df)) {
    return(ref_df)
  }
  
  if (!all(c("rsid", "strongest_snp_risk_allele") %in% names(ref_df))) {
    return(ref_df)
  }
  
  ref_df$risk_allele_match = extract_risk_allele_for_rsid(
    ref_df$strongest_snp_risk_allele,
    ref_df$rsid
  )
  
  ref_df
}

# Varianta rapida folosita dupa ce alela a fost deja extrasa in ref.
filter_matches_by_prepared_risk_allele = function(df) {
  if (!isTRUE(enable_risk_allele_filter)) {
    return(df)
  }
  
  if (is.null(df) || nrow(df) == 0 || !"risk_allele_match" %in% names(df)) {
    return(filter_matches_by_risk_allele(df))
  }
  
  if (!"genotype" %in% names(df)) {
    return(df)
  }
  
  genotype_clean = toupper(gsub("[^ACGT]", "", as.character(df$genotype)))
  allele = as.character(df$risk_allele_match)
  
  keep_rows = mapply(
    function(current_allele, current_genotype) {
      if (is.na(current_allele) || !nzchar(current_allele)) {
        return(TRUE)
      }
      
      grepl(current_allele, current_genotype, fixed = TRUE)
    },
    allele,
    genotype_clean,
    USE.NAMES = FALSE
  )
  
  df %>%
    filter(keep_rows) %>%
    select(-any_of(c("risk_allele_match", "genotype_clean")))
}

# Reduce rândurile repetate pentru același SNP, păstrând asocierea cu p-value cel mai mic.
collapse_matches_by_rsid = function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  rsid_col = if ("rsid_join" %in% names(df)) "rsid_join" else "rsid"
  group_cols = intersect(c("Sursa", "Sursă", "sample_id", rsid_col, "genotype"), names(df))
  
  if (length(group_cols) == 0 || !"p_value" %in% names(df)) {
    return(df %>% distinct())
  }
  
  df %>%
    mutate(
      p_value_numeric_sort = suppressWarnings(as.numeric(p_value))
    ) %>%
    arrange(
      is.na(p_value_numeric_sort),
      p_value_numeric_sort
    ) %>%
    group_by(across(all_of(group_cols))) %>%
    slice(1) %>%
    ungroup() %>%
    select(-any_of("p_value_numeric_sort"))
}



# BAZA DE DATE SQLite


get_db_connection = function() {
  dbConnect(SQLite(), dbname = db_file_path)
}

# Adaugă o coloană nouă doar dacă baza de date este mai veche și nu o are deja.
ensure_column_exists = function(conn, table_name, column_name, column_definition) {
  table_info = dbGetQuery(conn, paste0("PRAGMA table_info(", table_name, ")"))
  if (!(column_name %in% table_info$name)) {
    dbExecute(conn, paste0("ALTER TABLE ", table_name, " ADD COLUMN ", column_name, " ", column_definition))
  }
}

table_exists = function(conn, table_name) {
  result = dbGetQuery(
    conn,
    "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ?",
    params = list(table_name)
  )
  nrow(result) > 0
}

read_text_file_content = function(file_path) {
  paste(readLines(file_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

read_gwas_cache = function(variant_ids) {
  # Cache-ul evită apelurile repetate către API pentru aceleași rsID-uri.
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  variant_ids = unique(tolower(trimws(as.character(variant_ids))))
  variant_ids = variant_ids[grepl("^rs[0-9]+$", variant_ids)]
  
  if (length(variant_ids) == 0 || !table_exists(conn, "gwas_annotations_cache")) {
    return(empty_gwas_ref())
  }
  # Interogările sunt împărțite în bucăți pentru a evita eroarea SQLite "too many SQL variables".
  chunks = split_into_chunks(variant_ids, gwas_cache_lookup_chunk_size)
  
  results = lapply(chunks, function(chunk) {
    placeholders = paste(rep("?", length(chunk)), collapse = ",")
    query = paste0(
      "SELECT rsid, categoria, disease_trait, mapped_trait, strongest_snp_risk_allele, p_value, mapped_gene, study, pubmedid ",
      "FROM gwas_annotations_cache WHERE rsid IN (", placeholders, ")"
    )
    dbGetQuery(conn, query, params = as.list(chunk))
  })
  
  
  bind_rows(results)
}

read_gwas_cached_rsids = function(variant_ids) {
  # Reținem și rsID-urile căutate fără rezultat, ca să nu le trimitem repetat la API.
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  variant_ids = unique(tolower(trimws(as.character(variant_ids))))
  variant_ids = variant_ids[grepl("^rs[0-9]+$", variant_ids)]
  
  if (length(variant_ids) == 0 || !table_exists(conn, "gwas_cached_rsids")) {
    return(character())
  }
  chunks = split_into_chunks(variant_ids, gwas_cache_lookup_chunk_size)
  
  results = lapply(chunks, function(chunk) {
    placeholders = paste(rep("?", length(chunk)), collapse = ",")
    query = paste0("SELECT rsid FROM gwas_cached_rsids WHERE rsid IN (", placeholders, ")")
    dbGetQuery(conn, query, params = as.list(chunk))
  })
  
  unique(bind_rows(results)$rsid)
}

save_gwas_requested_rsids = function(requested_rsids) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  requested_rsids = unique(tolower(trimws(as.character(requested_rsids))))
  requested_rsids = requested_rsids[grepl("^rs[0-9]+$", requested_rsids)]
  
  if (length(requested_rsids) == 0 || !table_exists(conn, "gwas_cached_rsids")) {
    return(invisible(NULL))
  }
  
  marker_chunks = split_into_chunks(requested_rsids, gwas_cache_lookup_chunk_size)
  for (chunk in marker_chunks) {
    values_sql = paste(rep("(?)", length(chunk)), collapse = ",")
    query = paste0("INSERT OR IGNORE INTO gwas_cached_rsids (rsid) VALUES ", values_sql)
    dbExecute(conn, query, params = as.list(chunk))
  }
  
  invisible(NULL)
}

save_gwas_cache = function(ref_df, requested_rsids) {
  # Salvarea se face într-o tranzacție, ca datele să rămână coerente dacă apare o eroare.
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbBegin(conn)
  tryCatch({
    if (nrow(ref_df) > 0) {
      ref_df = unique(ref_df)
      row_chunks = split_into_chunks(seq_len(nrow(ref_df)), 1000)
      for (chunk_ids in row_chunks) {
        dbWriteTable(conn, "gwas_annotations_cache", ref_df[chunk_ids, , drop = FALSE], append = TRUE)
      }
    }
    
    requested_rsids = unique(tolower(trimws(as.character(requested_rsids))))
    requested_rsids = requested_rsids[grepl("^rs[0-9]+$", requested_rsids)]
    
    if (length(requested_rsids) > 0) {
      marker_chunks = split_into_chunks(requested_rsids, gwas_cache_lookup_chunk_size)
      for (chunk in marker_chunks) {
        values_sql = paste(rep("(?)", length(chunk)), collapse = ",")
        query = paste0("INSERT OR IGNORE INTO gwas_cached_rsids (rsid) VALUES ", values_sql)
        dbExecute(conn, query, params = as.list(chunk))
      }
    }
    
    dbCommit(conn)
  }, error = function(e) {
    dbRollback(conn)
    stop(e)
  })
}

get_gwas_data_cached = function(variant_ids, progress_callback = NULL) {
  # Varianta simplă cache + API, folosită intern pentru rsID-urile care chiar ajung la API.
  variant_ids = unique(tolower(trimws(as.character(variant_ids))))
  variant_ids = variant_ids[grepl("^rs[0-9]+$", variant_ids)]
  
  if (length(variant_ids) == 0) {
    return(empty_gwas_ref())
  }
  
  if (!is.null(progress_callback)) {
    progress_callback(1, 4, "Verific cache-ul local din SQLite")
  }
  
  cached_done = read_gwas_cached_rsids(variant_ids)
  missing_rsids = setdiff(variant_ids, cached_done)
  cached_ref = read_gwas_cache(variant_ids)
  
  fetched_ref = empty_gwas_ref()
  
  if (length(missing_rsids) > 0) {
    if (!is.null(progress_callback)) {
      progress_callback(2, 4, paste("Trimit către gwasrapidd", length(missing_rsids), "rsid-uri noi"))
    }
    
    fetched_ref = read_gwas_data_online(missing_rsids) %>% distinct()
    
    if (!is.null(progress_callback)) {
      progress_callback(3, 4, "Salvez rezultatele API în cache")
    }
    
    save_gwas_cache(fetched_ref, missing_rsids)
  }
  
  if (!is.null(progress_callback)) {
    progress_callback(4, 4, "Pregătesc rezultatele pentru interfață")
  }
  
  bind_rows(cached_ref, fetched_ref) %>%
    filter(rsid %in% variant_ids) %>%
    distinct()
}

# Obține date GWAS în ordinea corectă.
# Dacă TSV-ul local există: cache SQLite -> fișier local -> API doar pentru lipsuri.
# Dacă TSV-ul local lipsește: cache SQLite -> API limitat.
get_gwas_data_smart = function(rsids, progress_callback = NULL) {
  # Funcția principală a fluxului GWAS, folosită atât local, cât și online.
  notify_progress = function(step, total, message) {
    if (!is.null(progress_callback)) {
      progress_callback(step, total, message)
    }
  }
  
  # 1. Normalizare rsid-uri.
  rsids = unique(tolower(trimws(as.character(rsids))))
  rsids = rsids[grepl("^rs[0-9]+$", rsids)]
  
  if (length(rsids) == 0) {
    return(empty_gwas_ref())
  }
  
  notify_progress(1, 5, "Se caută în cache")
  cached_ref = read_gwas_cache(rsids)
  cached_done = read_gwas_cached_rsids(rsids)
  cached_with_results = unique(cached_ref$rsid)
  
  has_local_gwas = !is.na(resolve_gwas_file_path())
  local_ref = empty_gwas_ref()
  local_found = character()
  rsids_for_local = setdiff(rsids, cached_with_results)
  
  if (has_local_gwas) {
    # Doar rsID-urile care nu au rezultat în cache sunt căutate în TSV-ul local.
    notify_progress(2, 5, "Se verifică datele locale")
    local_ref = tryCatch(
      read_gwas_data_local(rsids_for_local),
      error = function(e) empty_gwas_ref()
    )
    local_found = unique(local_ref$rsid)
  } else {
    notify_progress(2, 5, "GWAS local indisponibil; se folosește API-ul ca fallback")
  }
  
  # API-ul este fallback: primește doar rsid-uri negăsite în cache sau local.
  missing_for_api = setdiff(rsids, union(cached_done, local_found))
  api_limit = gwas_api_max_new_rsids
  api_message = if (has_local_gwas) {
    "Se interoghează API doar pentru rsid-urile lipsă"
  } else {
    paste("Se interoghează API ca fallback, maximum", api_limit, "rsid-uri")
  }
  notify_progress(3, 5, api_message)
  
  # Limita protejează aplicația de interogări foarte mari.
  api_rsids = head(missing_for_api, api_limit)
  api_ref = empty_gwas_ref()
  
  if (length(api_rsids) > 0) {
    api_ref = get_gwas_data_cached(api_rsids)
  }
  
  notify_progress(4, 5, "Se combină rezultatele")
  # Datele din cache, local și API sunt aduse la aceeași structură înainte de combinare.
  final_ref = bind_rows(
    normalize_gwas_ref(cached_ref),
    normalize_gwas_ref(local_ref),
    normalize_gwas_ref(api_ref)
  ) %>%
    filter(rsid %in% rsids) %>%
    distinct()
  
  notify_progress(5, 5, "Analiza este finalizată")
  final_ref
}

migrate_legacy_tables = function(conn) {
  # Migrare pentru bazele de date create cu versiuni mai vechi ale aplicației.
  if (!(table_exists(conn, "uploads") && table_exists(conn, "sample_variants") && table_exists(conn, "matched_gwas"))) {
    return(invisible(NULL))
  }
  
  samples_count = dbGetQuery(conn, "SELECT COUNT(*) AS n FROM samples")$n[1]
  if (!is.na(samples_count) && samples_count > 0) {
    return(invisible(NULL))
  }
  
  dbExecute(conn, "
    INSERT INTO samples (
      sample_id,
      username,
      sample_name,
      sex,
      snp_count,
      original_file_name,
      original_file_content,
      created_at
    )
    SELECT
      upload_id,
      username,
      sample_name,
      sex,
      snp_count,
      original_file_name,
      original_file_content,
      created_at
    FROM uploads
  ")
  
  dbExecute(conn, "
    INSERT INTO variants (
      variant_id,
      sample_id,
      rsid,
      chrom,
      genotype
    )
    SELECT
      variant_id,
      upload_id,
      rsid,
      chrom,
      genotype
    FROM sample_variants
  ")
  
  dbExecute(conn, "
    INSERT INTO annotations (
      annotation_id,
      sample_id,
      rsid,
      chrom,
      genotype,
      categoria,
      disease_trait,
      mapped_trait,
      strongest_snp_risk_allele,
      p_value,
      mapped_gene,
      study,
      pubmedid
    )
    SELECT
      match_id,
      upload_id,
      rsid,
      chrom,
      genotype,
      categoria,
      disease_trait,
      mapped_trait,
      strongest_snp_risk_allele,
      p_value,
      mapped_gene,
      study,
      pubmedid
    FROM matched_gwas
  ")
}

# Creează tabelele necesare și contul admin implicit.
initialize_database = function() {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  # Tabela users păstrează conturile și rolurile; parolele sunt salvate ca hash sodium.
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS users (
      user_id INTEGER PRIMARY KEY AUTOINCREMENT,
      username TEXT NOT NULL UNIQUE,
      password TEXT NOT NULL,
      role TEXT NOT NULL DEFAULT 'user',
      created_at TEXT NOT NULL
    )
  ")
  
  ensure_column_exists(conn, "users", "role", "TEXT NOT NULL DEFAULT 'user'")
  
  # Tabela samples păstrează fiecare fișier încărcat și conținutul original.
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS samples (
      sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
      username TEXT NOT NULL,
      sample_name TEXT NOT NULL,
      sex TEXT,
      snp_count INTEGER,
      original_file_name TEXT,
      original_file_content TEXT,
      created_at TEXT NOT NULL
    )
  ")
  
  ensure_column_exists(conn, "samples", "original_file_name", "TEXT")
  ensure_column_exists(conn, "samples", "original_file_content", "TEXT")
  
  # Tabela variants conține SNP-urile extrase din fiecare fișier ADN.
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS variants (
      variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
      sample_id INTEGER NOT NULL,
      rsid TEXT NOT NULL,
      chrom TEXT,
      genotype TEXT,
      FOREIGN KEY(sample_id) REFERENCES samples(sample_id)
    )
  ")
  
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS annotations (
      annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
      sample_id INTEGER NOT NULL,
      rsid TEXT NOT NULL,
      chrom TEXT,
      genotype TEXT,
      categoria TEXT,
      disease_trait TEXT,
      mapped_trait TEXT,
      strongest_snp_risk_allele TEXT,
      p_value TEXT,
      mapped_gene TEXT,
      study TEXT,
      pubmedid TEXT,
      FOREIGN KEY(sample_id) REFERENCES samples(sample_id)
    )
  ")
  
  # Cache GWAS: rezultate API salvate local pentru a accelera executările următoare.
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS gwas_annotations_cache (
      cache_id INTEGER PRIMARY KEY AUTOINCREMENT,
      rsid TEXT NOT NULL,
      categoria TEXT,
      disease_trait TEXT,
      mapped_trait TEXT,
      strongest_snp_risk_allele TEXT,
      p_value TEXT,
      mapped_gene TEXT,
      study TEXT,
      pubmedid TEXT
    )
  ")
  
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS gwas_cached_rsids (
      rsid TEXT PRIMARY KEY
    )
  ")
  
  # Cereri prin care un utilizator poate solicita rol de admin.
  dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS admin_requests (
      request_id INTEGER PRIMARY KEY AUTOINCREMENT,
      username TEXT NOT NULL,
      status TEXT NOT NULL DEFAULT 'pending',
      requested_at TEXT NOT NULL,
      resolved_at TEXT,
      resolved_by TEXT
    )
  ")
  
  ensure_column_exists(conn, "admin_requests", "seen_by_user", "INTEGER NOT NULL DEFAULT 0")
  
  dbExecute(conn, "CREATE INDEX IF NOT EXISTS idx_variants_sample_id ON variants(sample_id)")
  dbExecute(conn, "CREATE INDEX IF NOT EXISTS idx_annotations_sample_id ON annotations(sample_id)")
  dbExecute(conn, "CREATE INDEX IF NOT EXISTS idx_gwas_annotations_cache_rsid ON gwas_annotations_cache(rsid)")
  dbExecute(conn, "CREATE INDEX IF NOT EXISTS idx_admin_requests_status ON admin_requests(status)")
  
  migrate_legacy_tables(conn)
  
  # La prima rulare se creează automat contul admin implicit.
  existing_user = dbGetQuery(
    conn,
    "SELECT username FROM users WHERE username = ?",
    params = list(default_app_username)
  )
  
  if (nrow(existing_user) == 0) {
    dbExecute(
      conn,
      "INSERT INTO users (username, password, role, created_at) VALUES (?, ?, ?, ?)",
      params = list(default_app_username, hash_password(default_app_password), "admin", as.character(Sys.time()))
    )
  } else {
    dbExecute(
      conn,
      "UPDATE users SET role = 'admin' WHERE username = ?",
      params = list(default_app_username)
    )
    
    existing_admin = dbGetQuery(
      conn,
      "SELECT password FROM users WHERE username = ?",
      params = list(default_app_username)
    )
    
    if (nrow(existing_admin) == 1 && !password_is_hashed(existing_admin$password[1])) {
      dbExecute(
        conn,
        "UPDATE users SET password = ? WHERE username = ?",
        params = list(hash_password(existing_admin$password[1]), default_app_username)
      )
    }
  }
}

# Verifică login-ul și întoarce datele utilizatorului.
get_user_account = function(username, password) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  user_row = dbGetQuery(
    conn,
    "SELECT username, role, password FROM users WHERE username = ?",
    params = list(username)
  )
  
  if (nrow(user_row) != 1) {
    return(user_row[0, c("username", "role"), drop = FALSE])
  }
  
  if (!verify_password(user_row$password[1], password)) {
    return(user_row[0, c("username", "role"), drop = FALSE])
  }
  
  if (!password_is_hashed(user_row$password[1])) {
    dbExecute(
      conn,
      "UPDATE users SET password = ? WHERE username = ?",
      params = list(hash_password(password), username)
    )
  }
  
  user_row[, c("username", "role"), drop = FALSE]
}

# Creează un cont nou pentru utilizator.
create_user_account = function(username, password) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  username = trimws(username)
  
  if (!nzchar(username)) {
    stop("Utilizatorul nu poate fi gol.")
  }
  
  if (nchar(password) < 6) {
    stop("Parola trebuie să aibă cel puțin 6 caractere.")
  }
  
  existing_user = dbGetQuery(
    conn,
    "SELECT username FROM users WHERE username = ?",
    params = list(username)
  )
  
  if (nrow(existing_user) > 0) {
    stop("Exista deja un cont cu acest utilizator.")
  }
  
  dbExecute(
    conn,
    "INSERT INTO users (username, password, role, created_at) VALUES (?, ?, ?, ?)",
    params = list(username, hash_password(password), "user", as.character(Sys.time()))
  )
  
  TRUE
}

# Trimite o cerere de promovare la rolul de admin.
request_admin_role = function(username) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  role_row = dbGetQuery(
    conn,
    "SELECT role FROM users WHERE username = ?",
    params = list(username)
  )
  
  if (nrow(role_row) != 1) {
    stop("Utilizatorul nu există.")
  }
  
  if (identical(role_row$role[1], "admin")) {
    stop("Acest cont are deja rol de admin.")
  }
  
  existing_request = dbGetQuery(
    conn,
    "SELECT request_id FROM admin_requests WHERE username = ? AND status = 'pending'",
    params = list(username)
  )
  
  if (nrow(existing_request) > 0) {
    stop("Există deja o cerere de admin în așteptare.")
  }
  
  dbExecute(
    conn,
    "INSERT INTO admin_requests (username, status, requested_at) VALUES (?, 'pending', ?)",
    params = list(username, as.character(Sys.time()))
  )
  
  TRUE
}

# Citește cererile de admin care așteaptă aprobare.
read_pending_admin_requests = function() {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbGetQuery(
    conn,
    "
      SELECT
        request_id,
        username,
        status,
        requested_at
      FROM admin_requests
      WHERE status = 'pending'
      ORDER BY datetime(requested_at) ASC
    "
  )
}

# Verifică dacă userul curent are deja o cerere de admin în așteptare.
get_admin_request_status = function(username) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  request_row = dbGetQuery(
    conn,
    "
      SELECT status
      FROM admin_requests
      WHERE username = ?
      ORDER BY datetime(requested_at) DESC
      LIMIT 1
    ",
    params = list(username)
  )
  
  if (nrow(request_row) == 0) {
    return(NA_character_)
  }
  
  request_row$status[1]
}

# Citește ultimul răspuns nevăzut al utilizatorului pentru cereri/rol admin.
get_unseen_admin_notice = function(username) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  notice_row = dbGetQuery(
    conn,
    "
      SELECT request_id, status
      FROM admin_requests
      WHERE username = ?
        AND status IN ('approved', 'rejected', 'revoked')
        AND COALESCE(seen_by_user, 0) = 0
      ORDER BY datetime(COALESCE(resolved_at, requested_at)) DESC
      LIMIT 1
    ",
    params = list(username)
  )
  
  if (nrow(notice_row) == 0) {
    return(NULL)
  }
  
  notice_row
}

# Marchează un răspuns admin ca văzut de utilizator.
mark_admin_notice_seen = function(request_id) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbExecute(
    conn,
    "UPDATE admin_requests SET seen_by_user = 1 WHERE request_id = ?",
    params = list(request_id)
  )
  
  TRUE
}

# Aprobă sau respinge o cerere de admin.
resolve_admin_request = function(request_id, admin_username, decision = c("approved", "rejected")) {
  decision = match.arg(decision)
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbBegin(conn)
  tryCatch({
    request_row = dbGetQuery(
      conn,
      "SELECT request_id, username, status FROM admin_requests WHERE request_id = ?",
      params = list(request_id)
    )
    
    if (nrow(request_row) != 1 || !identical(request_row$status[1], "pending")) {
      stop("Cererea selectată nu mai este disponibilă.")
    }
    
    if (identical(decision, "approved")) {
      dbExecute(
        conn,
        "UPDATE users SET role = 'admin' WHERE username = ?",
        params = list(request_row$username[1])
      )
    }
    
    dbExecute(
      conn,
      "
        UPDATE admin_requests
        SET status = ?, resolved_at = ?, resolved_by = ?, seen_by_user = 0
        WHERE request_id = ?
      ",
      params = list(decision, as.character(Sys.time()), admin_username, request_id)
    )
    
    dbCommit(conn)
    TRUE
  }, error = function(e) {
    dbRollback(conn)
    stop(e)
  })
}

# Citește adminii promovați care pot fi revocați de adminul principal.
read_promoted_admins = function() {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbGetQuery(
    conn,
    "
      SELECT
        username,
        role,
        created_at
      FROM users
      WHERE role = 'admin' AND username <> ?
      ORDER BY username ASC
    ",
    params = list(default_app_username)
  )
}

# Revocă rolul de admin pentru un utilizator promovat.
revoke_admin_role = function(username, admin_username) {
  if (!is_primary_admin(admin_username)) {
    stop("Doar adminul principal poate revoca roluri admin.")
  }
  
  if (is_primary_admin(username)) {
    stop("Rolul adminului principal nu poate fi revocat.")
  }
  
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  user_row = dbGetQuery(
    conn,
    "SELECT username, role FROM users WHERE username = ?",
    params = list(username)
  )
  
  if (nrow(user_row) != 1 || !identical(user_row$role[1], "admin")) {
    stop("Utilizatorul selectat nu este admin.")
  }
  
  dbExecute(
    conn,
    "UPDATE users SET role = 'user' WHERE username = ?",
    params = list(username)
  )
  
  dbExecute(
    conn,
    "
      INSERT INTO admin_requests (username, status, requested_at, resolved_at, resolved_by, seen_by_user)
      VALUES (?, 'revoked', ?, ?, ?, 0)
    ",
    params = list(username, as.character(Sys.time()), as.character(Sys.time()), admin_username)
  )
  
  TRUE
}

# Salvează upload-ul, fișierul original și rezultatele analizei în SQLite.
save_upload_to_db = function(username, sample_name, sample_df, sex_value, matched_df, original_file_name, original_file_content) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbBegin(conn)
  tryCatch({
    dbExecute(
      conn,
      "INSERT INTO samples (username, sample_name, sex, snp_count, original_file_name, original_file_content, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
      params = list(username, sample_name, sex_value, nrow(sample_df), original_file_name, original_file_content, as.character(Sys.time()))
    )
    
    sample_id = dbGetQuery(conn, "SELECT last_insert_rowid() AS sample_id")$sample_id[1]
    
    variants_to_save = sample_df %>%
      transmute(
        sample_id = sample_id,
        rsid = rsid,
        chrom = chrom,
        genotype = genotype
      )
    
    if (nrow(variants_to_save) > 0) {
      dbWriteTable(conn, "variants", variants_to_save, append = TRUE)
    }
    
    matched_to_save = matched_df %>%
      transmute(
        sample_id = sample_id,
        rsid = rsid_join,
        chrom = chrom,
        genotype = genotype,
        categoria = categoria,
        disease_trait = disease_trait,
        mapped_trait = mapped_trait,
        strongest_snp_risk_allele = strongest_snp_risk_allele,
        p_value = p_value,
        mapped_gene = mapped_gene,
        study = study,
        pubmedid = pubmedid
      )
    
    if (nrow(matched_to_save) > 0) {
      dbWriteTable(conn, "annotations", matched_to_save, append = TRUE)
    }
    
    dbCommit(conn)
    sample_id
  }, error = function(e) {
    dbRollback(conn)
    stop(e)
  })
}

# Verifică duplicatul în funcție de rol.
# Userul simplu verifică doar în propriul cont; adminul verifică global, deoarece vede toate upload-urile.
upload_exists_for_scope = function(username, role, original_file_name, original_file_content) {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    existing_upload = dbGetQuery(
      conn,
      "
        SELECT sample_id
        FROM samples
        WHERE original_file_name = ? AND original_file_content = ?
        LIMIT 1
      ",
      params = list(original_file_name, original_file_content)
    )
    
    return(nrow(existing_upload) > 0)
  }
  
  existing_upload = dbGetQuery(
    conn,
    "
      SELECT sample_id
      FROM samples
      WHERE username = ? AND original_file_name = ? AND original_file_content = ?
      LIMIT 1
    ",
    params = list(username, original_file_name, original_file_content)
  )
  
  nrow(existing_upload) > 0
}

# Citește istoricul upload-urilor. Adminul vede tot, userul doar ce a încărcat el.
read_user_upload_history = function(username, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          username,
          sample_name,
          sex,
          snp_count,
          original_file_name,
          created_at
        FROM samples
        ORDER BY datetime(created_at) DESC
      "
    )
  } else {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          sample_name,
          sex,
          snp_count,
          original_file_name,
          created_at
        FROM samples
        WHERE username = ?
        ORDER BY datetime(created_at) DESC
      ",
      params = list(username)
    )
  }
}

# Întoarce fișierul original salvat în baza de date.
get_upload_file_for_user = function(username, upload_id, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          username,
          sample_name,
          original_file_name,
          original_file_content
        FROM samples
        WHERE sample_id = ?
      ",
      params = list(upload_id)
    )
  } else {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          sample_name,
          original_file_name,
          original_file_content
        FROM samples
        WHERE username = ? AND sample_id = ?
      ",
      params = list(username, upload_id)
    )
  }
}

# Citește variantele ADN salvate pentru un upload.
read_upload_variants = function(username, upload_id, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    dbGetQuery(
      conn,
      "
        SELECT
          rsid,
          chrom,
          genotype
        FROM variants
        WHERE sample_id = ?
        ORDER BY variant_id
      ",
      params = list(upload_id)
    )
  } else {
    dbGetQuery(
      conn,
      "
        SELECT
          rsid,
          chrom,
          genotype
        FROM variants
        WHERE sample_id IN (
          SELECT sample_id
          FROM samples
          WHERE username = ? AND sample_id = ?
        )
        ORDER BY variant_id
      ",
      params = list(username, upload_id)
    )
  }
}

# Citește potrivirile dintre upload și GWAS pentru un upload salvat.
read_upload_matches = function(username, upload_id, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    dbGetQuery(
      conn,
      "
        SELECT
          rsid,
          chrom,
          genotype,
          categoria,
          disease_trait,
          mapped_trait,
          strongest_snp_risk_allele,
          p_value,
          mapped_gene,
          study,
          pubmedid
        FROM annotations
        WHERE sample_id = ?
        ORDER BY annotation_id
      ",
      params = list(upload_id)
    )
  } else {
    dbGetQuery(
      conn,
      "
        SELECT
          rsid,
          chrom,
          genotype,
          categoria,
          disease_trait,
          mapped_trait,
          strongest_snp_risk_allele,
          p_value,
          mapped_gene,
          study,
          pubmedid
        FROM annotations
        WHERE sample_id IN (
          SELECT sample_id
          FROM samples
          WHERE username = ? AND sample_id = ?
        )
        ORDER BY annotation_id
      ",
      params = list(username, upload_id)
    )
  }
}

# Citește metadatele unui upload din baza de date.
read_upload_metadata = function(username, upload_id, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  if (identical(role, "admin")) {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          username,
          sample_name,
          sex,
          snp_count,
          created_at
        FROM samples
        WHERE sample_id = ?
      ",
      params = list(upload_id)
    )
  } else {
    dbGetQuery(
      conn,
      "
        SELECT
          sample_id AS upload_id,
          sample_name,
          sex,
          snp_count,
          created_at
        FROM samples
        WHERE username = ? AND sample_id = ?
      ",
      params = list(username, upload_id)
    )
  }
}

# Șterge complet un upload și toate datele asociate lui.
delete_upload_for_user = function(username, upload_id, role = "user") {
  conn = get_db_connection()
  on.exit(dbDisconnect(conn), add = TRUE)
  
  dbBegin(conn)
  tryCatch({
    if (identical(role, "admin")) {
      upload_row = dbGetQuery(
        conn,
        "SELECT sample_id FROM samples WHERE sample_id = ?",
        params = list(upload_id)
      )
    } else {
      upload_row = dbGetQuery(
        conn,
        "SELECT sample_id FROM samples WHERE username = ? AND sample_id = ?",
        params = list(username, upload_id)
      )
    }
    
    if (nrow(upload_row) == 0) {
      stop("Upload-ul selectat nu a fost găsit.")
    }
    
    dbExecute(conn, "DELETE FROM variants WHERE sample_id = ?", params = list(upload_id))
    dbExecute(conn, "DELETE FROM annotations WHERE sample_id = ?", params = list(upload_id))
    dbExecute(conn, "DELETE FROM samples WHERE sample_id = ?", params = list(upload_id))
    
    dbCommit(conn)
    TRUE
  }, error = function(e) {
    dbRollback(conn)
    stop(e)
  })
}


# =========================
# Raport PDF
# =========================

# Împarte un text lung pe mai multe rânduri pentru PDF.
wrap_pdf_text = function(x, width = 70) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return("")
  }
  
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

# Scrie un bloc de text într-o poziție fixă în pagina PDF.
draw_pdf_block = function(text, x = 0.07, y = 0.95, cex = 0.9, font = 1, width = 70, lineheight = 1.25) {
  grid.text(
    wrap_pdf_text(text, width = width),
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    just = c("left", "top"),
    gp = gpar(cex = cex, fontface = font, lineheight = lineheight)
  )
}

# Generează o paletă cu nume pentru ggplot, astfel încât scale_fill_manual să nu mai dea warning-uri.
make_named_palette = function(levels_values) {
  levels_values = as.character(levels_values)
  levels_values = levels_values[!is.na(levels_values) & nzchar(levels_values)]
  levels_values = unique(levels_values)
  
  if (length(levels_values) == 0) {
    return(character(0))
  }
  
  setNames(rep(pastel_colors, length.out = length(levels_values)), levels_values)
}


# Curăță numele fișierelor generate, ca download-ul PDF să nu fie transformat în .htm
# din cauza caracterelor speciale din numele probei.
safe_download_name = function(x) {
  x = as.character(x)
  x = gsub("[^A-Za-z0-9_-]+", "_", x)
  x = gsub("_+", "_", x)
  x = gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "profil" else x
}

# Se asigură că un data.frame are coloanele cerute pentru raport/grafice.
ensure_columns = function(df, cols) {
  if (is.null(df) || !is.data.frame(df)) {
    df = data.frame(stringsAsFactors = FALSE)
  }
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] = NA_character_
    }
  }
  df
}

# Normalizează valorile cromozomilor, indiferent dacă vin ca 1, chr1, chromosome 1, 23, X, MT etc.
normalize_chrom_value = function(chrom_values) {
  chrom_values = toupper(trimws(as.character(chrom_values)))
  chrom_values = gsub("^CHR", "", chrom_values)
  chrom_values = gsub("^CHROMOSOME", "", chrom_values)
  chrom_values = gsub("\\s+", "", chrom_values)
  chrom_values = gsub("\\.0$", "", chrom_values)
  chrom_values = dplyr::case_when(
    chrom_values == "23" ~ "X",
    chrom_values == "24" ~ "Y",
    chrom_values == "25" ~ "MT",
    chrom_values %in% c("M", "MTDNA", "MITO", "MITOCHONDRIAL") ~ "MT",
    chrom_values %in% c(as.character(1:22), "X", "Y", "MT") ~ chrom_values,
    TRUE ~ NA_character_
  )
  chrom_values
}

# Normalizează cromozomii pentru graficele din profil individual și PDF.
make_chrom_summary = function(df_raw) {
  df_raw = ensure_columns(df_raw, c("chrom", "chromosome"))
  
  chrom_values = if ("chrom" %in% names(df_raw) && any(!is.na(df_raw$chrom) & nzchar(as.character(df_raw$chrom)) & !toupper(as.character(df_raw$chrom)) %in% c("NA", "NULL"))) {
    df_raw$chrom
  } else {
    df_raw$chromosome
  }
  
  data.frame(chrom = normalize_chrom_value(chrom_values), stringsAsFactors = FALSE) %>%
    filter(!is.na(chrom)) %>%
    group_by(chrom) %>%
    summarise(Count = as.integer(dplyr::n()), .groups = "drop") %>%
    mutate(chrom = factor(chrom, levels = c(as.character(1:22), "X", "Y", "MT"))) %>%
    filter(!is.na(chrom)) %>%
    arrange(chrom)
}

# Scurtează numele lungi doar pentru afișare; valorile interne rămân numele complete.
short_file_label = function(x, max_chars = 36) {
  x = as.character(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value) || nchar(value) <= max_chars) {
      return(value)
    }
    
    ext = tools::file_ext(value)
    suffix = if (nzchar(ext)) paste0(".", ext) else ""
    available = max_chars - nchar(suffix) - 3
    
    if (available < 10) {
      return(paste0(substr(value, 1, max_chars - 3), "..."))
    }
    
    paste0(substr(value, 1, available), "...", suffix)
  }, character(1), USE.NAMES = FALSE)
}

make_labeled_choices = function(values, max_chars = 42) {
  values = as.character(values)
  setNames(values, short_file_label(values, max_chars))
}

# Pentru grafice folosim o etichetă foarte scurtă, dar păstrăm numele complet în tooltip.
make_plot_sample_label = function(x, max_chars = 24) {
  x = tools::file_path_sans_ext(as.character(x))
  extracted = sub("^(user[0-9]+_file[0-9]+).*$", "\\1", x, ignore.case = TRUE)
  unchanged = identical(length(extracted), length(x)) & extracted == x
  extracted[unchanged] = short_file_label(x[unchanged], max_chars)
  extracted
}

# Desenează aceeași diagrama Venn într-un fișier PNG pentru descărcare.
draw_venn_png = function(sets, file) {
  selected = short_file_label(names(sets), 28)
  grDevices::png(file, width = 1200, height = 800, res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  op = par(mar = c(1, 1, 3, 1))
  on.exit(par(op), add = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 8), axes = FALSE, xlab = "", ylab = "", asp = 1, main = "Diagramă Venn - asocieri rsID")
  
  if (length(sets) == 2) {
    only_1 = length(setdiff(sets[[1]], sets[[2]]))
    only_2 = length(setdiff(sets[[2]], sets[[1]]))
    common_12 = length(intersect(sets[[1]], sets[[2]]))
    
    symbols(3.7, 4, circles = 2.2, inches = FALSE, add = TRUE, bg = adjustcolor("#B3CDE3", alpha.f = 0.45), fg = "#2c3e50")
    symbols(6.1, 4, circles = 2.2, inches = FALSE, add = TRUE, bg = adjustcolor("#FBB4AE", alpha.f = 0.45), fg = "#2c3e50")
    text(2.8, 6.5, selected[1], cex = 0.9)
    text(7.2, 6.5, selected[2], cex = 0.9)
    text(2.7, 4, only_1, cex = 1.5, font = 2)
    text(7.1, 4, only_2, cex = 1.5, font = 2)
    text(4.9, 4, common_12, cex = 1.5, font = 2)
  } else {
    only_1 = length(setdiff(sets[[1]], union(sets[[2]], sets[[3]])))
    only_2 = length(setdiff(sets[[2]], union(sets[[1]], sets[[3]])))
    only_3 = length(setdiff(sets[[3]], union(sets[[1]], sets[[2]])))
    common_12 = length(setdiff(intersect(sets[[1]], sets[[2]]), sets[[3]]))
    common_13 = length(setdiff(intersect(sets[[1]], sets[[3]]), sets[[2]]))
    common_23 = length(setdiff(intersect(sets[[2]], sets[[3]]), sets[[1]]))
    common_123 = length(Reduce(intersect, sets))
    
    symbols(3.8, 5.0, circles = 2.1, inches = FALSE, add = TRUE, bg = adjustcolor("#B3CDE3", alpha.f = 0.40), fg = "#2c3e50")
    symbols(6.2, 5.0, circles = 2.1, inches = FALSE, add = TRUE, bg = adjustcolor("#FBB4AE", alpha.f = 0.40), fg = "#2c3e50")
    symbols(5.0, 3.2, circles = 2.1, inches = FALSE, add = TRUE, bg = adjustcolor("#CCEBC5", alpha.f = 0.40), fg = "#2c3e50")
    text(2.4, 7.0, selected[1], cex = 0.85)
    text(7.6, 7.0, selected[2], cex = 0.85)
    text(5.0, 0.8, selected[3], cex = 0.85)
    text(3.0, 5.4, only_1, cex = 1.3, font = 2)
    text(7.0, 5.4, only_2, cex = 1.3, font = 2)
    text(5.0, 2.0, only_3, cex = 1.3, font = 2)
    text(5.0, 5.6, common_12, cex = 1.3, font = 2)
    text(4.0, 3.8, common_13, cex = 1.3, font = 2)
    text(6.0, 3.8, common_23, cex = 1.3, font = 2)
    text(5.0, 4.2, common_123, cex = 1.5, font = 2)
  }
}

ui = fluidPage(
  # Interfața are două zone: panoul de autentificare și aplicația principală.
  useShinyjs(),
  theme = shinytheme("flatly"),
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      html {
        font-size: 14px !important;
      }
      body {
        font-size: 14px !important;
        overflow-x: auto;
      }
      .container-fluid {
        width: 96%;
        max-width: 1920px;
      }
      #login_panel {
        max-width: 420px;
        margin: 40px auto;
        padding: 24px;
        background: #ffffff;
        border: 1px solid #d9e2ec;
        border-radius: 12px;
      }
      .form-control,
      .selectize-input,
      .btn,
      .nav-tabs > li > a,
      .dataTables_wrapper,
      .sidebar,
      table,
      .dataTable,
      input,
      label {
        font-size: 13px !important;
      }
      .btn-lg {
        min-height: 38px !important;
        font-size: 14px !important;
      }
      .shiny-input-container,
      input[type='file'] {
        max-width: 100%;
      }
      .btn-file {
        width: 100%;
        text-align: left;
      }
      .tab-content {
        padding-top: 8px;
      }
      #login_btn {
        margin-top: -6px;
      }
      .login-message {
        min-height: 18px;
        margin-top: 6px;
        color: #2c3e50;
      }
      .login-message:empty {
        min-height: 0;
        margin-top: 0;
      }
      .sidebar-status {
        max-height: 95px;
        overflow-y: auto;
        white-space: pre-wrap;
        word-break: break-word;
        margin-bottom: 6px;
      }
      body.sidebar-hidden .sidebar-column {
        display: none;
      }
      body.sidebar-hidden .main-column {
        width: 100%;
      }
      #hide_sidebar_btn,
      #show_sidebar_btn {
        display: none;
        position: fixed;
        left: 18px;
        bottom: 18px;
        z-index: 9998;
        box-shadow: 0 4px 12px rgba(0,0,0,0.18);
      }
      body.sidebar-hidden #show_sidebar_btn {
        display: block;
      }
      body.wide-tab:not(.sidebar-hidden) #hide_sidebar_btn {
        display: block;
      }
      .interpretation-table,
      .interpretation-table .dataTables_wrapper {
        width: 100%;
      }
      @media (max-width: 992px) {
        html,
        body {
          width: 100%;
          overflow-x: auto !important;
          -webkit-overflow-scrolling: touch;
        }
        html {
          font-size: 13px !important;
        }
        .container-fluid {
          width: max-content;
          min-width: 980px;
          padding-left: 10px;
          padding-right: 10px;
        }
        #app_content {
          min-width: 980px;
        }
        #login_panel {
          width: calc(100vw - 24px);
          max-width: none;
          margin: 12px auto;
          padding: 16px;
        }
        .nav-tabs {
          display: flex;
          flex-wrap: nowrap;
          overflow-x: auto;
          overflow-y: hidden;
          white-space: nowrap;
        }
        .nav-tabs > li {
          float: none;
          display: inline-block;
        }
        .sidebar-column,
        .main-column {
          float: left;
        }
        .sidebar-column {
          width: 300px !important;
        }
        .main-column {
          width: 650px !important;
        }
        body.sidebar-hidden .main-column {
          width: 960px !important;
        }
        #hide_sidebar_btn,
        #show_sidebar_btn {
          left: 10px;
          bottom: 10px;
          padding: 8px 10px;
          font-size: 12px !important;
        }
        .dataTables_wrapper {
          overflow-x: auto;
        }
      }
      @media (max-width: 992px) and (orientation: portrait) {
        .container-fluid,
        #app_content {
          min-width: 760px;
        }
        .sidebar-column {
          width: 260px !important;
        }
        .main-column {
          width: 480px !important;
        }
        body.sidebar-hidden .main-column {
          width: 740px !important;
        }
      }
      #loading-content {
        position: fixed;
        bottom: 20px;
        right: 20px;
        z-index: 9999;
        background: rgba(44,62,80,0.9);
        color: white;
        padding: 15px;
        border-radius: 10px;
        display: none;
      }
      #shiny-notification-panel {
        top: auto !important;
        right: 24px !important;
        bottom: 35px !important;
        left: auto !important;
        z-index: 10000 !important;
        transform: none;
      }
      .shiny-notification {
        max-width: 420px;
        white-space: normal;
        word-break: break-word;
      }
      .selectize-control.multi .selectize-input > div,
      .selectize-dropdown .option {
        max-width: 100%;
        white-space: normal;
        word-break: break-word;
      }
      .selectize-control.single .selectize-input {
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
        padding-right: 38px !important;
      }
      .selectize-control.single .selectize-input > div {
        max-width: calc(100% - 22px);
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      #startup_splash {
        position: fixed;
        inset: 0;
        z-index: 20000;
        display: flex;
        align-items: center;
        justify-content: center;
        background: linear-gradient(135deg, #f8fbfd 0%, #e9f1f5 100%);
        transition: opacity 0.7s ease, visibility 0.7s ease;
      }
      #startup_splash.splash-hidden {
        opacity: 0;
        visibility: hidden;
      }
      .startup-splash-card {
        text-align: center;
        padding: 34px 42px;
        background: rgba(255, 255, 255, 0.86);
        border: 1px solid #d9e2ec;
        border-radius: 22px;
        box-shadow: 0 18px 44px rgba(44, 62, 80, 0.16);
      }
      .startup-splash-logo {
        max-width: 360px;
        width: 58vw;
        height: auto;
        animation: splash-logo-in 0.8s ease both;
      }
      .startup-splash-text {
        margin-top: 18px;
        color: #2c3e50;
        font-size: 18px !important;
        font-weight: 600;
      }
      @keyframes splash-logo-in {
        from {
          opacity: 0;
          transform: translateY(10px) scale(0.96);
        }
        to {
          opacity: 1;
          transform: translateY(0) scale(1);
        }
      }
    "))
  ),
  tags$script(HTML("
    $(document).on('shiny:connected', function() {
      setTimeout(function() {
        $('#startup_splash').addClass('splash-hidden');
        setTimeout(function() {
          $('#startup_splash').remove();
        }, 800);
      }, 4200);
    });
  ")),
  tags$div(
    id = "startup_splash",
    tags$div(
      class = "startup-splash-card",
      tags$img(src = "logo_aplicatie.png", class = "startup-splash-logo", alt = "Logo aplicație"),
      #tags$div(class = "startup-splash-text", "Analiza profilului genetic")
    )
  ),
  titlePanel("🧬 Analiza profilului genetic"),
  tags$div(
    id = "login_panel",
    h3("🔐 Autentificare"),
    textInput("login_username", "Utilizator", value = ""),
    passwordInput("login_password", "Parolă"),
    actionButton("login_btn", "Intră în aplicație", class = "btn-primary", width = "100%"),
    tags$div(class = "login-message", textOutput("login_status")),
    hr(),
    h4("👤 Creează cont"),
    textInput("register_username", "Utilizator nou"),
    passwordInput("register_password", "Parolă nouă"),
    passwordInput("register_password_confirm", "Confirmă parola"),
    actionButton("register_btn", "Creează cont", width = "100%"),
    tags$div(class = "login-message", textOutput("register_status"))
  ),
  tags$div(
    id = "app_content",
    style = "display:none;",
    sidebarLayout(
      sidebarPanel(
        class = "sidebar-column",
        textOutput("current_user_text"),
        uiOutput("admin_request_status_ui"),
        hr(),
        fileInput(
          "adn_files",
          "⛓️ Încarcă fișiere ADN",
          multiple = TRUE
        ),
        tags$small("Acceptă 23andMe, Ancestry, FTDNA, CSV, VCF și fișiere .txt în format VCF."),
        hr(),
        tags$p("Fișierele se compară cu GWAS Catalog."),
        hr(),
        actionButton("process_btn", "🧠 ANALIZEAZĂ", class = "btn-primary btn-lg", width = "100%"),
        br(),
        div(
          class = "sidebar-status",
          textOutput("status_text")
        ),
        br(),
        actionButton("logout_btn", "🔓 Logout", width = "100%"),
        br(),
        br(),
        uiOutput("admin_request_action_ui")
      ),
      mainPanel(
        class = "main-column",
        tabsetPanel(
          id = "main_tabs",
          tabPanel(
            "📊 Statistici", value = "stats",
            br(),
            withSpinner(tableOutput("status_table"), type = 6, color = "#2c3e50"),
            hr(),
            withSpinner(plotlyOutput("stats_bar_plot", height = "380px"), type = 6, color = "#2c3e50")
          ),
          tabPanel(
            "🔗 Interpretare genetică", value = "interpretare",
            br(),
            downloadButton("download_main", "Descarcă tabel CSV"),
            hr(),
            tags$div(
              class = "interpretation-table",
              withSpinner(DTOutput("clinicalTable"), type = 6, color = "#2c3e50")
            )
          ),
          tabPanel(
            "🪪 Profil individual", value = "profil",
            br(),
            selectInput("selected_sample", "Selectează individul", choices = NULL),
            hr(),
            tableOutput("individual_summary"),
            downloadButton("download_individual", "Descarcă profil selectat (CSV)"),
            downloadButton("download_individual_pdf", "Descarcă raport PDF"),
            hr(),
            withSpinner(plotlyOutput("individual_pie_cat", height = "500px"), type = 6, color = "#2c3e50"),
            br(),
            withSpinner(plotlyOutput("risk_bar_plot", height = "400px"), type = 6, color = "#2c3e50"),
            hr(),
            withSpinner(plotlyOutput("chrom_plot", height = "400px"), type = 6, color = "#2c3e50"),
            br()
          ),
          tabPanel(
            "🔁 Comparare fișiere", value = "comparare",
            br(),
            selectizeInput(
              "compare_samples",
              "Selectează 2 sau 3 fișiere pentru comparație",
              choices = character(0),
              multiple = TRUE,
              options = list(
                placeholder = "Alege fișierele de comparat...",
                maxItems = 3
              )
            ),
            tags$small("Poți compara fișiere încărcate acum sau analize încărcate din baza de date, atâta timp cât sunt prezente în sesiunea curentă."),
            hr(),
            tableOutput("compare_summary"),
            br(),
            withSpinner(plotlyOutput("compare_venn_plot", height = "500px"), type = 6, color = "#2c3e50"),
            br(),
            downloadButton("download_compare_venn", "Descarcă diagrama Venn")
          ),
          tabPanel(
            "🗂️ Istoric", value = "istoric",
            br(),
            DTOutput("history_table"),
            uiOutput("admin_requests_panel"),
            uiOutput("revoke_admin_panel"),
            br(),
            fluidRow(
              column(8, uiOutput("history_file_ui")),
              column(
                4,
                tags$div(style = "margin-top: 25px;",
                         actionButton("load_history_upload", "Încarcă analiza din DB", class = "btn-primary", width = "100%"),
                         br(), br(),
                         actionButton("reset_session_results", "🔄 Resetează rezultatele afișate", class = "btn-warning", width = "100%"),
                         br(), br(),
                         downloadButton("download_history_file", "Descarcă fișierul original din DB", width = "100%"),
                         br(), br(),
                         uiOutput("admin_delete_ui")
                )
              )
            ),
            br(),
            
          )
        )
      )
    )
  ),
  tags$div(
    id = "loading-content",
    tags$i(class = "fa fa-spinner fa-spin", style = "font-size:24px; margin-right:10px;"),
    "Procesare fișiere"
  ),
  actionButton("show_sidebar_btn", "☰ Arată sidebar", class = "btn-info"),
  actionButton("hide_sidebar_btn", "☰ Ascunde sidebar", class = "btn-info")
)



# SERVER
server = function(input, output, session) {
  initialize_database()
  
  # Valori reactive principale: datele analizate, utilizatorul curent și istoricul upload-urilor.
  app_data <- reactiveVal(NULL)
  logged_in <- reactiveVal(FALSE)
  current_user <- reactiveVal(NULL)
  current_role <- reactiveVal("user")
  history_refresh <- reactiveVal(0)
  admin_request_refresh <- reactiveVal(0)
  login_status <- reactiveVal("")
  register_status <- reactiveVal("")
  status_text <- reactiveVal("Încarcă fișierele ADN, apoi apasă ANALIZEAZĂ.")
  
  # Afișare mai clară pentru sexul detectat, fără să modificăm valorile interne folosite la grafice/PDF.
  add_sex_icon = function(x) {
    x = as.character(x)
    dplyr::case_when(
      x == "Masculin" ~ "♂️ Masculin",
      x == "Feminin" ~ "♀️ Feminin",
      TRUE ~ x
    )
  }
  
  
  # Actualizează selectorul pentru compararea fișierelor pe baza fișierelor prezente în sesiune.
  refresh_compare_samples = function(selected = NULL) {
    current_data = app_data()
    
    if (is.null(current_data) || is.null(current_data$adn) || length(current_data$adn) == 0) {
      updateSelectizeInput(
        session,
        "compare_samples",
        choices = character(0),
        selected = character(0),
        server = FALSE
      )
      return(invisible(NULL))
    }
    
    choices = names(current_data$adn)
    selected = selected[selected %in% choices]
    
    updateSelectizeInput(
      session,
      "compare_samples",
      choices = make_labeled_choices(choices, 44),
      selected = selected,
      server = FALSE
    )
  }
  
  # Păstrează în selector atât fișierele încărcate acum, cât și analizele încărcate din baza de date.
  merge_app_data = function(new_data) {
    old_data = app_data()
    
    if (is.null(old_data) || is.null(old_data$adn) || length(old_data$adn) == 0) {
      return(new_data)
    }
    
    combined_adn = c(old_data$adn, new_data$adn)
    unique_names = make.unique(names(combined_adn), sep = " #")
    names(combined_adn) = unique_names
    
    old_stats = old_data$stats
    new_stats = new_data$stats
    combined_stats = bind_rows(old_stats, new_stats)
    
    if (nrow(combined_stats) == length(unique_names)) {
      combined_stats$Fișier = unique_names
    }
    
    old_ref = if (!is.null(old_data$ref)) old_data$ref else empty_gwas_ref()
    new_ref = if (!is.null(new_data$ref)) new_data$ref else empty_gwas_ref()
    
    old_ref = old_ref %>% mutate(across(everything(), as.character))
    new_ref = new_ref %>% mutate(across(everything(), as.character))
    
    combined_ref = bind_rows(old_ref, new_ref) %>% distinct()
    
    list(
      adn = combined_adn,
      ref = combined_ref,
      stats = combined_stats
    )
  }
  
  # Rezultatele NU se mai șterg automat când intri în Istoric sau când alegi alte fișiere.
  # Ștergerea rezultatelor afișate se face manual, din butonul de reset din tabul Istoric.
  observeEvent(input$reset_session_results, {
    app_data(NULL)
    shinyjs::reset("adn_files")
    updateSelectInput(session, "selected_sample", choices = character(0), selected = character(0))
    refresh_compare_samples(selected = character(0))
    updateSelectizeInput(session, "selected_history_upload", selected = character(0))
    status_text("Rezultatele afișate au fost resetate. Poți încărca fișiere noi sau poți alege alte analize din Istoric.")
  }, ignoreInit = TRUE)
  
  observe({
    # Afișează login-ul sau aplicația, în funcție de starea autentificării.
    if (isTRUE(logged_in())) {
      shinyjs::hide("login_panel")
      shinyjs::show("app_content")
    } else {
      shinyjs::show("login_panel")
      shinyjs::hide("app_content")
      shinyjs::runjs("$('body').removeClass('wide-tab sidebar-hidden');")
    }
  })
  
  observeEvent(input$main_tabs, {
    # Pe taburile cu tabele mari ascundem temporar sidebar-ul pentru mai mult spațiu.
    if (input$main_tabs %in% c("interpretare", "istoric")) {
      shinyjs::runjs("$('body').addClass('wide-tab').removeClass('sidebar-hidden');")
    } else {
      shinyjs::runjs("$('body').removeClass('wide-tab sidebar-hidden');")
    }
    shinyjs::runjs("setTimeout(function(){ $.fn.dataTable.tables({visible:true, api:true}).columns.adjust().draw(false); }, 250);")
  }, ignoreInit = FALSE)
  
  observeEvent(input$show_sidebar_btn, {
    # Buton manual pentru readucerea sidebar-ului fără schimbarea tabului curent.
    shinyjs::runjs("$('body').removeClass('sidebar-hidden');")
    shinyjs::runjs("setTimeout(function(){ $.fn.dataTable.tables({visible:true, api:true}).columns.adjust().draw(false); }, 250);")
  })
  
  observeEvent(input$hide_sidebar_btn, {
    # Buton manual pentru ascunderea sidebar-ului pe taburile cu tabele mari.
    shinyjs::runjs("$('body').addClass('sidebar-hidden');")
    shinyjs::runjs("setTimeout(function(){ $.fn.dataTable.tables({visible:true, api:true}).columns.adjust().draw(false); }, 250);")
  })
  
  observeEvent(input$login_btn, {
    # Verifică datele introduse și reține rolul utilizatorului în sesiunea curentă.
    req(input$login_username, input$login_password)
    
    user_row = get_user_account(trimws(input$login_username), input$login_password)
    
    if (nrow(user_row) == 1) {
      shinyjs::reset("adn_files")
      updateTabsetPanel(session, "main_tabs", selected = "stats")
      current_user(trimws(input$login_username))
      current_role(user_row$role[1])
      logged_in(TRUE)
      login_status("")
      
      admin_notice = tryCatch(
        get_unseen_admin_notice(trimws(input$login_username)),
        error = function(e) NULL
      )
      
      if (!is.null(admin_notice)) {
        notice_message = dplyr::case_when(
          identical(admin_notice$status[1], "approved") ~ "Cererea pentru rol admin a fost aprobată.",
          identical(admin_notice$status[1], "rejected") ~ "Cererea pentru rol admin a fost respinsă.",
          identical(admin_notice$status[1], "revoked") ~ "Rolul admin a fost revocat.",
          TRUE ~ ""
        )
        
        if (nzchar(notice_message)) {
          status_text(notice_message)
        }
        
        mark_admin_notice_seen(admin_notice$request_id[1])
      }
    } else {
      login_status("Utilizator sau parolă invalidă.")
      current_role("user")
      logged_in(FALSE)
    }
  })
  
  observeEvent(input$logout_btn, {
    # Logout-ul curăță doar sesiunea afișată; datele salvate în SQLite rămân intacte.
    logged_in(FALSE)
    current_user(NULL)
    current_role("user")
    app_data(NULL)
    shinyjs::reset("adn_files")
    shinyjs::runjs("$('body').removeClass('wide-tab sidebar-hidden');")
    updateTabsetPanel(session, "main_tabs", selected = "stats")
    updateSelectInput(session, "selected_sample", choices = character(0), selected = character(0))
    refresh_compare_samples(selected = character(0))
    login_status("")
    status_text("Încarcă fișierele ADN, apoi apasă ANALIZEAZĂ.")
  })
  
  output$login_status = renderText({
    login_status()
  })
  
  observeEvent(input$register_btn, {
    # Conturile noi sunt salvate în SQLite, cu parola hash-uită.
    req(input$register_username, input$register_password, input$register_password_confirm)
    
    tryCatch({
      if (input$register_password != input$register_password_confirm) {
        stop("Parolele introduse nu coincid.")
      }
      
      create_user_account(
        username = input$register_username,
        password = input$register_password
      )
      
      updateTextInput(session, "login_username", value = trimws(input$register_username))
      updateTextInput(session, "register_username", value = "")
      updateTextInput(session, "register_password", value = "")
      updateTextInput(session, "register_password_confirm", value = "")
      
      register_status("Cont creat cu succes. Acum puteți completa câmpurile de logare.")
    }, error = function(e) {
      register_status(paste("Eroare creare cont:", conditionMessage(e)))
    })
  })
  
  output$register_status = renderText({
    register_status()
  })
  
  output$current_user_text = renderText({
    req(current_user())
    paste("Utilizator:", current_user(), "| Rol:", current_role())
  })
  
  output$admin_request_status_ui = renderUI({
    req(logged_in(), current_user())
    admin_request_refresh()
    
    if (identical(current_role(), "admin")) {
      admin_label = if (is_primary_admin(current_user())) {
        "Cont admin principal activ."
      } else {
        "Cont admin promovat activ."
      }
      
      return(tags$div(
        style = "margin: 8px 0 4px 0; color: #2c3e50;",
        tags$small(admin_label)
      ))
    }
    
    request_status = tryCatch(
      get_admin_request_status(current_user()),
      error = function(e) NA_character_
    )
    
    status_text_admin = dplyr::case_when(
      identical(request_status, "pending") ~ "Cererea pentru rol admin este în așteptare.",
      identical(request_status, "approved") ~ "Cererea pentru rol admin a fost aprobată.",
      identical(request_status, "rejected") ~ "Cererea pentru rol admin a fost respinsă.",
      TRUE ~ "Nu există cerere pentru rol admin."
    )
    
    tags$div(
      style = "margin: 8px 0 4px 0; color: #2c3e50;",
      tags$small(status_text_admin)
    )
  })
  
  output$admin_request_action_ui = renderUI({
    req(logged_in(), current_user())
    admin_request_refresh()
    
    if (identical(current_role(), "admin")) {
      return(NULL)
    }
    
    request_status = tryCatch(
      get_admin_request_status(current_user()),
      error = function(e) NA_character_
    )
    
    if (identical(request_status, "pending") || identical(request_status, "approved")) {
      return(NULL)
    }
    
    tags$div(
      style = "margin: 8px 0 4px 0;",
      actionButton("request_admin_btn", "Solicită rol admin", class = "btn-info", width = "100%")
    )
  })
  
  observeEvent(input$request_admin_btn, {
    req(logged_in(), current_user())
    
    tryCatch({
      request_admin_role(current_user())
      admin_request_refresh(admin_request_refresh() + 1)
      status_text("Cererea pentru rol admin a fost trimisă. Așteaptă aprobarea unui admin.")
    }, error = function(e) {
      status_text(paste("Cerere admin:", conditionMessage(e)))
    })
  })
  
  # Procesează fișierele încărcate, interoghează GWAS prin API și salvează rezultatele.
  observeEvent(input$process_btn, {
    # Fluxul principal: citire fișiere ADN -> căutare GWAS -> salvare în DB -> afișare rezultate.
    req(logged_in())
    req(input$adn_files)
    shinyjs::show("loading-content")
    status_text("⏳ Procesare în curs. Aplicația citește fișierele ADN și pregătește interogarea GWAS.")
    
    evaluation_metrics = start_evaluation_metrics(input$adn_files$name)
    
    result = tryCatch({
      adn_list = list()
      stats = data.frame(
        Fișier = character(),
        SNP = numeric(),
        Sex = character(),
        stringsAsFactors = FALSE
      )
      original_contents = list()
      invalid_files = character()
      
      for (i in seq_len(nrow(input$adn_files))) {
        # Fiecare fișier încărcat este convertit la structură standard rsid/chrom/genotype.
        fname = input$adn_files$name[i]
        res = read_adn_data(
          input$adn_files$datapath[i],
          original_file_name = input$adn_files$name[i]
        )
        
        if (!is.null(res)) {
          sample_data = res$data
          adn_list[[fname]] = sample_data
          original_contents[[fname]] = read_text_file_content(input$adn_files$datapath[i])
          
          sex_val = if (is.null(res$sex) || length(res$sex) == 0) "Necunoscut (markeri sexuali lipsă sau ambigui)" else res$sex
          
          stats = rbind(
            stats,
            data.frame(
              Fișier = fname,
              SNP = nrow(sample_data),
              Sex = sex_val,
              stringsAsFactors = FALSE
            )
          )
        } else {
          invalid_files = c(invalid_files, fname)
        }
      }
      
      if (length(invalid_files) > 0) {
        showNotification(
          paste(
            "Următoarele fișiere nu au putut fi citite ca fișiere ADN valide:",
            paste(short_file_label(invalid_files, 36), collapse = ", ")
          ),
          type = "warning",
          duration = 15
        )
      }
      
      if (length(adn_list) == 0) {
        stop("Niciun fișier ADN valid nu a putut fi procesat.")
      }
      
      # Userul simplu verifică duplicatele doar în contul lui; adminul verifică global.
      # Dacă utilizatorul încarcă simultan un fișier duplicat și unul nou,
      # se procesează mai departe doar fișierele noi.
      duplicate_files = character()
      for (fname in names(original_contents)) {
        if (upload_exists_for_scope(current_user(), current_role(), fname, original_contents[[fname]])) {
          duplicate_files = c(duplicate_files, fname)
        }
      }
      
      if (length(duplicate_files) > 0) {
        
        showNotification(
          paste(
            "Următoarele fișiere nu au fost încărcate deoarece există deja:",
            paste(duplicate_files, collapse = ", ")
          ),
          type = "warning",
          duration = 15
        )
        
        keep_files = setdiff(names(adn_list), duplicate_files)
        
        if (length(keep_files) == 0) {
          
          stop(paste(
            "Toate fișierele selectate există deja pentru nivelul de acces al contului curent:",
            paste(duplicate_files, collapse = ", ")
          ))
          
        }
        
        
        
        adn_list = adn_list[keep_files]
        original_contents = original_contents[keep_files]
        stats = stats %>% filter(Fișier %in% keep_files)
        
        status_text(paste(
          "ℹ️ Am ignorat fișierele deja existente:",
          paste(duplicate_files, collapse = ", "),
          ". Continui procesarea pentru fișierele noi:",
          paste(keep_files, collapse = ", ")
        ))
      }
      
      rsid_sets = lapply(adn_list, function(x) unique(tolower(trimws(as.character(x$rsid)))))
      rsid_scope_label = if (length(rsid_sets) == 1) "rsid-uri din fișier" else "rsid-uri comune între fișiere"
      # Pentru un singur fișier folosim toate rsID-urile; pentru mai multe fișiere folosim reuniunea.
      all_rsids = if (length(rsid_sets) == 1) {
        rsid_sets[[1]]
      } else {
        Reduce(union, rsid_sets)
      }
      if (length(rsid_sets) > 1) {
        rsid_scope_label = "rsid-uri totale din fisiere"
      }
      all_rsids = unique(all_rsids[grepl("^rs[0-9]+$", all_rsids)])
      
      if (length(all_rsids) == 0) {
        stop("Nu au fost găsite rsid-uri valide în fișierele încărcate.")
      }
      
      shinyjs::hide("loading-content")
      
      ref = withProgress(
        message = "🔎 Analiză GWAS",
        detail = "Se verifică datele locale",
        value = 0,
        session = session,
        {
          get_gwas_data_smart(
            all_rsids,
            progress_callback = function(i, n, detail = NULL) {
              # Progresul GWAS este afișat în interfață, ca utilizatorul să vadă etapa curentă.
              incProgress(1 / max(1, n), detail = detail, session = session)
              status_text(if (is.null(detail)) "⏳ Analiza GWAS este în curs." else paste("⏳", detail))
            }
          )
        }
      )
      
      status_text(paste("✅ Analiza GWAS este finalizată. Pregătesc alelele de risc pentru", nrow(ref), "asocieri."))
      
      ref = withProgress(
        message = "🧬 Pregătire alele de risc",
        detail = "Extrag alela de risc din rezultatele GWAS",
        value = 0,
        session = session,
        {
          incProgress(0.2, detail = "Normalizez câmpul STRONGEST SNP-RISK ALLELE", session = session)
          prepared_ref = prepare_gwas_ref_for_join(ref)
          incProgress(0.8, detail = "Alelele de risc au fost pregătite", session = session)
          prepared_ref
        }
      )
      ref_for_join = ref
      matched_by_file = list()
      
      withProgress(
        message = "🧬 Filtrare după alela de risc",
        detail = "Se verifică genotipurile față de alelele raportate în GWAS",
        value = 0,
        session = session,
        {
          for (i in seq_along(adn_list)) {
            # Salvăm separat variantele brute și potrivirile GWAS pentru fiecare fișier.
            fname = names(adn_list)[i]
            status_text(paste("⏳ Se filtrează după alela de risc pentru", fname))
            incProgress(1 / max(1, length(adn_list)), detail = paste("Se procesează", fname), session = session)
            
            sample_data = adn_list[[i]]
            sample_data$rsid_join = tolower(trimws(sample_data$rsid))
            matched_data = inner_join(
              sample_data,
              ref_for_join,
              by = c("rsid_join" = "rsid"),
              relationship = "many-to-many"
            )
            
            matched_before_filter = nrow(matched_data)
            matched_data = filter_matches_by_prepared_risk_allele(matched_data)
            matched_after_filter = nrow(matched_data)
            matched_data = collapse_matches_by_rsid(matched_data)
            adn_list[[fname]] = sample_data
            matched_by_file[[fname]] = matched_data
            
            status_text(paste(
              "⏳ Filtrare alelă de risc:", fname,
              "- păstrate", matched_after_filter,
              "din", matched_before_filter,
              "asocieri. Se salvează în SQLite."
            ))
            
            sex_val = stats$Sex[stats[[1]] == fname][1]
            
            save_upload_to_db(
              username = current_user(),
              sample_name = fname,
              sample_df = sample_data,
              sex_value = sex_val,
              matched_df = matched_data,
              original_file_name = fname,
              original_file_content = original_contents[[fname]]
            )
          }
        }
      )
      
      new_session_data = list(adn = adn_list, ref = ref, stats = stats)
      app_data(merge_app_data(new_session_data))
      session_names = names(app_data()$adn)
      updateSelectInput(session, "selected_sample", choices = make_labeled_choices(session_names, 46), selected = tail(session_names, 1))
      refresh_compare_samples(selected = head(session_names, 2))
      history_refresh(history_refresh() + 1)
      final_message = "✅ Analiza este finalizată. Rezultatele GWAS au fost combinate și salvate în SQLite."
      if (length(invalid_files) > 0) {
        final_message = paste(
          final_message,
          "Fișiere ignorate deoarece nu au fost recunoscute:",
          paste(short_file_label(invalid_files, 36), collapse = ", ")
        )
      }
      status_text(final_message)
      TRUE
    }, error = function(e) {
      status_text(paste("❌ Eroare:", conditionMessage(e)))
      FALSE
    })
    
    shinyjs::hide("loading-content")
    finish_evaluation_metrics(evaluation_metrics, result)
    invisible(result)
  })
  
  output$status_text = renderText({
    status_text()
  })
  
  output$status_table = renderTable({
    req(app_data())
    df_stats = app_data()$stats
    if ("Fișier" %in% names(df_stats)) {
      df_stats$Fișier = short_file_label(df_stats$Fișier, 44)
    }
    if ("Sex" %in% names(df_stats)) {
      df_stats$Sex = add_sex_icon(df_stats$Sex)
    }
    df_stats
  })
  
  output$stats_bar_plot = renderPlotly({
    req(app_data())
    
    df_stats = app_data()$stats
    df_stats$Eticheta = make_plot_sample_label(df_stats[[names(df_stats)[1]]])
    df_stats$Tooltip = paste0(
      "Fișier: ", df_stats[[names(df_stats)[1]]],
      "<br>SNP-uri: ", df_stats$SNP,
      "<br>Status profil: ", df_stats$Sex
    )
    
    p = ggplot(df_stats, aes(x = Eticheta, y = SNP, fill = Sex, text = Tooltip)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c(
        "Masculin" = "#B3CDE3",
        "Feminin" = "#FBB4AE",
        "Necunoscut (Lipsă markeri sexuali)" = "#FFFFCC",
        "Necunoscut (markeri sexuali lipsă sau ambigui)" = "#FFFFCC"
      )) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
        plot.margin = margin(10, 20, 45, 10)
      ) +
      labs(
        title = "Distribuția variantelor genetice și sexul detectat",
        x = "Eșantion",
        y = "Număr total SNP-uri",
        fill = "Status profil"
      )
    
    ggplotly(p, tooltip = "text") %>%
      layout(margin = list(b = 70))
  })
  
  output$dynamic_filtres = renderUI({
    NULL
  })
  
  observe({
    # Filtrul pe categorie a fost eliminat din interfa??.
    # Evit?m construirea unei liste mari de categorii GWAS dup? fiecare analiz?.
    NULL
  })
  
  filtred_report = reactive({
    # Combină variantele ADN cu adnotările GWAS și pregătește linkurile pentru tabel.
    req(app_data())
    combined = bind_rows(app_data()$adn, .id = "Sursa")
    combined$rsid_join = tolower(combined$rsid)
    ref_for_join = prepare_gwas_ref_for_join(app_data()$ref)
    report = inner_join(
      combined,
      ref_for_join,
      by = c("rsid_join" = "rsid"),
      relationship = "many-to-many"
    )
    report = filter_matches_by_prepared_risk_allele(report)
    report = collapse_matches_by_rsid(report)
    
    report$rsid = sprintf(
      '<a href="https://www.ebi.ac.uk/gwas/labs/variants/%s" target="_blank">%s</a>',
      report$rsid_join,
      report$rsid
    )
    
    report$pubmedid = ifelse(
      is.na(report$pubmedid) | report$pubmedid == "",
      "",
      sprintf(
        '<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank">%s</a>',
        report$pubmedid,
        report$pubmedid
      )
    )
    
    return(report)
  })
  
  output$clinicalTable = renderDT({
    # Tabelul principal afișează interpretările GWAS, cu filtrare și scroll orizontal.
    display_df = filtred_report() %>%
      select(-any_of(c("rsid_join", "categoria")))
    
    datatable(
      display_df,
      escape = FALSE,
      filter = "top",
      options = list(pageLength = 10, scrollX = TRUE, searchHighlight = TRUE, search = list(regex = TRUE, caseInsensitive = TRUE))
    )
  })
  
  indiv_data = reactive({
    # Datele pentru profilul individual sunt calculate doar pentru proba selectată.
    req(app_data(), input$selected_sample)
    adn = app_data()$adn[[input$selected_sample]]
    adn$rsid_join = tolower(trimws(as.character(adn$rsid)))
    ref_for_join = prepare_gwas_ref_for_join(app_data()$ref)
    df_merge = inner_join(adn, ref_for_join, by = c("rsid_join" = "rsid"))
    df_merge = filter_matches_by_prepared_risk_allele(df_merge)
    df_merge = collapse_matches_by_rsid(df_merge)
    return(df_merge)
  })
  
  output$individual_summary = renderTable({
    df = indiv_data()
    req(df)
    
    total_variante = nrow(df)
    total_gene = dplyr::n_distinct(df$mapped_gene)
    
    if (total_variante == 0) return(data.frame(Status = "Lipsă date"))
    
    data.frame(
      "Asocieri" = total_variante,
      "Gene" = total_gene
    )
  })
  
  output$risk_bar_plot = renderPlotly({
    df = indiv_data()
    req(df)
    
    df_stats = df %>%
      count(categoria, sort = TRUE) %>%
      slice_head(n = 15)
    
    p = ggplot(df_stats, aes(x = reorder(categoria, n), y = n, fill = n)) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "#CCEBC5", high = "#FBB4AE") +
      theme_minimal() +
      labs(
        title = "<b>Asocieri GWAS pe categorii</b>",
        x = "Categorie",
        y = "Număr asocieri"
      )
    
    ggplotly(p)
  })
  
  output$individual_pie_cat = renderPlotly({
    df_plot = indiv_data() %>%
      filter(!is.na(categoria), nzchar(as.character(categoria))) %>%
      count(categoria, sort = TRUE) %>%
      slice_head(n = 10)
    
    if (nrow(df_plot) == 0) {
      return(plot_ly() %>% layout(title = "Nu există categorii disponibile pentru acest profil."))
    }
    
    plot_ly(
      df_plot,
      labels = ~categoria,
      values = ~n,
      type = "pie",
      hole = 0.4,
      marker = list(colors = pastel_colors[seq_len(min(length(pastel_colors), nrow(df_plot)))])
    ) %>%
      layout(title = list(text = "<b>Repartizarea pe categorii</b>", y = 0.95), margin = list(t = 100))
  })
  
  output$chrom_plot = renderPlotly({
    req(app_data(), input$selected_sample)
    
    df_raw = app_data()$adn[[input$selected_sample]]
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      return(plot_ly() %>% layout(title = "Nu există date brute pentru cromozomi."))
    }
    
    chrom_summary = make_chrom_summary(df_raw)
    
    if (nrow(chrom_summary) == 0) {
      return(plot_ly() %>% layout(title = "Nu există cromozomi valizi pentru acest profil."))
    }
    
    chrom_summary$chrom_label = as.character(chrom_summary$chrom)
    
    plot_ly(
      chrom_summary,
      x = ~chrom_label,
      y = ~Count,
      type = "bar",
      marker = list(color = rep(pastel_colors, length.out = nrow(chrom_summary)))
    ) %>%
      layout(
        title = "Distribuția variantelor pe cromozomi",
        xaxis = list(title = "Cromozom", categoryorder = "array", categoryarray = c(as.character(1:22), "X", "Y", "MT")),
        yaxis = list(title = "Număr variante")
      )
  })
  
  output$download_main = downloadHandler(
    filename = function() {
      paste0("Analiza_genetica_", Sys.Date(), ".csv")
    },
    content = function(file) {
      clean_html_text = function(x) {
        x = as.character(x)
        gsub("<[^>]+>", "", x)
      }
      
      data_to_save = tryCatch({
        filtred_report()
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(data_to_save)) return(NULL)
      
      # Tabelul din interfață conține linkuri HTML; la export păstrăm doar textul curat.
      data_to_save$rsid = clean_html_text(data_to_save$rsid)
      if ("pubmedid" %in% names(data_to_save)) {
        data_to_save$pubmedid = clean_html_text(data_to_save$pubmedid)
      }
      
      if ("rsid_join" %in% names(data_to_save)) {
        data_to_save$rsid_join <- NULL
      }
      
      names(data_to_save)[names(data_to_save) %in% c("Sursă", "SursÄƒ", "Surs?")] = "Sursa"
      
      temp_csv = tempfile(fileext = ".csv")
      write.table(
        data_to_save,
        temp_csv,
        sep = ";",
        row.names = FALSE,
        na = "",
        quote = TRUE,
        qmethod = "double",
        fileEncoding = "UTF-8"
      )
      
      con = file(file, open = "wb")
      on.exit(close(con), add = TRUE)
      writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
      writeBin(readBin(temp_csv, what = "raw", n = file.info(temp_csv)$size), con)
    },
    contentType = "text/csv"
  )
  
  output$download_individual = downloadHandler(
    filename = function() { paste0("Profil_genetic_", input$selected_sample, "_", Sys.Date(), ".csv") },
    content = function(file) {
      d_indiv = indiv_data()
      d_indiv = d_indiv[, !(names(d_indiv) %in% "rsid_join")]
      con = file(file, open = "wb")
      on.exit(close(con), add = TRUE)
      writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
      write.table(d_indiv, con, sep = ";", row.names = FALSE, na = "", fileEncoding = "UTF-8")
    }
  )
  
  output$download_individual_pdf = downloadHandler(
    filename = function() {
      paste0("Raport_profil_", safe_download_name(input$selected_sample), "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # Raportul PDF folosește aceleași date ca profilul individual din interfață.
      req(app_data(), input$selected_sample)
      
      df = indiv_data()
      df = ensure_columns(df, c("rsid", "genotype", "categoria", "disease_trait", "mapped_trait", "mapped_gene", "p_value", "study"))
      
      sample_name = input$selected_sample
      sample_raw = ensure_columns(app_data()$adn[[sample_name]], c("chrom", "chromosome", "rsid", "genotype"))
      total_matches = nrow(df)
      total_genes = dplyr::n_distinct(df$mapped_gene)
      total_traits = dplyr::n_distinct(df$mapped_trait)
      sex_value = app_data()$stats %>%
        filter(Fișier == sample_name) %>%
        pull(Sex)
      
      if (length(sex_value) == 0) {
        sex_value = "Necunoscut"
      }
      
      top_traits = df %>%
        filter(!is.na(categoria), nzchar(as.character(categoria))) %>%
        count(categoria, sort = TRUE) %>%
        slice_head(n = 10)
      
      chrom_summary = make_chrom_summary(sample_raw)
      
      top_rows = df %>%
        select(rsid, genotype, disease_trait, mapped_trait, mapped_gene, p_value, study) %>%
        distinct() %>%
        slice_head(n = 12)
      
      grDevices::pdf(file, width = 11.69, height = 8.27)
      on.exit(grDevices::dev.off(), add = TRUE)
      
      # Pagina 1: rezumat
      grid.newpage()
      draw_pdf_block("Raport profil genetic", x = 0.07, y = 0.93, cex = 1.6, font = 2)
      draw_pdf_block(paste("Data raportului:", Sys.Date()), x = 0.07, y = 0.87, cex = 0.95)
      
      grid.lines(x = unit(c(0.07, 0.93), "npc"), y = unit(c(0.84, 0.84), "npc"), gp = gpar(col = "#b8c4d0", lwd = 1))
      
      draw_pdf_block("Rezumat profil", x = 0.07, y = 0.80, cex = 1.15, font = 2)
      draw_pdf_block(paste("Profil selectat:", sample_name), x = 0.07, y = 0.74, cex = 1.0, width = 55)
      draw_pdf_block(paste("Sex detectat:", sex_value[1]), x = 0.07, y = 0.68, cex = 1.0)
      draw_pdf_block(paste("Asocieri totale cu GWAS:", total_matches), x = 0.07, y = 0.62, cex = 1.0)
      draw_pdf_block(paste("Gene unice:", total_genes), x = 0.07, y = 0.56, cex = 1.0)
      draw_pdf_block(paste("Mapped traits unice:", total_traits), x = 0.07, y = 0.50, cex = 1.0)
      
      # Pagina 2: top categorii
      grid.newpage()
      draw_pdf_block("Top categorii GWAS", x = 0.07, y = 0.93, cex = 1.4, font = 2)
      draw_pdf_block(paste("Profil:", sample_name), x = 0.07, y = 0.87, cex = 0.95, width = 70)
      grid.lines(x = unit(c(0.07, 0.93), "npc"), y = unit(c(0.84, 0.84), "npc"), gp = gpar(col = "#b8c4d0", lwd = 1))
      
      if (nrow(top_traits) > 0) {
        top_traits_page <- top_traits %>% slice_head(n = 8)
        y_positions_traits = seq(0.78, 0.20, length.out = nrow(top_traits_page))
        
        for (i in seq_len(nrow(top_traits_page))) {
          trait_text = paste0(
            i, ". ", top_traits_page$categoria[i],
            "\nNumăr asocieri: ", top_traits_page$n[i]
          )
          draw_pdf_block(trait_text, x = 0.09, y = y_positions_traits[i], cex = 0.98, width = 90, lineheight = 1.45)
        }
      } else {
        draw_pdf_block("Nu există categorii disponibile.", x = 0.07, y = 0.78, cex = 0.95, width = 70)
      }
      
      # Pagina 3: asocieri
      grid.newpage()
      draw_pdf_block("Asocieri identificate", x = 0.07, y = 0.93, cex = 1.4, font = 2)
      draw_pdf_block(paste("Profil:", sample_name), x = 0.07, y = 0.87, cex = 0.95, width = 70)
      grid.lines(x = unit(c(0.07, 0.93), "npc"), y = unit(c(0.84, 0.84), "npc"), gp = gpar(col = "#b8c4d0", lwd = 1))
      
      if (nrow(top_rows) > 0) {
        y_positions = seq(0.79, 0.12, length.out = nrow(top_rows))
        
        for (i in seq_len(nrow(top_rows))) {
          row_text = paste0(
            i, ". ",
            top_rows$rsid[i], " | Genotip: ", top_rows$genotype[i],
            "\nGenă: ", ifelse(is.na(top_rows$mapped_gene[i]) || top_rows$mapped_gene[i] == "", "NA", top_rows$mapped_gene[i]),
            "\nTrait: ", ifelse(is.na(top_rows$mapped_trait[i]) || top_rows$mapped_trait[i] == "", "NA", top_rows$mapped_trait[i]),
            "\nP-value: ", ifelse(is.na(top_rows$p_value[i]) || top_rows$p_value[i] == "", "NA", top_rows$p_value[i])
          )
          draw_pdf_block(row_text, x = 0.07, y = y_positions[i], cex = 0.82, width = 95, lineheight = 1.3)
        }
      } else {
        draw_pdf_block("Nu există rânduri disponibile.", x = 0.07, y = 0.79, cex = 0.9, width = 70)
      }
      
      # Pagina 4: grafic categorii
      if (nrow(top_traits) > 0) {
        p1 = ggplot(top_traits, aes(x = reorder(categoria, n), y = n, fill = n)) +
          geom_col() +
          coord_flip() +
          scale_fill_gradient(low = "#CCEBC5", high = "#FBB4AE") +
          theme_minimal() +
          labs(
            title = "Top categorii GWAS",
            x = "Categorie",
            y = "Număr variante"
          )
        print(p1)
      }
      
      # Pagina 5: grafic cromozomi
      if (nrow(chrom_summary) > 0) {
        p2 = ggplot(chrom_summary, aes(x = chrom, y = Count, fill = chrom)) +
          geom_col(show.legend = FALSE) +
          scale_fill_manual(values = make_named_palette(chrom_summary$chrom)) +
          theme_minimal() +
          labs(
            title = "Distributia variantelor pe cromozomi",
            x = "Cromozom",
            y = "Număr variante"
          )
        print(p2)
      }
    },
    contentType = "application/pdf"
  )
  
  # Pregătește datele pentru compararea rsid-urilor între două sau trei fișiere din sesiunea curentă.
  compare_sets = reactive({
    # Pregătește seturile de rsID pentru diagrama Venn și tabelul de comparație.
    req(app_data())
    selected = input$compare_samples
    selected = selected[selected %in% names(app_data()$adn)]
    
    if (length(selected) < 2) {
      return(NULL)
    }
    
    if (length(selected) > 3) {
      return(list(error = "Pentru diagrama de tip Venn selectează maximum trei fișiere."))
    }
    
    sets = lapply(selected, function(sample_name) {
      df = app_data()$adn[[sample_name]]
      validate(need("rsid" %in% names(df), paste("Fișierul", sample_name, "nu are coloană rsid.")))
      rsids = unique(tolower(trimws(as.character(df$rsid))))
      rsids[grepl("^rs[0-9]+$", rsids)]
    })
    
    names(sets) = selected
    sets
  })
  
  output$compare_summary = renderTable({
    sets = compare_sets()
    validate(need(!is.null(sets), "Selectează cel puțin două fișiere pentru comparație."))
    validate(need(is.null(sets$error), sets$error))
    selected = names(sets)
    selected_display = short_file_label(selected, 34)
    common_all = Reduce(intersect, sets)
    union_all = Reduce(union, sets)
    
    summary_df = data.frame(
      Indicator = c(
        "Fișiere comparate",
        "rsID-uri comune în toate fișierele",
        "rsID-uri totale unice"
      ),
      Valoare = as.character(c(
        paste(selected_display, collapse = " | "),
        length(common_all),
        length(union_all)
      )),
      stringsAsFactors = FALSE
    )
    
    if (length(sets) == 2) {
      only_1 = length(setdiff(sets[[1]], sets[[2]]))
      only_2 = length(setdiff(sets[[2]], sets[[1]]))
      
      extra_df = data.frame(
        Indicator = c(
          paste("Doar", selected_display[1]),
          paste("Doar", selected_display[2])
        ),
        Valoare = as.character(c(only_1, only_2)),
        stringsAsFactors = FALSE
      )
      summary_df = bind_rows(summary_df, extra_df)
    }
    
    if (length(sets) == 3) {
      extra_df = data.frame(
        Indicator = c(
          paste("Doar", selected_display[1]),
          paste("Doar", selected_display[2]),
          paste("Doar", selected_display[3]),
          paste("Comune", selected_display[1], "+", selected_display[2]),
          paste("Comune", selected_display[1], "+", selected_display[3]),
          paste("Comune", selected_display[2], "+", selected_display[3])
        ),
        Valoare = as.character(c(
          length(setdiff(sets[[1]], union(sets[[2]], sets[[3]]))),
          length(setdiff(sets[[2]], union(sets[[1]], sets[[3]]))),
          length(setdiff(sets[[3]], union(sets[[1]], sets[[2]]))),
          length(setdiff(intersect(sets[[1]], sets[[2]]), sets[[3]])),
          length(setdiff(intersect(sets[[1]], sets[[3]]), sets[[2]])),
          length(setdiff(intersect(sets[[2]], sets[[3]]), sets[[1]]))
        )),
        stringsAsFactors = FALSE
      )
      summary_df = bind_rows(summary_df, extra_df)
    }
    
    summary_df %>% mutate(across(everything(), as.character))
  })
  
  output$compare_venn_plot = renderPlotly({
    sets = compare_sets()
    req(!is.null(sets), is.null(sets$error))
    selected = names(sets)
    selected_display = short_file_label(selected, 28)
    
    if (length(sets) == 2) {
      only_1 = length(setdiff(sets[[1]], sets[[2]]))
      only_2 = length(setdiff(sets[[2]], sets[[1]]))
      common_12 = length(intersect(sets[[1]], sets[[2]]))
      
      plot_ly(type = "scatter", mode = "text") %>%
        layout(
          title = "Diagramă Venn - asocieri rsID",
          xaxis = list(visible = FALSE, range = c(0, 10)),
          yaxis = list(visible = FALSE, range = c(0, 7), scaleanchor = "x"),
          shapes = list(
            list(type = "circle", x0 = 1.5, x1 = 6.0, y0 = 1.2, y1 = 5.7, fillcolor = "rgba(179,205,227,0.45)", line = list(color = "#2c3e50")),
            list(type = "circle", x0 = 4.0, x1 = 8.5, y0 = 1.2, y1 = 5.7, fillcolor = "rgba(251,180,174,0.45)", line = list(color = "#2c3e50"))
          ),
          annotations = list(
            list(x = 2.6, y = 5.9, text = selected_display[1], showarrow = FALSE),
            list(x = 7.4, y = 5.9, text = selected_display[2], showarrow = FALSE),
            list(x = 2.8, y = 3.4, text = as.character(only_1), showarrow = FALSE, font = list(size = 18)),
            list(x = 7.2, y = 3.4, text = as.character(only_2), showarrow = FALSE, font = list(size = 18)),
            list(x = 5.0, y = 3.4, text = as.character(common_12), showarrow = FALSE, font = list(size = 18))
          )
        )
    } else {
      only_1 = length(setdiff(sets[[1]], union(sets[[2]], sets[[3]])))
      only_2 = length(setdiff(sets[[2]], union(sets[[1]], sets[[3]])))
      only_3 = length(setdiff(sets[[3]], union(sets[[1]], sets[[2]])))
      common_12 = length(setdiff(intersect(sets[[1]], sets[[2]]), sets[[3]]))
      common_13 = length(setdiff(intersect(sets[[1]], sets[[3]]), sets[[2]]))
      common_23 = length(setdiff(intersect(sets[[2]], sets[[3]]), sets[[1]]))
      common_123 = length(Reduce(intersect, sets))
      
      plot_ly(type = "scatter", mode = "text") %>%
        layout(
          title = "Diagramă Venn - asocieri rsID",
          xaxis = list(visible = FALSE, range = c(0, 10)),
          yaxis = list(visible = FALSE, range = c(0, 8), scaleanchor = "x"),
          shapes = list(
            list(type = "circle", x0 = 1.6, x1 = 6.2, y0 = 2.4, y1 = 7.0, fillcolor = "rgba(179,205,227,0.40)", line = list(color = "#2c3e50")),
            list(type = "circle", x0 = 3.8, x1 = 8.4, y0 = 2.4, y1 = 7.0, fillcolor = "rgba(251,180,174,0.40)", line = list(color = "#2c3e50")),
            list(type = "circle", x0 = 2.7, x1 = 7.3, y0 = 0.7, y1 = 5.3, fillcolor = "rgba(204,235,197,0.40)", line = list(color = "#2c3e50"))
          ),
          annotations = list(
            list(x = 2.5, y = 7.2, text = selected_display[1], showarrow = FALSE),
            list(x = 7.5, y = 7.2, text = selected_display[2], showarrow = FALSE),
            list(x = 5.0, y = 0.5, text = selected_display[3], showarrow = FALSE),
            list(x = 2.8, y = 5.3, text = as.character(only_1), showarrow = FALSE, font = list(size = 17)),
            list(x = 7.2, y = 5.3, text = as.character(only_2), showarrow = FALSE, font = list(size = 17)),
            list(x = 5.0, y = 2.1, text = as.character(only_3), showarrow = FALSE, font = list(size = 17)),
            list(x = 5.0, y = 5.4, text = as.character(common_12), showarrow = FALSE, font = list(size = 17)),
            list(x = 3.8, y = 3.7, text = as.character(common_13), showarrow = FALSE, font = list(size = 17)),
            list(x = 6.2, y = 3.7, text = as.character(common_23), showarrow = FALSE, font = list(size = 17)),
            list(x = 5.0, y = 4.1, text = as.character(common_123), showarrow = FALSE, font = list(size = 19))
          )
        )
    }
  })
  
  output$download_compare_venn = downloadHandler(
    filename = function() {
      paste0("Diagrama_Venn_", Sys.Date(), ".png")
    },
    content = function(file) {
      sets = compare_sets()
      draw_venn_png(sets, file)
    },
    contentType = "image/png"
  )
  
  output$history_table = renderDT({
    # Istoricul citește upload-urile din SQLite; adminul vede toate fișierele.
    req(logged_in(), current_user())
    history_refresh()
    history_df = read_user_upload_history(current_user(), current_role())
    datatable(history_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$admin_requests_panel = renderUI({
    req(logged_in(), current_user())
    history_refresh()
    
    if (!is_primary_admin(current_user())) {
      return(NULL)
    }
    
    pending_requests = read_pending_admin_requests()
    
    if (nrow(pending_requests) == 0) {
      return(tags$div(
        style = "margin: 12px 0; color: #6c757d;",
        "Nu există cereri noi pentru rol admin."
      ))
    }
    
    tagList(
      hr(),
      h4("Cereri pentru rol admin"),
      DTOutput("admin_requests_table"),
      br(),
      fluidRow(
        column(
          8,
          selectInput(
            "selected_admin_request",
            "Selectează cererea",
            choices = setNames(
              pending_requests$request_id,
              paste0("#", pending_requests$request_id, " - ", pending_requests$username, " - ", pending_requests$requested_at)
            )
          )
        ),
        column(
          2,
          tags$div(
            style = "margin-top: 25px;",
            actionButton("approve_admin_request", "Aprobă", class = "btn-success", width = "100%")
          )
        ),
        column(
          2,
          tags$div(
            style = "margin-top: 25px;",
            actionButton("reject_admin_request", "Respinge", class = "btn-danger", width = "100%")
          )
        )
      )
    )
  })
  
  output$admin_requests_table = renderDT({
    req(logged_in(), current_user())
    req(is_primary_admin(current_user()))
    history_refresh()
    
    pending_requests = read_pending_admin_requests()
    datatable(pending_requests, options = list(pageLength = 5, scrollX = TRUE))
  })
  
  observeEvent(input$approve_admin_request, {
    req(logged_in(), current_user(), input$selected_admin_request)
    req(is_primary_admin(current_user()))
    
    tryCatch({
      resolve_admin_request(as.integer(input$selected_admin_request), current_user(), "approved")
      admin_request_refresh(admin_request_refresh() + 1)
      history_refresh(history_refresh() + 1)
      status_text("Cererea a fost aprobată. Utilizatorul are acum rol de admin.")
    }, error = function(e) {
      status_text(paste("Eroare aprobare admin:", conditionMessage(e)))
    })
  })
  
  observeEvent(input$reject_admin_request, {
    req(logged_in(), current_user(), input$selected_admin_request)
    req(is_primary_admin(current_user()))
    
    tryCatch({
      resolve_admin_request(as.integer(input$selected_admin_request), current_user(), "rejected")
      admin_request_refresh(admin_request_refresh() + 1)
      history_refresh(history_refresh() + 1)
      status_text("Cererea pentru rol admin a fost respinsă.")
    }, error = function(e) {
      status_text(paste("Eroare respingere admin:", conditionMessage(e)))
    })
  })
  
  output$revoke_admin_panel = renderUI({
    req(logged_in(), current_user())
    history_refresh()
    
    if (!is_primary_admin(current_user())) {
      return(NULL)
    }
    
    promoted_admins = read_promoted_admins()
    
    if (nrow(promoted_admins) == 0) {
      return(tags$div(
        style = "margin: 12px 0; color: #6c757d;",
        "Nu există admini promovați care pot fi revocați."
      ))
    }
    
    tagList(
      hr(),
      h4("Revocare rol admin"),
      DTOutput("promoted_admins_table"),
      br(),
      fluidRow(
        column(
          8,
          selectInput(
            "selected_revoke_admin",
            "Selectează adminul",
            choices = promoted_admins$username
          )
        ),
        column(
          4,
          tags$div(
            style = "margin-top: 25px;",
            actionButton("revoke_admin_btn", "Revocă rol admin", class = "btn-danger", width = "100%")
          )
        )
      )
    )
  })
  
  output$promoted_admins_table = renderDT({
    req(logged_in(), current_user())
    req(is_primary_admin(current_user()))
    history_refresh()
    
    promoted_admins = read_promoted_admins()
    datatable(promoted_admins, options = list(pageLength = 5, scrollX = TRUE))
  })
  
  observeEvent(input$revoke_admin_btn, {
    req(logged_in(), current_user(), input$selected_revoke_admin)
    req(is_primary_admin(current_user()))
    
    tryCatch({
      revoke_admin_role(input$selected_revoke_admin, current_user())
      admin_request_refresh(admin_request_refresh() + 1)
      history_refresh(history_refresh() + 1)
      status_text(paste("Rolul admin a fost revocat pentru utilizatorul", input$selected_revoke_admin))
    }, error = function(e) {
      status_text(paste("Eroare revocare admin:", conditionMessage(e)))
    })
  })
  
  output$history_file_ui = renderUI({
    req(logged_in(), current_user())
    history_refresh()
    history_df = read_user_upload_history(current_user(), current_role())
    req(nrow(history_df) > 0)
    
    if ("username" %in% names(history_df)) {
      choice_labels = paste0(history_df$upload_id, " - ", history_df$username, " - ", history_df$sample_name, " - ", history_df$created_at)
    } else {
      choice_labels = paste0(history_df$upload_id, " - ", history_df$sample_name, " - ", history_df$created_at)
    }
    
    choices = setNames(history_df$upload_id, short_file_label(choice_labels, 64))
    selectizeInput(
      "selected_history_upload",
      "Alege fișierele salvate",
      choices = choices,
      selected = character(0),
      multiple = TRUE,
      options = list(
        placeholder = "Alege unul sau mai multe fișiere...",
        plugins = list("remove_button")
      )
    )
  })
  
  output$admin_delete_ui = renderUI({
    req(logged_in(), current_user())
    history_refresh()
    
    if (!identical(current_role(), "admin")) {
      return(NULL)
    }
    
    actionButton("delete_history_upload", "🗑️ Șterge upload-ul selectat din DB", class = "btn-danger")
  })
  
  observeEvent(input$load_history_upload, {
    # Încarcă una sau mai multe analize salvate și le readuce în sesiunea curentă.
    req(logged_in(), current_user(), input$selected_history_upload)
    
    upload_ids = as.integer(input$selected_history_upload)
    upload_ids = upload_ids[!is.na(upload_ids)]
    req(length(upload_ids) > 0)
    
    adn_list = list()
    stats_list = list()
    ref_list = list()
    
    for (upload_id in upload_ids) {
      variants_df = read_upload_variants(current_user(), upload_id, current_role())
      matches_df = read_upload_matches(current_user(), upload_id, current_role())
      metadata_df = read_upload_metadata(current_user(), upload_id, current_role())
      
      if (nrow(metadata_df) == 1 && nrow(variants_df) > 0) {
        sample_name = metadata_df$sample_name[1]
        display_name = if (sample_name %in% names(adn_list)) paste0(sample_name, " #", upload_id) else sample_name
        
        adn_list[[display_name]] = variants_df
        stats_list[[display_name]] = data.frame(
          Fișier = display_name,
          SNP = nrow(variants_df),
          Sex = metadata_df$sex[1],
          stringsAsFactors = FALSE
        )
        
        if (nrow(matches_df) > 0) {
          ref_list[[display_name]] = matches_df %>%
            select(
              rsid,
              categoria,
              disease_trait,
              mapped_trait,
              strongest_snp_risk_allele,
              p_value,
              mapped_gene,
              study,
              pubmedid
            )
        }
      }
    }
    
    req(length(adn_list) > 0)
    
    stats = bind_rows(stats_list)
    ref_df = if (length(ref_list) > 0) bind_rows(ref_list) %>% distinct() else empty_gwas_ref()
    selected_sample = names(adn_list)[1]
    
    new_session_data = list(adn = adn_list, ref = ref_df, stats = stats)
    app_data(new_session_data)
    session_names = names(app_data()$adn)
    if (!selected_sample %in% session_names) {
      selected_sample = tail(session_names, 1)
    }
    updateSelectInput(session, "selected_sample", choices = make_labeled_choices(session_names, 46), selected = selected_sample)
    refresh_compare_samples(selected = character(0))
    updateSelectizeInput(session, "selected_history_upload", selected = character(0))
    
    if (length(upload_ids) == 1) {
      status_text(paste("Analiza pentru", selected_sample, "a fost încărcată din baza de date."))
    } else {
      status_text(paste(length(upload_ids), "analize au fost încărcate din baza de date."))
    }
  })
  
  observeEvent(input$delete_history_upload, {
    # Ștergerea este permisă aici doar pentru admin și elimină datele asociate din toate tabelele.
    req(logged_in(), current_user(), input$selected_history_upload)
    req(identical(current_role(), "admin"))
    
    upload_ids = as.integer(input$selected_history_upload)
    upload_ids = upload_ids[!is.na(upload_ids)]
    req(length(upload_ids) > 0)
    
    tryCatch({
      for (upload_id in upload_ids) {
        delete_upload_for_user(current_user(), upload_id, current_role())
      }
      
      app_data(NULL)
      history_refresh(history_refresh() + 1)
      updateSelectInput(session, "selected_sample", choices = character(0), selected = character(0))
      refresh_compare_samples(selected = character(0))
      status_text(paste(length(upload_ids), "upload-uri au fost șterse din baza de date."))
    }, error = function(e) {
      status_text(paste("Eroare la ștergere:", conditionMessage(e)))
    })
  })
  
  output$download_history_file = downloadHandler(
    filename = function() {
      req(logged_in(), current_user(), input$selected_history_upload)
      selected_upload_id = as.integer(input$selected_history_upload[1])
      file_row = get_upload_file_for_user(current_user(), selected_upload_id, current_role())
      req(nrow(file_row) == 1)
      
      if (!is.na(file_row$original_file_name[1]) && nzchar(file_row$original_file_name[1])) {
        file_row$original_file_name[1]
      } else {
        paste0("upload_", file_row$upload_id[1], ".txt")
      }
    },
    content = function(file) {
      req(logged_in(), current_user(), input$selected_history_upload)
      selected_upload_id = as.integer(input$selected_history_upload[1])
      file_row = get_upload_file_for_user(current_user(), selected_upload_id, current_role())
      req(nrow(file_row) == 1)
      
      if (is.na(file_row$original_file_content[1]) || !nzchar(file_row$original_file_content[1])) {
        stop("Fișierul original nu există în baza de date pentru acest upload.")
      }
      
      writeLines(file_row$original_file_content[1], con = file, useBytes = TRUE)
    },
    contentType = "text/plain"
  )
  
  output$history_matches_table = renderDT({
    req(logged_in(), current_user(), input$selected_history_upload)
    history_refresh()
    
    upload_ids = as.integer(input$selected_history_upload)
    upload_ids = upload_ids[!is.na(upload_ids)]
    req(length(upload_ids) > 0)
    
    matches_list = lapply(upload_ids, function(upload_id) {
      matches_df = read_upload_matches(current_user(), upload_id, current_role())
      metadata = read_upload_metadata(current_user(), upload_id, current_role())
      
      if (nrow(matches_df) == 0 || nrow(metadata) == 0) {
        return(NULL)
      }
      
      sample_label = metadata$sample_name[1]
      if ("username" %in% names(metadata) && identical(current_role(), "admin")) {
        sample_label = paste0(metadata$username[1], " - ", sample_label)
      }
      
      matches_df$Sursă = sample_label
      matches_df
    })
    
    matches_df = bind_rows(matches_list)
    req(nrow(matches_df) > 0)
    
    matches_df$rsid = sprintf(
      '<a href="https://www.ebi.ac.uk/gwas/labs/variants/%s" target="_blank">%s</a>',
      matches_df$rsid,
      matches_df$rsid
    )
    
    datatable(matches_df, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
}

shinyApp(ui, server)




