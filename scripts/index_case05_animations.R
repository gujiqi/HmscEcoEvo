#!/usr/bin/env Rscript

# Verify the Case05 videos and build a local, filterable gallery.
suppressPackageStartupMessages(library(av))

default_root <- file.path(getwd(), "..", "outputs", "HmscEcoEvo")
root <- normalizePath(Sys.getenv("HMSCEE_OUTPUT_ROOT", default_root),
                      winslash = "/", mustWork = TRUE)
output <- file.path(root, "case05_plant200_time_animations_20260928")
index <- read.csv(file.path(output, "animation_index.csv"),
                  stringsAsFactors = FALSE)
if (anyDuplicated(index$group_id)) stop("Duplicate animation IDs")

validation <- lapply(seq_len(nrow(index)), function(i) {
  row <- index[i, ]
  if (!file.exists(row$video) || !file.exists(row$poster)) {
    stop("Missing video or poster: ", row$group_id)
  }
  info <- av::av_media_info(row$video)
  duration <- as.numeric(info$duration)
  width <- as.integer(info$video$width[1])
  height <- as.integer(info$video$height[1])
  codec <- as.character(info$video$codec[1])
  if (!is.finite(duration) || duration < 5 ||
      width != 1600L || height != 900L || codec != "h264") {
    stop("Invalid MP4 metadata: ", row$group_id)
  }
  data.frame(group_id = row$group_id, category = row$category,
             scheme = row$scheme, metric = row$metric,
             n_time_slices = row$n_time_slices,
             duration_seconds = duration,
             mp4_bytes = file.info(row$video)$size,
             width = width, height = height, codec = codec)
})
validation <- do.call(rbind, validation)
write.csv(validation, file.path(output, "animation_validation.csv"),
          row.names = FALSE)

titles <- c(
  "01_density" = "Lineage location density",
  "02_log_density" = "Low-density lineage support (log scale)",
  "03_ancestral_eta" = "Ancestral environmental response",
  "04_land_fraction" = "Historical land fraction",
  "05_support_intensity" = "Lineage-support intensity",
  "06_support_shannon" = "Support-mixture Shannon index",
  "07_support_effective" = "Effective support lineages",
  "effective_support" = "Corrected effective support diversity",
  "active_lineage_environmental_support_count" = "Environmentally supported lineages",
  "active_lineage_environmental_support_fraction" = "Fraction of lineages supported",
  "modern_training_noanalog_flag" = "No-analog environment flag",
  "ancestral_environmental_support" = "Ancestral environmental support",
  "candidate_ecological_refugia_evidence_class" = "Candidate ecological refugia",
  "location_density_difference" = "Location-density scheme difference"
)
schemes <- c(
  spherical_no_ldd = "Spherical / no rare jumps",
  spherical_land = "Spherical / rare jumps",
  terrain_resistance = "Terrain resistance",
  particle_topography = "Finite propagules / terrain",
  manuscript_corrected = "Manuscript-corrected result",
  all_lineages = "All active lineages"
)
categories <- c(
  corrected_diversity_66 = "Corrected diversity",
  full_66 = "Four dispersal schemes / 66 slices",
  sampled_3_or_6 = "Selected output slices",
  sampled_scheme_difference = "Selected scheme differences"
)
label <- function(values, key) {
  value <- unname(values[key])
  if (is.na(value)) gsub("_", " ", key) else value
}
escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

priority <- match(index$category,
                  c("corrected_diversity_66", "full_66",
                    "sampled_3_or_6", "sampled_scheme_difference"))
index <- index[order(priority, index$scheme, index$metric), ]
cards <- vapply(seq_len(nrow(index)), function(i) {
  row <- index[i, ]
  title <- label(titles, row$metric)
  scheme <- label(schemes, row$scheme)
  category <- label(categories, row$category)
  coverage <- sprintf("%s to %s Ma | %s time slices",
                      row$oldest_ma, row$youngest_ma, row$n_time_slices)
  if (row$n_time_slices < 66) {
    coverage <- paste(coverage, "| sampled snapshots only")
  }
  sprintf(paste0(
    '<article class="item" data-category="%s" data-scheme="%s" ',
    'data-search="%s"><div class="kind">%s</div><h2>%s</h2>',
    '<p>%s</p><video controls preload="none" poster="posters/%s">',
    '<source src="videos/%s" type="video/mp4"></video>',
    '<div class="foot"><span>%s</span><a href="videos/%s">Open MP4</a></div>',
    '</article>'),
    escape(row$category), escape(row$scheme),
    escape(tolower(paste(title, scheme, category))),
    escape(category), escape(title), escape(scheme),
    escape(basename(row$poster)), escape(basename(row$video)),
    escape(coverage), escape(basename(row$video)))
}, character(1))

scheme_ids <- sort(unique(index$scheme))
scheme_options <- vapply(scheme_ids, function(x) {
  sprintf('<option value="%s">%s</option>', escape(x),
          escape(label(schemes, x)))
}, character(1))

html <- c(
  '<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>HmscEE Case05 time animations</title><style>',
  ':root{font-family:Segoe UI,Arial,sans-serif;color:#18364a;background:#f7f9f9}',
  '*{box-sizing:border-box}body{margin:0}header{background:white;border-bottom:1px solid #d9e3e4}',
  '.wrap{max-width:1450px;margin:auto;padding:22px 28px}',
  'h1{font-size:30px;margin:0 0 6px}header p{margin:0;color:#526b79;line-height:1.5}',
  '.tools{display:flex;gap:12px;flex-wrap:wrap;padding:18px 0 8px}',
  'select,input{background:white;color:#18364a;border:1px solid #bfcfd3;',
  'border-radius:5px;padding:10px 12px;font:inherit;min-width:220px}',
  'input{flex:1;max-width:430px}.count{color:#607583;font-size:14px;padding:8px 0}',
  '.note{color:#526b79;font-size:14px;line-height:1.5;margin:4px 0 22px}',
  '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}',
  '.item{background:white;border:1px solid #d9e3e4;border-radius:6px;',
  'padding:16px;min-width:0}.item[hidden]{display:none}',
  '.kind{font-size:12px;text-transform:uppercase;color:#177b70;font-weight:700}',
  'h2{font-size:19px;margin:6px 0;letter-spacing:0}.item p{color:#58707e;margin:0 0 12px}',
  'video{width:100%;aspect-ratio:16/9;background:#e8eeed;display:block}',
  '.foot{display:flex;justify-content:space-between;gap:12px;font-size:13px;',
  'color:#5d7280;padding-top:10px}.foot a{color:#126e68;text-decoration:none;font-weight:700}',
  '@media(max-width:900px){.grid{grid-template-columns:1fr}.wrap{padding:18px}',
  'h1{font-size:25px}}',
  '</style></head><body>',
  '<header><div class="wrap"><h1>Case05 | Plant200 through time</h1>',
  '<p>Maps from the ancestral-response and movement workflows, with a matched',
  ' global trajectory below each map.</p></div></header><main class="wrap">',
  '<div class="tools"><select id="category"><option value="">All result groups</option>',
  '<option value="corrected_diversity_66">Corrected diversity</option>',
  '<option value="full_66">Four schemes / 66 slices</option>',
  '<option value="sampled_3_or_6">Selected outputs</option>',
  '<option value="sampled_scheme_difference">Scheme differences</option></select>',
  '<select id="scheme"><option value="">All schemes / lineages</option>',
  scheme_options, '</select>',
  '<input id="search" type="search" placeholder="Search metrics or lineages"></div>',
  '<div id="count" class="count"></div>',
  '<p class="note">Each metric has a fixed colour scale; ocean is white.',
  ' The lower curve is the land-area-weighted spatial mean, except scheme',
  ' comparisons show mean absolute difference and class maps show land share.',
  ' Sampled snapshots do not interpolate',
  ' missing times. Model support is not observed occupancy. The 0 Ma endpoint',
  ' uses modern conditioning and should be interpreted separately.</p>',
  '<section class="grid">', cards, '</section></main>',
  '<script>',
  'const cat=document.getElementById("category"),',
  'scheme=document.getElementById("scheme"),',
  'query=document.getElementById("search"),',
  'cards=[...document.querySelectorAll(".item")];',
  'function filter(){let n=0;for(const card of cards){',
  'const show=(!cat.value||card.dataset.category===cat.value)&&',
  '(!scheme.value||card.dataset.scheme===scheme.value)&&',
  'card.dataset.search.includes(query.value.toLowerCase().trim());',
  'card.hidden=!show;if(show)n++;else card.querySelector("video").pause();}',
  'document.getElementById("count").textContent=n+" animations shown";}',
  'for(const el of [cat,scheme,query])el.addEventListener("input",filter);filter();',
  '</script></body></html>'
)
writeLines(html, file.path(output, "index.html"), useBytes = TRUE)
message("Validated ", nrow(index), " videos and wrote ", file.path(output, "index.html"))
