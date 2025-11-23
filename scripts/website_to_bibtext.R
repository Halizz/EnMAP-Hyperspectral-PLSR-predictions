# Install required packages if needed
# install.packages(c("rvest", "stringr", "lubridate"))

library(rvest)
library(stringr)
library(lubridate)

website_to_bibtex <- function(url, output_file = "sample.bib") {
  # Read the website
  cat("Fetching website:", url, "\n")
  webpage <- tryCatch({
    read_html(url)
  }, error = function(e) {
    stop("Error reading website: ", e$message)
  })

  # Extract basic information
  title <- webpage %>%
    html_nodes("title") %>%
    html_text() %>%
    trimws()

  # Try to find author information - this will vary by website
  authors <- webpage %>%
    html_nodes("meta[name='author'], meta[property='og:author'], .author, [rel='author']") %>%
    html_text() %>%
    paste(collapse = " and ")

  if (authors == "") {
    authors <- "Unknown Author"
  }

  # Try to find publication date
  date_nodes <- webpage %>%
    html_nodes("time, meta[property='article:published_time'], [itemprop='datePublished']")

  if (length(date_nodes) > 0) {
    date_text <- date_nodes[1] %>% html_text() %>% trimws()
    # Try to extract just the year
    year <- str_extract(date_text, "\\d{4}")
    if (is.na(year)) {
      year <- format(Sys.Date(), "%Y") # Current year as fallback
    }
  } else {
    year <- format(Sys.Date(), "%Y") # Current year as fallback
  }

  # Generate a key for the BibTeX entry
  key <- str_replace_all(tolower(substr(title, 1, 20)), "[^a-z0-9]", "") %>%
         paste0(year)

  # Create BibTeX entry
  bibtex_entry <- sprintf(
    "@online{%s,
  title = {%s},
  author = {%s},
  year = {%s},
  url = {%s},
  urldate = {%s}
}",
    key,
    title,
    authors,
    year,
    url,
    format(Sys.Date(), "%Y-%m-%d")
  )

  # Write to file
  cat("Writing BibTeX entry to", output_file, "\n")
  writeLines(bibtex_entry, output_file)

  cat("BibTeX entry:\n\n")
  cat(bibtex_entry, "\n")

  return(invisible(bibtex_entry))
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur manuell aufgerufen
# # Example usage
# website_to_bibtex("https://www.usgs.gov/landsat-missions/ground-control-points")
