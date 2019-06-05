library(dplyr)
library(tidyr)
library(purrr)
library(MMWRweek)
source("R/delphi_epidata.R")

# Function to fetch and combined data
fetch_delphi_data_multi_issue <- function(
    source,
    regions = "nat",
    issues = 1,
    epiweeks_range) {
    
    epiweeks <- list(Epidata$range(epiweeks_range[1], epiweeks_range[2]))
    
    all_obs <- map_dfr(regions, function(region) {
        map_dfr(issues, function(issue) {
            if(identical(source, "fluview")) {
                obs_one_issue <- Epidata$fluview(
                    regions = list(region),
                    epiweeks = epiweeks,
                    issue = list(issue))
            } else if(identical(source, "twitter")) {
                obs_one_issue <- Epidata$twitter(
                    locations = list(region),
                    epiweeks = epiweeks
                )
            } else if(identical(source, "wiki")) {
                obs_one_issue <- Epidata$wiki(
                    articles = list("influenza", "common_cold", "cough"),
                    language = "en",
                    epiweeks = epiweeks
                )
            } else {
                stop("Unsupported Epidata source")
            }

            temp <- map_dfr(obs_one_issue$epidata,
                function(x) {
                    x[sapply(x, function(comp) is.null(comp))] <- NA
                    return(as.data.frame(x))
                })
            
            if(identical(source, "wiki")) {
                temp <- temp %>%
                    select(article, epiweek, value) %>%
                    spread(article, value) %>%
                    rename(
                        wiki_influenza = influenza,
                        wiki_common_cold = common_cold,
                        wiki_cough = cough
                    )
            }
            
            return(temp)
        })
    })
    
    all_obs <- all_obs %>%
        separate(epiweek, c("year", "week"), sep=4, remove=FALSE) %>%
        mutate(
            year = as.integer(year),
            week = as.integer(week))
#            release_date = as.Date(release_date))
    
    return(all_obs)
}


# Fetch ILI data
all_issues <- expand.grid(
    year = 2000:2019,
    week = sprintf("%02d", 1:53)
  ) %>%
  apply(1, function(x) paste(x, collapse = "")) %>%
  as.integer() %>%
  sort()
#all_issues <- all_issues[all_issues >= 201510 & all_issues <= 201530]
#all_issues <- all_issues[all_issues >= 200040 & all_issues <= 201530]

fluview_data_nyc <- fetch_delphi_data_multi_issue(
    source = "fluview",
#    regions = c("ny", "jfk", "ny_minus_jfk"),
    regions = c("jfk"),
    issues = all_issues,
    epiweeks_range = c(201040, 201530))

fluview_data_ny_upstate <- fetch_delphi_data_multi_issue(
    source = "fluview",
#    regions = c("ny", "jfk", "ny_minus_jfk"),
    regions = c("ny_minus_jfk"),
    issues = all_issues,
    epiweeks_range = c(201040, 201530))

#twitter_data_ny <- fetch_delphi_data_multi_issue(
#    source = "twitter",
#    regions = c("ny"),
#    epiweeks_range = c(201040, 201530))

wiki_data_ny <- fetch_delphi_data_multi_issue(
    source = "wiki",
    epiweeks_range = c(201040, 201530))

all_obs <- fluview_data_nyc %>%
    filter(issue == 201740) %>%
    transmute(
        epiweek = epiweek,
        wili_jfk = wili
    ) %>%
    left_join(
        fluview_data_ny_upstate %>%
        filter(issue == 201740) %>%
        transmute(
            epiweek = epiweek,
            wili_ny_upstate = wili
        )
    ) %>%
    left_join(
        wiki_data_ny
    ) %>%
    select(
        epiweek,
        year,
        week,
        wili_jfk,
        wili_ny_upstate,
        wiki_common_cold,
        wiki_cough,
        wiki_influenza
    )

saveRDS(all_obs, file = "data/ny_flu_wiki.rds")
