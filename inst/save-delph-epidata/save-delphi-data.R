library(dplyr)
library(tidyr)
library(purrr)
library(MMWRweek)
library(lubridate)
source("R/delphi_epidata.R")

#' For each week specified by a vector of dates giving the week start dates,
#' determine whether the week contains a specified date.
#'
#' @param week_start_date A vector of Date objects specifying the date of the first day in
#' the weeks of interest
#' @param year_to_pick integer or character giving the year to pick, e.g. "2010"
#' @param month_to_pick integer or character giving the month to pick, e.g. "12"
#' @param day_to_pick integer or character giving the day to pick, e.g. "22"
#'
#' @return a logical vector of the same length as time.  Entry i is TRUE if the
#' week beginning on week_start_date[i] contains the date specified by
#' year_to_pick, month_to_pick, and day_to_pick; FALSE otherwise.
#'
#' @export
pick_week <- function(
  week_start_date,
  year_to_pick,
  month_to_pick,
  day_to_pick) {
  date_to_pick <- lubridate::ymd(paste(year_to_pick, month_to_pick, day_to_pick, sep = "-"))
  selTF <- (week_start_date >= date_to_pick - 6) & (week_start_date <= date_to_pick)
  return(selTF)
}

# Function to fetch and combine data
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

location_codes <- c(
  'al', 'ak', 'az', 'ar', 'ca', 'co', 'ct', 'de', 'fl', 'ga', 'hi', 'id', 'il',
  'in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt',
  'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or',
  'pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy',
  'as', 'mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')

fluview_data <- fetch_delphi_data_multi_issue(
  source = "fluview",
  #    regions = c("ny", "jfk", "ny_minus_jfk"),
  regions = location_codes,
  issues = all_issues,
  epiweeks_range = c(201040, 201930))

fluview_data_wide <- fluview_data %>%
  group_by(region, epiweek) %>%
  filter(issue == max(issue)) %>%
  ungroup() %>%
  select(region, epiweek, wili) %>%
  spread(key = "region", value = "wili")
names(fluview_data_wide) <- c("epiweek", paste0("wili_", names(fluview_data_wide)[-1]))

wiki_data <- fetch_delphi_data_multi_issue(
    source = "wiki",
    epiweeks_range = c(201040, 201930))

all_obs <- fluview_data_wide %>%
    left_join(
        wiki_data
    ) %>%
    select(
        epiweek,
        year,
        week,
        starts_with("wili_"),
        starts_with("wiki_")
    ) %>% dplyr::mutate(
      time = MMWRweek::MMWRweek2Date(year, week),
      christmas_week = pick_week(
        week_start_date = time,
        year_to_pick = lubridate::year(time),
        month_to_pick = "12",
        day_to_pick = "25"),
      postchristmas_week = pick_week(
        week_start_date = time,
        year_to_pick = lubridate::year(time),
        month_to_pick = "1",
        day_to_pick = "1") |
        pick_week(
          week_start_date = time,
          year_to_pick = lubridate::year(time) + 1,
          month_to_pick = "1",
          day_to_pick = "1"),
      thanksgiving_day = lubridate::year(time) %>%
        splusTimeDate::holiday.Thanksgiving() %>%
        as.character() %>%
        strsplit(split = "/") %>%
        `[[`(1) %>%
        `[`(2) %>%
        as.numeric(),
      postthanksgiving_day = lubridate::year(time) %>%
        splusTimeDate::holiday.Thanksgiving() %>%
        `+`(7) %>%
        as.character() %>%
        strsplit(split = "/") %>%
        `[[`(1) %>%
        `[`(2) %>%
        as.numeric(),
      postthanksgiving_month = lubridate::year(time) %>%
        splusTimeDate::holiday.Thanksgiving() %>%
        `+`(7) %>%
        as.character() %>%
        strsplit(split = "/") %>%
        `[[`(1) %>%
        `[`(1) %>%
        as.numeric(),
      thanksgiving_week = pick_week(
        week_start_date = time,
        year_to_pick = lubridate::year(time),
        month_to_pick = "11",
        day_to_pick = thanksgiving_day),
      postthanksgiving_week = pick_week(
        week_start_date = time,
        year_to_pick = lubridate::year(time),
        month_to_pick = postthanksgiving_month,
        day_to_pick = postthanksgiving_day)
    ) %>% dplyr::select(
      -thanksgiving_day,
      -postthanksgiving_day,
      -postthanksgiving_month
    )


saveRDS(all_obs, file = "data/all_locations_flu_wiki.rds")
