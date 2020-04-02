make_names <- function(names) {
  make.names(names) %>% make.unique()
}

new_events_tbl <- function(.id = integer(0),
                           .initial = sample_int(integer(0), integer(0)),
                           .final = sample_int(integer(0), integer(0)),
                           .channel = character(0),
                           descriptions_dt = data.table::data.table(),
                           sampling_rate = NULL) {
  if (length(c(.id, .initial, .final, .channel, descriptions_dt)) == 0) {
    events <- data.table::data.table(
      .id = .id,
      .initial = .initial,
      .final = .final,
      .channel = .channel
    )
  } else {
    if (length(.channel) == 0) .channel <- NA_character_
    
    if (length(descriptions_dt) == 0) {
      events <- data.table::data.table(
        .id = .id,
        .initial = .initial,
        .final = .final,
        .channel = .channel
      )
    } else {
      events <- data.table::data.table(
        .id = .id,
        descriptions_dt,
        .initial = .initial,
        .final = .final,
        .channel = .channel
      )
    }
  }
  if (!is.null(sampling_rate)) {
    events[, .initial := sample_int(as.integer(.initial),
                                    sampling_rate = sampling_rate)
           ]
    events[, .final := sample_int(as.integer(.final),
                                  sampling_rate = sampling_rate)
           ]
  }
  data.table::setattr(events, "class", c("events_tbl", class(events)))
  events[]
}