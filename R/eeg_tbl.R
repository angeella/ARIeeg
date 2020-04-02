#################Create eeg_lst object##################################
rm(list=ls())
source("library.R")
source("utils.R")
source("utils1.R")
source("new_events_tbl.R")
source("validate_eeg_lst.R")
source("signal_tbl.R")

#load data
load("data_eeg_emotion.RData")

#Builds an eeg_lst object composed of two data.table`objects and one tibble.
# All three are linked by a unique identifier `.id`. Amplitude values and timestamps appear in the `signal`
#table. Triggers, blinks, artifact rejection markings, and other
#events logged by the EEG recording software appear in the events table.
#Segment information and recording IDs appear in the segments tibble.

#signal tbl
#events tbl
#segments_tbl

####signal tbl

#' The `signal` table is organised into columns representing timestamps
#' (`.sample`) and individual electrodes. Each `.sample` corresponds to
#' 1 sample in the original recording, i.e. if the sampling rate of the EEG
#' recording is 500 Hz, then each `.sample` corresponds to 2 milliseconds.
#' These timestamps correspond to `.initial` in the `events` table, which
#' displays only the timestamps where logged events began.

###.id = dati$timings$epoch
## .sample = dati$timings$sample with attribute sampling_rate (dati$srate)
## .channels_dbl

.sample <- sample_int(dati$timings$sample,sampling_rate = dati$srate)

.id <- dati$timings$epoch

ch_n <- sapply(c(1:32), function(x) paste0("Ch",x))

nam <- dati$chan_info$electrode
for(i in 1:length(dati$chan_info$electrode)){
  assign(nam[i],channel_dbl(
    values = as.numeric(unlist(dati$signals[,i])),
    .x = dati$chan_info$cart_x[i],
    .y = dati$chan_info$cart_y[i],
    .z = dati$chan_info$cart_z[i],
    reference = "",
    resolution = 1,
    number = ch_n[i],
    unit = "microvolt",
    radius = dati$chan_info$radius[i],
    theta = dati$chan_info$theta[i],
    phi = dati$chan_info$phi[i]
  ))
}


.signal <- data.table::data.table(.id,.sample,Fp1,Fp2,F3,F4,F7,F8,FC1,FC2,C3,C4,T7,T8,
               CP1, CP2 ,P7  ,  P8  ,  P3  ,  P4,   PO7 ,  PO8 , O1 , O2  ,  Fz  ,  FCz , 
                Cz, CPz,Pz,RM,EOGvo,EOGvu,EOGhl,EOGhr)
#class(.signal) <- c("signal_tbl", "data.table", "data.frame")
#data.table::setkey(.signal, .id, .sample)
.signal[, .id := .id][, .sample := .sample]
data.table::setnames(.signal, make_names(colnames(.signal)))
data.table::setcolorder(.signal, c(".id", ".sample"))
data.table::setattr(.signal, "class", c("signal_tbl", class(.signal)))
data.table::setkey(.signal, .id, .sample)
#.segments, .id, .recording, segment
####.events tbl
###.id Integers indicating to which group the row of the signal matrix belongs.
## .initial = dati$events$event_onset
## .final =  as_sample_int(dati$events$event_onset + dati$events$event_time, sampling_rate =  dati$srate, unit = "s") - 1L,
## .channel = dati$chan_info$electrode
## descriptions_dt = data.table::data.table()

#init_events = sample_int(round(dati$events$event_onset * dati$srate) + 1L, sampling_rate = dati$srate)
.events <- new_events_tbl (
  .id = dati$events$epoch,
  #.initial = dati$events$event_onset %>%
  #  as_sample_int(sampling_rate = dati$srate, unit = "s"),
  .initial = dati$events$event_onset,
  #.initial =  init_events,
  #.final =  init_events,
  #.final =  round(dati$events$event_time* dati$srate) %>% as.integer() +init_events,
  .final = dati$events$event_onset + dati$srate - 1L,
  #.final = as_sample_int(dati$events$event_onset, 
  #                      sampling_rate =  dati$srate, unit = "s"),
  #.final = rep(1,length(dati$events$epoch)) %>%
  #  as_sample_int(sampling_rate = dati$srate, unit = "s"), 
  .channel = NA_character_,
  descriptions_dt = data.table::data.table(
 # .recording = dati$epochs$recording, 
                                           .type = rep("Stimulus",length(dati$events$event_type)),
                                           .description = paste0(dati$events$event_type)))
 #                                          .subj = dati$events$subj))

attr(.events$.initial,"sampling_rate")=dati$srate
attr(.events$.initial,"class")="sample_int"
attr(.events$.final,"sampling_rate")=dati$srate
attr(.events$.final,"class")="sample_int"
data.table::setkey(.events, .id)
#segments_tbl <- dplyr::tibble(.id =  dati$epochs$epoch, .recording = dati$epochs$recording)
segments_tbl <- dplyr::tibble(.id = dati$epochs$epoch, 
                              .recording = dati$epochs$recording,
                              #description = dati$events$event_type,
                              segment = rep(1,length(dati$epochs$recording)),
                              .subj = dati$events$subj)
                              #type = rep("Stimulus", length(dati$epochs$recording)))

segments_tbl <- validate_segments(segments_tbl)

data <-   eeg_lst(
  signal_tbl = .signal,
  events_tbl = .events,
  segments_tbl = segments_tbl)

is_eeg_lst(data)
is_events_tbl(data$.events)
is_signal_tbl(data$.signal)

#plot(data) 

summary(data)


channels_tbl(data)
signal_tbl(data)
events_tbl(data)
#Drop off final 5 channels 
chan_to_rm <- c("RM"  ,  "EOGvo" ,"EOGvu"
                , "EOGhl", "EOGhr")
data <- 
  data %>%
  select(-one_of(chan_to_rm))

data_seg <- data %>%
  eeg_segment(.description %in% c(1,2,3,4,5),
              lim = c(min(dati$timings$time), max(dati$timings$time))
  ) %>% eeg_baseline()

data_segs_some <- data_seg %>%
  mutate(
    condition =
      description
  ) %>%
  select(-type)

save(data_segs_some, file = "data_seg_some_all_cond.RData")
#Filter by condtions happy and neutral faces

data_seg <- data %>%
  eeg_segment(.description %in% c(3,5),
              lim = c(min(dati$timings$time), max(dati$timings$time))
  ) %>% eeg_baseline()

data_segs_some <- data_seg %>%
  mutate(
    condition =
      if_else(description == "3", "disgust", "object")
  ) %>%
  select(-type)

save(data_segs_some, file = "data_segs_some.RData")



#Some plots

#Plot of all the ERP of the O1 electrode

data_segs_some %>%
  select(O1) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line()

#Plot ERP of each subject (average across condition):

data_segs_some %>%
  select(O1) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(aes(group = .subj))  +
stat_summary(
  fun = "mean", geom = "line", alpha = 1, size = 1.5,
  aes(color = "red"),show.legend = FALSE
) 

##### average by condition
Dav=D0_segs%>%group_by(.sample,condition)%>%summarize_all(mean, na.rm = TRUE)
dim(Dav$.signal)
plot(Dav%>% select(c(30:40)) )


##### average POTENTIAL in range
AvePot=D0_segs%>%filter(between(as_time(.sample, unit = "s"), .1, .2)) %>%
  group_by(condition)%>%summarize_all(mean, na.rm = TRUE)

#Plot ERP of condition (average across subject)

data_segs_some %>%
  select(O1) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = 0.1,aes(group = condition))  +
  stat_summary(
    fun = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition),show.legend = FALSE
  ) 

###all condition
data %>%
  eeg_segment(.description %in% c(1,2,3,4,5),
              lim = c(min(dati$timings$time), max(dati$timings$time))
  ) %>% eeg_baseline() %>%
  mutate(
    condition =
      description
  ) %>%
  select(-type) %>%
  select("T7", "T8") %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(aes(group = condition))  +
  stat_summary(
    fun = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition),show.legend = TRUE
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")


data_segs_some %>%
  select(O1, O2, PO7, PO8) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")

data_segs_some %>%
  transmute(
    Occipital = chs_mean(O1, O2, na.rm = TRUE),
    Parietal = chs_mean(P3, P4, P7, P8, Pz, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun.y = "mean", geom = "line", alpha = 1, size = 1,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  theme(legend.position = "bottom")


ERP_data <- data_segs_some %>%
  group_by(.sample, condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE)



ERP_data %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(aes(color = condition)) +
  facet_wrap(~.key) +
  theme(legend.position = "bottom") +
  ggtitle("ERPs for happy vs neutral") +
  theme_eeguana()

ERP_plot %>% plot_in_layout()

data_segs_some %>%
  filter(between(as_time(.sample, unit = "s"), .1, .2)) %>%
  group_by(condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE) %>%
  plot_topo() +
  annotate_head() +
  geom_contour() +
  geom_text(colour = "black") +
  facet_grid(~condition)

df <- data_segs_some %>%
  select(O1, O2, P7, P8) %>%
  as_tibble() %>%
  # We can use regular dplyr functions now
  group_by(.key, .time) %>%
  summarize(
    `t-value` = t.test(
      .value[condition == "happy"],
      .value[condition == "neutral"]
    )$statistic
  )


ggplot(df, aes(x = .time, y = `t-value`)) + geom_line() +
  facet_wrap(~.key)


faces_seg_t <-
  data_segs_some %>%
  select(O1, O2, P7, P8) %>%
  group_by(.sample) %>%
  summarize_at(channel_names(.), list(t =  ~t.test(
    .[condition == "happy"],
    .[condition == "neutral"]
  )$statistic))

faces_seg_t %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id)) +
  stat_summary(fun.y = "mean", geom = "line", alpha = 1, size = 1) +
  facet_wrap(~.key) +
  theme(legend.position = "bottom")

