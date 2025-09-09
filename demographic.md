demographic
================

``` r
demographic_df = read_csv("./data/wtc_demographic.csv") |> 
  janitor::clean_names() |> 
  mutate(sid = sid %% 1000)
```

    ## New names:
    ## Rows: 329 Columns: 23
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," dbl
    ## (22): ...1, sid, sid_new, hosp, DaysSinceWTC, GA_days, gender, c_sectio... date
    ## (1): child_DOB
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
view(demographic_df)
```
