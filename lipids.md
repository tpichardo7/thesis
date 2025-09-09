lipids
================

``` r
lipids_df = read_excel("./data/lipids.xls",
            sheet = "Final Data",
            skip = 1, 
            col_names = TRUE) |> 
  janitor::clean_names() |> 
  mutate(sample_id = as.numeric(sub("^T", "", sample_id))) |> 
  mutate(sample_id = sample_id %% 1000) |> 
  rename(sid = sample_id) |> 
  drop_na()

view(lipids_df)
```
