# EFA mock-data project

This project runs SVD-based evolving factor analysis (EFA) on mock q-by-time data.

Expected data location:

```text
C:\Users\ming0\Desktop\MATLAB\mock data\q.csv
C:\Users\ming0\Desktop\MATLAB\mock data\t_81points.csv
C:\Users\ming0\Desktop\MATLAB\mock data\q_delta_S_noise_10%_81points.csv
```

Run from MATLAB:

```matlab
cd("C:\Users\ming0\Documents\MATLAB\Projects\EFA")
run_efa_mockdata
```

Main outputs are written to:

```text
C:\Users\ming0\Documents\MATLAB\Projects\EFA\results
```

The script first looks for data in `C:\Users\ming0\Desktop\MATLAB\mock data`.
If that folder does not exist, it falls back to
`C:\Users\ming0\Documents\MATLAB\mock data`.

The script produces:

- q-time signal map with a log-scaled time axis.
- Forward EFA singular-value traces from SVD of columns `1:k`.
- Backward EFA singular-value traces from SVD of columns `k:end`.
- Transition-time graph showing forward appearance and backward disappearance markers on the log time axis.
- q-spectrum similarity heatmap and matched overlay plots for separating whether front/back components refer to the same spectral species.
- Component-count table and figures.
- Summary table showing how many columns are needed before each component is above threshold in forward and backward scans.

The practical component threshold is controlled near the top of `run_efa_mockdata.m`:

```matlab
opts.relativeSingularValueThreshold = 1e-2;
opts.energyCutoff = 0.995;
```

By default, the main component threshold is estimated from the largest
spectral gap in the first few full-data singular values:

```matlab
opts.thresholdMode = "spectralGap";
opts.gapSearchStartComponent = 2;
opts.gapSearchEndComponent = opts.maxComponents;
```

Use `thresholdMode = "relative"` if you want the older fixed-fraction
threshold based on the largest full-data singular value.
