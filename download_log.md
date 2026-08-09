# Download / audit log - EQ-1 replication

## Protocol freeze record

Frozen 2026-07-20 23:01:01 PDT, **before** any regional event query was issued.

SHA-256 at freeze time:

```
23355ee57356b69c92570d98d5fe2cd6f58570c620215ed474748b6ca665c12e  PROTOCOL.md
87aa64e6d071e487ca7b317ba2039985cc192aa046f0c2299c29b9f27a4794d7  protocol_params.json
b0affd004e7a739762f09952ca7843139a3452fec0a44d9feb25727390b7ad07  download_catalogs.py
```

For true pre-registration these hashes should be pinned somewhere third-party-timestamped
(git commit pushed to GitHub, or OSF registration) - repo is not yet a git repo; recommend
`git init` + push before analysis unblinding. (Left to Jim per repo conventions.)

## Amendments to non-analytic tooling

- 2026-07-20 ~23:05 PDT: first download run failed with `SSL: CERTIFICATE_VERIFY_FAILED`
  (miniconda Python lacks a system CA bundle). `download_catalogs.py` amended to build an SSL
  context from `certifi`. **No analytic parameter changed**; no event data had been received.
  New hash: `f211ec37dabb5cc303db69bbc68621686e5dd45c8c8e75b4334557cf08e51954  download_catalogs.py`

## Amendment 1 hash record

2026-07-20 23:19:42 PDT - original Kosmos r89 tidal-window definitions adopted (PROTOCOL.md
Amendment 1; params v1.1.0). Raw catalogs downloaded but NOT analyzed at amendment time.

```
76404c44df02665645a26f343114f478f8f8944b2bd03723f3c36e30d305c35d  PROTOCOL.md
ef4d0a37e25ed67aadc245141836276a0b064f289416ee4da563d3ed9f4271cd  protocol_params.json
```

## EQ-18 cross-test protocol freeze record

2026-07-21 00:14:24 PDT - XUE_LU_PROTOCOL.md frozen BEFORE any analysis of the Lu/Xue Zenodo
catalogs (only head-rows/row-counts had been inspected, per INVENTORY.md).

```
3e36181a112b758ab7858ce244f2c5baf753fd63bf731d026e32a50a60011e09  XUE_LU_PROTOCOL.md
cd2b1e3c8d96a37a36a47f6db95d2a928b49c99b01a1a44d0e30f09dc6b4a59d  xue_lu_params.json
```

## EQ-18 tooling amendments (non-analytic)

- 2026-07-21 ~00:25 PDT: first H1 run produced invalid associations because the RAW Lu/Xue
  catalogs order columns `... sec EID lat lon depth mag` while the declustered files use
  `... sec lat lon depth mag EID`. `xue_lu_crosstest.py::load_catalog` now auto-detects order and
  asserts lat/lon validity. No results from the misparsed run were interpreted or retained.
  Environment note: astropy install had upgraded numpy (2.5.1), requiring pandas/scipy upgrades
  (pandas 3.0.3, scipy 1.18.0) to fix a binary-incompatibility crash before any analysis ran.

## Retrieval runs

See `data/retrieval_log_*.tsv` written by `download_catalogs.py`.

## Coso Fig 4c reproduction protocol freeze record

aa97685b1000e7f5da0affc05ff81483a0d168a01191da385bb7faccd053bae0 *COSO_FIG4C_PROTOCOL.md
Thu Aug 6 2026 (frozen BEFORE any analysis in Weifan Lu's Fig 4c bin — lat [36.2,36.6],
lon [−118,−117.6]; parameters received from Weifan by email 2026-08-06. Prior Coso runs
used a different, mostly non-overlapping box; no analysis of this bin had been run.)

## EQ-22 Long Valley protocol freeze record

826ddc072275789e7f85276560f5ad9592033cdc6ea23bcb048e83a74920cffc *LONG_VALLEY_PROTOCOL.md
Tue Jul 21 09:30:40 PDT 2026
(frozen BEFORE any Long Valley data download)

## Overnight prediction-experiments protocol freeze record

62b259286d568a8ef0a446cf9749a8b1850f4ece0a3040fb22b8b02908755adf *OVERNIGHT_PREDICTION_PROTOCOL.md
Sun Aug 9 2026 (frozen BEFORE any test-window (2010+) analysis under these definitions;
train/test split 2010-01-01 chosen from span/power only)

## Pattern-experiments protocol freeze record (round 2)

165527d14b28bd1a0ea1cf5340e6b0252548d0e5ea8c2d33e6f1edca0ad2aa16 *PATTERN_PROTOCOL.md
Sun Aug 9 2026 (frozen BEFORE any test-window analysis under these definitions)

## EXP-J addendum hash
1e126abc8a42a1e6ca7e0a5b9a874d83a9084ff9b5496c2bb5d57047ed119a6d *PATTERN_PROTOCOL.md
Sun Aug 9 2026 (EXP-J frozen before computation)

## EXP-K/L addendum hash
4b347599113aca2dd7ae6313c178f2e142d3fc31e0632e8e82b81a33fa581e54 *PATTERN_PROTOCOL.md
Sun Aug 9 2026

## EXP-M addendum hash
aca4b729277762fe1ca9f9fdf561291e3527710f5240394b7a7a12e06d6995b2 *PATTERN_PROTOCOL.md
Sun Aug 9 2026 (frozen before any global download)
