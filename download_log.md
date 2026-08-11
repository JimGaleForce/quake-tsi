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

## Amendment: LONG_VALLEY_PROTOCOL.md hash discrepancy (recorded 2026-08-09)

Faraday's audit: the recorded freeze hash 826ddc07... does not match the committed (fd5b2978...)
or working-copy (d3141cd1...) versions. PROTOCOL.md and XUE_LU_PROTOCOL.md verify EXACTLY at
pre-copyedit commit 78f2227, so the 0d8f897 copyedit commit is the likely cause: the record is
shifted by a documented copyedit, not broken. Continuity is via git history (see also the
0d8f897 note in HANDOFF). Do not recompute or "fix" the recorded hash; this amendment is the
record.

## K-034 pre-registration freeze record (positive control, Landers 1992)

Frozen 2026-08-11 05:51:00 UTC, **before** any K-034 download was issued and before any K-034
statistic was computed. Sealed literature values (Popper R2-2 mandate 1) and the pre-registered
ranked cell list / windows / statistics / PASS rule (mandates 2-4).

```
eae95839fcdea8b9ca62097fed25a02b74559c03b59c1913bcb59cac4f8c320e *K034_SEALED_LITERATURE.md
01e41f971a73f449a09b6faf09f569cd47f3395e1ca4d584c76ade699ac1c226 *K034_PREREGISTERED_CELLS.md
```

DEVIATION (flagged in both files and in results_k034.json): R2-2 assigns authorship of the seal to
the supervisor. Executed by a single worker agent, it is analyst-authored; the hash pins the
comparison target against later editing but not against analyst foreknowledge.

## K-034 data freeze record (ComCat FDSN, target cells)

Frozen 2026-08-11 05:52:28 UTC, immediately after retrieval and **before** any K-034 statistic was
computed. ComCat FDSN event/1/query, format=csv, M>=1.5, 1985-01-01 -> 2023-01-01, orderby=time-asc,
limit=20000 with recursive halving on cap-hit (none of the 14 cells hit the cap after halving; see
data/k034/download.log). Retrieval script: download_k034.py. Boxes are exactly those pre-registered
in K034_PREREGISTERED_CELLS.md (hashed above, before this download).

```
5b2026110c9f5e4805137e94327698c9f83e06af3c137a635c30aef510e014af  data/k034/long_valley.csv  rows=29241  1985-01-01 .. 2022-12-26
0a913da71aa0a9e18ae475edf45749a5c99a550d913ad464c1b9e0d87d0430cc  data/k034/coso.csv  rows=22320  1985-01-01 .. 2022-12-28
a09ecf92a452b907b874c19496ad8f6898609d20dee597cc3b7ca086d6052b45  data/k034/geysers.csv  rows=34546  1985-01-01 .. 2022-12-29
cd3e38a45f760e85b50972f7a63e8a484a1f59487a453b08726748872c636cc5  data/k034/yellowstone.csv  rows=6862  1985-01-16 .. 2022-12-18
72b9ac5eed92efd916273be8dc8f2df19a1e6d87dd713f781ceff8718340330f  data/k034/salton_brawley.csv  rows=13998  1985-01-09 .. 2022-12-29
6e27d9c2257799dcc0132325cefd20e7e3c9316281b954d0ab5e5e4c155bf809  data/k034/lassen.csv  rows=944  1985-02-17 .. 2022-11-10
3ba9f2de000a9cdcf375454dccd7a078790df614d62df5539230e8d6a86248da  data/k034/mono_west_nv_mina.csv  rows=992  1985-06-12 .. 2022-12-31
1cc296bf953e323f9468926c772a88bf0d9af9bf44bcaf82c85883723e8eeb89  data/k034/little_skull_mtn.csv  rows=374  1985-01-20 .. 2022-12-01
4e96f4f535bb4600eef07e60a327892074c0bc81edcd4f996b1985283ff1c9e0  data/k034/cedar_city_ut.csv  rows=789  1985-02-15 .. 2022-11-15
4a46be7da5e13c9292cb44c80f1b5f3de5efc50bdc35533a3c016a34edae1eef  data/k034/smith_valley_nv.csv  rows=662  1985-01-17 .. 2022-12-22
46d86d0d7ca63c29112e8a6590895cf822655e167acc7497482e032159d949db  data/k034/parkfield.csv  rows=3308  1985-01-04 .. 2022-12-11
c7cc4b08ba2c6cc978dc5213379899e1d72b00a156c6e0dad8e647424219b83e  data/k034/mendocino.csv  rows=11118  1985-01-02 .. 2022-12-31
32bbb12f48d16de7a8f49d34ca62718f56ac04ed4e70c418ad7a586777c8bc5a  data/k034/wasatch_slc.csv  rows=1016  1985-04-21 .. 2022-11-19
14869a7de15c8ac7a2143503a7343c84ba13202b6def6cc9e56c1df8a8bb2d51  data/k034/san_jacinto.csv  rows=11310  1985-01-01 .. 2022-12-30
7dcf6ec28b2d4259709bc8854cefbb0cb55467035ed2acb611b972a2ce24646f *download_k034.py
```
