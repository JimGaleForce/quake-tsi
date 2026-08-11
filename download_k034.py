"""K-034 target-cell catalogue download (ComCat FDSN). Freeze+hash before analysis."""
import os, sys, time, json, hashlib, urllib.request, ssl, io
import pandas as pd, certifi

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'data', 'k034')
os.makedirs(OUT, exist_ok=True)
CTX = ssl.create_default_context(cafile=certifi.where())
URL = 'https://earthquake.usgs.gov/fdsnws/event/1/query'

CELLS = {
 'long_valley':      (37.40, 37.90, -119.10, -118.60),
 'coso':             (35.80, 36.20, -118.00, -117.60),
 'geysers':          (38.70, 38.92, -122.95, -122.70),
 'yellowstone':      (44.30, 44.90, -111.20, -110.40),
 'salton_brawley':   (32.90, 33.30, -115.80, -115.40),
 'lassen':           (40.30, 40.70, -121.70, -121.30),
 'mono_west_nv_mina':(38.20, 38.60, -118.40, -117.90),
 'little_skull_mtn': (36.60, 37.00, -116.50, -116.00),
 'cedar_city_ut':    (37.40, 37.90, -113.40, -112.80),
 'smith_valley_nv':  (38.60, 39.00, -119.60, -119.20),
 'parkfield':        (35.70, 36.10, -120.70, -120.30),
 'mendocino':        (40.20, 40.70, -124.70, -124.10),
 'wasatch_slc':      (40.40, 41.00, -112.20, -111.60),
 'san_jacinto':      (33.20, 33.70, -116.90, -116.40),
}
T0, T1, MINMAG, CAP = '1985-01-01', '2023-01-01', 1.5, 20000
log = []

def fetch(box, t0, t1, depth=0):
    la0, la1, lo0, lo1 = box
    q = (f'{URL}?format=csv&starttime={t0}&endtime={t1}&minlatitude={la0}&maxlatitude={la1}'
         f'&minlongitude={lo0}&maxlongitude={lo1}&minmagnitude={MINMAG}&orderby=time-asc&limit={CAP}')
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(q, context=CTX, timeout=300) as r:
                body = r.read().decode('utf-8')
            break
        except Exception as e:
            log.append(f'  retry {attempt} {t0}->{t1}: {e}')
            time.sleep(4 * attempt)
    else:
        raise RuntimeError(f'failed {t0}->{t1}')
    df = pd.read_csv(io.StringIO(body)) if body.strip() else pd.DataFrame()
    n = len(df)
    log.append(f'  {t0} -> {t1}  rows={n}')
    if n >= CAP:
        if depth > 6:
            raise RuntimeError('cap hit too deep')
        mid = (pd.Timestamp(t0) + (pd.Timestamp(t1) - pd.Timestamp(t0)) / 2).strftime('%Y-%m-%dT%H:%M:%S')
        log.append(f'  CAP HIT -> split at {mid}')
        return pd.concat([fetch(box, t0, mid, depth + 1), fetch(box, mid, t1, depth + 1)], ignore_index=True)
    return df

manifest = {}
for name, box in CELLS.items():
    p = os.path.join(OUT, f'{name}.csv')
    if os.path.exists(p):
        log.append(f'{name}: cached');
    else:
        log.append(f'{name}: box={box}')
        df = fetch(box, T0, T1)
        if len(df):
            df = df.drop_duplicates(subset=['id']).sort_values('time')
        df.to_csv(p, index=False)
        log.append(f'{name}: FINAL rows={len(df)}')
    h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    d = pd.read_csv(p)
    manifest[name] = dict(box=box, rows=int(len(d)), sha256=h,
                          first=str(d['time'].iloc[0]) if len(d) else None,
                          last=str(d['time'].iloc[-1]) if len(d) else None)
    print(name, manifest[name]['rows'], h[:16], flush=True)

json.dump(dict(minmag=MINMAG, t0=T0, t1=T1, cells=manifest), open(os.path.join(OUT, 'manifest.json'), 'w'), indent=1)
open(os.path.join(OUT, 'download.log'), 'w').write('\n'.join(log))
print('DONE')
