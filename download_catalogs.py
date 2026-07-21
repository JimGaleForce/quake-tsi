"""EQ-1 catalog downloader. Reads protocol_params.json; writes raw FDSN responses
under data/raw/<region>/<year>.txt unmodified, plus a retrieval log.

Stdlib only. Chunked per calendar year to stay under per-query event caps.
"""
import json
import ssl
import time
import urllib.request
import urllib.parse
from datetime import date
from pathlib import Path

try:
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    SSL_CTX = ssl.create_default_context()

HERE = Path(__file__).parent
PARAMS = json.loads((HERE / "protocol_params.json").read_text(encoding="utf-8"))

ENDPOINTS = {
    "INGV": ("https://webservices.ingv.it/fdsnws/event/1/query", "text"),
    "GEONET": ("https://service.geonet.org.nz/fdsnws/event/1/query", "text"),
    "USGS": ("https://earthquake.usgs.gov/fdsnws/event/1/query", "csv"),
}

UA = {"User-Agent": "TSI-replication/1.0 (jim@jimgale.net; EQ-1 frozen protocol)"}


def fetch(url: str, tries: int = 4) -> tuple[int, bytes]:
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120, context=SSL_CTX) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return 204, b""
            if e.code in (429, 503) and attempt < tries - 1:
                time.sleep(10 * (attempt + 1))
                continue
            return e.code, e.read() if hasattr(e, "read") else b""
        except Exception:
            if attempt < tries - 1:
                time.sleep(5 * (attempt + 1))
                continue
            raise
    return -1, b""


def main() -> None:
    span_start = int(PARAMS["time_span"]["start"][:4])
    span_end_date = date.fromisoformat(PARAMS["time_span"]["end"])
    minmag = PARAMS["magnitude"]["primary_floor"]
    log_lines = []

    for key, region in PARAMS["regions"].items():
        base, fmt = ENDPOINTS[region["primary_source"]]
        outdir = HERE / "data" / "raw" / key
        outdir.mkdir(parents=True, exist_ok=True)
        total = 0
        for year in range(span_start, span_end_date.year + 1):
            start = f"{year}-01-01"
            end = f"{year + 1}-01-01" if year < span_end_date.year else PARAMS["time_span"]["end"]
            q = {
                "starttime": start,
                "endtime": end,
                "minlatitude": region["min_lat"],
                "maxlatitude": region["max_lat"],
                "minlongitude": region["min_lon"],
                "maxlongitude": region["max_lon"],
                "minmagnitude": minmag,
                "format": fmt,
                "orderby": "time-asc" if region["primary_source"] != "USGS" else "time-asc",
            }
            url = f"{base}?{urllib.parse.urlencode(q)}"
            status, body = fetch(url)
            outfile = outdir / f"{year}.{ 'csv' if fmt == 'csv' else 'txt' }"
            outfile.write_bytes(body)
            # count data rows (skip headers/comments)
            n = sum(
                1
                for ln in body.decode("utf-8", "replace").splitlines()
                if ln and not ln.startswith("#") and not ln.lower().startswith(("time", "eventid", "publicid"))
            )
            total += n
            log_lines.append(f"{key}\t{year}\tHTTP {status}\t{n} events\t{url}")
            time.sleep(1.0)  # politeness
        log_lines.append(f"{key}\tTOTAL\t\t{total} events")
        print(f"{key}: {total} events")

    ts = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    (HERE / "data" / f"retrieval_log_{ts.replace(':', '')}.tsv").write_text(
        "\n".join(log_lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
