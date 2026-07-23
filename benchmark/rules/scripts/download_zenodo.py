import requests
from pathlib import Path
from multiprocessing import Pool
import gzip
import shutil

odir = Path(snakemake.params.odir)
odir.mkdir(exist_ok=True, parents=True)

record_id = snakemake.params.zenodo_id
api_url = f"https://zenodo.org/api/records/{record_id}"

r = requests.get(api_url, timeout=30)
if r.status_code != 200:
    raise RuntimeError(f"Zenodo API returned {r.status_code} for {api_url}:\n{r.text[:500]}")
r = r.json()
targets = [(f["key"], f["links"]["self"]) for f in r["files"] if 'triticum' in f["key"] or 'human' in f["key"]]
targets += [(f["key"], f["links"]["self"]) for f in r["files"] if 'benchmark' in f["key"]]


def download_file(tup):
    key, url = tup
    print(f"Downloading {key}")
    of = odir / key.replace('benchmark_', '')
    with requests.get(url, stream=True) as cramfile:
        with open(of, "wb") as out:
            for chunk in cramfile.iter_content(chunk_size=8192):
                out.write(chunk)
    if key.endswith(".gz"):
        with gzip.open(of, 'rb') as f_in:
            with open(odir / key.replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        of.unlink()

with Pool(snakemake.threads) as p:
    p.map(download_file, targets)
