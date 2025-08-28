import requests
from pathlib import Path
from multiprocessing import Pool
import gzip
import shutil

odir = Path(snakemake.params.odir)
odir.mkdir(exist_ok=True, parents=True)

record_id = snakemake.params.zenodo_id
api_url = f"https://zenodo.org/api/records/{record_id}"

r = requests.get(api_url).json()
targets = [(f["key"], f["links"]["self"]) for f in r["files"] if 'example' in f["key"] or 'mouse.fna.gz' in f["key"]]

def download_file(tup):
    key, url = tup
    print(f"Downloading {key}")
    of = odir / key.replace('example_', 'NPC_')
    with requests.get(url, stream=True) as cramfile:
        with open(of, "wb") as out:
            for chunk in cramfile.iter_content(chunk_size=8192):
                out.write(chunk)
    if key == 'mouse.fna.gz':
        with gzip.open(of, 'rb') as f_in:
            with open(odir / 'mouse.fna', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        of.unlink()

with Pool(snakemake.threads) as p:
    p.map(download_file, targets)
