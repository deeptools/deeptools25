import requests
from pathlib import Path
from multiprocessing import Pool
import gzip
import shutil
import hashlib

odir = Path(snakemake.params.odir)
odir.mkdir(exist_ok=True, parents=True)

record_id = snakemake.params.zenodo_id
api_url = f"https://zenodo.org/api/records/{record_id}"

r = requests.get(api_url, timeout=30)
if r.status_code != 200:
    raise RuntimeError(f"Zenodo API returned {r.status_code} for {api_url}:\n{r.text[:500]}")
r = r.json()

def download_file(r, targetfile):
    for f in r["files"]:
        if f["key"].replace('.gz', '').replace('benchmark_', '') == Path(targetfile).name:
            key = f["key"]
            url = f["links"]["self"]
            checksum = f["checksum"]

            # download
            of = odir / key.replace('benchmark_', '')
            with requests.get(url, stream=True, timeout=60) as targetfile:
                with open(of, "wb") as out:
                    for chunk in targetfile.iter_content(chunk_size=4096):
                        out.write(chunk)

            # Verify checksum
            hash_md5 = hashlib.md5()
            with open(of, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            got = hash_md5.hexdigest()
            assert got == checksum.replace('md5:', ''), f"MD5 mismatch for {key}: expected {checksum}, got {got}"
            (Path(of).parent / (Path(of).name.replace('.gz', '') + '.valid')).touch()

            if key.endswith(".gz"):
                with gzip.open(of, 'rb') as f_in:
                    with open(odir / key.replace('.gz', ''), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                of.unlink()
            return
    raise ValueError(f"Target file {targetfile} not found")

download_file(r, snakemake.output.ofile)
