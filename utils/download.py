from io import BytesIO
from typing import List
from zipfile import ZipFile
from urllib.request import urlopen
import os


def download_and_extract(remote_model_url: str, local_model_dir) -> List[str]:
    resp = urlopen(remote_model_url)
    os.makedirs(local_model_dir, exist_ok=True)
    with ZipFile(BytesIO(resp.read())) as zip_file:
        all_files_and_dirs = zip_file.namelist()
        zip_file.extractall(local_model_dir)
    return all_files_and_dirs
