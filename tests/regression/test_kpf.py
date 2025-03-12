import requests
import os
from core.models.level2 import RV2
from instruments.kpf.level2 import KPFRV2
from tests.regression.compliance import check_l2_extensions, check_l2_header

file_urls = {'KPF': ['https://zenodo.org/records/15009170/files/KP.20240715.46231.07.fits?download=1',
                     'https://zenodo.org/records/15009170/files/KP.20240715.46231.07_L1.fits?download=1']}

def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, 'wb') as file:
        file.write(response.content)

def test_kpf():
    l0file = 'kpf_L0.fits'
    l1file = 'kpf_L1.fits'
    if not os.path.exists(l0file):
        download_file(file_urls['KPF'][0], l0file)
    if not os.path.exists(l1file):
        download_file(file_urls['KPF'][1], l1file)

    kpf = RV2.from_fits(l1file, l0file=l0file, instrument="KPF")
    standard_out = './kpf_L2_standard.fits'
    kpf.to_fits(standard_out)
    l2_standard = RV2.from_fits(standard_out)

    check_l2_extensions(standard_out)
    # check_l2_header(l2_standard.headers['PRIMARY'])

if __name__ == "__main__":
    test_kpf()
