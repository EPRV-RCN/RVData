from rvdata.core.models.level2 import RV2
from tqdm.auto import tqdm
from os import makedirs
from os.path import basename, dirname, join, exists
from glob import glob

def main():
    BASE_DIR = '/Volumes/expres_solar/extracted/'
    ORIGINAL_DIR = join(BASE_DIR, 'fitspec')
    TRANSLATED_DIR = join(BASE_DIR, 'translated')
    dates = [basename(_) for _ in glob(join(ORIGINAL_DIR, '2025', '2503*_solar'))]
    for date in tqdm(dates):
        files = glob(join(ORIGINAL_DIR, f"20{date[:2]}", date,  'Sun*.fits'))
        for fh in tqdm(files):
            out_fh = join(TRANSLATED_DIR, f"20{date[:2]}", date, basename(fh).replace('.fits', '_expres.fits'))
            # if not exists(out_fh):
            makedirs(dirname(out_fh), exist_ok=True)
            try:
                rv2_obj = RV2.from_fits(fh, instrument='EXPRES')
                print(out_fh)
                rv2_obj.to_fits(out_fh)
            except:
                with open('log.txt', 'a+') as logfh:
                    logfh.write(f"Failed to convert: {fh}\n")
                pass
        else:
            continue

if __name__=='__main__':
    main()