from rvdata.core.models.level2 import RV2
from tqdm.auto import tqdm
from os import makedirs, listdir
from os.path import basename, dirname, join, exists
from glob import glob
from shutil import copyfile
import pandas as pd
import sqlalchemy as db
from astropy.io import fits 
import numpy as np

ENGINE_RV = db.create_engine("postgresql://solar:lowell@10.10.115.134/solar")

BASE_DIR = '/Volumes/expres_solar/extracted/'
PYRH_DIR = '/Volumes/expres_solar/pyrheliometer/'
ORIGINAL_DIR = join(BASE_DIR, 'fitspec')
TRANSLATED_DIR = join(BASE_DIR, 'translated')

def main():
    dates = [basename(_) for _ in glob(join(ORIGINAL_DIR, '2025', '2503*_solar'))]
    for date in tqdm(dates):
        pyrh_file = join(PYRH_DIR, f"20{date[0:2]}", f"expres_ljpyrheliometer_20{date.split('_solar')[0]}.csv")
        files = glob(join(ORIGINAL_DIR, f"20{date[:2]}", date,  'Sun*.fits'))
        if len(files) > 0:
            out_dir = join(TRANSLATED_DIR, f"20{date[0:2]}", f"20{date.split('_solar')[0]}")
            makedirs(out_dir, exist_ok=True)
            if exists(pyrh_file):
                copyfile(pyrh_file, join(out_dir, basename(pyrh_file)))
            for fh in tqdm(files):
                # out_fh = join(TRANSLATED_DIR, f"20{date[:2]}", date, basename(fh).replace('.fits', '_expres_lvl2.fits'))
                out_fh = join(out_dir, basename(fh).replace('Sun_', 'Sun_20').replace('.fits', '_expres_lvl2.fits') )
                try:
                    rv2_obj = RV2.from_fits(fh, instrument='EXPRES')
                    tqdm.write(f"Translating: {basename(fh)} to {out_fh}")
                    rv2_obj.to_fits(out_fh)
                except:
                    with open('log.txt', 'a+') as logfh:
                        logfh.write(f"Failed to convert: {fh}\n")
                    pass
            else:
                continue
        
def checkout_list():
    summary_lines = []

    summary_lines.append("Solar Data Checkout Summary")
    summary_lines.append(f"Prepared on: {__import__('datetime').date.today()}\n")

    total_fits = 0
    total_csv = 0
    dates = glob(join(TRANSLATED_DIR,  '2*'))
    dates = glob('../expres_standard_test/2*')
    obs_ids = []
    for day_dir in sorted(dates):
        fits_files = glob(join(day_dir, "*.fits"))
        obs_ids.append([basename(fh).split("Sun_20")[1][:11] for fh in fits_files])
        csv_files = glob(join(day_dir, "*.csv"))
        total_fits += len(fits_files)
        total_csv += len(csv_files)

        summary_lines.append(f"{basename(day_dir)}/")
        summary_lines.append(f"    Number of FITS files: {len(fits_files)}:")
        for fh in fits_files:
            summary_lines.append(f"        {basename(fh)}")
        summary_lines.append(f"    Pyrheliometer file: {basename(csv_files[0]) if csv_files else 'Missing'}")
        summary_lines.append(f"    Notes:\n")

    summary_lines.append(f"Total Dates Included: {len(dates)}")
    summary_lines.append(f"Total Solar FITS Files: {total_fits}")
    summary_lines.append(f"Total Pyrheliometer CSV Files: {total_csv} \n\n")
    summary_lines.append("Data format:\n")
    summary_lines.append("- FITS files: EPRV Standard Level 2 EXPRES spectra")
    summary_lines.append("- CSV files: Time [MST], Irradiance")
    summary_fh = join(TRANSLATED_DIR, 
                      f"expres_summary_{min([basename(d) for d in dates])}-{max([basename(d) for d in dates])}.txt")
    
    flat_list = [float(item) for sublist in obs_ids for item in sublist]

    query = """
    SELECT
        obsinfo.obs_id,
        obsinfo.tobs,
        obsinfo.doy,
        obsinfo.airmass,
        obsinfo.quality,
        obsinfo.exptime,
        ccf_vel.ccf_mnvel / 100 AS ccf_mnvel,
        ccf_vel.ccf_errvel / 100 AS ccf_errvel
    FROM
        obsinfo
        RIGHT JOIN ccf_vel ON ccf_vel.obs_id = obsinfo.obs_id
    WHERE
        obsinfo.obs_id = ANY (%s)
    """

    df = pd.read_sql(query, con=ENGINE_RV, params=(flat_list,))
    if len(df) != len(flat_list):
        df = df[df['obs_id'].isin(flat_list)]
        found_ids = set(df['obs_id'])
        missing_ids = [obs_id for obs_id in flat_list if obs_id not in found_ids]

    quality = np.zeros(len(df))
    file_names = []
    for (jj, x) in df.iterrows():
        quality[jj] = fits.getval(join(BASE_DIR, 'fitspec',
                                       f"20{str(x['obs_id'])[:2]}",
                                       f"{str(x['obs_id'])[:6]}_solar",
                                       f"Sun_{x['obs_id']:.4f}.fits"), 'QUALITY', ext=1)
        file_names.append(f"Sun_20{x['obs_id']:.4f}_expres_lvl2.fits")
    quality[quality < 0] = 0
    df['quality'] = quality
    df['file'] = file_names
    df[['file', 'tobs','airmass', 'quality', 'exptime', 'ccf_mnvel',
       'ccf_errvel']].to_csv(summary_fh.replace("summary", "rv").replace('.txt','.csv'))
    
    
    with open(summary_fh, "w") as f:
        f.write("\n".join(summary_lines))
   

if __name__=='__main__':
    main()