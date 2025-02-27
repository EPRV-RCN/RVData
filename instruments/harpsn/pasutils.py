'''
RVData/instruments/harpsn/utils.py

UNIGE-ESO - EPRV
Author: Loris JACQUES
Created: Mon Jan 20 2025
Last Modified: Mon Jan 20 2025
Version: 1.0.0
Description: 
'''

'''
---------------------
external libraries
---------------------
'''
from astropy.io import fits
from astropy.constants import c
import os
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import numpy as np
import pandas as pd
from astropy.time import Time
import math
from datetime import datetime
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body, get_body_barycentric_posvel
from astropy import units as u
'''
---------------------
internal libraries
---------------------
'''
from core.models.level2 import RV2
import instruments.harpsn.config.config as config

# TODO get the file adaptativ to every configuration of fibers number
def get_files_names(full_path:str) -> dict:
    """
    This function retrieves the names of related FITS files based on a given file's path
    and constructs a dictionary containing paths to these files.

    Args:
        full_path (str): The full path to the raw file.

    Returns:
        dict: A dictionary with keys representing the file types and values containing their respective paths.
    """
    # Get the directory path and base file name from the full path
    repo_path = os.path.dirname(full_path)
    base_file_name = os.path.basename(full_path)
    
    # Construct paths for the S2D and BLAZE FITS files (both A and B versions)
    s2d_blaze_file_A = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_A.fits')
    s2d_blaze_file_B = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_B.fits')
    drift_file_B = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_DRIFT_MATRIX_B.fits')


    if not os.path.isfile(drift_file_B):
        with fits.open(full_path) as hdu_raw:
            dpr_type = hdu_raw['PRIMARY'].header['HIERARCH TNG DPR TYPE']
            if dpr_type.split(",")[1] == 'SKY':
                print('SKY type doesn\'t have any DRIFT correction')
                drift_file_B = None
            elif dpr_type.split(",")[1] == 'DARK':
                print('DARK type doesn\'t have any DRIFT correction')
                drift_file_B = None
            else:
                print('ERROR: NO DRIFT FILE FOUND')
                return

    # Open the S2D BLAZE FITS file (A version) to retrieve the BLAZE file names
    # These names are stored in specific header fields: HIERARCH ESO PRO REC1 CALn NAME
    with fits.open(s2d_blaze_file_A) as hdul:
        for i in hdul['PRIMARY'].header['ESO PRO REC1 CAL* CATG']:
            if 'BLAZE_A' == hdul['PRIMARY'].header[i]:
                blaze_file_A = os.path.join(
                    repo_path, 
                    hdul["PRIMARY"].header[i[:-4]+'NAME']
                )
            if 'BLAZE_B' == hdul['PRIMARY'].header[i]:
                blaze_file_B = os.path.join(
                    repo_path, 
                    hdul["PRIMARY"].header[i[:-4]+'NAME']
                )
    
    # Construct a dictionary of all the file paths
    # `full_path` corresponds to the raw file path (the input to this function)
    names = {
        "raw_file": full_path,
        "s2d_blaze_file_A": s2d_blaze_file_A, 
        "s2d_blaze_file_B": s2d_blaze_file_B,
        "blaze_file_A": blaze_file_A,
        "blaze_file_B": blaze_file_B,
        "drift_file_B": drift_file_B
    }
    return names


def convert_S2D_BLAZE(RV2: RV2, file_path: str, trace_ind_start: int, slice_nb: int) -> None:
    """_summary_
    get usufull data from S2D_BLAZE.fits files
    """

    with fits.open(file_path) as hdul:
        # A ne faire qu'une fois donc à la premiere itération
        if(trace_ind_start == 1):
            # On fixe le header de l'extension INSTRUMENT_HEADER avec le header de S2D_BLAZE_A
            RV2.set_header("INSTRUMENT_HEADER", hdul["PRIMARY"].header)

            # On récupère la valeur de BARYCORR_KMS
            barycorr_kms_data = hdul["PRIMARY"].header['HIERARCH TNG QC BERV']

            # On calcul la valeur de BARYCORR_Z en faisant attention à ce que la vitesse de la lumière soit de mêmes unitées que les km/s 
            barycorr_z_data = (barycorr_kms_data/(c.to('km/s'))).value
            
            # On récupère la valeur de BJD_TBD
            bjd_tdb_data = hdul["PRIMARY"].header['HIERARCH TNG QC BJD']

            # TODO Create a function to concat this :
            # On crée les 3 extensions 
            barycorr_kms_hdu = fits.ImageHDU(data = np.ones(1)*barycorr_kms_data)
            barycorr_kms_hdu.header['EXTNAME'] = 'BARYCORR_KMS'
            barycorr_kms_hdu.header['CTYPE1'] = ('BARYCORR_KMS', 'Name of axis 1')
            barycorr_kms_hdu.add_datasum()

            barycorr_z_hdu = fits.ImageHDU(data = np.ones(1)*barycorr_z_data)
            barycorr_z_hdu.header['EXTNAME'] = 'BARYCORR_Z'
            barycorr_z_hdu.header['CTYPE1'] = ('BARYCORR_Z', 'Name of axis 1')
            barycorr_z_hdu.add_datasum()

            bjd_tdb_hdu = fits.ImageHDU(data = np.ones(1)*bjd_tdb_data)
            bjd_tdb_hdu.header['EXTNAME'] = 'BJD_TDB'
            bjd_tdb_hdu.header['CTYPE1'] = ('BJD_TDB', 'Name of axis 1')
            bjd_tdb_hdu.add_datasum()
            
            # On modifie les extensions qui avaient déja été crée par le constructeur de classe
            RV2.set_header("BARYCORR_KMS", barycorr_kms_hdu.header)
            RV2.set_data("BARYCORR_KMS", barycorr_kms_hdu.data)

            RV2.set_header("BARYCORR_Z", barycorr_z_hdu.header)
            RV2.set_data("BARYCORR_Z", barycorr_z_hdu.data)
            
            RV2.set_header("BJD_TDB", bjd_tdb_hdu.header)
            RV2.set_data("BJD_TDB", bjd_tdb_hdu.data)
        
        for field in config.extnames.keys():
            for slice in range(1, slice_nb+1):
                #We extract the values of the specific slice
                single_cam_values = hdul[field].data[slice-1::slice_nb,:]

                # We need to remove barycentric correction from the datas that have it because format need to be without correction
                if('BARY' in field):
                    single_cam_values = doppler_shift(single_cam_values, RV2.data["BARYCORR_KMS"][0])

                #We create a new HDU with the values of the slice and the original header
                hdu_l2 = fits.ImageHDU(data = single_cam_values.copy(), header = hdul[field].header)

                # And we replace or add new elements
                hdu_l2.header['EXTNAME'] = 'TRACE'+str(trace_ind_start+slice-1)+config.extnames[field]
                hdu_l2.header['CTYPE1'] = (config.extnames[field][1:], 'Name of axis 1')
                hdu_l2.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

                # We add the headers and the datas to the extensions
                if(hdu_l2.header['EXTNAME'] not in RV2.extensions):
                    RV2.create_extension(ext_name = hdu_l2.header['EXTNAME'], ext_type = 'ImageHDU', header = hdu_l2.header, data = hdu_l2.data)
                else:
                    RV2.set_data(ext_name = hdu_l2.header['EXTNAME'], data = hdu_l2.data)
                    RV2.set_header(ext_name = hdu_l2.header['EXTNAME'], header = hdu_l2.header)


def create_PRIMARY(RV2: RV2, names: list[str], nb_trace: int, nb_fiber: int):
    """Creates the L2 header by copying the information from the raw file and adding the necessary information for the L2 file.
    """
    #We create an empty HDU to store the L2 Primary header
    l2_hdu = fits.PrimaryHDU(data = None)
    l2_hdu.header['EXTNAME'] = 'PRIMARY'
    #We load the header map to convert between raw file headers and L2 header
    header_map = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/config/header_map.csv')
    
    for index, values in header_map.iterrows():
        if(header_map['skip'].iloc[index] == True):
            continue
        #Add the HIERARCH keyword to the header if the keyword is longer than 8 characters
        if(len(values.iloc[0]) > 8):
            values.iloc[0] = 'HIERARCH ' + values.iloc[0]
        values.iloc[0] = values.iloc[0].strip()

        try:
            #If there is a fixed value to set, we set it
            if(pd.notna(header_map['value'].iloc[index])):
                l2_hdu.header[values.iloc[0]] = (
                    header_map['value'].iloc[index],
                    header_map['Description'].iloc[index]
                )

            #Otherwise, we copy the value from the good file
            elif(pd.notna(header_map['TNG_keyword'].iloc[index])):
                if(header_map['from'].iloc[index]=='S2D_BLAZE_A'):
                    l2_hdu.header[values.iloc[0]] = (
                        RV2.headers['INSTRUMENT_HEADER'][header_map['TNG_keyword'].iloc[index]], 
                        header_map['Description'].iloc[index]
                    )
                elif(header_map['from'].iloc[index]=='RAW'):
                    with fits.open(names["raw_file"]) as hdu_raw:
                        l2_hdu.header[values.iloc[0]] = (
                            hdu_raw['PRIMARY'].header[header_map['TNG_keyword'].iloc[index]], 
                            header_map['Description'].iloc[index]
                        )
                elif(header_map['from'].iloc[index]=='CONFIG'):
                    l2_hdu.header[values.iloc[0]] = (
                        getattr(config, header_map['TNG_keyword'].iloc[index], None), 
                        header_map['Description'].iloc[index]
                    )

            #If the value is not present in the raw file, we set it to Null 
            else:
                l2_hdu.header[values.iloc[0]] = (
                'Null',
                header_map['Description'].iloc[index]
            )
                
        except Exception as e:
            l2_hdu.header[values.iloc[0]] = (
                'Null',
                header_map['Description'].iloc[index]
            )
            key = header_map['Keyword'].iloc[index]
            print(f'{e} Also named {key}.')
            
    
    # FILENAME KEYWORD
    l2_hdu.header['FILENAME'] = (
        'inst' 
        + config.data_format 
        + '_'
        + RV2.filename.split('.')[1].replace("-", "").replace("_", "")
        + '.fits', 
        header_map[header_map['Keyword'] == 'FILENAME']['Description'].iloc[0]
    )

    # Getting SIMBAD/GAIA Catalog datas
    catalog_data = get_simbad_data(RV2.headers['INSTRUMENT_HEADER'][header_map[header_map['Keyword'] == 'OBJECT']['TNG_keyword'].iloc[0]])
    if(catalog_data['CID']!='Null'):
        try:
            catalog_data['CCLR'] = get_gaia_data(catalog_data['CID'])
        except: 
            print('Gaia request failed.')
    
    catalog_data['CRA'] = RV2.headers['INSTRUMENT_HEADER'][header_map[header_map['Keyword']=='CRA']['TNG_keyword'].iloc[0]].replace("h", ":").replace("m", ":")

    add_keyword_cat = ['CDEC', 'CEQNX', 'CEPCH', 'CPMR', 'CPMD', 'CRV']
    for keyword in add_keyword_cat:
        catalog_data[keyword] = RV2.headers['INSTRUMENT_HEADER'][header_map[header_map['Keyword'] == keyword]['TNG_keyword'].iloc[0]]

    rv = catalog_data['CRV']
    rv_z = round(rv/(c/1e3).value, 8)
    catalog_data['CZ'] = rv_z

    # Keywords qui dependent du numéro de la TRACE
    keyword_list = ['CSRC', 'CID', 'CRA', 'CDEC', 'CEQNX', 'CEPCH', 'CPLX', 'CPMR', 'CPMD', 'CRV', 'CZ', 'CCLR']

    with fits.open(names["raw_file"]) as hdu_raw:
        dpr_type = hdu_raw['PRIMARY'].header['HIERARCH TNG DPR TYPE'].split(",")
        for i in range(1, nb_trace+1):
            if(dpr_type[math.ceil(i*nb_fiber/nb_trace)-1] == 'STAR'):
                l2_hdu.header['TRACE'+str(i)] = (
                    'SCI',
                    header_map[header_map['Keyword'] == 'TRACE']['Description'].iloc[0]
                )
            elif(dpr_type[math.ceil(i*nb_fiber/nb_trace)-1] == 'WAVE'):
                l2_hdu.header['TRACE'+str(i)] = (
                    'CAL',
                    header_map[header_map['Keyword'] == 'TRACE']['Description'].iloc[0]
                )
            elif(dpr_type[math.ceil(i*nb_fiber/nb_trace)-1] == 'DARK'):
                l2_hdu.header['TRACE'+str(i)] = (
                    'DARK',
                    header_map[header_map['Keyword'] == 'TRACE']['Description'].iloc[0]
                )
            elif(dpr_type[math.ceil(i*nb_fiber/nb_trace)-1] == 'SKY'):
                l2_hdu.header['TRACE'+str(i)] = (
                    dpr_type[math.ceil(i*nb_fiber/nb_trace)-1],
                    header_map[header_map['Keyword'] == 'TRACE']['Description'].iloc[0]
                )
            else:
                l2_hdu.header['TRACE'+str(i)] = (
                    'UNKNOWN',
                    header_map[header_map['Keyword'] == 'TRACE']['Description'].iloc[0]
                )
            
            # CALIBRATION SOURCE KEYWORD
            if(l2_hdu.header['TRACE'+str(i)]=='CAL'):
                l2_hdu.header['CLSRC'+str(i)] = (
                    RV2.headers['INSTRUMENT_HEADER'][header_map[header_map['Keyword']=='CLSRC']['TNG_keyword'].iloc[0]].split('_')[i-1], 
                    header_map[header_map['Keyword'] == 'CLSRC']['Description'].iloc[0]
                )
            else:
                l2_hdu.header['CLSRC'+str(i)] = (
                    'Null', 
                    header_map[header_map['Keyword'] == 'CLSRC']['Description'].iloc[0]
                )

            # CATALOG KEYWORDS
            if(l2_hdu.header['TRACE'+str(i)]=='SCI'):
                for keyword in keyword_list:
                    l2_hdu.header[keyword+str(i)] = (
                        catalog_data[keyword], 
                        header_map[header_map['Keyword'] == keyword]['Description'].iloc[0]
                    )
            else:
                for keyword in keyword_list:
                    l2_hdu.header[keyword+str(i)] = (
                        'Null', 
                        header_map[header_map['Keyword'] == keyword]['Description'].iloc[0]
                    )


    # OBSLON KEYWORD
    l2_hdu.header['OBSLON']  = (
        parse_geo_coord(RV2.headers['INSTRUMENT_HEADER']['GEOLON']), 
        header_map[header_map['Keyword'] == 'OBSLON']['Description'].iloc[0]
    )

    # OBSLAT KEYWORD
    l2_hdu.header['OBSLAT']  = (
        parse_geo_coord(RV2.headers['INSTRUMENT_HEADER']['GEOLAT']), 
        header_map[header_map['Keyword'] == 'OBSLAT']['Description'].iloc[0]
    )

    # BINNING KEYWORD
    # binx = str(RV2.headers['INSTRUMENT_HEADER']['HIERARCH TNG DET WIN1 BINX'])
    # biny = str(RV2.headers['INSTRUMENT_HEADER']['HIERARCH TNG DET WIN1 BINY'])
    binx = str(1)
    biny = str(1)
    l2_hdu.header['BINNING']  = (
        f"{binx}x{biny}", 
        header_map[header_map['Keyword'] == 'BINNING']['Description'].iloc[0]
    )

    # NUMTRACE KEYWORD
    l2_hdu.header['NUMTRACE']  = (
        nb_trace, 
        header_map[header_map['Keyword'] == 'NUMTRACE']['Description'].iloc[0]
    )

    # DATE KEYWORD
    current_time = Time.now()
    l2_hdu.header['DATE']  = (
        current_time.iso, 
        header_map[header_map['Keyword'] == 'DATE']['Description'].iloc[0]
    )

    # JD_UTC KEYWORD 
    l2_hdu.header['JD_UTC'] = (
        RV2.headers['INSTRUMENT_HEADER']['MJD-OBS'] + 2400000.5,
        header_map[header_map['Keyword'] == 'JD_UTC']['Description'].iloc[0]
    )
    
    # TLST KEYWORD
    l2_hdu.header['TLST'] = (
        convert_lst(RV2.headers['INSTRUMENT_HEADER']['LST']), 
        header_map[header_map['Keyword'] == 'TLST']['Description'].iloc[0]
    )

    # TRA KEYWORD
    l2_hdu.header['TRA'] = (
        deg_to_sexagesimal(RV2.headers['INSTRUMENT_HEADER']['RA-DEG'],True), 
        header_map[header_map['Keyword'] == 'TRA']['Description'].iloc[0]
    )

    # TDEC KEYWORD
    l2_hdu.header['TDEC'] = (
        deg_to_sexagesimal(RV2.headers['INSTRUMENT_HEADER']['DEC-DEG'],False), 
        header_map[header_map['Keyword'] == 'TDEC']['Description'].iloc[0]
    )

    # TEL KEYWORD
    l2_hdu.header['TEL'] = (
        RV2.headers['INSTRUMENT_HEADER']['EL']*180/np.pi, 
        header_map[header_map['Keyword'] == 'TEL']['Description'].iloc[0]
    )
    
    # TZA KEYWORD
    l2_hdu.header['TZA'] = (
        np.round(90 - l2_hdu.header['TEL'], 3), 
        header_map[header_map['Keyword'] == 'TZA']['Description'].iloc[0]
    )

    # THA KEYWORD
    l2_hdu.header['THA'] = (
        compute_hour_angle(l2_hdu.header['TLST'], l2_hdu.header['TRA']), 
        header_map[header_map['Keyword'] == 'THA']['Description'].iloc[0]
    )

    # SUNEL/MOONANG/MOONEL/MOONILLU/MOONRV KEYWORDS
    moon_sun_params = get_moon_sun_info(
        RV2.headers['INSTRUMENT_HEADER']['RA-DEG'],
        RV2.headers['INSTRUMENT_HEADER']['DEC-DEG'],
        l2_hdu.header['OBSLAT'], 
        l2_hdu.header['OBSLON'],
        l2_hdu.header['OBSALT'],
        l2_hdu.header['DATE-OBS'],
        l2_hdu.header['JD_UTC']
    )

    # List of corresponding keywords
    moon_sun_keywords = ['SUNEL', 'MOONANG', 'MOONEL', 'MOONILLU', 'MOONRV']

    # Assign values to headers dynamically
    for key, value in zip(moon_sun_keywords, moon_sun_params):
        l2_hdu.header[key] = (
            value, 
            header_map[header_map['Keyword'] == key]['Description'].iloc[0]
        )

    # INSTERA KEYWORD
    l2_hdu.header['INSTERA'] = (
        get_instrument_version(l2_hdu.header['DATE-OBS']), 
        header_map[header_map['Keyword'] == 'INSTERA']['Description'].iloc[0]
    )

    # EXSNR-N KEYWORD
    for i in range(1, int(l2_hdu.header['NUMORDER'])+1):
        l2_hdu.header[f'EXSNR{str(i)}'] = (
            RV2.headers['INSTRUMENT_HEADER'][f"HIERARCH TNG QC ORDER{str(i)} SNR"], 
            header_map[header_map['Keyword'] == 'EXSNR']['Description'].iloc[0]
        )

    # EXSNRW-N KEYWORD
    for i in range(int(l2_hdu.header['NUMORDER'])):
        l2_hdu.header[f'EXSNRW{str(i+1)}'] = (
            round(RV2.data["TRACE1_WAVE"][i,0]+(RV2.data["TRACE1_WAVE"][i,-1] - RV2.data["TRACE1_WAVE"][i,0])/2), 
            header_map[header_map['Keyword'] == 'EXSNRW']['Description'].iloc[0]
        )

    # COLOFLAG KEYWORD
    try:
        color_flag = RV2.headers['INSTRUMENT_HEADER'][header_map[header_map['Keyword'] == 'COLOFLAG']['TNG_keyword'].iloc[0]]
        if(color_flag == 1):
            l2_hdu.header['COLOFLAG'] = (
                'Pass', 
                header_map[header_map['Keyword'] == 'COLOFLAG']['Description'].iloc[0]
            )
        else:
            l2_hdu.header['COLOFLAG'] = (
                'Fail', 
                header_map[header_map['Keyword'] == 'COLOFLAG']['Description'].iloc[0]
            )
    except:
        l2_hdu.header['COLOFLAG'] = (
            'Fail', 
            header_map[header_map['Keyword'] == 'COLOFLAG']['Description'].iloc[0]
        )

    # SUMMFLAG KEYWORD
    flags = ["TELFLAG", "INSTFLAG", "DRPFLAG", "COLOFLAG", "OBSFLAG"]

    # Récupérer toutes les valeurs des flags
    flag_values = [l2_hdu.header.get(flag, "Pass") for flag in flags]

    # Priorité des états : Fail > Warn > Pass
    if "Fail" in flag_values:
        l2_hdu.header['SUMMFLAG'] = (
            "Fail", 
            header_map[header_map['Keyword'] == 'SUMMFLAG']['Description'].iloc[0]
        )
    elif "Warn" in flag_values:
        l2_hdu.header['SUMMFLAG'] = (
            "Warn", 
            header_map[header_map['Keyword'] == 'SUMMFLAG']['Description'].iloc[0]
        )
    else:
        l2_hdu.header['SUMMFLAG'] = (
            "Pass", 
            header_map[header_map['Keyword'] == 'SUMMFLAG']['Description'].iloc[0]
        )

    if('PRIMARY' not in RV2.extensions):
        RV2.create_extension(ext_name = 'PRIMARY', ext_type = 'PrimaryHDU', header = l2_hdu.header)
    else:
        RV2.set_header(ext_name = 'PRIMARY', header = l2_hdu.header)
    return
        

def convert_BLAZE(RV2: RV2, file_path: str, trace_ind_start: int, slice_nb: int) -> None:
    """
    This function processes a FITS file and converts its data into a single or multiple 'BLAZE' extensions, 
    which are then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'BLAZE' extensions will be added or updated.
        file_path (str): The path to the FITS file containing the data to process.
        trace_ind_start (int): The starting index for naming the extensions (e.g., TRACE<trace_ind_start>_BLAZE).
        slice_nb (int): The number of slices to extract from the FITS file, each corresponding to one extension.

    Returns:
        None: Modifies the RV2 object in place by adding or updating the 'BLAZE' extensions.
    """

    with fits.open(file_path) as hdul:
        # Loop through each slice from 1 to slice_nb
        for slice in range(1, slice_nb+1):
            # Extract the corresponding data for this slice
            # Each slice takes every slice_nb-th row starting from (slice-1)
            blaze_hdu = fits.ImageHDU(
                data = hdul[1].data[slice-1::slice_nb,:],
                header = hdul[1].header
            )

            # Update the header of the new HDU with relevant metadata
            blaze_hdu.header['EXTNAME'] = 'TRACE'+str(trace_ind_start+slice-1)+'_BLAZE'
            blaze_hdu.header['CTYPE1'] = ('Pixels', 'Name of axis 1')
            blaze_hdu.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

            # Check if the extension already exists in the RV2 object
            if(blaze_hdu.header['EXTNAME'] not in RV2.extensions):
                # If the extension does not exist, create it
                RV2.create_extension(
                    ext_name = blaze_hdu.header['EXTNAME'], 
                    ext_type = 'ImageHDU', 
                    header = blaze_hdu.header, 
                    data = blaze_hdu.data
                )
            else:
                # If the extension exists, update its data and header
                RV2.set_header(blaze_hdu.header['EXTNAME'], blaze_hdu.header)
                RV2.set_data(blaze_hdu.header['EXTNAME'], blaze_hdu.data)


def convert_DRIFT(RV2: RV2, file_path: str) -> None:
    """
    Processes a FITS file and converts its data into a 'DRIFT' extension, 
    which is then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'DRIFT' extension will be added or updated.
        file_path (str): The path to the FITS file containing the drift data.

    Returns:
        None: The function modifies the RV2 object in place by adding or updating the 'DRIFT' extension.
    """
    if(file_path != None):
        with fits.open(file_path) as hdul:
            # Extract drift data from the FITS file (2nd HDU)
            drift_data = hdul[1].data

            drift_hdu = fits.ImageHDU(
                data = drift_data,
                header = hdul[1].header
            )
    else:
        # If no file is provided, create an empty ImageHDU with default dimensions
        # This case occurs when Fiber B is SKY or DARK.
        drift_hdu = fits.ImageHDU(
            data = np.zeros((config.NUMORDER, config.num_pixel))
        )

    # Update the header with relevant metadata
    drift_hdu.header['EXTNAME'] = 'DRIFT'
    drift_hdu.header['CTYPE1'] = ('Pixels', 'Name of axis 1')
    drift_hdu.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

    # Check if the extension already exists in the RV2 object
    if(drift_hdu.header['EXTNAME'] not in RV2.extensions):
        # If the extension does not exist, create it
        RV2.create_extension(
            ext_name = drift_hdu.header['EXTNAME'], 
            ext_type = 'ImageHDU', 
            header = drift_hdu.header, 
            data = drift_hdu.data
        )
    else:
        # If the extension exists, update its data and header
        RV2.set_header(drift_hdu.header['EXTNAME'], drift_hdu.header)
        RV2.set_data(drift_hdu.header['EXTNAME'], drift_hdu.data)
        

def get_simbad_data(obj: str) -> dict:
    """
    Retrieves astrometric data for a given object from the Simbad database.

    Args:
        obj (str): The name of the astronomical object to search in Simbad.

    Returns:
        dict: A dictionary containing the following keys:
        - 'CSRC': The catalog source (Gaia DR3 or DR2).
        - 'CID': The Gaia identifier of the object.
        - 'CPLX': The parallax value (in milliarcseconds).
        - If no data is found, default values ('Null') are assigned.

    Raises:
        Exception: If the Simbad query fails or if the object is not found.
    """

    data = {}

    try:
        # Configure Simbad with custom settings
        custom_simbad = Simbad()
        custom_simbad.TIMEOUT = config.timeout  # Increase timeout if needed
        custom_simbad.add_votable_fields('ids', 'plx')

        # Query Simbad for the object
        result = custom_simbad.query_object(obj)

        # Extract Gaia DR3 or DR2 identifiers
        for name in result['IDS'][0].split('|'):
            if(name.lower().startswith('gaia dr3')):
                gaia_dr3_source = name[:8]
                gaia_dr3_name = name[5:]
            elif(name.lower().startswith('gaia dr2')):
                gaia_dr2_source = name[:8]
                gaia_dr2_name = name[5:]

        # Prioritize DR3 over DR2
        if gaia_dr3_name:
            data['CSRC'] = gaia_dr3_source
            data['CID'] = gaia_dr3_name
        elif gaia_dr2_name:
            data['CSRC'] = gaia_dr2_source
            data['CID'] = gaia_dr2_source

        # Retrieve parallax value
        if not np.ma.is_masked(result['PLX_VALUE'][0]):
            data['CPLX'] = result['PLX_VALUE'][0]
        else:
            data['CPLX'] = 'Null'

        return data
    
    except Exception as e:
        print(f"Catalog not found for the name of {obj}, err:{e}")

        # Return default values if the object is not found
        cat_list = ['CSRC', 'CID', 'CPLX', 'CCLR']
        for key in cat_list:
            data[key] = 'Null'
        return data


def get_gaia_data(gaia_name: str) -> float:
    """
    Retrieves Gaia photometric data and computes the color index (BP - RP).

    Args:
        gaia_name (str): The Gaia source identifier in the format "DR3 <source_id>".

    Returns:
        color_br (float): The computed color index (BP - RP) for the specified Gaia source.
    """

    # Extract the numerical part of the Gaia source ID
    name = gaia_name[4:]

    # Construct the ADQL query to fetch Gaia DR3 photometric data
    query = f"""
    SELECT phot_bp_mean_mag, phot_rp_mean_mag 
    FROM gaiadr3.gaia_source
    WHERE source_id = '{name}'
    """

    # Execute the query using the Gaia TAP service
    job = Gaia.launch_job(query)
    result = job.get_results()

    # Compute and return the color index (BP - RP)
    color_br = result["phot_bp_mean_mag"][0] - result["phot_rp_mean_mag"][0]
    return color_br


def parse_geo_coord(coord_str: str) -> float:
    """
    Convert a coordinate string (latitude or longitude) in DMS format to decimal degrees.
    
    Args:
        coord_str (str): Coordinate string in the format 'DD MM SS.s D', where D is the direction (N, S, E, W).
    
    Returns:
        decimal_degrees (float): Coordinate in decimal degrees.
    """
    
    # Séparer les composants de la chaîne
    parts = coord_str.split()
    degrees = int(parts[0])
    minutes = int(parts[1])
    seconds = float(parts[2].rstrip("WENS"))  # Retirer la direction si attachée
    direction = parts[3].upper()  # Récupérer la direction (N, S, E, W)

    # Conversion en degrés décimaux
    decimal_degrees = degrees + minutes / 60 + seconds / 3600

    # Appliquer le signe pour S (latitude) ou W (longitude)
    if direction in ('S', 'W'):
        decimal_degrees = -decimal_degrees

    return decimal_degrees


def check_and_remove_empty_extensions(RV2: RV2) -> None:
    """
    Checks and removes empty extensions from the RV2 object.

    This function iterates through the extensions in the RV2 object 
    and removes those that are completely empty (`None`). 
    It does not remove extensions filled with zeros, such as 'drift',
    to avoid checksum errors.

    Args:
        RV2 (RV2): The RV2 object containing extensions to be checked.

    Returns:
        None: The function modifies the RV2 object in place.
    """

    # Ensure RV2 has an 'extensions' attribute and it's a dictionary
    if not hasattr(RV2, 'extensions') or not isinstance(RV2.extensions, dict):
        raise ValueError("RV2 object does not have accessible extensions.")
    
    # Iterate through RV2 extensions
    for ext_name, ext_value in RV2.data.items():
        # Check if the extension is empty
        if ext_value is None:
            print(f"Empty extension found, removing: {ext_name}")

            # Call the method to remove the extension
            RV2.del_extension(ext_name)
        else:
            print(f"Extension {ext_name} contains data.")


def doppler_shift(wave: np.ndarray, rv: float) -> np.ndarray:
    """
    Performs the doppler shift on the wavelength values.

    Args:
        wave (np.ndarray): the original wavelength values
        rv (float): the radial velocity of the object in km/s

    Returns:
        wave_shifted (np.ndarray): the doppler shifted wavelength values
    """

    wave_shifted = wave + wave*rv/(c/1e3).value
    return wave_shifted


def convert_lst(lst: float) -> str:
    """
    Converts the local sidereal time from seconds to hh:mm:ss.ms

    Args:
        lst (foat): the local sidereal time in seconds

    Returns:
        lst_sexa (str): the local sidereal time in the format hh:mm:ss.ms
    """

    if(lst < 0):
        res = convert_lst(-lst)
        return '-' + res
    
    hours = int(lst//3600)
    minutes = int((lst %3600)//60)
    seconds = int(lst %60)
    ms = round(lst % 1 * 1000)

    lst_sexa = f"{hours:02}:{minutes:02}:{seconds:02}.{ms:03}"
    return lst_sexa


def deg_to_sexagesimal(value_deg: float, is_ra: bool=False) -> str:
    """
    Converts a value in degrees to sexagesimal format (HH:MM:SS.sss or ±DD:MM:SS.sss).
    
    Args:
        value_deg (float): The value in degrees.
        is_ra (bool): If True, the value is treated as right ascension (RA).
                      If False, the value is treated as declination (Dec).
    
    Returns:
        sexagesimal_str (str): The converted value in sexagesimal format.
    """

    # Handle the signs for Dec (and ignore for RA)
    sign = '-' if value_deg < 0 and not is_ra else ''
    value_deg = abs(value_deg)
    
    # Conversion for RA (hours) or Dec (degrees)
    if is_ra:
        value_h = value_deg / 15
        hours = int(value_h)
        minutes = int((value_h - hours) * 60)
        seconds = ((value_h - hours) * 60 - minutes) * 60
        return f"{hours:02}:{minutes:02}:{seconds:06.3f}"
    else:
        degrees = int(value_deg)
        minutes = int((value_deg - degrees) * 60)
        seconds = ((value_deg - degrees) * 60 - minutes) * 60
        return f"{sign}{degrees:02}:{minutes:02}:{seconds:06.3f}"


def compute_hour_angle(lst: str, ra: str) -> str:
    """
    Computes the hour angle from the local sidereal time and the right ascension of the object.

    Args:
        lst (str): local sideral time in the format hh:mm:ss.ms
        ra (str): right ascension of the object in the format hh:mm:ss.ms

    Returns:
        str: the hour angle in the format hh:mm:ss.ms
    """

    h_ra = int(ra[:2])
    m_ra = int(ra[3:5]) + h_ra*60
    s_ra = float(ra[6:12]) + m_ra*60

    h_lst = int(lst[:2])
    m_lst = int(lst[3:5]) + h_lst*60
    s_lst = float(lst[6:12]) + m_lst*60

    s_ha = s_lst - s_ra
    return convert_lst(s_ha)


def get_instrument_version(date_obs_str: str) -> str:
    """
    Determines the version of the instrument based on the observation date.

    Args:
        date_obs_str (str): The observation date in ISO 8601 format ("YYYY-MM-DDTHH:MM:SS.sss").

    Returns:
        int: The version of the instrument corresponding to the observation date.

    Raises:
        ValueError: If the observation date does not fall within any defined version range.
    """
    # Convert the observation date string into a datetime object
    date_obs = datetime.fromisoformat(date_obs_str)

    # Iterate through the version ranges defined in the configuration
    for version_info in config.INSTRUMENT_VERSIONS:
        # Convert the start date to a datetime object
        start_date = datetime.fromisoformat(version_info["start_date"])

        # Convert the end date to a datetime object, or use datetime.max for open-ended ranges
        end_date = (
            datetime.fromisoformat(version_info["end_date"])
            if version_info["end_date"] is not None
            else datetime.max
        )

        # Check if the observation date falls within the current version's date range
        if start_date <= date_obs <= end_date:
            return version_info["version"]

    # If no version matches, raise an exception
    raise ValueError(f"No instrument version corresponds to the date {date_obs_str}")


def get_moon_sun_info(target_ra :float, target_dec: float, 
                  obs_lat: float, obs_lon: float, 
                  obs_alt: float, obs_time: str, 
                  jd_utc: float) -> list:
    """
    Calculates information about the Moon's position and its relationship 
    with a given target object based on input coordinates and observation 
    parameters. Calculates the elevation of the Sun above the horizon for 
    a given location and time.

    Args:
        target_ra (float): Right Ascension of the target object in degrees.
        target_dec (float): Declination of the target object in degrees.
        obs_lat (float): Latitude of the observation location in degrees.
        obs_lon (float): Longitude of the observation location in degrees.
        obs_alt (float): Altitude of the observation location in meters.
        obs_time (str): Observation time in ISO format (e.g., 'YYYY-MM-DD HH:MM:SS').
        jd_utc (float): Julian Date (UTC) for the observation.

    Returns:
        list: A list containing the following values:
            - sun_el (float): The elevation of the Sun at the observation site (in degrees).
            - moon_ang (float): The angular separation between the Moon and the target (in degrees).
            - moon_el (float): The elevation of the Moon at the observation site (in degrees).
            - moon_illu (float): The illumination of the Moon in percentage.
            - moon_rv (float): The radial velocity of reflected sunlight off moon (in km/s).
    """

    # Observer's location and observation time
    location = EarthLocation(lat=obs_lat, lon=obs_lon, height=obs_alt)
    time = Time(obs_time)

    # Get the Moon's position at the given time and location (GCRS)
    moon_coord = get_body("moon", time, location=location)
    
    # We need to convert moon_coord in ICRS to match target_coord ref
    moon_coord_icrs = moon_coord.transform_to("icrs")
    
    # Target's position
    target_coord = SkyCoord(ra=target_ra, dec=target_dec, frame='icrs', unit="deg", obstime=time)

    # Calculate the angular separation between the target and the Moon
    moon_ang = round(moon_coord_icrs.separation(target_coord).deg, 4)

    # Calculate the Moon's elevation above the horizon (GCRS here but could be ICRS, no impact)
    moon_altaz = moon_coord.transform_to(AltAz(obstime=time, location=location))
    moon_el = round(moon_altaz.alt.deg, 4)

    # Get the Sun's position and transform it to the observer's AltAz frame
    sun = get_body("sun", time, location=location)

    sun_altaz = sun.transform_to(AltAz(obstime=time, location=location))
    sun_el = sun_altaz.alt.deg

    # Calculate the Moon's illumination
    elongation = moon_coord.separation(sun)
    moon_phase_angle = np.arctan2(sun.distance*np.sin(elongation), moon_coord.distance - sun.distance*np.cos(elongation))
    moon_illu = round((1 + np.cos(moon_phase_angle)) / 2 * 100, 4).value

    # Calculate the RV of reflected sunlight off moon
    moon_rv = round(get_moon_velocity_in_target_direction(target_ra, target_dec, jd_utc), 4)
    print('sun_el', sun_el)
    print('moon_ang', moon_ang)
    print('moon_el', moon_el)
    print('moon_illu', moon_illu)
    print('moon_rv', moon_rv)

    return [sun_el, moon_ang, moon_el, moon_illu, moon_rv]


def get_moon_velocity_in_target_direction(alpha_deg: float, delta_deg: float, julian_day: float) -> float:
    """
    Compute the velocity of the Moon projected in the direction of a given target in the sky.

    Parameters:
        alpha_deg (float): Right Ascension (RA) of the target in degrees.
        delta_deg (float): Declination (Dec) of the target in degrees.
        julian_day (float): Julian date for the calculation.

    Returns:
        projected_velocity_km_s (float): Radial velocity of the Moon in km/s in the target direction.
    """

    # Convert Julian Day to Astropy Time object
    t = Time(julian_day, format='jd')

    # Get the Moon's barycentric position & velocity (in AU and AU/day) using JPL ephemeris
    moon_pos, moon_vel = get_body_barycentric_posvel('moon', t, ephemeris='jpl')

    # Extract velocity components (AU/day)
    x_vel, y_vel, z_vel = moon_vel.xyz.to_value(u.AU/u.day)

    # Convert target coordinates (RA, Dec) to radians
    alpha_rad = np.deg2rad(alpha_deg)
    delta_rad = np.deg2rad(delta_deg)

    # Compute projected radial velocity
    projected_velocity_au_per_day = (
        x_vel * np.cos(alpha_rad) * np.cos(delta_rad) +
        y_vel * np.sin(alpha_rad) * np.cos(delta_rad) +
        z_vel * np.sin(delta_rad)
    )

    # Convert AU/day to km/s
    au_per_day_to_km_s = 149597870.691 / 86400.0
    projected_velocity_km_s = projected_velocity_au_per_day * au_per_day_to_km_s

    return projected_velocity_km_s