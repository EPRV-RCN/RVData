'''
RVData/rvdata/instruments/nirps/utils/create_PRIMARY.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wen Mar 07 2025
Last Modified: Wen Mar 07 2025
Version: 1.0.0
Description:
Extracts and processes the necessary data to create the PRIMARY header. It uses
the header_map.csv file to match ESO keywords with those of the new format.
The process works in two phases: first, automatically retrieving keywords that
are already in the correct format, and then modifying the remaining ones on a
case-by-case basis. Stores all the keywords in an `RV2` PRIMARY object.

---------------------
Libraries
---------------------
'''
from astropy.io import fits
from astropy.time import Time
from astropy.constants import c
from astropy.coordinates import (
    SkyCoord,
    EarthLocation,
    AltAz,
    get_body,
    get_body_barycentric_posvel
)
from astropy import units as u
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from datetime import datetime

import os
import pandas as pd
import math
import numpy as np

import rvdata.instruments.nirps.config.config as config
from rvdata.core.models.level2 import RV2


def create_PRIMARY(
        RV2: RV2, names: list[str], nb_trace: int, nb_slice: int
) -> None:
    """
    Create the PRIMARY HDU for the L2 FITS file by copying relevant metadata
    from different files and applying necessary transformations.

    This function reads a mapping file (`header_map.csv`) that specifies how
    to translate, copy, or compute header keywords for the L2 file. It then
    constructs a `PrimaryHDU` with the appropriate metadata.

    Parameters:
        RV2 (RV2): An instance of the RV2 class containing metadata and headers
            required for processing.
        names (list[str]): A dictionary mapping different file types (e.g., raw file)
            to their corresponding file paths.
        nb_trace (int): Number of traces in the dataset.
        nb_slice (int): Number of slices in the dataset.

    Returns:
        None : The function modifies the RV2 object in place by adding or
            updating the extnames_raw extensions.

    Notes:
    - If a keyword has `skip = True` in `header_map.csv`, it is not copied
      automatically but requires a specific computation.
    - If a keyword value is missing in the source file, it is set to `Null`.
    - Special handling is applied for some keyword.
    """

    # We create an empty HDU to store the L2 Primary header
    l2_hdu = fits.PrimaryHDU(data=None)
    l2_hdu.header['EXTNAME'] = 'PRIMARY'

    # Get the parent directory of the "utils" folder
    base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    # Properly construct the file path
    header_map_path = os.path.join(base_dir, "config", "header_map.csv")

    # Load the CSV file
    header_map = pd.read_csv(header_map_path)

    # We iterate through the header_map file to translate each keyword.
    for index, values in header_map.iterrows():
        # If the keyword has its skip value set to True, it is not copied
        # automatically but requires a specific calculation or conversion.
        if (bool(header_map['skip'].iloc[index]) is True):
            continue

        # Add the HIERARCH keyword to the header if the keyword is longer than
        # 8 characters
        if (len(values.iloc[0]) > 8):
            values.iloc[0] = 'HIERARCH ' + values.iloc[0]
        values.iloc[0] = values.iloc[0].strip()

        try:
            # If there is a fixed value to set, we set it
            if (pd.notna(header_map['value'].iloc[index])):
                l2_hdu.header[values.iloc[0]] = (
                    header_map['value'].iloc[index],
                    header_map['Description'].iloc[index]
                )

            # Otherwise, we copy the value from the good file
            elif (pd.notna(header_map['ESO_keyword'].iloc[index])):
                if (header_map['from'].iloc[index] == 'S2D_BLAZE_A'):
                    l2_hdu.header[values.iloc[0]] = (
                        RV2.headers['INSTRUMENT_HEADER'][
                            header_map['ESO_keyword'].iloc[index]
                        ],
                        header_map['Description'].iloc[index]
                    )
                elif (header_map['from'].iloc[index] == 'RAW'):
                    with fits.open(names["raw_file"]) as hdu_raw:
                        l2_hdu.header[values.iloc[0]] = (
                            hdu_raw['PRIMARY'].header[
                                header_map['ESO_keyword'].iloc[index]
                            ],
                            header_map['Description'].iloc[index]
                        )
                elif (header_map['from'].iloc[index] == 'CONFIG'):
                    l2_hdu.header[values.iloc[0]] = (
                        getattr(
                            config, header_map['ESO_keyword'].iloc[index], None
                        ),
                        header_map['Description'].iloc[index]
                    )

            # If the value is not present in the specified file , we set it to the
            # default value define in the header map
            else:
                l2_hdu.header[values.iloc[0]] = (
                    header_map['default_value'].iloc[index],
                    header_map['Description'].iloc[index]
                )
        except Exception as e:
            # If an error occurs (mostly due to the absence of the keyword in the
            # specified file), the value is set to the default value define in the
            # header map.
            l2_hdu.header[values.iloc[0]] = (
                header_map['default_value'].iloc[index],
                header_map['Description'].iloc[index]
            )
            key = header_map['Keyword'].iloc[index]
            print(f'{e} Also named {key}.')

    # Here, we handle each skipped keyword by applying specific
    # translations/conversions.

    # FILENAME KEYWORD
    if (os.name == 'nt'):
        l2_hdu.header['FILENAME'] = (
            'inst'
            + config.data_format
            + '_'
            + RV2.filename.split('.')[1].replace("-", "").replace("_", "")
            + '.fits',
            header_map[
                header_map['Keyword'] == 'FILENAME'
            ]['Description'].iloc[0]
        )
    else:
        l2_hdu.header['FILENAME'] = (
            'inst'
            + config.data_format
            + '_'
            + RV2.filename.split('.')[1].replace("-", "").replace(":", "")
            + '.fits',
            header_map[
                header_map['Keyword'] == 'FILENAME'
            ]['Description'].iloc[0]
        )

    # Getting SIMBAD/GAIA Catalog datas
    catalog_data = get_simbad_data(
        RV2.headers['INSTRUMENT_HEADER'][header_map[
            header_map['Keyword'] == 'OBJECT'
        ]['ESO_keyword'].iloc[0]]
    )
    cid_default = header_map[
        header_map['Keyword'] == 'CID'
    ]['default_value'].iloc[0]

    if (catalog_data['CID'] != cid_default):
        try:
            catalog_data['CCLR'] = get_gaia_data(catalog_data['CID'])
        except Exception:
            print('Gaia request failed.')
    else:
        print("Gaia request can't be done because SIMBAD request failed")

    # Conversion coord format CRA and CDEC
    cra_raw = RV2.headers['INSTRUMENT_HEADER'][
        header_map[header_map['Keyword'] == 'CRA']['ESO_keyword'].iloc[0]
    ]
    cdec_raw = RV2.headers['INSTRUMENT_HEADER'][
        header_map[header_map['Keyword'] == 'CDEC']['ESO_keyword'].iloc[0]
    ]
    catalog_data['CRA'] = convert_to_sexagesimal(cra_raw)
    catalog_data['CDEC'] = convert_to_sexagesimal(cdec_raw)

    add_keyword_cat = ['CEQNX', 'CEPCH', 'CPMR', 'CPMD', 'CRV']
    for keyword in add_keyword_cat:
        catalog_data[keyword] = (
            RV2.headers['INSTRUMENT_HEADER']
            [header_map[
                header_map['Keyword'] == keyword
            ]['ESO_keyword'].iloc[0]]
        )

    rv = catalog_data['CRV']
    rv_z = round(rv/(c/1e3).value, 8)
    catalog_data['CZ'] = rv_z

    # Keywords that depend on the TRACE number.
    keyword_list = [
        'CSRC', 'CID', 'CRA', 'CDEC', 'CEQNX', 'CEPCH',
        'CPLX', 'CPMR', 'CPMD', 'CRV', 'CZ', 'CCLR'
    ]

    with fits.open(names["raw_file"]) as hdu_raw:
        dpr_type = (
            hdu_raw['PRIMARY']
            .header['HIERARCH ESO DPR TYPE'].split(",")
        )
        for i in range(1, nb_trace+1):
            if (dpr_type[math.ceil(i/nb_slice)-1] == 'OBJECT'):
                l2_hdu.header['TRACE'+str(i)] = (
                    'SCI',
                    header_map[
                        header_map['Keyword'] == 'TRACE'
                    ]['Description'].iloc[0]
                )
            elif (dpr_type[math.ceil(i/nb_slice)-1] == 'FP'):
                l2_hdu.header['TRACE'+str(i)] = (
                    'CAL',
                    header_map[
                        header_map['Keyword'] == 'TRACE'
                    ]['Description'].iloc[0]
                )
            elif (dpr_type[math.ceil(i/nb_slice)-1] == 'SKY'):
                l2_hdu.header['TRACE'+str(i)] = (
                    dpr_type[math.ceil(i/nb_slice)-1],
                    header_map[
                        header_map['Keyword'] == 'TRACE'
                    ]['Description'].iloc[0]
                )
            else:
                l2_hdu.header['TRACE'+str(i)] = (
                    header_map[
                        header_map['Keyword'] == 'TRACE'
                    ]['default_value'].iloc[0],
                    header_map[
                        header_map['Keyword'] == 'TRACE'
                    ]['Description'].iloc[0]
                )

            # CALIBRATION SOURCE KEYWORD
            if (l2_hdu.header['TRACE'+str(i)] == 'CAL'):
                clsrc_value = RV2.headers['INSTRUMENT_HEADER'][
                    header_map[
                        header_map['Keyword'] == 'CLSRC'
                    ]['ESO_keyword'].iloc[0]
                ]
                if clsrc_value == 'HEADER':
                    l2_hdu.header['CLSRC'+str(i)] = (
                        RV2.headers['INSTRUMENT_HEADER'][
                            'HIERARCH ESO PRO REC1 RAW2 CATG'
                        ].split('_')[math.ceil(i/nb_slice)-1],
                        header_map[
                            header_map['Keyword'] == 'CLSRC'
                        ]['Description'].iloc[0]
                    )
                else:
                    l2_hdu.header['CLSRC'+str(i)] = (
                        RV2.headers['INSTRUMENT_HEADER'][
                            header_map[
                                header_map['Keyword'] == 'CLSRC'
                            ]['ESO_keyword'].iloc[0]
                        ].split('_')[math.ceil(i/nb_slice)-1],
                        header_map[
                            header_map['Keyword'] == 'CLSRC'
                        ]['Description'].iloc[0]
                    )
            else:
                l2_hdu.header['CLSRC'+str(i)] = (
                    header_map[
                        header_map['Keyword'] == 'CLSRC'
                    ]['default_value'].iloc[0],
                    header_map[
                        header_map['Keyword'] == 'CLSRC'
                    ]['Description'].iloc[0]
                )

            # CATALOG KEYWORDS
            if (l2_hdu.header['TRACE'+str(i)] == 'SCI'):
                for keyword in keyword_list:
                    l2_hdu.header[keyword+str(i)] = (
                        catalog_data[keyword],
                        header_map[
                            header_map['Keyword'] == keyword
                        ]['Description'].iloc[0]
                    )
            else:
                for keyword in keyword_list:
                    l2_hdu.header[keyword+str(i)] = (
                        header_map[
                            header_map['Keyword'] == keyword
                        ]['default_value'].iloc[0],
                        header_map[
                            header_map['Keyword'] == keyword
                        ]['Description'].iloc[0]
                    )

    # BINNING KEYWORD
    binx = str(RV2.headers['INSTRUMENT_HEADER']['HIERARCH ESO DET WIN1 BINX'])
    biny = str(RV2.headers['INSTRUMENT_HEADER']['HIERARCH ESO DET WIN1 BINY'])
    l2_hdu.header['BINNING'] = (
        f"{binx}x{biny}",
        header_map[header_map['Keyword'] == 'BINNING']['Description'].iloc[0]
    )

    # NUMTRACE KEYWORD
    l2_hdu.header['NUMTRACE'] = (
        nb_trace,
        header_map[header_map['Keyword'] == 'NUMTRACE']['Description'].iloc[0]
    )

    # DATE KEYWORD
    current_time = Time.now()
    l2_hdu.header['DATE'] = (
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
        deg_to_sexagesimal(RV2.headers['INSTRUMENT_HEADER']['RA'], True),
        header_map[header_map['Keyword'] == 'TRA']['Description'].iloc[0]
    )

    # TDEC KEYWORD
    l2_hdu.header['TDEC'] = (
        deg_to_sexagesimal(RV2.headers['INSTRUMENT_HEADER']['DEC'], False),
        header_map[header_map['Keyword'] == 'TDEC']['Description'].iloc[0]
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

    # MOONANG/MOONEL/MOONILLU/MOONRV/SUNEL KEYWORDS
    moon_sun_params = get_moon_sun_info(
        RV2.headers['INSTRUMENT_HEADER']['RA'],
        RV2.headers['INSTRUMENT_HEADER']['DEC'],
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
            RV2.headers['INSTRUMENT_HEADER'][
                f"HIERARCH ESO QC ORDER{str(i)} SNR"
            ],
            header_map[header_map['Keyword'] == 'EXSNR']['Description'].iloc[0]
        )

    # EXSNRW-N KEYWORD
    for i in range(int(l2_hdu.header['NUMORDER'])):
        l2_hdu.header[f'EXSNRW{str(i+1)}'] = (
            round(
                RV2.data["TRACE1_WAVE"][i, 0]
                + (
                    RV2.data["TRACE1_WAVE"][i, -1]
                    - RV2.data["TRACE1_WAVE"][i, 0]
                )/2
            ),
            header_map[
                header_map['Keyword'] == 'EXSNRW'
            ]['Description'].iloc[0]
        )

    # DRPFLAG KEYWORD
    drp_flag = RV2.headers['INSTRUMENT_HEADER'][
        header_map[header_map['Keyword'] == 'DRPFLAG']['ESO_keyword'].iloc[0]
    ]
    if drp_flag == 1:
        drpflag = 'Pass'
    else:
        drpflag = 'Fail'

    l2_hdu.header['DRPFLAG'] = (
        drpflag,
        header_map[
            header_map['Keyword'] == 'DRPFLAG'
        ]['Description'].iloc[0]
    )

    # COLOFLAG KEYWORD
    try:
        color_flag = RV2.headers['INSTRUMENT_HEADER'][
            header_map[
                header_map['Keyword'] == 'COLOFLAG'
            ]['ESO_keyword'].iloc[0]
        ]
        if (color_flag == 1):
            coloflag = 'Pass'
        else:
            coloflag = 'Fail'
    except Exception:
        coloflag = 'Fail'

    l2_hdu.header['COLOFLAG'] = (
        coloflag,
        header_map[
            header_map['Keyword'] == 'COLOFLAG'
            ]['Description'].iloc[0]
    )

    # SUMMFLAG KEYWORD
    flags = ["COLOFLAG", "TELFLAG", "INSTFLAG", "DRPFLAG", "OBSFLAG"]

    # Retrieve all flag values.
    flag_values = [l2_hdu.header.get(flag, "Pass") for flag in flags]

    # Priority of states: Fail > Warn > Pass.
    if "Fail" in flag_values:
        if "Fail" == flag_values[0] and "Fail" not in flag_values[1:]:
            summflag = 'Warn'
        else:
            summflag = 'Fail'
    elif "Warn" in flag_values:
        summflag = 'Warn'
    else:
        summflag = 'Pass'

    l2_hdu.header['SUMMFLAG'] = (
        summflag,
        header_map[
            header_map['Keyword'] == 'SUMMFLAG'
        ]['Description'].iloc[0]
    )
    l2_hdu.header['GEOSYS'] = ('WGS84', 'Coordinate system for observatory location')
    if ('PRIMARY' not in RV2.extensions):
        RV2.create_extension(
            ext_name='PRIMARY', ext_type='PrimaryHDU', header=l2_hdu.header
        )
    else:
        RV2.set_header(ext_name='PRIMARY', header=l2_hdu.header)
    return


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
        for name in result['ids'][0].split('|'):
            if (name.lower().startswith('gaia dr3')):
                gaia_dr3_source = name[:8]
                gaia_dr3_name = name[5:]
            elif (name.lower().startswith('gaia dr2')):
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
        if not np.ma.is_masked(result['plx_value'][0]):
            data['CPLX'] = result['plx_value'][0]
        else:
            data['CPLX'] = 'Null'

        return data

    except Exception as e:
        print(f"Catalog not found for the name of {obj}, err:{e}")

        # Return default values if the object is not found
        cat_list = ['CSRC', 'CID', 'CPLX', 'CCLR']

        base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        header_map_path = os.path.join(base_dir, "config", "header_map.csv")
        header_map = pd.read_csv(header_map_path)

        for key in cat_list:
            data[key] = header_map[
                header_map['Keyword'] == key
            ]['default_value'].iloc[0]
        return data


def get_gaia_data(gaia_name: str) -> float:
    """
    Retrieves Gaia photometric data and computes the color index (BP - RP).

    Args:
        gaia_name (str): The Gaia source identifier in the format
            "DR3 <source_id>".

    Returns:
        color_br (float): The computed color index (BP - RP) for the
            specified Gaia source.
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


def convert_to_sexagesimal(value: float) -> str:
    """
    Converts a numerical value in HHMMSS.SSS or DDMMSS.SSS format
    into a properly formatted sexagesimal string (HH:MM:SS.SSS or DD:MM:SS.SSS).

    Args:
        value (float): The numerical value to convert, where
            HHMMSS.SSS represents hours, minutes, and seconds
            or DDMMSS.SSS represents degrees, minutes, and seconds.

    Returns:
        str: The formatted sexagesimal string in HH:MM:SS.SSS or DD:MM:SS.SSS format.
    """

    # Preserve the negative sign if present
    sign = "-" if value < 0 else ""
    value = abs(value)

    hours_or_degrees = int(value // 10000)  # Extract HH or DD
    minutes = int((value % 10000) // 100)   # Extract MM
    seconds = (value % 100)                 # Extract SS.SSS

    # Return formatted string with leading zeros and three decimal places for seconds
    return f"{sign}{hours_or_degrees:02}:{minutes:02}:{seconds:06.3f}"


def convert_lst(lst: float) -> str:
    """
    Converts the local sidereal time from seconds to hh:mm:ss.ms

    Args:
        lst (foat): the local sidereal time in seconds

    Returns:
        lst_sexa (str): the local sidereal time in the format hh:mm:ss.ms
    """

    if (lst < 0):
        res = convert_lst(-lst)
        return '-' + res

    hours = int(lst//3600)
    minutes = int((lst % 3600)//60)
    seconds = int(lst % 60)
    ms = round(lst % 1 * 1000)

    lst_sexa = f"{hours:02}:{minutes:02}:{seconds:02}.{ms:03}"
    return lst_sexa


def deg_to_sexagesimal(value_deg: float, is_ra: bool = False) -> str:
    """
    Converts a value in degrees to sexagesimal format
    (HH:MM:SS.sss or ±DD:MM:SS.sss).

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
    Computes the hour angle from the local sidereal time and the right
    ascension of the object.

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
        date_obs_str (str): The observation date in ISO 8601 format
            ("YYYY-MM-DDTHH:MM:SS.sss").

    Returns:
        int: The version of the instrument corresponding to the observation
            date.

    Raises:
        ValueError: If the observation date does not fall within any defined
            version range.
    """
    # Convert the observation date string into a datetime object
    date_obs = datetime.fromisoformat(date_obs_str)

    # Iterate through the version ranges defined in the configuration
    for version_info in config.INSTRUMENT_VERSIONS:
        # Convert the start date to a datetime object
        start_date = datetime.fromisoformat(version_info["start_date"])

        # Convert the end date to a datetime object, or use datetime.max for
        # open-ended ranges
        end_date = (
            datetime.fromisoformat(version_info["end_date"])
            if version_info["end_date"] is not None
            else datetime.max
        )

        # Check if the observation date falls within the current version's
        # date range
        if start_date <= date_obs <= end_date:
            return version_info["version"]

    # If no version matches, raise an exception
    raise ValueError(
        f"No instrument version corresponds to the date {date_obs_str}"
    )


def get_moon_sun_info(
        target_ra: float, target_dec: float,
        obs_lat: float, obs_lon: float,
        obs_alt: float, obs_time: str,
        jd_utc: float
) -> list:
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
        obs_time (str): Observation time in ISO format
            (e.g., 'YYYY-MM-DD HH:MM:SS').
        jd_utc (float): Julian Date (UTC) for the observation.

    Returns:
        list: A list containing the following values:
            - sun_el (float): The elevation of the Sun at the observation site
                (in degrees).
            - moon_ang (float): The angular separation between the Moon and
                the target (in degrees).
            - moon_el (float): The elevation of the Moon at the observation
                site (in degrees).
            - moon_illu (float): The illumination of the Moon in percentage.
            - moon_rv (float): The radial velocity of reflected sunlight off
                moon (in km/s).
    """

    # Observer's location and observation time
    location = EarthLocation(lat=obs_lat, lon=obs_lon, height=obs_alt)
    time = Time(obs_time)

    # Get the Moon's position at the given time and location (GCRS)
    moon_coord = get_body("moon", time, location=location)

    # We need to convert moon_coord in ICRS to match target_coord ref
    moon_coord_icrs = moon_coord.transform_to("icrs")

    # Target's position
    target_coord = SkyCoord(
        ra=target_ra, dec=target_dec, frame='icrs', unit="deg", obstime=time
    )

    # Calculate the angular separation between the target and the Moon
    moon_ang = round(moon_coord_icrs.separation(target_coord).deg, 4)

    # Calculate the Moon's elevation above the horizon (GCRS here but could be
    # ICRS, no impact)
    moon_altaz = moon_coord.transform_to(
        AltAz(obstime=time, location=location)
    )
    moon_el = round(moon_altaz.alt.deg, 4)

    # Get the Sun's position and transform it to the observer's AltAz frame
    sun = get_body("sun", time, location=location)

    sun_altaz = sun.transform_to(AltAz(obstime=time, location=location))
    sun_el = round(sun_altaz.alt.deg, 4)

    # Calculate the Moon's illumination
    elongation = moon_coord.separation(sun)
    moon_phase_angle = np.arctan2(
        sun.distance*np.sin(elongation),
        moon_coord.distance - sun.distance*np.cos(elongation)
    )
    moon_illu = round((1 + np.cos(moon_phase_angle).value) / 2 * 100, 4)

    # Calculate the RV of reflected sunlight off moon
    moon_rv = round(get_moon_velocity_in_target_direction(
        target_ra, target_dec, jd_utc), 4
    )

    return [sun_el, moon_ang, moon_el, moon_illu, moon_rv]


def get_moon_velocity_in_target_direction(
        alpha_deg: float, delta_deg: float, julian_day: float
) -> float:
    """
    Compute the velocity of the Moon projected in the direction of a given
    target in the sky.

    Parameters:
        alpha_deg (float): Right Ascension (RA) of the target in degrees.
        delta_deg (float): Declination (Dec) of the target in degrees.
        julian_day (float): Julian date for the calculation.

    Returns:
        projected_velocity_km_s (float): Radial velocity of the Moon in km/s
            in the target direction.
    """

    # Convert Julian Day to Astropy Time object
    t = Time(julian_day, format='jd')

    # Get the Moon's barycentric position & velocity (in AU and AU/day)
    # using JPL ephemeris
    moon_pos, moon_vel = get_body_barycentric_posvel(
        'moon', t, ephemeris='jpl'
    )

    # Extract velocity components (AU/day)
    x_vel, y_vel, z_vel = moon_vel.xyz.to_value(u.km/u.s)

    # Convert target coordinates (RA, Dec) to radians
    alpha_rad = np.deg2rad(alpha_deg)
    delta_rad = np.deg2rad(delta_deg)

    # Compute projected radial velocity
    projected_velocity_km_s = (
        x_vel * np.cos(alpha_rad) * np.cos(delta_rad) +
        y_vel * np.sin(alpha_rad) * np.cos(delta_rad) +
        z_vel * np.sin(delta_rad)
    )

    return projected_velocity_km_s
