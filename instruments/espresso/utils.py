from astroquery.simbad import Simbad
from astropy.io import fits
from astropy.time import Time
import numpy as np
import pandas as pd
import os
import instruments.espresso.config.config as config

def do_conversion(RV2):
    """Complete all the necessary steps to convert the raw file to the level 2 format.

    Args:
        RV2 (RV2 instance): the RV2 object created with the from_fits() method

    """
    RV2.blaze_file = []
    convert_S2D(RV2)
    convert_raw_file(RV2)
    create_primary_header(RV2)
    convert_blaze(RV2)
    return


def create_primary_header(RV2):
    """Create.

        Args:
            fn (str): file path (relative to the repository)
            instrument (str): instrument name. None implies FITS file is in EPRV standard format.
            overwrite (bool): if this instance is not empty, specifies whether to overwrite

        Raises:
            IOError: when a invalid file is presented

        Note:
            This is not a @classmethod so initialization is
            required before calling this function

        """
    #We create an empty HDU to store the L2 header
    primary_HDU = fits.PrimaryHDU(data = None)
    primary_HDU.header['EXTNAME'] = 'PRIMARY'
    #We load the header map to convert between raw file headers and L2 header
    header_map = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/config/header_map.csv')
    #We replace the %UT% in the header map with the front end ID to get the correct keywords depending on which UT was used
    header_map['ESO_keyword'] = header_map['ESO_keyword'].str.replace('%UT%', RV2.UT)
    #Special case for UT1
    if RV2.UT == '1':
        header_map['ESO_keyword'] = header_map['ESO_keyword'].str.replace('INS1', 'INS')
    #We iterate over all the keys of the header map
    for index, values in header_map.iterrows():
        if(header_map.skip.iloc[index] == 1):
            continue
        #Add the HIERARCH keyword to the header if the keyword is longer than 8 characters
        if(len(values[0]) > 8):
            values[0] = 'HIERARCH ' + values[0]
        values[0] = values[0].strip()

        try:
            #If there is a fixed value to set, we set it
            if(pd.notna(header_map.value.iloc[index])):
                primary_HDU.header[values[0]] = header_map.value.iloc[index]
            #Otherwise, we copy the value from the raw file
            elif(pd.notna(header_map['ESO_keyword'].iloc[index])):
                primary_HDU.header[values[0]] = (RV2.headers['INSTRUMENT_HEADER'][header_map['ESO_keyword'].iloc[index]], RV2.headers['INSTRUMENT_HEADER'].comments[header_map['ESO_keyword'].iloc[index]].strip())
            #If the value is not present in the raw file, we set it to None
            else:
                primary_HDU.header[values[0]] = None
        except Exception as e:
            print(e, values[0],RV2.headers['INSTRUMENT_HEADER'][header_map['ESO_keyword'].iloc[index]]) 
    

    primary_HDU.header['FILENAME'] = ('r.'+RV2.filename[:-5] + '_L2.fits', 'Name of the file')

    #Getting Gaia identifier
    try:
        gaia_name = get_gaia_name(RV2.headers['INSTRUMENT_HEADER']['OBJECT'])
        primary_HDU.header['CATSRC'] = (gaia_name[:9],'Catalog name')
        primary_HDU.header['CATID'] = (gaia_name[9:], 'Catalog ID')

    except Exception as e:
        print('Could not find Gaia name: ' + str(e))

    current_time = Time.now()
    primary_HDU.header['DATE']  = (current_time.iso, 'Date of file creation')
    primary_HDU.header['TCSZA'] = (np.round(90 - primary_HDU.header['TCSEL'], 3), '[deg] Zenith angle')

    #Converting LST in seconds to hh:mm:ss.ms
    lst = convert_lst(primary_HDU.header['TCSLST'])
    primary_HDU.header['TCSLST'] = (lst, 'Local Sidereal Time')
    ra, dec = convert_to_sexagesimal(primary_HDU.header['CATRA'], primary_HDU.header['CATDEC'])
    primary_HDU.header['CATRA'] = (ra, 'Catalog Right Ascension')
    primary_HDU.header['CATDEC'] = (dec, 'Catalog Declination')
    tcs_ra, tcs_dec = convert_to_sexagesimal(primary_HDU.header['TCSRA'], primary_HDU.header['TCSDEC'])
    primary_HDU.header['TCSRA'] = (tcs_ra, 'Telescope Right Ascension')
    primary_HDU.header['TCSDEC'] = (tcs_dec, 'Telescope Declination')
    primary_HDU.header['TCSHA'] = (compute_hour_angle(lst, tcs_ra), 'Hour angle')
    primary_HDU.header['CATPMRA'] = (primary_HDU.header['CATPMRA']*1000, '[marcsec/yr] Proper Motion Alpha')
    primary_HDU.header['CATPMDEC'] = (primary_HDU.header['CATPMDEC']*1000, '[marcsec/yr] Proper Motion Delta')    
    if('PRIMARY' not in RV2.extensions):
        RV2.create_extension(ext_name = 'PRIMARY', ext_type = 'PrimaryHDU', header = primary_HDU.header)
    else:
        RV2.set_header(ext_name = 'PRIMARY', header = primary_HDU.header)
    return

def convert_raw_file(RV2):
    """Transfer the raw file information to the level 2 hdul file by keeping the same primary header
    and copying the exposure meter, pupil image and guiding frame information.

    """
    
    with fits.open(RV2.dirname +'/'+ RV2.filename) as hdul:        
        ins_header = fits.ImageHDU(header=hdul['PRIMARY'].header.copy(), data = None)
        ins_header.header['EXTNAME'] = 'INSTRUMENT_HEADER'

        if(ins_header.header['EXTNAME'] not in RV2.extensions):
            RV2.create_extension(ext_name = ins_header.header['EXTNAME'], ext_type = 'ImageHDU', header = ins_header.header)
        else:
            RV2.set_header(ext_name = ins_header.header['EXTNAME'], header = ins_header.header)
        try:
            if('HIERARCH ESO OCS ENABLED FE' in hdul['PRIMARY'].header):
                front_end_id  = hdul['PRIMARY'].header['HIERARCH ESO OCS ENABLED FE']
            elif('TELESCOP' in hdul['PRIMARY'].header):
                front_end_id = hdul['PRIMARY'].header['TELESCOP'][-1]
            RV2.UT = front_end_id
            pupil_key = 'PS' + str(front_end_id)
            guiding_key = 'FS' + str(front_end_id) + 'INT'
        except Exception as e: 
            print('Could not find front end ID ' + str(e))
            
        try:
            exp_meter_hdu = hdul['Exp Meter bin table'].copy()
            exp_meter_hdu.header['EXTNAME'] = 'EXPOMETER'
        except Exception as e:
            print('Error when converting exposure meter: ' + str(e))
            exp_meter_hdu = fits.ImageHDU(data = None)
            exp_meter_hdu.header['EXTNAME'] = 'EXPOMETER'
        try:
            pupil_img_hdu = hdul[pupil_key].copy()
            pupil_img_hdu.header['EXTNAME'] = 'PUPILIMAGE'
        except Exception as e:
            print('Error when converting pupil image: ' + str(e))
            pupil_img_hdu = fits.ImageHDU(data = None)
            pupil_img_hdu.header['EXTNAME'] = 'PUPILIMAGE'
        try:
            guiding_frame_hdu = hdul[guiding_key].copy()
            guiding_frame_hdu.header['EXTNAME'] = 'GUIDINGIMAGE'
        except Exception as e: 
            print('Error when converting guiding image: ' + str(e))
            guiding_frame_hdu = fits.ImageHDU(data = None)
            guiding_frame_hdu.header['EXTNAME'] = 'GUIDINGIMAGE'
        if('Exp Meter bin table' not in RV2.extensions):
            RV2.create_extension(ext_name = 'Exp Meter bin table', ext_type = 'BinTableHDU', header = exp_meter_hdu.header, data = exp_meter_hdu.data)
        else:
            RV2.set_data(ext_name = 'Exp Meter bin table', data = exp_meter_hdu.data)
            RV2.set_header(ext_name = 'Exp Meter bin table', header = exp_meter_hdu.header)
        if('PUPIL_IMAGE' not in RV2.extensions):
            RV2.create_extension(ext_name = 'PUPIL_IMAGE', ext_type = 'ImageHDU', header = pupil_img_hdu.header, data = pupil_img_hdu.data)
        else:
            RV2.set_data(ext_name = 'Pupil image', data = pupil_img_hdu.data)
            RV2.set_header(ext_name = 'Pupil image', header = pupil_img_hdu.header)
        if('GUIDING FRAME' not in RV2.extensions):
            RV2.create_extension(ext_name = 'GUIDING_FRAME', ext_type = 'ImageHDU', header = guiding_frame_hdu.header, data = guiding_frame_hdu.data)
        else:
            RV2.set_data(ext_name = 'GUIDING_FRAME', data = guiding_frame_hdu.data)
            RV2.set_header(ext_name = 'GUIDING_FRAME', header = guiding_frame_hdu)
    return 

def convert_S2D(RV2):
    """Converts the e2ds files to the level 2 format by extracting the values of the e2ds files and adding them to the level 2 file.
    """
    #Loop over the two fibers (science and FP)
    trace_number = 1
    for i, fiber in enumerate(config.fibers.keys()):
        #Load the s2d file
        fname = RV2.dirname + '/r.' +RV2.filename[:-5] + '_S2D_BLAZE_' + fiber + '.fits'

        with fits.open(fname) as hdul:
            #
            prim_header = hdul['PRIMARY'].header.copy()            
            if(fiber == 'SCI'):
                blaze_file = prim_header['HIERARCH ESO PRO REC1 CAL28 NAME']
            else:
                blaze_file = prim_header['HIERARCH ESO PRO REC1 CAL29 NAME']
            RV2.blaze_file.append(blaze_file)
            #We iterate over the cameras, slices and fields to extract the values from the e2ds file
            for slice in config.slices:
                for field in config.extnames.keys():
                    #We extract the values from one invividual camera
                    values = hdul[field].data
                    if(field == 'WAVEDATA_VAC_BARY'):
                        #We doppler shift the wavelength values
                        single_slice_bary = values[slice::2,:].copy()
                        values = doppler_shift(values, prim_header['HIERARCH ESO QC BERV'])
                        bary_data = np.ones((1,1))*prim_header['HIERARCH ESO QC BERV']
                        single_slice_berv = bary_data#[slice::2,:]
                        berv_hdu = fits.ImageHDU(data = single_slice_berv.copy())
                        berv_hdu.header['EXTNAME'] = 'TRACE'+str(trace_number)  + '_BERV'

                    #We extract the values from one individual slice (every other row in the camera values)
                    single_slice_values = values[slice::2,:]
                    #We create a new HDU with the values of the slice and the original header
                    hdu_l2 = fits.ImageHDU(data = single_slice_values.copy(), header = hdul[field].header)

                    hdu_l2.header['EXTNAME'] = 'TRACE'+str(trace_number)  + config.extnames[field]
                    if(hdu_l2.header['EXTNAME'] not in RV2.extensions):
                        print('Creating extension: ' + hdu_l2.header['EXTNAME'] + ' with data of shape ', hdu_l2.data.shape)
                        RV2.create_extension(ext_name = hdu_l2.header['EXTNAME'], ext_type = 'ImageHDU', header = hdu_l2.header, data = hdu_l2.data)
                    else:
                        print('Setting data for extension: ' + hdu_l2.header['EXTNAME'] + ' with data of shape ', hdu_l2.data.shape)
                        RV2.set_data(ext_name = hdu_l2.header['EXTNAME'], data =hdu_l2.data)
                        RV2.set_header(ext_name = hdu_l2.header['EXTNAME'], header = hdu_l2.header)
                trace_number += 1
            
                    
    return
def convert_blaze(RV2):
    """Copies the blaze calibration files to the level 2 file by extracting the values from the blaze files and adding them to the level 2 file.
    """
    #Loop over the two fibers (science and FP)
    trace_number = 1
    for it, fiber in enumerate(config.fibers.keys()):
        #If file format is L2, we split by cameras and slices
    
        try:
            with fits.open(RV2.dirname +'/'+ RV2.blaze_file[it]) as hdul:
                for slice in config.slices:
                    blaze = hdul[1].data#[config.cam_range[camera][0]:config.cam_range[camera][1], :]
                    single_slice_blaze = blaze[slice::2,:]
                    blaze_hdu_l2 = fits.ImageHDU(data = single_slice_blaze.copy(), header = hdul[1].header)

                    blaze_hdu_l2.header['EXTNAME'] =  'TRACE'+str(trace_number)  +'_BLAZE'
                    if(blaze_hdu_l2.header['EXTNAME'] not in RV2.extensions):
                        print('Creating extension: ' + blaze_hdu_l2.header['EXTNAME'] + ' with data of shape ', blaze_hdu_l2.data.shape)
                        RV2.create_extension(ext_name = blaze_hdu_l2.header['EXTNAME'], ext_type = 'ImageHDU', header = blaze_hdu_l2.header, data = blaze_hdu_l2.data)
                    else:
                        print('Setting data for extension: ' + blaze_hdu_l2.header['EXTNAME'] + ' with data of shape ', blaze_hdu_l2.data.shape)
                        RV2.set_data(ext_name = blaze_hdu_l2.header['EXTNAME'], data =blaze_hdu_l2.data)
                        RV2.set_header(ext_name = blaze_hdu_l2.header['EXTNAME'], header = blaze_hdu_l2.header)
                    trace_number += 1
        except Exception as e:
            print('Error when converting blaze file: ' + str(e))
    #If file format is original, we extract the values from the blaze file without separating by cameras and slices
        
    return

def get_gaia_name(obj):
        """Query the Simbad database to get the Gaia name of the object.

        Args:
            obj (str): the name of the object

        Returns:
            str: the Gaia identifier of the object in the format "Gaia DR2 1234567890"
        """
        try:
            #We query the object in the Simbad database
            result_table = Simbad.query_objectids(obj)
            #We iterate over the results to find the Gaia identifier
            for name in result_table:
                if(name[0].lower().startswith('gaia dr3')):
                    
                    return name[0]
        except:
            return 'None'
        return

def convert_to_sexagesimal(ra, dec):
        """This function makes the conversion between the ESPRESSO drs format for the right ascension and declination to the sexagesimal format.

        Args:
            ra (float): The right ascension obtained from the TARG ALPHA keyword of the header
            dec (float): The declination obtained from the TARG DELTA keyword of the header

        Returns:
            tuple of string: The right ascension and declination in the sexagesimal format
        """
        ra = f'{ra:10.3f}'
        ra = ra[:2] + 'h' + ra[2:4] + 'm' + ra[4:] + 's'
        dec = f'{dec:011.3f}' if dec < 0 else f'{dec:10.3f}'
        if dec.startswith('-'):

            dec = dec[:3] + 'd' + dec[3:5] + 'm' + dec[5:] + 's'

        else:

            dec = dec[:2] + 'd' + dec[2:4] + 'm' + dec[4:] + 's'

        return ra, dec

def compute_hour_angle(lst, ra):
        """Computes the hour angle from the local sidereal time and the right ascension of the object.

        Args:
            lst (str): local sideral time in the format hh:mm:ss.ms
            ra (ra): right ascension of the object in the format hh:mm:ss.ms

        Returns:
            str: the hour angle in the format hh:mm:ss.ms
        """
        h_ra = int(ra[:ra.find('h')])
        m_ra = int(ra[ra.find('h')+1:ra.find('m')])
        s_ra = float(ra[ra.find('m')+1:ra.find('s')])

        h_lst = int(lst[:lst.find('h')])
        m_lst = int(lst[lst.find('h')+1:lst.find('m')])
        s_lst = float(lst[lst.find('m')+1:lst.find('s')])
        ha_s = 3600*(h_lst - h_ra) + 60*(m_lst - m_ra) + (s_lst - s_ra)
        return convert_lst(ha_s)

def convert_lst(lst):
        """Converts the local sidereal time from seconds to hh:mm:ss.ms

        Args:
            lst (foat): the local sidereal time in seconds

        Returns:
            str: the local sidereal time in the format hh:mm:ss.ms
        """
        if(lst < 0):
            res = convert_lst(-lst)
            return '-' + res
        hours = int(lst//3600)
        minutes = int((lst %3600)//60)
        seconds = int(lst %60)
        ms = round(lst % 1 * 1000)
        
        return f"{hours:02}h{minutes:02}m{seconds:02}.{ms:03}s"

def doppler_shift(wave, rv):
        """Performs the doppler shift on the wavelength values.

        Args:
            wave (numpy array): the original wavelength values
            rv (radial velocity): the radial velocity of the object in km/s

        Returns:
            numpy array: the doppler shifted wavelength values
        """
        c = 299792.458 #speed of light in km/s
        return wave + wave *rv/c   