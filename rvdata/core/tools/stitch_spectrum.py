import numpy as np
from astropy import units as u
from astropy import constants
from specutils import SpectrumCollection
from specutils.manipulation import FluxConservingResampler
from scipy.interpolate import interp1d


def calculate_spectral_envelope(
    wave, flux, wavepoints=None, boxsize=100, percentile=99, fsr=None
):
    """Calculate a envelope function across spectral orders.

    Parameters
    ----------
    wave : wavelength array (norders x npixels)

    flux : flux array (norders x npixels)


    Keywords
    --------
    wavepoints : array of wavelengths to use for calculating the envelope
                 flux.  If None, then the wavelengths are derived based on the
                 peak of the spectrum

    boxsize : int - size of box centered on each value in wavepoints for which
              to calculate fluxout

    percentile : int - threshold in percent to use for deriving the peak
                flux.  Not used if wave is specified

    fsr : Free spectral range mask array, corresponding in size to the format
          of nd flux and wavelength extensions.  If specified, restrict the
          peak search to within these bounds

    Returns
    -------
    waveout, fluxout : arrays of length equal to the number of orders in nd.


    2021-04-30, CFB
    UPDATED 2025-12-04, LPA
    """

    if wave.shape == flux.shape:
        norders, npixels = wave.shape
    else:
        print("wave and flux array not equal shape")
        return (None, None)

    # Parse the FSR mask
    pixstart = np.zeros(norders, dtype=float)
    pixend = np.full(norders, npixels - 1, dtype=float)
    if fsr is not None:
        if fsr.shape == wave.shape:
            fsr_mask = ~fsr
            has_valid_pixels = np.any(fsr_mask, axis=1)
            first_valid = np.argmax(fsr_mask, axis=1)
            last_valid = np.argmax(fsr_mask[:, ::-1], axis=1)
            pixstart = np.where(has_valid_pixels, first_valid, np.nan)
            pixend = np.where(has_valid_pixels, npixels - 1 - last_valid, np.nan)
        else:
            print("FSR and wavelength array are not equal shape. Not using FSR.")

    if wavepoints is not None:
        # Wavelengths are specified
        # Calculate fluxes at these points
        wavepoints = np.array(wavepoints)  # Ensure this is a numpy array
        # Remove any nans and zeros
        wavepoints = wavepoints[~np.isnan(wavepoints)]
        wavepoints = wavepoints[wavepoints != 0]
        if wavepoints.size == 0:
            return (wavepoints, np.array([], dtype=float))
        wave_midpoints = wave[:, npixels // 2]
        # Identify the order that best matches each wavepoint using broadcasted distances
        worder = np.nanargmin(
            np.abs(wave_midpoints[:, None] - wavepoints[None, :]), axis=0
        )
        selected_wave = wave[worder, :]
        # Locate the nearest pixel index within each selected order
        windex = np.nanargmin(np.abs(selected_wave - wavepoints[:, None]), axis=1)
        hbox = boxsize // 2
        start = np.clip(windex - hbox, 0, npixels - 1)
        end = np.clip(windex + hbox, 0, npixels - 1)
        pixel_indices = np.arange(npixels)
        window_mask = (pixel_indices >= start[:, None]) & (pixel_indices < end[:, None])
        window_flux = np.where(window_mask, flux[worder], np.nan)
        fluxout = np.nanpercentile(window_flux, percentile, axis=1)
        waveout = wavepoints
    else:
        # Calculate the envelope
        valid_orders = ~np.isnan(pixstart)
        start = np.where(valid_orders, pixstart, 0).astype(int)
        end = np.where(valid_orders, pixend, 0).astype(int)
        pixel_indices = np.arange(npixels)
        order_window_mask = (pixel_indices >= start[:, None]) & (
            pixel_indices < end[:, None]
        )
        masked_flux = np.where(order_window_mask, flux, np.nan)
        masked_wave = np.where(order_window_mask, wave, np.nan)
        fluxout = np.full(norders, np.nan)
        if np.any(valid_orders):
            fluxout[valid_orders] = np.nanpercentile(
                masked_flux[valid_orders], percentile, axis=1
            )
        waveout = np.nanmean(
            np.where(masked_flux > fluxout[:, None], masked_wave, np.nan), axis=1
        )

    return (waveout, fluxout)


def get_wavelength_grid_with_constant_velocity(wavegrid_start, wavegrid_end, velpix):
    sl = constants.c.value  # 299792458. m/s
    avg_dlwavegrid = velpix / sl
    lwavegrid = np.arange(
        np.log(wavegrid_start), np.log(wavegrid_end) + avg_dlwavegrid, avg_dlwavegrid
    )
    wavegrid = np.exp(lwavegrid)
    return wavegrid


def resample_flux_conserving(sci_wav, sci_dflx, spec_mask, inst_stitch_config):
    """flux-conserving rebinning and stitching of spectral orders

    Parameters
    ----------
    sci_wav : 2D numpy array
        Wavelengths of the science spectrum, shape (norders, npixels).
    sci_dflx : 2D numpy array
        Fluxes of the science spectrum, shape (norders, npixels).
    spec_mask : 2D boolean numpy array
        Mask for the science spectrum, shape (norders, npixels).
    inst_stitch_config : dict
        Instrument-specific stitching configuration parameters.

    Returns
    -------
    st_wave : 1D numpy array
        Wavelengths of the stitched spectrum.
    st_flux : 1D numpy array
        Fluxes of the stitched spectrum.

    2025-06-21, LPA
    UPDATED 2025-12-04, LPA
    """

    nbins = inst_stitch_config["nbins"]

    # Create SpectrumCollection from the science data

    flux = u.Quantity(sci_dflx, unit="photon")
    wave = u.Quantity(sci_wav, unit="AA")
    spec_coll = SpectrumCollection(flux=flux, spectral_axis=wave, mask=spec_mask)

    # Define a common output grid

    min_wave = spec_coll.spectral_axis.min().value
    max_wave = spec_coll.spectral_axis.max().value
    new_wave = np.linspace(min_wave, max_wave, nbins) * spec_coll.spectral_axis.unit

    # Rebin each order (flux-conserving)

    resampler = FluxConservingResampler()

    rebinned_orders = []
    for iorder in range(spec_coll.shape[0]):
        rebinned = resampler(spec_coll[iorder], new_wave)
        rebinned_orders.append(rebinned)

    # Combine overlapping orders

    # Stack the flux arrays and take a nan-mean
    flux_stack = np.vstack(
        [specrb_order.flux.value for specrb_order in rebinned_orders]
    )
    combined_flux = np.nanmean(flux_stack, axis=0) * rebinned_orders[0].flux.unit

    st_wave = new_wave.value
    st_flux = combined_flux.value

    return st_wave, st_flux


def stitch_orders(sci_wav, sci_flx, sci_blz, inst_stitch_config=None):
    """Stitch the spectral orders of a science spectrum.

    Parameters
    ----------
    sci_wav : 2D numpy array
        Wavelengths of the science spectrum, shape (norders, npixels).
    sci_flx : 2D numpy array
        Fluxes of the science spectrum, shape (norders, npixels).
    sci_blz : 2D numpy array
        Blaze function of the science spectrum, shape (norders, npixels).
    inst_stitch_config : dict, optional
        Instrument-specific stitching configuration parameters.
        Default is for NEID.

    Returns
    -------
    st_wave : 1D numpy array
        Wavelengths of the stitched spectrum.
    st_flux : 1D numpy array
        Fluxes of the stitched spectrum.

    2025-06-21, LPA
    """
    # instrument configuration parameters
    iordermin = inst_stitch_config["iordermin"]
    iorderflatbreak = inst_stitch_config["iorderflatbreak"]
    iordermax = inst_stitch_config["iordermax"]
    nbins = inst_stitch_config["nbins"]

    # Prepare the science data
    sci_dflx = sci_flx / (
        sci_blz / np.nanmax(sci_blz, axis=1, keepdims=True)
    )  # Normalize by blaze function

    sciflx_mask = np.isnan(sci_flx)
    sciwav_mask = (sci_wav == 0.0) | np.isnan(sci_wav)
    sciblz_mask = np.isnan(sci_blz)
    spec_mask = np.logical_or(sciflx_mask, sciwav_mask, sciblz_mask)

    # Mask data where the spectrum is bad or missing
    # sci_flxm = np.where(spec_mask, float("nan"), sci_flx)
    sci_wavm = np.where(spec_mask, float("nan"), sci_wav)
    sci_blzm = np.where(spec_mask, float("nan"), sci_blz)
    sci_dflxm = np.where(spec_mask, float("nan"), sci_dflx)

    # Calculate the spectral envelope for the blaze extension
    sci_bzewavb, sci_bzbeenvb = calculate_spectral_envelope(
        sci_wavm[iordermin:iorderflatbreak, :], sci_blzm[iordermin:iorderflatbreak, :]
    )
    sci_bzewavr, sci_bzbeenvr = calculate_spectral_envelope(
        sci_wavm[iorderflatbreak:iordermax, :], sci_blzm[iorderflatbreak:iordermax, :]
    )
    # sci_bzewav = np.concatenate((sci_bzewavb, sci_bzewavr), axis=0)
    # sci_bzbeenv = np.concatenate((sci_bzbeenvb, sci_bzbeenvr), axis=0)

    # Interpolate the spectral envelope to the wavelength grid of the science data
    sci_bzenvb = interp1d(
        sci_bzewavb, sci_bzbeenvb, bounds_error=False, fill_value="extrapolate"
    )(sci_wavm[iordermin:iorderflatbreak, :])
    sci_bzenvr = interp1d(
        sci_bzewavr, sci_bzbeenvr, bounds_error=False, fill_value="extrapolate"
    )(sci_wavm[iorderflatbreak:iordermax, :])
    sci_bzenv = np.concatenate((sci_bzenvb, sci_bzenvr), axis=0)
    sci_bzenvn = sci_bzenv / np.nanmax(sci_bzenv, axis=1, keepdims=True)

    # Normalize the science data by the spectral envelope, correcting the lamp SED
    sci_dblzed_flx = sci_dflxm[iordermin:iordermax, :] * sci_bzenvn

    # Resample the flux-conserving and stitch the orders
    st_wave, st_flux = resample_flux_conserving(
        sci_wav[iordermin:iordermax, :],
        sci_dblzed_flx,
        spec_mask[iordermin:iordermax, :],
        inst_stitch_config,
    )

    return st_wave, st_flux
