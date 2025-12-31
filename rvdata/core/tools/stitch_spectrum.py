import numpy as np
from astropy import constants
from scipy.signal import savgol_filter
from scipy.interpolate import PchipInterpolator
import bindensity


def echelle_order_to_order_index(echelle_order, order_table):
    """
    Convert echelle order number to order index using the order table DataFrame.
    Parameters:
    -----------
    echelle_order : int
        The echelle order number to convert.
    order_table : pandas.DataFrame
        DataFrame containing 'echelle_order' and 'order_index' columns.
    Returns:
    --------
    int or None
        The corresponding order index, or None if not found.
    """

    try:
        return order_table.loc[order_table["echelle_order"] == echelle_order][
            "order_index"
        ].values[0]
    except IndexError:
        return None


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
    """
    Generate a wavelength grid with constant velocity spacing.

    Parameters
    ----------
    wavegrid_start : float
        Starting wavelength of the grid.
    wavegrid_end : float
        Ending wavelength of the grid.
    velpix : float
        Velocity spacing per pixel (in m/s).

    Returns
    -------
    wavegrid : 1D numpy array
        Wavelength grid with constant velocity spacing.

    2025-12-04, LPA
    """

    sl = constants.c.value  # 299792458. m/s
    avg_dlwavegrid = velpix / sl
    lwavegrid = np.arange(
        np.log(wavegrid_start), np.log(wavegrid_end) + avg_dlwavegrid, avg_dlwavegrid
    )
    wavegrid = np.exp(lwavegrid)
    return wavegrid


def _clean_sort_unique(wave, flux, mask=None, use_logwave=True):
    """Filter finite points, apply optional mask, sort by wavelength, and merge duplicate x."""
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)

    if wave.ndim != 1 or flux.ndim != 1 or wave.shape != flux.shape:
        raise ValueError(
            "wave_order and flux_order must be 1D arrays with the same shape"
        )

    good = np.isfinite(wave) & np.isfinite(flux)
    if mask is not None:
        good &= np.asarray(mask, dtype=bool)

    wave = wave[good]
    flux = flux[good]

    if wave.size < 2:
        raise ValueError(
            "Need at least 2 valid (wave, flux) points to build an envelope"
        )

    # Sort by wavelength
    s = np.argsort(wave)
    wave = wave[s]
    flux = flux[s]

    x = np.log(wave) if use_logwave else wave

    # Merge duplicate x (can happen if orders share same representative wavelength)
    xu, inv = np.unique(x, return_inverse=True)
    if xu.size != x.size:
        # average flux for duplicates
        fsum = np.zeros_like(xu, dtype=float)
        cnt = np.zeros_like(xu, dtype=int)
        np.add.at(fsum, inv, flux)
        np.add.at(cnt, inv, 1)
        flux = fsum / cnt
        x = xu

        # reconstruct wave for range checks (only needed for meta; keep consistent)
        wave = np.exp(x) if use_logwave else x

    return wave, x, flux


def sed_envelope_pchip(
    wave_order,
    flux_order,
    use_logwave=True,
    use_logflux=False,
    extrapolate=False,
    mask=None,
):
    """
    Shape-preserving PCHIP envelope through one (wave, flux) point per order.

    Returns
    -------
    evaluate : callable
        evaluate(wave_grid) -> envelope flux on wave_grid
    meta : dict
        Basic info about the fit.
    """
    wave, x, flux = _clean_sort_unique(
        wave_order, flux_order, mask=mask, use_logwave=use_logwave
    )

    if use_logflux:
        if np.any(flux <= 0):
            raise ValueError("use_logflux=True requires all used flux points to be > 0")
        y = np.log(flux)
    else:
        y = flux

    interp = PchipInterpolator(x, y, extrapolate=extrapolate)
    x_min, x_max = x.min(), x.max()

    def evaluate(wave_grid):
        lg = np.asarray(wave_grid, dtype=float)
        xg = np.log(lg) if use_logwave else lg
        yg = interp(xg)

        if not extrapolate:
            inside = np.isfinite(xg) & (xg >= x_min) & (xg <= x_max)
            yg = np.where(inside, yg, np.nan)

        return np.exp(yg) if use_logflux else yg

    meta = {
        "method": "pchip",
        "use_logwave": bool(use_logwave),
        "use_logflux": bool(use_logflux),
        "extrapolate": bool(extrapolate),
        "n_points": int(x.size),
    }
    return evaluate, meta


def mask_bad_spectrum_data(sci_flx, sci_wav, sci_blz, sci_var):
    """
    Mask bad data in the science flux, wavelength, blaze, and variance arrays.

    Parameters:
    -----------
    sci_flx : np.ndarray
        Science flux array.
    sci_wav : np.ndarray
        Science wavelength array.
    sci_blz : np.ndarray
        Science blaze array.
    sci_var : np.ndarray
        Science variance array.

    Returns:
    --------
    tuple
        Masked science flux, wavelength, blaze, variance arrays, and the combined bad data mask
    """

    # Create masks for bad data
    sciflx_mask = (sci_flx <= 1.0) | np.isnan(sci_flx)
    sciwav_mask = (sci_wav <= 0.0) | np.isnan(sci_wav)
    sciblz_mask = (sci_blz <= 1.0) | np.isnan(sci_blz)
    spec_mask = sciflx_mask | sciwav_mask | sciblz_mask

    # Mask data where the spectrum is bad or missing
    sci_flxm = np.where(spec_mask, float("nan"), sci_flx)
    sci_wavm = np.where(spec_mask, float("nan"), sci_wav)
    sci_blzm = np.where(spec_mask, float("nan"), sci_blz)
    sci_varm = np.where(spec_mask, float("nan"), sci_var)

    return sci_flxm, sci_wavm, sci_blzm, sci_varm, spec_mask


def calculate_blaze_envelope(sci_wav, sci_blz, smooth_blaze_envelope=True):
    """
    Calculate the blaze envelope for the given wavelength and blaze arrays.

    Parameters:
    -----------
    sci_wav : np.ndarray
        Science wavelength array.
    sci_blz : np.ndarray
        Science blaze array.
    smooth_blaze_envelope : bool, optional
        Whether to apply Savitzky-Golay smoothing to the blaze envelope. Default is True.

    Returns:
    --------
    np.ndarray
        The calculated blaze envelope array.
    """

    # Calculate the spectral envelope for the blaze extension
    sci_wavenv, sci_blzenv = calculate_spectral_envelope(sci_wav, sci_blz)

    # (Optional) smooth the blaze envelope using Savitzky-Golay filter
    if smooth_blaze_envelope:
        sci_blzenvs = savgol_filter(sci_blzenv, window_length=5, polyorder=2, axis=0)
    else:
        sci_blzenvs = sci_blzenv.copy()

    # Resample the blaze envelope into each order wavelength axis using PCHIP interpolation
    eval_sedenv, meta_sedenv = sed_envelope_pchip(
        sci_wavenv,
        sci_blzenvs,
        use_logwave=True,
        use_logflux=True,  # True -> if envelope varies over orders of magnitude (and flux>0)
        extrapolate=True,
    )
    return eval_sedenv(sci_wav)


def calculate_normalized_blaze_function(
    sci_wav, sci_blz, inst_stitch_config, order_table
):
    """
    Calculate the normalized blaze function for the science data.
    Handles blue and red segments separately if specified.
    Blue and red segments are defined in the instrument stitching configuration.
    Blue and red are defined by the flat calibration sources used to create the blaze.
    Ex. LDLS for blue, Quartz lamp for red.

    Parameters:
    -----------
    sci_wav : np.ndarray
        Science wavelength array.
    sci_blz : np.ndarray
        Science blaze array.
    inst_stitch_config : dict
        Instrument stitching configuration dictionary.
    order_table : pandas.DataFrame
        DataFrame containing 'echelle_order' and 'order_index' columns.

    Returns:
    --------
    np.ndarray
        The normalized blaze function array.
    """
    # Initialize normalized blaze array
    sci_blze = np.full_like(sci_blz, fill_value=np.nan)

    # Determine order index ranges for blue and red segments
    # if any of the echorder start/end values are None, pass

    orderib_start = echelle_order_to_order_index(
        inst_stitch_config["echorder_blue_start"], order_table
    )
    orderib_end = echelle_order_to_order_index(
        inst_stitch_config["echorder_blue_end"], order_table
    )
    orderir_start = echelle_order_to_order_index(
        inst_stitch_config["echorder_red_start"], order_table
    )
    orderir_end = echelle_order_to_order_index(
        inst_stitch_config["echorder_red_end"], order_table
    )

    if None in [orderib_end, orderir_start]:
        # Calculate the blaze envelope for the entire range
        sci_blzenv = calculate_blaze_envelope(
            sci_wav[orderib_start : orderir_end + 1, :],
            sci_blz[orderib_start : orderir_end + 1, :],
            smooth_blaze_envelope=True,
        )
        # Calculate the normalized blaze function for the entire range
        sci_blze[orderib_start : orderir_end + 1, :] = (
            sci_blz[orderib_start : orderir_end + 1, :] / sci_blzenv
        )
        return sci_blze
    else:
        # Calculate the blaze envelope for blue and red segments
        sci_blzenvb = calculate_blaze_envelope(
            sci_wav[orderib_start : orderib_end + 1, :],
            sci_blz[orderib_start : orderib_end + 1, :],
            smooth_blaze_envelope=True,
        )
        sci_blzenvr = calculate_blaze_envelope(
            sci_wav[orderir_start : orderir_end + 1, :],
            sci_blz[orderir_start : orderir_end + 1, :],
            smooth_blaze_envelope=True,
        )

        # Calculate the normalized blaze function and join blue and red segments
        sci_blze[orderib_start : orderib_end + 1, :] = (
            sci_blz[orderib_start : orderib_end + 1, :] / sci_blzenvb
        )
        sci_blze[orderir_start : orderir_end + 1, :] = (
            sci_blz[orderir_start : orderir_end + 1, :] / sci_blzenvr
        )
        return sci_blze


def stitch_deblazed_spectrum(wavegrid, sci_wav, sci_dflx, sci_dcov, min_orders=1):
    """
    Stitch deblazed spectrum from multiple echelle orders onto a common wavelength grid.

    Parameters:
    -----------
    wavegrid : np.ndarray
        Common wavelength grid to stitch onto.
    sci_wav : np.ndarray
        Science wavelength array for each order.
    sci_dflx : np.ndarray
        Deblazed science flux array for each order.
    sci_dcov : np.ndarray
        Deblazed science covariance array for each order.

    Returns:
    --------
    tuple
        Stitched wavelength array, stitched flux array, stitched variance array.
    """
    # Rebin each order (flux-conserving with bindensity package)
    flx_stack = []
    var_stack = []
    for iord in range(sci_dflx.shape[0]):
        flx_var_grid = bindensity.resampling(
            wavegrid,
            sci_wav[iord, :],
            sci_dflx[
                iord, :-1
            ],  # bindensity expects flux array pixel dimension to be length N-1
            cov=sci_dcov[
                :, iord, :-1
            ],  # bindensity expects cov array pixel dimension to be length N-1
            kind="cubic",
        )
        flx_stack.append(flx_var_grid[0])
        var_stack.append(flx_var_grid[1][0])

    flx_stack = np.ma.array(flx_stack, mask=~np.isfinite(flx_stack), fill_value=np.nan)
    var_stack = np.ma.array(var_stack, mask=~np.isfinite(var_stack), fill_value=np.nan)

    # Combine overlapping orders using inverse-variance weighting
    # Note: bindensity.resampling returns flux arrays of length len(wavegrid) - 1,
    # so slice wavegrid[:-1] to ensure st_wav and st_flx have matching lengths.

    st_wav = wavegrid[:-1]
    valid = (
        np.isfinite(flx_stack.filled())
        & ~flx_stack.mask
        & np.isfinite(var_stack.filled())
        & ~var_stack.mask
        & (var_stack.filled() > 0)
    )
    w = np.divide(
        1.0, var_stack.filled(), where=valid, out=np.zeros_like(var_stack.filled())
    )

    wsum = w.sum(axis=0)
    fwsum = np.nansum(w * flx_stack.filled(), axis=0, where=~flx_stack.mask)
    st_flx = np.divide(fwsum, wsum, where=wsum > 0, out=np.full_like(fwsum, np.nan))
    st_var = np.divide(1.0, wsum, where=wsum > 0, out=np.full_like(wsum, np.nan))

    # Mask regions with fewer than min_orders contributing orders
    n_contrib = valid.sum(axis=0)
    if min_orders > 1:
        bad = n_contrib < min_orders
        st_flx[bad] = np.nan
        st_var[bad] = np.nan

    return st_wav, st_flx, st_var
