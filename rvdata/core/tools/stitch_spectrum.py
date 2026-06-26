import numpy as np
from astropy import constants
from scipy.signal import savgol_filter
from scipy.interpolate import PchipInterpolator
from bindensity import libbindensity


def bindensity_resampling_fixed(new_x, x, y, cov=None, kind="cubic"):
    r"""
    Resample data on new bins using a linear or cubic
    interpolation of the cumulative.

    The linear interpolation of the cumulative corresponds to
    a uniform density on the bins.
    The only rule necessary to compute the interpolation is
    the conservation of the integral over each original bin
    (the value of the cumulative is fixed at each edge of a bin)
    This rule allows to determine the 2 coefficients of the interpolation.

    The cubic interpolation of the cumulative corresponds to
    a quadratic density on the bins.
    In order to minimize the correlation between new bins,
    the cubic interpolation is much simplified compared to
    a cubic spline interpolation.
    In particular, the interpolation is not C2 but only C1,
    which means that the density is continuous but not derivable.
    The rules to compute the interpolation are:
    - Conservation of the integral over each original bin
    (the value of the cumulative is fixed at each edge of a bin)
    - The derivative at an edge between two bins is fixed at
    the mean between the densities over each of the two bins.
    Those two rules allow to determine the 4 coefficients of the interpolation.

    Parameters
    ----------
    new_x : (new_n+1,) ndarray
      Edges of the new bins.
    x : (n+1,) ndarray
      Edges of the original bins.
    y : (n,) ndarray
      Density over each original bin.
    cov : (nd+1, n) ndarray, optional
      Covariance matrix of the density over original bins in lower banded form.
    kind : str, optional
      Kind of interpolation to use ('linear' or 'cubic').
      Default is 'cubic'.

    Returns
    -------
    new_y : (new_n,) ndarray
      The density resampled on the new bins.
    new_cov : (new_nd+1, new_n) ndarray
      The covariance of the new bins in lower banded form (if cov is provided).

    Notes
    -----
    bindensity_resampling_fixed is adapted from bindensity.resampling to prevent
    libbindensity.resampling_covariance function getting stuck in an inifinite
    loop for certain inputs.
    2025-01, Leonardo A. Paredes
    """

    # Check shapes
    new_n = new_x.size - 1
    n = y.size
    if x.size != n + 1:
        raise ValueError(
            "Incompatible sizes. "
            + "For n bins, x should be of size n+1 and y of size n."
        )
    if cov is not None:
        if cov.shape[1] != n:
            raise ValueError("The shape of cov should be (nd+1, n).")
        nd = cov.shape[0] - 1
    dx = x[1:] - x[:-1]
    if np.any(dx <= 0):
        raise ValueError(
            "The bin edges should be provided " + "in strictly increasing order."
        )
    # Check nans
    isdef = (y == y).astype(int)
    # Find the original bins in which the new bin edges are.
    ix = np.searchsorted(x, new_x, "right").astype(int) - 1
    # Added condition to handle case where size of y is larger than max of ix
    if n == 0 or ix.min() < 0 or ix.max() > n:
        n = int(ix.max())
    # Restrict to new bins that fall inside the original range.
    kstart = np.searchsorted(ix, 0, "left").astype(int)
    kend = np.searchsorted(ix, n, "left").astype(int)
    new_n_in = kend - kstart - 1
    new_x_in = new_x[kstart:kend]
    ix_in = ix[kstart:kend]
    new_dx_in = new_x_in[1:] - new_x_in[:-1]
    if np.any(new_dx_in <= 0):
        raise ValueError(
            "The bin edges should be provided " + "in strictly increasing order."
        )
    # Position on each original bin
    delta = new_x_in - x[ix_in]
    # Range of original bins on which explicitly depends each new bin
    if kind == "linear":
        istart = ix_in[:-1]
        iend = ix_in[1:] + 1
    elif kind == "cubic":
        isdefleft = np.insert(isdef[:-1], 0, 0)
        isdefright = np.append(isdef[1:], 0)
        dl = isdefleft[ix_in]
        dr = isdefright[ix_in]
        istart = (ix_in - dl)[:-1]
        iend = (ix_in + dr + 1)[1:]
        # Precompute useful quantities for cubic weights
        t = delta / dx[ix_in]
        t2 = t * t
        t3 = t2 * t
        # t = 0 on the left border of the original bin and 1 on the right
        # Cubic Hermite basis
        # For a cubic polynomial f
        # f(t) = f(0)*h00(t) + f(1)*h01(t) + f'(0)*h10(t) + f'(1)*h11(t)
        h01 = 3 * t2 - 2 * t3
        h10 = t3 - 2 * t2 + t
        h11 = t3 - t2
        # Integral between the original bin left edge
        # and the position of the new edge,
        # as a function of the density over the 3 consecutive original bins.
        Fkcenter = (h01 + 0.5 * (h10 + h11)) * dx[ix_in]
        Fkleft = 0.5 * h10 * dx[ix_in]
        Fkright = 0.5 * h11 * dx[ix_in]
    else:
        raise Exception("The interpolation kind must be 'linear' or 'cubic'.")

    # Avoid to compute undefined bins
    libbindensity.resampling_check_def(new_n_in, isdef, istart, iend)
    isize = iend - istart

    # Weight of each original bin density to compute the new bins
    w = np.empty(np.sum(isize))
    if kind == "linear":
        libbindensity.resampling_linear_weights(
            new_n_in, dx, new_dx_in, delta, istart, isize, w
        )
    else:
        libbindensity.resampling_cubic_weights(
            new_n_in, dl, dr, dx, new_dx_in, Fkleft, Fkcenter, Fkright, istart, isize, w
        )

    # Compute new bins density
    new_y = np.full(new_n, np.nan)
    libbindensity.resampling_y(new_n_in, kstart, istart, iend, isize, y, w, new_y)

    if cov is None:
        return new_y

    # Compute new covariance shape
    new_nd = np.empty(1, dtype=int)
    libbindensity.resampling_covariance_nd(nd, new_n_in, istart, iend, new_nd)
    new_nd = new_nd[0]

    # Compute new covariance
    new_cov = np.zeros((new_nd + 1, new_n))
    libbindensity.resampling_covariance(
        n, nd, new_n, kstart, new_n_in, cov, istart, iend, isize, w, new_cov
    )

    return (new_y, new_cov)


def echelle_order_to_order_index(echelle_order, order_table):
    """
    Convert echelle order number to order index using the order table.
    Parameters:
    -----------
    echelle_order : int
        The echelle order number to convert.
    order_table : astropy.table.Table or pandas.DataFrame
        Table containing 'ECHELLE_ORDER' and 'ORDER_INDEX' columns.
    Returns:
    --------
    int or None
        The corresponding order index, or None if not found.
    """

    try:
        mask = np.array(order_table["ECHELLE_ORDER"]) == echelle_order
        matches = np.array(order_table["ORDER_INDEX"])[mask]
        return int(matches[0])
    except (IndexError, KeyError):
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
    # Note on blaze threshold: Set to 0 (not 1.0) to handle both raw blaze
    # functions (counts, can be >1e6) and normalized blaze functions (0-1).
    # This global threshold works for all supported instruments (NEID, KPF,
    # ESPRESSO, HARPS, HARPSN, EXPRES) since valid blaze values are always > 0.
    sciflx_mask = (sci_flx <= 1.0) | np.isnan(sci_flx)
    sciwav_mask = (sci_wav <= 0.0) | np.isnan(sci_wav)
    sciblz_mask = (sci_blz <= 0.0) | np.isnan(sci_blz)
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
        Science blaze array. Can be raw (counts) or pre-normalized.
        The function normalizes by dividing by an envelope fit to the
        blaze peaks, so the output will have maximum values near 1.0
        regardless of the input scale.
    inst_stitch_config : dict
        Instrument stitching configuration dictionary.
    order_table : pandas.DataFrame
        DataFrame containing 'ECHELLE_ORDER' and 'ORDER_INDEX' columns.

    Returns:
    --------
    np.ndarray
        The normalized blaze function array with maximum values near 1.0.
        Each order's blaze is divided by an envelope fit through the
        blaze peaks across orders.
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
        Science wavelength array for each order. Can be in increasing or
        decreasing order along the pixel axis.
    sci_dflx : np.ndarray
        Deblazed science flux array for each order.
    sci_dcov : np.ndarray
        Deblazed science covariance array for each order.

    Returns:
    --------
    tuple
        Stitched wavelength array, stitched flux array, stitched variance array.
    """
    # Check wavelength direction and flip if decreasing
    # bindensity requires strictly increasing wavelengths
    # Try multiple orders to find one with sufficient valid data for direction check
    n_orders = sci_wav.shape[0]
    wav_diff = np.nan
    for offset in [0, 1, -1, 2, -2]:
        sample_order = n_orders // 2 + offset
        if 0 <= sample_order < n_orders:
            order_wav = sci_wav[sample_order, :]
            valid_mask = np.isfinite(order_wav)
            if np.sum(valid_mask) > 10:  # Need at least 10 valid points
                wav_diff = np.nanmean(np.diff(order_wav))
                if np.isfinite(wav_diff):
                    break
    # Fallback: if no order has enough valid data, assume increasing
    if not np.isfinite(wav_diff):
        wav_diff = 1.0  # Assume increasing wavelengths
    if wav_diff < 0:
        # Wavelengths are decreasing, flip all arrays along pixel axis
        sci_wav = sci_wav[:, ::-1]
        sci_dflx = sci_dflx[:, ::-1]
        sci_dcov = sci_dcov[:, :, ::-1]

    # Rebin each order (flux-conserving with bindensity package)
    flx_stack = []
    var_stack = []
    for iord in range(sci_dflx.shape[0]):
        # temporary replacement of bindensity.resampling until fixed in bindensity package
        flx_var_grid = bindensity_resampling_fixed(
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
