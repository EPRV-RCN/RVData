"""
Standard models for RV data
"""

import datetime
import functools
import hashlib
import importlib
import inspect
import os
import re
import warnings
from collections import OrderedDict

import git
from git.exc import InvalidGitRepositoryError

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.definitions import (
    FITS_TYPE_MAP,
    INSTRUMENT_READERS,
    LEVEL2_PRIMARY_KEYWORDS,
    LEVEL3_PRIMARY_KEYWORDS,
    LEVEL4_PRIMARY_KEYWORDS,
    BASE_RECEIPT_COLUMNS,
)
from rvdata.core.tools.git import get_git_branch, get_git_revision_hash, get_git_tag
from rvdata.core.tools.headers import parse_value_to_datatype


def _format_receipt_args(bound_args):
    """Render bound call arguments as a key=value string for the RECEIPT ARGS column.

    Skips ``self`` and ``cls`` (so the receipt only shows real parameters).
    Long or multi-line values are replaced with ``<TypeName>`` to keep the
    string readable in the U256 column. The final string is truncated to
    240 characters to stay safely under the column width.
    """
    parts = []
    for name, value in bound_args.arguments.items():
        if name in ("self", "cls"):
            continue
        try:
            rendered = str(value)
        except Exception:
            rendered = f"<{type(value).__name__}>"
        # Default object repr (``<module.Class object at 0x...>``) is noise;
        # collapse it to ``<Class>``. Same for anything too long/multi-line.
        if rendered.startswith("<") and " object at 0x" in rendered:
            rendered = f"<{type(value).__name__}>"
        if "\n" in rendered or len(rendered) > 80:
            rendered = f"<{type(value).__name__}>"
        parts.append(f"{name}={rendered}")
    out = ", ".join(parts)
    if len(out) > 240:
        out = out[:237] + "..."
    return out


def receipt_logged(func):
    """Wrap an ``RVDataModel`` instance method so each call adds a RECEIPT entry.

    The wrapped method's parameters are rendered key=value into the ARGS
    column. A successful return logs ``STATUS=PASS``; an exception logs
    ``STATUS=FAIL`` and is re-raised unchanged.

    Intended for high-level public methods (e.g. ``read``,
    ``convert_level2_to_level3``). Avoid decorating low-level data ops like
    ``set_data``/``set_header`` — they are called heavily from internal code
    paths and would flood the receipt.
    """
    sig = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            bound = sig.bind(self, *args, **kwargs)
            bound.apply_defaults()
            args_str = _format_receipt_args(bound)
        except TypeError:
            # If signature binding fails, fall back to a generic marker so
            # the underlying function's TypeError still surfaces below.
            args_str = "<unbindable>"
        try:
            result = func(self, *args, **kwargs)
        except Exception:
            self.receipt_add_entry(func.__name__, args_str, "FAIL")
            raise
        self.receipt_add_entry(func.__name__, args_str, "PASS")
        return result

    return wrapper


class RVDataModel(object):
    """
    The base class for all RV data models.

    Warning:
        This class (RVDataModel) should not be used directly.
        Based on the data level of your .fits file, used the appropriate
        level specific data model.

    This is the base model for all data models. Level specific data inherit
    from this class, so any attribute and method listed here applies to all data
    models.

    Attributes
    ----------
    extensions : dict
        A dictionary of extensions. This maps extension name to their FITS data type,
        e.g. PrimaryHDU, ImageHDU, BinTableHDU.
    headers : dict
        A dictionary of headers of each extension (HDU). This stores all header information
        from the FITS file as a dictionary with extension name as the keys and the header
        content as the values. Headers are stored as OrderedDict types.
    data : dict
        A dictionary of data of each extension (HDU). This stores all extension data from
        the FITS file as a dictionary with extension name as the keys and the data content
        as the values. Data type is translated from the FITS type to an appropriate Python
        data type by core.model.definitions.FITS_TYPE_MAP.
    receipt : pandas.DataFrame
        A table that records the history of this data. The receipt keeps track of the data
        process history, so that the information stored by this instance can be reproduced
        from the original data. It is structured as a pandas.DataFrame table, with each row
        as an entry. Anything that modifies the content of a data product are expected to
        also write to the receipt. Three string inputs from the primitive are required: name,
        any relevant parameters, and a status. The receipt will also automatically fill in
        additional information, such as the time of execution, code release version, current
        branch, ect. It is not recommended to modify the receipt Dataframe directly. Use the
        provided methods to make any adjustments, such as:
            >>> from core.models.level1 import RV1
            >>> data = RV1()
            >>> data.receipt_add_entry('primitive1', 'param1=value', 'PASS')
    """

    def __init__(self):
        """
        Constructor
        """
        self.filename: str = None
        self.level = None  # level of data model is set in each derived class
        self.read_methods = INSTRUMENT_READERS

        self.extensions = OrderedDict()  # map name to FITS type
        self.headers = (
            OrderedDict()
        )  # map name to extension header; the dict is keyword to (value, comment)
        self.data = OrderedDict()  # map name to extension data
        self.receipt = pd.DataFrame([])

    # =============================================================================
    # I/O related methods
    @classmethod
    def from_fits(cls, fn, instrument=None, **kwargs):
        """
        Create a data instance from a file

        This method implys the ``read`` method for reading the file. Refer to
        it for more detail. It is assumed that the input FITS file is in RVData standard format

        Args:
            fn (str): file path (relative to the repository)
            instrument (str): name of instrument. None implies FITS file is in EPRV standard format.

        Returns:
            cls (data model class): the data instance containing the file content

        """

        # Check whether the file exists and is a FITS file
        if not os.path.isfile(fn):
            raise IOError(f"{fn} does not exist.")

        if not fn.endswith(".fits") and not fn.endswith(".fit"):
            # Can only read .fits files
            raise IOError("input files must be FITS files")

        # create an empty instance
        this_data = cls()
        # populate it with self.read()
        with fits.open(fn, memmap=False) as hdul:
            this_data.filename = os.path.basename(fn)
            this_data.dirname = os.path.dirname(fn)
            this_data.read(hdul, instrument, **kwargs)

        # compute MD5 sum of source file and write it into a receipt entry for tracking.
        # Note that MD5 sum has known security vulnerabilities, but we are only using
        # this to ensure data integrity, and there is no known reason for someone to try
        # to hack astronomical data files.  If something more secure is is needed,
        # substitute hashlib.sha256 for hashlib.md5
        md5 = hashlib.md5()
        with open(fn, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5.update(chunk)
        this_data.receipt_add_entry(
            "from_fits", f"fn={fn}, instrument={instrument}", "PASS"
        )

        # Keep the RECEIPT extension in sync with the live receipt DataFrame so
        # the read/from_fits entries added above are visible without a
        # write/read cycle.
        this_data._sync_receipt_to_extension()

        # Return this instance
        return this_data

    # TODO: overwrite is not yet used
    @receipt_logged
    def read(self, hdu_list, instrument=None, overwrite=False, **kwargs):
        """
        Read the content of a RVData standard .fits file and populate this
        data structure.

        Args:
            hdu_list (fits.HDUList): list of HDUs read from a FITS file
            instrument (str): instrument name. None implies FITS file is in EPRV standard format.
            overwrite (bool): if this instance is not empty, specifies whether to overwrite

        Raises:
            IOError: when a invalid file is presented

        Note:
            This is not a @classmethod so initialization is
            required before calling this function

        """

        # Handles the Receipt and the auxilary HDUs
        for hdu in hdu_list:
            if isinstance(hdu, fits.PrimaryHDU):
                self.headers[hdu.name] = hdu.header
            elif isinstance(hdu, fits.BinTableHDU):
                if "RECEIPT" in hdu.name:
                    t = Table.read(hdu)
                    # Table contains the RECEIPT
                    df: pd.DataFrame = t.to_pandas()
                    receipt_columns = BASE_RECEIPT_COLUMNS["Name"].tolist()
                    if df.empty:
                        df = pd.DataFrame(columns=receipt_columns)
                    else:
                        # Reindex to include all receipt_columns
                        # Avoid using fill_value in reindex() due to pandas bug
                        # with string dtype when there are multiple columns
                        all_cols = df.columns.union(receipt_columns, sort=False)
                        df = df.reindex(columns=all_cols)
                        # Fill missing columns (NaN values) with empty strings
                        df = df.fillna("")
                    self.receipt = df
                # else:
                #     t = Table.read(hdu)
                #     self.headers[hdu.name] = hdu.header
                #     self.data[hdu.name] = t.to_pandas()

            # TODO: we can fill in all the object headers and data, not just tables
            #       It's confusing because here we're using Astropy Tables,
            #       but the FITS_TYPE_MAP is calling for Pandas DataFrames.
            # TODO: really we should be using create_extension()

        # Leave the rest of HDUs to level specific readers
        # assume reader method and class names folow RV2 conventions
        lvl = self.level
        if instrument is None:
            if lvl == 2:
                import rvdata.core.models.level2

                method = rvdata.core.models.level2.RV2._read
                method(self, hdu_list)
            elif lvl == 3:
                import rvdata.core.models.level3

                method = rvdata.core.models.level3.RV3._read
                method(self, hdu_list)
            elif lvl == 4:
                import rvdata.core.models.level4

                method = rvdata.core.models.level4.RV4._read
                method(self, hdu_list)
        elif instrument in self.read_methods.keys():
            clsname = self.read_methods[instrument]["class"]
            methname = self.read_methods[instrument]["method"]
            modname = self.read_methods[instrument]["module"]
            if lvl != 2:
                modname = modname.replace("level2", "level{}".format(lvl))
                clsname = clsname.replace("RV2", "RV{}".format(lvl))
                methname = methname.replace("level2", "level{}".format(lvl))

            module = importlib.import_module(modname)

            cls = getattr(module, clsname)
            method = getattr(cls, methname)
            method(self, hdu_list, **kwargs)
        else:
            # the provided data_type is not recognized, ie.
            # not in the self.read_methods list
            raise IOError("cannot recognize data type {}".format(instrument))

        # check and recast the headers into appropriate types
        for _, row in pd.concat(
            [LEVEL2_PRIMARY_KEYWORDS, LEVEL3_PRIMARY_KEYWORDS, LEVEL4_PRIMARY_KEYWORDS]
        ).iterrows():
            key = row["Keyword"]
            header = self.headers["PRIMARY"]
            if key in header:
                # A fits.Header keeps the comment in a parallel store, so
                # indexing returns only the scalar. Pass the comment through
                # as a (value, comment) tuple so the recast preserves it
                # rather than overwriting it with an empty string.
                if isinstance(header, fits.Header):
                    value = (header[key], header.comments[key])
                else:
                    value = header[key]
                parsed_value = parse_value_to_datatype(key, row["DataType"], value)
                header[key] = parsed_value

    def to_fits(self, out_filedir=None, out_filename=None):
        """
        Collect the content of this instance into a monolithic FITS file

        Args:
            out_filedir (str, optional): file path. If not provided,
                automatically set to the current working directory.
            out_filename (str, optional): file base name. If not
                provided, automatically generates a filename following
                the EPRV standard naming convention.

        Returns:
            out_filepath (str): The full path to the output file.

        Note:
            Filename should follow the EPRV naming convention:
            inst_SL#_YYYYMMDDThhmmss.fits

        """
        # Set the output file path to working directory if not given one
        if out_filedir is None:
            out_filedir = os.getcwd()

        # Auto-generate filename if not provided
        if out_filename is None:
            out_filename = self.generate_standard_filename()
        # If filename is provided, see if it has a path attached
        else:
            if os.path.dirname(out_filename):
                out_filedir = os.path.dirname(out_filename)
                out_filename = os.path.basename(out_filename)

        # Ensure that the output file name is a fits file!
        if not out_filename.endswith(".fits"):
            # we only want to write to a '.fits file
            raise NameError("filename must end with .fits")

        # Construct full output path
        out_filepath = os.path.join(out_filedir, out_filename)

        # Add receipt entry before writing (so it's included in the file).
        # Kept manual (not @receipt_logged) because the entry must be added
        # *before* _sync_receipt_to_extension below, so the to_fits row
        # actually lands in the serialized RECEIPT extension.
        self.receipt_add_entry("to_fits", f"out_filepath={out_filepath}", "PASS")
        # Check filename convention and warn if it doesn't match
        self.check_filename_convention(out_filepath)

        # Update FILENAME header to match the actual filename being written
        basename = os.path.basename(out_filepath)
        if "PRIMARY" in self.headers:
            self.headers["PRIMARY"]["FILENAME"] = (basename, "Name of the FITS file")

        # Materialize the accumulated receipt entries into the RECEIPT
        # extension so they get serialized below.
        self._sync_receipt_to_extension()

        hdu_list = self._create_hdul()
        # finish up writing
        hdul = fits.HDUList(hdu_list)
        dirname = os.path.dirname(out_filepath)
        if dirname and not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
        hdul.writeto(out_filepath, overwrite=True, output_verify="silentfix")
        hdul.close()

        return out_filepath

    # =============================================================================
    # Filename convention methods

    # Pattern for standard filename: inst_SL#_YYYYMMDDThhmmss.fits
    FILENAME_PATTERN = re.compile(
        r"^([a-zA-Z]+)_SL([234])_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})\.fits$"
    )

    def generate_standard_filename(self):
        """
        Generate a standard filename based on the EPRV naming convention.

        The convention is: inst_SL#_YYYYMMDDThhmmss.fits
        where:
        - inst: instrument name (lowercase)
        - SL#: standard level (SL2, SL3, or SL4)
        - YYYYMMDDThhmmss: observation date/time

        Returns:
            str: The generated standard filename

        Raises:
            ValueError: If required header values (INSTRUME, DATE-OBS) are missing
                or if the data level is not set
        """
        # Get instrument name from PRIMARY header
        if "PRIMARY" not in self.headers:
            raise ValueError("PRIMARY header not found")

        instrume = self.headers["PRIMARY"].get("INSTRUME")
        if instrume is None:
            raise ValueError("INSTRUME keyword not found in PRIMARY header")
        # Handle tuple format (value, comment)
        if isinstance(instrume, tuple):
            instrume = instrume[0]
        instrume = str(instrume).lower()

        # Get level
        if self.level is None:
            raise ValueError("Data level not set")
        level = self.level

        # Get DATE-OBS from PRIMARY header
        date_obs = self.headers["PRIMARY"].get("DATE-OBS")
        if date_obs is None:
            raise ValueError("DATE-OBS keyword not found in PRIMARY header")
        # Handle tuple format (value, comment)
        if isinstance(date_obs, tuple):
            date_obs = date_obs[0]

        # Parse DATE-OBS (format: YYYY-MM-DDTHH:MM:SS.sss or similar)
        # Remove any fractional seconds and parse
        date_str = str(date_obs).split(".")[0]  # Remove fractional seconds
        try:
            dt = datetime.datetime.fromisoformat(date_str)
        except ValueError:
            raise ValueError(f"Cannot parse DATE-OBS value: {date_obs}")

        # Format as YYYYMMDDThhmmss
        datetime_str = dt.strftime("%Y%m%dT%H%M%S")

        return f"{instrume}_SL{level}_{datetime_str}.fits"

    def validate_filename(self, filename):
        """
        Validate if a filename matches the EPRV naming convention.

        The convention is: inst_SL#_YYYYMMDDThhmmss.fits

        Args:
            filename (str): The filename to validate (basename only, no path)

        Returns:
            bool: True if the filename matches the convention, False otherwise
        """
        # Get just the basename if a path was provided
        basename = os.path.basename(filename)
        return bool(self.FILENAME_PATTERN.match(basename))

    def check_filename_convention(self, filename):
        """
        Check if the filename follows the EPRV naming convention and issue
        a warning if it doesn't.

        Args:
            filename (str): The filename to check

        Returns:
            bool: True if the filename matches the convention, False otherwise
        """
        basename = os.path.basename(filename)
        if not self.validate_filename(basename):
            try:
                suggested = self.generate_standard_filename()
                warnings.warn(
                    f"Filename '{basename}' does not follow the EPRV naming convention. "
                    f"Suggested filename: '{suggested}'"
                )
            except ValueError:
                warnings.warn(
                    f"Filename '{basename}' does not follow the EPRV naming convention "
                    f"(inst_SL#_YYYYMMDDThhmmss.fits)"
                )
            return False
        return True

    # =============================================================================
    # Receipt related members
    def receipt_info(self):
        """
        Print the short version of the receipt to stdout.
        """
        rec = self.receipt
        if rec.empty:
            print("(receipt is empty)")
            return
        cols = [c for c in ("TIME", "FUNCTION", "STATUS") if c in rec.columns]
        print(rec[cols].to_string(index=False))

    def receipt_add_entry(self, function, args, status):
        """
        Add an entry to the receipt

        Args:
            function (str): Name of the function/primitive making this entry
            args (str): Arguments or parameters relevant to this entry
                (use ``""`` if not applicable)
            status (str): Status to be recorded
        """

        # time of execution in ISO format
        time = datetime.datetime.now().isoformat()

        # get version control info (git)
        repo = git.Repo(search_parent_directories=True)
        try:
            git_commit_hash = repo.head.object.hexsha
            git_branch = repo.active_branch.name
            git_tag = str(repo.tags[-1])
        except (
            TypeError,
            IndexError,
            InvalidGitRepositoryError,
        ):  # expected if running in testing env or using pip-installed package
            from packaging.version import parse
            from importlib.metadata import PackageNotFoundError, version

            try:
                # setuptools_scm is now used to manage versioning based on git tags and hashes.
                # If the version lies beyond the last tag, the version will be the last tag with
                # the minor version incremented and a string appended of the form:
                # .dev{N}+g{hash} where {N} is the number of commits beyond the last tag and {hash}
                # is the abbreviated commit hash. importlib and packaging provide standard ways
                # of obtaining the version of a package and parsing the version string.
                translator_version = parse(version("rvdata"))
                git_commit_hash = str(
                    translator_version.local
                )  # this can be None so cast to str
                git_tag = translator_version.public
                if "dev" in git_tag:
                    git_branch = "develop"
                else:
                    git_branch = "main"
            except PackageNotFoundError:
                git_commit_hash = ""
                git_branch = ""
                git_tag = ""
        except (
            ValueError,
            BrokenPipeError,
        ):  # 1/10/23 behavior under Docker uncovered by hour-long testing
            git_commit_hash = get_git_revision_hash()
            git_branch = get_git_branch()
            git_tag = get_git_tag()

        # add the row to the bottom of the table. Column names match
        # BASE-RECEIPT-columns.csv (the published RECEIPT schema).
        row = {
            "TIME": time,
            "CODE_RELEASE": git_tag,
            "BRANCH_NAME": git_branch,
            "COMMIT_HASH": git_commit_hash,
            "FUNCTION": function,
            "ARGS": args,
            "STATUS": status,
        }

        self.receipt = pd.concat([self.receipt, pd.DataFrame([row])], ignore_index=True)

    def _sync_receipt_to_extension(self):
        """
        Copy ``self.receipt`` (the live DataFrame) into ``self.data["RECEIPT"]``
        so the RECEIPT extension serializes with the correct rows and string-
        typed columns.

        Values are coerced to strings (with NaN/None becoming ``""``) so
        ``Table.from_pandas`` produces string Astropy columns rather than the
        ``float64`` columns it would infer from an empty/numeric DataFrame.
        """
        df = self.receipt.fillna("").astype(str)
        self.set_data("RECEIPT", df)

    # =============================================================================
    # Extension methods

    def create_extension(self, ext_name: str, ext_type: str, header=None, data=None):
        """
        Create a new empty extension

        Args:
            ext_name (str): name of extension, will be forced uppercase
            ext_type: FITS data type as string (e.g. BinTableHDU)
            header (OrderedDict): optional header to initialize extension header
            data: optional data to initialize extension data
        """

        # enforce upper-case
        ext_name = ext_name.upper()

        # check whether the extension already exist
        if ext_name in self.extensions:
            raise NameError("Name {} already exists as extension".format(ext_name))

        self.extensions[ext_name] = ext_type

        # NOTE: OrderDict(None) and np.array(None) don't work, so use [] to init.
        if header is None:
            self.headers[ext_name] = OrderedDict([])
        else:
            self.headers[ext_name] = header
        if data is None:
            self.data[ext_name] = FITS_TYPE_MAP[ext_type]([])
        else:
            if ext_type == "BinTableHDU" and isinstance(data, pd.DataFrame):
                self.data[ext_name] = Table.from_pandas(data)
            else:
                self.data[ext_name] = FITS_TYPE_MAP[ext_type](data)

    def del_extension(self, ext_name):
        """
        Delete an existing extension

        Args:
            ext_name (str): extension name
        """

        if ext_name in self.extensions.keys():
            del self.headers[ext_name]
            del self.data[ext_name]
            del self.extensions[ext_name]
        else:
            warnings.warn(f"Cannot delete nonexistent extension {ext_name}")

    def set_header(self, ext_name, header):
        """
        Set extension header

        Args:
            ext_name (str): name of extension, will be forced uppercase
            header (OrderedDict): header to set extension header
        """
        # check whether the extension already exist
        if ext_name in self.extensions.keys():
            self.headers[ext_name] = header
        else:
            raise NameError("Name {} does not exist as extension".format(ext_name))

    def _get_min_bit_depth(self, ext_name):
        """Return the MinBitDepth requirement for an extension, or None.

        Overridden in level-specific models to look up requirements
        from the extensions config CSV.
        """
        return None

    def set_data(self, ext_name, data):
        """
        Set extension data

        Args:
            ext_name (str): name of extension, will be forced uppercase
            data: data to set extension data
        """
        # check whether the extension already exist
        if ext_name in self.extensions.keys():
            ext_type = self.extensions[ext_name]
            if ext_type == "BinTableHDU" and isinstance(data, pd.DataFrame):
                data = Table.from_pandas(data)
            if isinstance(data, type(FITS_TYPE_MAP[ext_type]([]))):
                # Enforce MinBitDepth for ImageHDU arrays
                if ext_type == "ImageHDU" and isinstance(data, np.ndarray):
                    min_depth = self._get_min_bit_depth(ext_name)
                    if (min_depth is not None
                            and data.size > 0
                            and data.dtype.itemsize * 8 < min_depth):
                        if np.issubdtype(data.dtype, np.floating):
                            target = {32: np.float32, 64: np.float64}.get(
                                min_depth, np.float64
                            )
                        elif np.issubdtype(data.dtype, np.unsignedinteger):
                            target = {8: np.uint8, 16: np.uint16,
                                      32: np.uint32, 64: np.uint64}.get(
                                min_depth, np.uint64
                            )
                        else:
                            target = {8: np.int8, 16: np.int16,
                                      32: np.int32, 64: np.int64}.get(
                                min_depth, np.int64
                            )
                        warnings.warn(
                            f"Extension '{ext_name}' has dtype {data.dtype} "
                            f"({data.dtype.itemsize * 8}-bit) but "
                            f"MinBitDepth={min_depth}. "
                            f"Upcasting to {target.__name__}.",
                            UserWarning,
                            stacklevel=2,
                        )
                        data = data.astype(target)
                self.data[ext_name] = data
            else:
                raise TypeError(
                    f"Data type does not match extension {ext_name} data type"
                )
        else:
            raise NameError("Name {} does not exist as extension".format(ext_name))

    @staticmethod
    def _restore_column_metadata(hdu, stored_header):
        """Restore TUNIT, TDISP, and TNULL cards that BinTableHDU may overwrite.

        The BinTableHDU constructor can normalize or drop column metadata
        (e.g. ``m/s`` becomes ``m s-1``).  This method copies the original
        values back from the stored header into both the HDU header and
        the Column objects so they survive ``writeto()``.
        """
        for i in range(1, len(hdu.columns) + 1):
            for kw in ("TUNIT", "TDISP", "TNULL"):
                card = f"{kw}{i}"
                if card in stored_header:
                    hdu.header[card] = stored_header[card]
            if f"TUNIT{i}" in stored_header:
                hdu.columns[i - 1].unit = stored_header[f"TUNIT{i}"]

    def _create_hdul(self):
        """
        Create an hdul in FITS format.
        This is used by the base model for writing data context to file
        """
        hdu_list = []
        hdu_definitions = self.extensions.items()
        for key, value in hdu_definitions:
            hduname = key
            if value == "PrimaryHDU":
                # Preserve comments regardless of the in-memory header form.
                # A fits.Header (e.g. set by from_fits) keeps comments in a
                # parallel store, so iterating .items() would yield bare
                # scalars and silently drop them; copy it directly instead.
                # A plain dict holds (value, comment) tuples, which the
                # element-wise assignment expands correctly.
                if isinstance(self.headers[key], fits.Header):
                    head = self.headers[key].copy()
                else:
                    head = fits.Header()
                    for keyword, content in self.headers[key].items():
                        head[keyword] = content
                hdu = fits.PrimaryHDU(header=head)
                hdu_list.insert(0, hdu)
            elif value == "ImageHDU":
                data = self.data[key]
                if data is None:
                    ndim = 0
                else:
                    ndim = len(data.shape)
                self.headers[key]["NAXIS"] = ndim
                if ndim == 0:
                    self.headers[key]["NAXIS1"] = 0
                else:
                    for d in range(ndim):
                        self.headers[key]["NAXIS{}".format(d + 1)] = data.shape[d]
                head = fits.Header(self.headers[key])
                try:
                    hdu = fits.ImageHDU(data=data, header=head)
                    hdu.name = hduname
                    hdu_list.append(hdu)
                except KeyError as ke:
                    print("KeyError exception raised: -->ke=" + str(ke))
                    print("Attempting to handle it...")
                    if str(ke) == "'bool'":
                        data = data.astype(float)
                        print("------>SHAPE=" + str(data.shape))
                        hdu = fits.ImageHDU(data=data, header=head)
                        hdu_list.append(hdu)
                    else:
                        raise KeyError("A different error...")
            elif value == "BinTableHDU":
                table = self.data[key]
                self.headers[key]["NAXIS2"] = len(table)
                head = fits.Header(self.headers[key])
                hdu = fits.BinTableHDU(data=table, header=head)
                hdu.name = hduname
                self._restore_column_metadata(hdu, self.headers[key])
                hdu_list.append(hdu)
            else:
                print(
                    "Can't translate {} into a valid FITS format.".format(
                        type(self.data[key])
                    )
                )
                continue

        return hdu_list
