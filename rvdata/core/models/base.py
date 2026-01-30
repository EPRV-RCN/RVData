"""
Standard models for RV data
"""

import datetime
import hashlib
import importlib
import os
import warnings
from collections import OrderedDict

import git
from git.exc import InvalidGitRepositoryError

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
            >>> data.receipt_add_entry('primitive1', 'param1', 'PASS')
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
        with fits.open(fn) as hdul:
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
        this_data.receipt_add_entry("from_fits", "PASS")

        # Return this instance
        return this_data

    # TODO: overwrite is not yet used
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
            if key in self.headers["PRIMARY"]:
                value = self.headers["PRIMARY"][key]
                parsed_value = parse_value_to_datatype(key, row["DataType"], value)
                self.headers["PRIMARY"][key] = parsed_value

    def to_fits(self, fn):
        """
        Collect the content of this instance into a monolithic FITS file

        Args:
            fn (str): file path

        Note:
            Can only write to KPF formatted FITS

        """
        if not fn.endswith(".fits"):
            # we only want to write to a '.fits file
            raise NameError("filename must end with .fits")

        hdu_list = self._create_hdul()
        # finish up writing
        hdul = fits.HDUList(hdu_list)
        dirname = os.path.dirname(fn)
        if dirname and not os.path.isdir(dirname):
            os.makedirs(dirname, exist_ok=True)
        hdul.writeto(fn, overwrite=True, output_verify="silentfix")
        hdul.close()

    # =============================================================================
    # Receipt related members
    def receipt_add_entry(self, module, status):
        """
        Add an entry to the receipt

        Args:
            module (str): Name of the module making this entry
            status (str): status to be recorded
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

        # add the row to the bottom of the table
        row = {
            "Time": time,
            "Code_Release": git_tag,
            "Commit_Hash": git_commit_hash,
            "Branch_Name": git_branch,
            "Module_Name": module,
            "Status": status,
        }

        self.receipt = pd.concat([self.receipt, pd.DataFrame([row])], ignore_index=True)

    def receipt_info(self):
        """
        Print the short version of the receipt

        Args:
            receipt_name (string): name of the receipt
        """
        rec = self.receipt
        msg = rec["Time", "Module_Name", "Status"]
        print(msg)

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

    def set_data(self, ext_name, data):
        """
        Set extension data

        Args:
            ext_name (str): name of extension, will be forced uppercase
            data: data to set extension data
        """
        # check whether the extension already exist
        if ext_name in self.extensions.keys():
            if isinstance(data, type(FITS_TYPE_MAP[self.extensions[ext_name]]([]))):
                self.data[ext_name] = data
            else:
                raise TypeError(
                    f"Data type does not match extension {ext_name} data type"
                )
        else:
            raise NameError("Name {} does not exist as extension".format(ext_name))

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
                table = Table.from_pandas(self.data[key])
                self.headers[key]["NAXIS1"] = len(table)
                head = fits.Header(self.headers[key])
                hdu = fits.BinTableHDU(data=table, header=head)
                hdu.name = hduname
                hdu_list.append(hdu)
            else:
                print(
                    "Can't translate {} into a valid FITS format.".format(
                        type(self.data[key])
                    )
                )
                continue

        return hdu_list
