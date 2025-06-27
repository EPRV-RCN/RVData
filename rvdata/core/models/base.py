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
import numpy as np
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.definitions import (
    FITS_TYPE_MAP,
    INSTRUMENT_READERS,
    LEVEL2_PRIMARY_KEYWORDS,
    LEVEL3_PRIMARY_KEYWORDS,
    LEVEL4_PRIMARY_KEYWORDS,
)
from rvdata.core.models.receipt_columns import RECEIPT_COL
from rvdata.core.tools.git import get_git_branch, get_git_revision_hash, get_git_tag


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
        A dictionary of extensions. This maps extension name to their FITS data type, e.g. PrimaryHDU, ImageHDU, BinTableHDU.
    headers : dict
        A dictionary of headers of each extension (HDU). This stores all header information from the FITS file as a dictionary with extension name as the keys and the header content as the values. Headers are stored as OrderedDict types.
    data : dict
        A dictionary of data of each extension (HDU). This stores all extension data from the FITS file as a dictionary with extension name as the keys and the data content as the values. Data type is translated from the FITS type to an appropriate Python data type by core.model.definitions.FITS_TYPE_MAP.
    receipt : pandas.DataFrame
        A table that records the history of this data. The receipt keeps track of the data process history, so that the information stored by this instance can be reproduced from the original data. It is structured as a pandas.DataFrame table, with each row as an entry. Anything that modifies the content of a data product are expected to also write to the receipt. Three string inputs from the primitive are required: name, any relevant parameters, and a status. The receipt will also automatically fill in additional information, such as the time of execution, code release version, current branch, ect. It is not recommended to modify the receipt Dataframe directly. Use the provided methods to make any adjustments, such as:
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
        self.headers = OrderedDict()  # map name to extension header
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
        this_data = cls()
        if not os.path.isfile(fn):
            raise IOError(f"{fn} does not exist.")

        # populate it with self.read()
        this_data.read(fn, instrument, **kwargs)
        # Return this instance
        return this_data

    def read(self, fn, instrument=None, overwrite=False, **kwargs):
        """
        Read the content of a RVData standard .fits file and populate this
        data structure.

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

        if not fn.endswith(".fits"):
            # Can only read .fits files
            raise IOError("input files must be FITS files")

        self.filename = os.path.basename(fn)
        self.dirname = os.path.dirname(fn)

        with fits.open(fn) as hdu_list:

            # Handles the Receipt and the auxilary HDUs
            for hdu in hdu_list:
                if isinstance(hdu, fits.PrimaryHDU):
                    self.headers[hdu.name] = hdu.header
                elif isinstance(hdu, fits.BinTableHDU):
                    t = Table.read(hdu)
                    if "RECEIPT" in hdu.name:
                        # Table contains the RECEIPT
                        df: pd.DataFrame = t.to_pandas()
                        # TODO: get receipt columns from core.models.config.BASE-RECEIPT-columns.csv
                        df = df.reindex(
                            df.columns.union(RECEIPT_COL, sort=False),
                            axis=1,
                            fill_value="",
                        )
                        setattr(self, hdu.name, df)
                        setattr(self, hdu.name.lower(), getattr(self, hdu.name))
                        self.headers[hdu.name] = hdu.header
                        setattr(self, hdu.name, t.to_pandas())

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

                    method = rvdata.core.models.level4.RV3._read
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
        for i, row in pd.concat(
            [LEVEL2_PRIMARY_KEYWORDS, LEVEL3_PRIMARY_KEYWORDS, LEVEL4_PRIMARY_KEYWORDS]
        ).iterrows():
            key = row["Keyword"]
            if key in self.headers["PRIMARY"]:
                value = self.headers["PRIMARY"][key]
                if value is None:
                    continue
                try:
                    if row["Data type"].lower() == "uint":
                        self.headers["PRIMARY"][key] = int(value)
                    elif row["Data type"].lower() == "float":
                        self.headers["PRIMARY"][key] = float(value)
                    elif row["Data type"].lower() == "string":
                        self.headers["PRIMARY"][key] = str(value)
                    elif row["Data type"].lower() == "double":
                        self.headers["PRIMARY"][key] = np.float64(value)
                    else:
                        warnings.warn(f"Unknown type {row['Type']} for keyword {key}")
                except (TypeError, AttributeError, ValueError):
                    warnings.warn(
                        f"Cannot convert value {value} for keyword {key} to type {row['Data type']}"
                    )

        # compute MD5 sum of source file and write it into a receipt entry for tracking.
        # Note that MD5 sum has known security vulnerabilities, but we are only using
        # this to ensure data integrity, and there is no known reason for someone to try
        # to hack astronomical data files.  If something more secure is is needed,
        # substitute hashlib.sha256 for hashlib.md5
        md5 = hashlib.md5()
        with open(fn, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5.update(chunk)
        self.receipt_add_entry("from_fits", "PASS")

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

        gen_hdul = getattr(self, "_create_hdul", None)
        if gen_hdul is None:
            raise TypeError("Write method not found. Is this the base class?")
        else:
            hdu_list = gen_hdul()
        # finish up writing
        hdul = fits.HDUList(hdu_list)
        if not os.path.isdir(os.path.dirname(fn)):
            os.makedirs(os.path.dirname(fn), exist_ok=True)
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
