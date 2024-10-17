"""
Standard models for RV data
"""

import datetime
import hashlib
import importlib

# Standard dependencies
import os
import warnings
from collections import OrderedDict

import git
import pandas as pd

# External dependencies
from astropy.io import fits
from astropy.table import Table

from core.models.definitions import (
    FITS_TYPE_MAP,
    INSTRUMENT_READERS,
    LEVEL2_EXTENSIONS,
)
from core.models.receipt_columns import RECEIPT_COL

# Pipeline dependencies
from core.tools.git import get_git_branch, get_git_revision_hash, get_git_tag


class RVDataModel(object):
    """The base class for all RV data models.

    Warning:
        This class (RVDataModel) should not be used directly.
        Based on the data level of your .fits file, used the appropriate
        level specific data model.

    This is the base model for all data models. Level specific data inherit
    from this class, so any attribute and method listed here applies to all data
    models.

    Attributes:
        header (dict): a dictionary of headers of each extension (HDU)

            Header stores all header information from the FITS file. Since Each file is
            organized into extensions (HDUs), and Astropy parses each extension's
            header cards into a dictionary, this attribute is structured as a dictionary
            of Astropy header objects. The first layer is the name of the header, and the second layer
            is the name of the key.

        receipt (pandas.DataFrame): a table that records the history of this data

            The receipt keeps track of the data process history, so that the information
            stored by this instance can be reproduced from the original data. It is
            structured as a pandas.DataFrame table, with each row as an entry

            Anything that modifies the content of a data product are expected to also
            write to the receipt. Three string inputs from the primitive are required: name,
            any relevant parameters, and a status. The receipt will also automatically fill
            in additional information, such as the time of execution, code release version,
            current branch, ect.

            Note:
                It is not recommended to modify the receipt Dataframe directly. Use the provided methods
                to make any adjustments.

            Examples:
                >>> from core.models.level1 import RV1
                >>> data = RV1()
                # Add an entry into the receipt
                # Three args are required: name_of_primitive, param, status
                >>> data.receipt_add_entry('primitive1', 'param1', 'PASS')
                >>> data.receipt
                                        Time     ...  Module_Param Status
                0  2020-06-22T15:42:18.360409     ...        input1   PASS

        extensions (dict): a dictionary of extensions.

            This attribute stores any additional information that any primitive may wish to
            record to FITS. Creating an extension creates an empty extension of the given type and
            one may modify it directly. Creating an extension will also create a new key-value
            pair in header, so that one can write header keywords to the extension. When writing to
            FITS extensions are stored in the FITS data type as specified in
            core.models.definitions.FITS_TYPE_MAP (image or binary table). Whitespace or
            any symbols that may be interpreted by Python as an operator (e.g. -) are not
            allowed in extension names.

            Examples:
                >>> from core.models.level1 import RV1
                >>> data = RV1()
                # Add an extension
                # A unique name is required
                >>> data.create_extension('extension1', pd.DataFrame)
                # Access the extension by using its name as an attribute
                # Add a column called 'col1' to the Dataframe
                >>> data.extension1['col1'] = [1, 2, 3]
                >>> data.extension1['extension1']
                col1
                0     1
                1     2
                2     3
                # add a key-value pair to the header
                >>> data.header['extension1']['key'] = 'value'
                # delete the extension we just made
                >>> data.del_extension['extension1']
        config (DataFrame): two-column dataframe that stores each line of the input configuration file
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

        for i, row in LEVEL2_EXTENSIONS.iterrows():
            if row["Required"]:
                # TODO: set description and comment
                self.create_extension(row["Name"], row["DataType"])

    # =============================================================================
    # I/O related methods
    @classmethod
    def from_fits(cls, fn, instrument=None):
        """Create a data instance from a file

        This method emplys the ``read`` method for reading the file. Refer to
        it for more detail. It is assume that the input FITS file is in RVData standard format

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
        this_data.read(fn, instrument)
        # Return this instance
        return this_data

    def read(self, fn, instrument=None, overwrite=False):
        """Read the content of a RVData standard .fits file and populate this
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
            if instrument is None:
                import core.models.level2

                method = core.models.level2.RV2._read
                method(self, hdu_list)
            elif instrument in self.read_methods.keys():
                module = importlib.import_module(
                    self.read_methods[instrument]["module"]
                )
                cls = getattr(module, self.read_methods[instrument]["class"])
                method = getattr(cls, self.read_methods[instrument]["method"])
                method(self, hdu_list)
            else:
                # the provided data_type is not recognized, ie.
                # not in the self.read_methods list
                raise IOError("cannot recognize data type {}".format(instrument))
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
        except TypeError:  # expected if running in testing env
            git_commit_hash = ""
            git_branch = ""
            git_tag = ""
        except ValueError:  # 12/22/22 new behavior under Docker
            git_commit_hash = get_git_revision_hash()
            git_branch = get_git_branch()
            git_tag = get_git_tag()
        except (
            BrokenPipeError
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
        self.RECEIPT = self.receipt

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

        # check whether the extension already exist
        if ext_name in self.extensions:
            raise NameError("Name {} already exists as extension".format(ext_name))

        ext_name = ext_name.upper()
        self.extensions[ext_name] = ext_type
        self.headers[ext_name] = header
        # NOTE: can't init OrderDict(None), so use OrderedDict([])
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
            if FITS_TYPE_MAP[self.extensions[ext_name]] == type(data):
                self.data[ext_name] = data
            else:
                raise TypeError(
                    f"Data type does not match extension {ext_name} data type"
                )
        else:
            raise NameError("Name {} does not exist as extension".format(ext_name))
