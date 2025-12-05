import unicodedata
import warnings

import numpy as np


def to_ascii_safe(input_string):
    # Normalize the string to NFKD form
    normalized_string = unicodedata.normalize("NFKD", input_string)

    # Encode to ASCII bytes, ignoring errors
    ascii_bytes = normalized_string.encode("ascii", "ignore")

    # Convert bytes back to string
    ascii_string = ascii_bytes.decode("ascii")

    return ascii_string


def parse_value_to_datatype(keyword: str, datatype: str, value):
    """
    Parse a value into a target datatype used for FITS header fields.

    Parameters
    - keyword (str): Header keyword name (used only for warning messages).
    - datatype (str): Target datatype. Supported types:
        'uint', 'float', 'double', 'string', 'boolean'.
    - value: The value to convert

    Returns
    - Converted value of the requested type
    """

    try:
        if value is None:
            return None
        if isinstance(value, str) and value.lower() == "undefined":
            return None
        if datatype.lower() == "uint":
            return int(value)
        elif datatype.lower() == "float":
            return float(value)
        elif datatype.lower() == "string":
            return str(value)
        elif datatype.lower() == "double":
            return np.float64(value)
        elif datatype.lower() == "boolean":
            if isinstance(value, bool):
                return value
            elif isinstance(value, str):
                return value[0].lower() == "t"
        else:
            warnings.warn(f"Unknown type {datatype} for keyword {keyword}")
    except (TypeError, AttributeError, ValueError):
        warnings.warn(
            f"Cannot convert value {value} for keyword {keyword} to type {datatype}"
        )
