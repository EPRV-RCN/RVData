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


def parse_value_to_datatype(keyword: str, datatype: str, in_value):
    """
    Parse a value into a target datatype used for FITS header fields.

    Parameters
    - keyword (str): Header keyword name (used only for warning messages).
    - datatype (str): Target datatype. Supported types:
        'uint', 'float', 'double', 'string', 'boolean'.
    - value: The value to convert; if it's a tuple, handle it

    Returns
    - Converted value of the requested type as (value, comment) tuple.
    """

    out_value = None
    comment = ""

    if isinstance(in_value, tuple):
        value = in_value[0]
        comment = in_value[1]
    else:
        value = in_value

    try:
        if value is None:
            out_value = None
        elif isinstance(value, str) and value.lower() == "undefined":
            out_value = None
        elif datatype.lower() == "uint":
            out_value = int(value)
        elif datatype.lower() == "float":
            out_value = float(value)
        elif datatype.lower() == "string":
            out_value = str(value)
        elif datatype.lower() == "double":
            out_value = np.float64(value)
        elif datatype.lower() == "boolean":
            if isinstance(value, bool):
                out_value = value
            elif isinstance(value, str):
                out_value = value[0].lower() == "t"
        else:
            warnings.warn(f"Unknown type {datatype} for keyword {keyword}")
    except (TypeError, AttributeError, ValueError):
        warnings.warn(
            f"Cannot convert value {value} for keyword {keyword} to type {datatype}"
        )

    return (out_value, comment)
