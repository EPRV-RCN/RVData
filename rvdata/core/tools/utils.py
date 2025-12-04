import configparser

# --- Utility functions for configuration file parsing ---


def parse_str_to_types(string):
    """Converts string to different object types they represent.
    Supported formats: True,Flase,None,int,float,list,tuple"""
    if string == "True":
        return True
    elif string == "False":
        return False
    elif string == "None":
        return None
    elif string.lstrip("-+ ").isdigit():
        return int(string)
    elif string == "":
        return ""
    elif (string[0] in "[(") and (
        string[-1] in ")]"
    ):  # Recursively parse a list/tuple into a list
        if len(string.strip("()[]")) == 0:
            return []
        else:
            return [parse_str_to_types(s) for s in string.strip("()[]").split(",")]
    else:
        try:
            return float(string)
        except ValueError:
            return string


def create_configdict_from_file(
    configFilename, listOfConfigSections=None, flattenSections=True
):
    """Returns a configuration object as a dictionary by loading the config file.
        Values in the config files are parsed appropriately to python objects.
    Parameters
    ----------
    configFilename : str
                    File name of the config file to load
    listOfConfigSections : list (default:None)
                    Only the sections in the listOfConfigSections will be loaded to the dictionary.
                    if listOfConfigSections is None (default), all the sections in config file will be loaded.
    flattenSections: (bool, default=True)
                    True: Flattens the sections in the config file into a single level Config dictionary.
                    False: will return a dictionary of dictionaries for each section.
    Returns
    -------
    configDictionary : dictionary
                    if `flattenSections` is True (default): Flattened {key:value,..} dictionary is returned
                    else: A dictionary of dictionary is returned {Section:{key:Value},..}
    """
    configLoader = configparser.ConfigParser()
    configLoader.optionxform = str  # preserve the Case sensitivity of keys
    with open(configFilename) as cfgFile:
        configLoader.read_file(cfgFile)

    # Create a Config Dictionary
    config = {}
    if isinstance(listOfConfigSections, str):  # Convert the single string to a list
        listOfConfigSections = [listOfConfigSections]
    elif listOfConfigSections is None:
        listOfConfigSections = configLoader.sections()

    for configSection in listOfConfigSections:
        if flattenSections:  # Flatten the sections
            for key, value in configLoader.items(configSection):
                config[key] = parse_str_to_types(value)
        else:
            subConfig = {}
            for key, value in configLoader.items(configSection):
                subConfig[key] = parse_str_to_types(value)
            config[configSection] = subConfig
    return config


# ---
