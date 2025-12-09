

.. |missing| replace:: **TBD**

.. _new-reader:
Adding a New Instrument Reader
==============================

#. Make a new directory for your instrument under the instruments direcotry (e.g. ``RVData/instruments/kpf``).
#. Create an empty file name ``__init__.py`` inside your new instrument directory.
#. If useful you can add a config directory to hold static files that are needed for your reader. For example, under ``RVData/instruments/kpf/config`` we have a file named header_map.csv which defines the mapping between the KPF header keywords and the RV standard header keywords.
#. Create a new python file for the instrument reader named according to the data level you will be reading. (e.g. ``RVData/instruments/kpf/level2.py``)
#. Create a new class for this data level. This class should inherit from the appropriate base model. For example, if you are building a level2 data reader you should inherit from the ``core.models.level2.RV2`` class. See an example of the KPF level2 reader at the bottom of this page.
#. Define a method named ``_read`` which should take a FITS ``HDUlist`` object as an argument. This function should be able to take the ``HDUlist`` from your instrument and convert it into the standard Python data model as defined |missing|.
#. Add your instrument to the ``INSTRUMENT_READERS`` dictionary in ``core/models/definitions.py``. 

.. code-block:: python

    {
        'KPF': # name of your instrument { 
            'class': 'KPFRV2', # name of the class containing the instrument reader
            'method': '_read', # name of the method that reads your instrument's data
            'module': 'instruments.kpf.level2' # module path to the file that contains your instrument's class
        }
    }

.. literalinclude:: ../../rvdata/instruments/kpf/level2.py
   :language: python
   :start-after: # import base class
