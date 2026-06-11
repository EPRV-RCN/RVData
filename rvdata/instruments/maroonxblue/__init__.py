"""MAROON-X blue camera package.

Exists so the level 3 stitching config resolves at
``instruments/maroonxblue/config/maroonxblue_level3.config``. The base
``RV3.convert_level2_to_level3`` locates the config as
``instruments/{INSTRUME}/config/{INSTRUME}_level3.config`` and MAROON-X sets
``INSTRUME`` to the per-camera name (``MAROONXBLUE``/``MAROONXRED``). Making
this a package ensures the config ships as package data in the built wheel.
The reader/translator code lives in the shared ``maroonx`` package.
"""
