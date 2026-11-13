"""Shared test setup."""

from astroquery.simbad import conf as simbad_conf

# SIMBAD at CDS Strasbourg has been offline since 2026-11-09 after a
# cybersecurity incident (notice at https://simbad.cds.unistra.fr/), which
# fails every test that queries it. Use the Harvard mirror until CDS is back,
# then delete this file.
simbad_conf.server = "simbad.harvard.edu"
