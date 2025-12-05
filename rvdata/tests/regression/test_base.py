from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4

from rvdata.tests.regression.compliance import (
    check_l2_extensions,
    check_l2_header,
    check_l4_extensions,
    check_l4_header,
)


def test_base():
    l2 = RV2()
    l2_standard = "./base_L2_standard.fits"
    l2.to_fits(l2_standard)
    l2_obj = RV2.from_fits(l2_standard)

    check_l2_extensions(l2_standard)
    check_l2_header(l2_obj.headers["PRIMARY"])

    l4 = RV4()
    l4_standard = "./base_L4_standard.fits"
    l4.to_fits(l4_standard)
    l4_obj = RV4.from_fits(l4_standard)

    check_l4_extensions(l4_standard)
    check_l4_header(l4_obj.headers["PRIMARY"])


if __name__ == "__main__":
    test_base()
