from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4

from rvdata.tests.regression.compliance import (
    check_l2_compliance, check_l3_compliance, check_l4_compliance,
)


def test_base():
    l2 = RV2()
    l2_standard = "./base_L2_standard.fits"
    l2.to_fits(l2_standard)
    check_l2_compliance(l2_standard)

    l3 = RV3()
    l3_standard = "./base_L3_standard.fits"
    l3.to_fits(l3_standard)
    check_l3_compliance(l3_standard)

    l4 = RV4()
    l4_standard = "./base_L4_standard.fits"
    l4.to_fits(l4_standard)
    check_l4_compliance(l4_standard)


if __name__ == "__main__":
    test_base()
