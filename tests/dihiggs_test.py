import pytest

from rosetta import HiggsBasis as HB
from rosetta.interfaces.dihiggs import dihiggs
from itertools import combinations_with_replacement as comb


@pytest.mark.xfail(reason="Requires eHDECAY to be installed")
def test_dihiggs_utilities(request):
    h_channels = {
        "bb": (5, -5),
        "mumu": (13, -13),
        "tautau": (15, -15),
        "gammagamma": (22, 22),
        "ZZ": (23, 23),
        "WW": (24, -24),
    }
    hh_channels = {}

    for ch1, ch2 in comb(list(h_channels.keys()), 2):
        if ch2 == "bb":
            ch1, ch2 = ch2, ch1

        id1, id2 = h_channels[ch1], h_channels[ch2]

        if ch1 == ch2:
            ch = "4" + ch1[: len(ch1) // 2]
        else:
            ch = "2" + ch1[: len(ch1) // 2] + "2" + ch2[: len(ch2) // 2]

        hh_channels[ch] = (h_channels[ch1], h_channels[ch2])

    cards_dir = request.path.parent / 'Cards'

    instance = HB.HiggsBasis(
        flavor="universal", param_card=cards_dir / "HiggsBasis_universal_1e-3.dat"
    )

    xs, err, decays = dihiggs.get_xs_and_br(instance)

    print(xs)
    # print channels

    print("  channel       BR             xs x BR (fb)")
    print("-------------------------------------------")

    for ch, (id1, id2) in sorted(hh_channels.items()):
        BR = decays[id1] * decays[id2]
        if id1 != id2:
            BR *= 2.0

        mustr = "  {:<10}    {:.3e}      {:.3e}".format(ch, BR, xs * BR)

        print(mustr)
