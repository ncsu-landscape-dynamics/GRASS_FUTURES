#!/usr/bin/env python3
#
##############################################################################
#
# MODULE:       r.futures.validation
#
# AUTHOR(S):    Anna Petrasova (kratochanna gmail.com)
#
# PURPOSE:      FUTURES Validation (Quantity/Allocation Disagreement, Kappa Simulation)
#
# COPYRIGHT:    (C) 2016-2021 by the GRASS Development Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
##############################################################################

#%module
#% description: Module for gloabal validation of FUTURES simulation
#% keyword: raster
#% keyword: statistics
#% keyword: validation
#%end
#%option G_OPT_R_INPUT
#% key: simulated
#% description: Simulated land use raster
#% required: yes
#%end
#%option G_OPT_R_INPUT
#% key: reference
#% description: Reference land use raster
#% required: yes
#%end
#%option G_OPT_R_INPUT
#% key: original
#% label: Original land use raster
#% description: Required for kappa simulation
#% required: no
#%end
#%option
#% key: format
#% type: string
#% description: Output format
#% options: plain,shell,json
#% required: no
#% answer: plain
#%end


import sys
import atexit
import json
import grass.script as gs

try:
    from grass.script.utils import append_random
except ImportError:
    import random
    import string

    def append_random(name, suffix_length=8):
        """Add a random part to of a specified length to a name (string)
        ..note::
            This function is simplified copy from grass79.
        """
        allowed_chars = string.ascii_lowercase + string.digits
        suffix = "".join(random.choice(allowed_chars) for _ in range(suffix_length))
        return "{name}_{suffix}".format(**locals())


TMP = []


def cleanup():
    if TMP:
        gs.run_command("g.remove", flags="f", type="raster", name=TMP, quiet=True)


def main():
    simulated = options["simulated"]
    original = options["original"]
    reference = options["reference"]
    oformat = options["format"]
    # reclassify FUTURES output to binary (developed/undeveloped)
    simulated_name = simulated.split('@')[0]
    tmp_simulated = append_random(simulated_name, 8)
    TMP.append(tmp_simulated)
    gs.write_command(
        "r.reclass",
        input=simulated,
        output=tmp_simulated,
        rules="-",
        stdin="-1 = 0\n0 thru 1000 = 1",
    )
    # input_maps = [tmp_simulated, reference]
    input_maps = [simulated, reference]
    val_col = 2
    if original:
        input_maps.append(original)
        val_col = 3
    data = gs.read_command("r.stats", flags="cn", input=input_maps).strip()
    O0 = O1 = 0
    S0 = R0 = S1 = R1 = 0
    S1O0 = S0O0 = S0O1 = S1O1 = 0
    R1O0 = R0O0 = R0O1 = R1O1 = 0
    TP = FP = FN = TN = 0
    all = 0
    for line in data.splitlines():
        line = [int(n) for n in line.strip().split()]
        if line[0] == 0:
            S0 += line[val_col]
        if line[0] == 1:
            S1 += line[val_col]
        if line[1] == 0:
            R0 += line[val_col]
        if line[1] == 1:
            R1 += line[val_col]
        if line[0] == 0 and line[1] == 0:
            TN += line[val_col]
        if line[0] == 0 and line[1] == 1:
            FN += line[val_col]
        if line[0] == 1 and line[1] == 0:
            FP += line[val_col]
        if line[0] == 1 and line[1] == 1:
            TP += line[val_col]
        if original and line[0] == 1 and line[2] == 0:
            S1O0 += line[val_col]
        if original and line[0] == 0 and line[2] == 0:
            S0O0 += line[val_col]
        if original and line[0] == 0 and line[2] == 1:
            S0O1 += line[val_col]
        if original and line[0] == 1 and line[2] == 1:
            S1O1 += line[val_col]
        if original and line[2] == 1:
            O1 += line[val_col]
        if original and line[2] == 0:
            O0 += line[val_col]
        if original and line[1] == 1 and line[2] == 0:
            R1O0 += line[val_col]
        if original and line[1] == 0 and line[2] == 0:
            R0O0 += line[val_col]
        if original and line[1] == 0 and line[2] == 1:
            R0O1 += line[val_col]
        if original and line[1] == 1 and line[2] == 1:
            R1O1 += line[val_col]
        all += line[val_col]

    # compute
    all = float(all)
    quantity_disagreement = abs((TP + FP) - (TP + FN)) / all
    allocation_disagreement = 2 * min(FP, FN) / all
    p0 = (TP / all) + (TN / all)
    pe = (R0 / all) * (S0 / all) + (R1 / all) * (S1 / all)
    kappa = (p0 - pe) / (1 - pe)
    if original:
        a = (R0O0 / all) * (S0O0 / all) + (R1O0 / all) * (S1O0 / all)
        b = (R0O1 / all) * (S0O1 / all) + (R1O1 / all) * (S1O1 / all)
        pe_tr = (O0 / all) * a + (O1 / all) * b
        kappasim = (p0 - pe_tr) / (1 - pe_tr)

    # print
    if oformat == "plain":
        print(f"Quantity disagreement: {100 * quantity_disagreement:.3f} %")
        print(f"Allocation disagreement: {100 * allocation_disagreement:.3f} %")
        print(f"Kappa: {kappa:.3f}")
        if original:
            print(f"Kappa simulation: {kappasim:.5f}")
    elif oformat == "shell":
        print(f"quantity={quantity_disagreement:.5f}")
        print(f"allocation={allocation_disagreement:.5f}")
        print(f"kappa={kappa:.3f}")
        if original:
            print(f"kappasimulation={kappasim:.5f}")
    elif oformat == "json":
        out = {
            "quantity": quantity_disagreement,
            "allocation": allocation_disagreement,
            "kappa": kappa,
        }
        if original:
            out.update({"kappasimulation": kappasim})
        print(
            json.dumps(
                json.loads(json.dumps(out), parse_float=lambda x: round(float(x), 5))
            )
        )


if __name__ == "__main__":
    options, flags = gs.parser()
    atexit.register(cleanup)
    sys.exit(main())
