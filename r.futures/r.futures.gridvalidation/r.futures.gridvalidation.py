#!/usr/bin/env python3
#
##############################################################################
#
# MODULE:       r.futures.gridvalidation
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
#% description: Module for validating FUTURES simulation on a grid
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
#%option G_OPT_R_OUTPUT
#% key: allocation_disagreement
#% required: no
#% description: Output quantity disagreement raster
#%end
#%option G_OPT_R_OUTPUT
#% key: quantity_disagreement
#% description: Output quantity disagreement raster
#% required: no
#%end
#%option G_OPT_R_OUTPUT
#% key: kappa
#% description: Output Cohen's kappa raster
#% required: no
#%end
#%option G_OPT_R_OUTPUT
#% key: kappasimulation
#% description: Output kappa simulation raster
#% required: no
#%end
#%rules
#% required: allocation_disagreement, quantity_disagreement, kappasimulation
#%end
#%rules
#% requires: kappasimulation, original
#%end


import sys
import atexit
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
    kappa = options["kappa"]
    kappasimulation = options["kappasimulation"]
    quantity_disagreement = options["quantity_disagreement"]
    allocation_disagreement = options["allocation_disagreement"]

    maps = {
        "R1S1": (
            append_random("R1S1", 8),
            append_random("R1S1", 8),
            f"{reference} == 1 && {simulated} >= 0",
        ),
        "R0S0": (
            append_random("R0S0", 8),
            append_random("R0S0", 8),
            f"{reference} == 0 && {simulated} < 0",
        ),
        "R1S0": (
            append_random("R1S0", 8),
            append_random("R1S0", 8),
            f"{reference} == 1 && {simulated} < 0",
        ),
        "R0S1": (
            append_random("R0S1", 8),
            append_random("R0S1", 8),
            f"{reference} == 0 && {simulated} >= 0",
        ),
        "R1O0": (
            append_random("R1O0", 8),
            append_random("R1O0", 8),
            f"{reference} == 1 && {original} == 0",
        ),
        "S1O0": (
            append_random("S1O0", 8),
            append_random("S1O0", 8),
            f"{simulated} >= 0 && {original} == 0",
        ),
        "R1O1": (
            append_random("R1O0", 8),
            append_random("R1O0", 8),
            f"{reference} == 1 && {original} == 1",
        ),
        "S1O1": (
            append_random("S1O0", 8),
            append_random("S1O0", 8),
            f"{simulated} >= 0 && {original} == 1",
        ),
        "R0O0": (
            append_random("R0O0", 8),
            append_random("R0O0", 8),
            f"{reference} == 0 && {original} == 0",
        ),
        "S0O0": (
            append_random("S0O0", 8),
            append_random("S0O0", 8),
            f"{simulated} < 0 && {original} == 0",
        ),
        "R0O1": (
            append_random("R0O1", 8),
            append_random("R0O1", 8),
            f"{reference} == 0 && {original} == 1",
        ),
        "S0O1": (
            append_random("S0O1", 8),
            append_random("S0O1", 8),
            f"{simulated} < 0 && {original} == 1",
        ),
        "O0": (append_random("O0", 8), append_random("O0", 8), f"{original} == 0"),
        "O1": (append_random("O1", 8), append_random("O1", 8), f"{original} == 1"),
        "R1": (
            append_random("R1", 8),
            append_random("R1", 8),
            f"{reference} == 1",
        ),
        "S1": (
            append_random("S1", 8),
            append_random("S1", 8),
            f"{simulated} >= 0",
        ),
        "R0": (
            append_random("R0", 8),
            append_random("R0", 8),
            f"{reference} == 0",
        ),
        "S0": (
            append_random("S0", 8),
            append_random("S0", 8),
            f"{simulated} < 0",
        ),
        "all": (
            append_random("all", 8),
            append_random("all", 8),
            f" ! isnull({simulated})",
        ),
    }
    kappasim = [
        "O0",
        "O1",
        "all",
        "R1O0",
        "R0S0",
        "R1S1",
        "R0O0",
        "S0O0",
        "R1O0",
        "S1O0",
        "R0O1",
        "S0O1",
        "R1O1",
        "S1O1",
    ]
    kapp = ["R1S1", "R0S0", "R1", "S1", "R0", "S0", "all"]
    quantity = ["R1", "S1", "all"]
    allocation = ["R1S0", "R0S1", "all"]

    compute_maps = []
    if kappa:
        compute_maps.extend(kapp)
    if kappasimulation:
        compute_maps.extend(kappasim)
    if quantity_disagreement:
        compute_maps.extend(quantity)
    if allocation_disagreement:
        compute_maps.extend(allocation)
    compute_maps = list(set(compute_maps))

    for compute_map in compute_maps:
        cmap = maps[compute_map]
        gs.use_temp_region()
        gs.run_command("g.region", align=simulated)
        gs.mapcalc(f"{cmap[0]} = if ({cmap[2]}, 1, null())")
        gs.del_temp_region()
        gs.run_command("r.resamp.stats", input=cmap[0], output=cmap[1], method="count")
        gs.run_command("g.remove", type="raster", name=cmap[0], flags="f", quiet=True)
        TMP.append(cmap[1])

    maps = {k: v[1] for k, v in maps.items()}
    if kappa:
        expr = f"""eval(R1S1 = float({maps['R1S1']}) / {maps['all']}, \
                        R0S0 = float({maps['R0S0']}) / {maps['all']}, \
                        R0 = float({maps['R0']}) / {maps['all']}, \
                        S0 = float({maps['S0']}) / {maps['all']}, \
                        R1 = float({maps['R1']}) / {maps['all']}, \
                        S1 = float({maps['S1']}) / {maps['all']}, \
                        p0 = R1S1 + R0S0, \
                        pe = R0 * S0 + R1 * S1)
                        {kappa} = (p0 - pe) / (1 - pe)
                        """
        gs.mapcalc(expr)
    if kappasimulation:
        expr = f"""eval(R0O0 = float({maps['R0O0']}) / {maps['all']}, \
                        S0O0 = float({maps['S0O0']}) / {maps['all']}, \
                        R1O0 = float({maps['R1O0']}) / {maps['all']}, \
                        S1O0 = float({maps['S1O0']}) / {maps['all']}, \
                        R0O1 = float({maps['R0O1']}) / {maps['all']}, \
                        S0O1 = float({maps['S0O1']}) / {maps['all']}, \
                        R1O1 = float({maps['R1O1']}) / {maps['all']}, \
                        S1O1 = float({maps['S1O1']}) / {maps['all']}, \
                        R0S0 = float({maps['R0S0']}) / {maps['all']}, \
                        R1S1 = float({maps['R1S1']}) / {maps['all']}, \
                        O0 = float({maps['O0']}) / {maps['all']}, \
                        O1 = float({maps['O1']}) / {maps['all']}, \
                        p0 = R1S1 + R0S0, \
                        pe = O0 * (R0O0 * S0O0 + R1O0 * S1O0) + O1 * (R0O1 * S0O1 + R1O1 * S1O1))
                        {kappasimulation} = (p0 - pe) / (1 - pe)
            """
        gs.mapcalc(expr)
    if quantity_disagreement:
        expr = f"""{quantity_disagreement} = abs({maps['R1']} - {maps['S1']}) / float({maps['all']})"""
        gs.mapcalc(expr)
    if allocation_disagreement:
        expr = f"""{allocation_disagreement} = 2 * min({maps['R1S0']}, {maps['R0S1']}) / float({maps['all']})"""
        gs.mapcalc(expr)


if __name__ == "__main__":
    options, flags = gs.parser()
    atexit.register(cleanup)
    sys.exit(main())
