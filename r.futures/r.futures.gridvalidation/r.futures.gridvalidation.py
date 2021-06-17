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


def main():
    simulated = options["simulated"]
    original = options["original"]
    reference = options["reference"]
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
        "O0": (append_random("O0", 8), append_random("O0", 8), f"{original} == 0"),
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
        "all": (
            append_random("all", 8),
            append_random("all", 8),
            f" ! isnull({simulated})",
        ),
    }
    kappa = ["R1S1", "O0", "all", "R1O0", "S1O0", "R0S0"]
    quantity = ["R1", "S1", "all"]
    allocation = ["R1S0", "R0S1", "all"]

    compute_maps = []
    if kappasimulation:
        compute_maps.extend(kappa)
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
        gs.run_command("g.remove", type="raster", name=cmap[0], flags="f")

    if kappasimulation:
        expr = f"""eval(pe = (float({maps['O0'][1]}) / {maps['all'][1]}) * (float({maps['R1O0'][1]}) / {maps['all'][1]}) * (float({maps['S1O0'][1]}) / {maps['all'][1]}), \
             p0 = float({maps['R1S1'][1]}) / {maps['all'][1]} + float({maps['R0S0'][1]}) / {maps['all'][1]})
        {kappasimulation} = (p0 - pe) / (1 - pe)
        """
        gs.mapcalc(expr)
    if quantity_disagreement:
        expr = f"""{quantity_disagreement} = abs({maps['R1'][1]} - {maps['S1'][1]}) / float({maps['all'][1]})"""
        gs.mapcalc(expr)
    if allocation_disagreement:
        expr = f"""{allocation_disagreement} = 2 * min({maps['R1S0'][1]}, {maps['R0S1'][1]}) / float({maps['all'][1]})"""
        gs.mapcalc(expr)
    gs.run_command(
         "g.remove", type="raster", name=[maps[r][1] for r in compute_maps], flags="f"
    )


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
