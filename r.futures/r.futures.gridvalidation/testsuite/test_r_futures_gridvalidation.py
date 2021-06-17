#!/usr/bin/env python3

from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestGridValidation(TestCase):

    kappasim = "kappa_output"
    allocation = "allocation"
    quantity = "quantity"
    actual = "actual"

    @classmethod
    def setUpClass(cls):
        cls.runModule("g.region", raster="lsat7_2002_30@PERMANENT")
        cls.runModule("r.unpack", input="data/result.pack", output=cls.actual)
        cls.runModule(
            "r.mapcalc",
            expression=f"{cls.actual} = if ({cls.actual} >= 0, 1, 0)",
            overwrite=True,
        )
        cls.runModule(
            "r.mapcalc",
            expression="ndvi_2002 = double(lsat7_2002_40@PERMANENT - lsat7_2002_30@PERMANENT) / double(lsat7_2002_40@PERMANENT + lsat7_2002_30@PERMANENT)",
        )
        cls.runModule(
            "r.mapcalc",
            expression="ndvi_1987 = double(lsat5_1987_40@landsat - lsat5_1987_30@landsat) / double(lsat5_1987_40@landsat + lsat5_1987_30@landsat)",
        )
        cls.runModule(
            "r.mapcalc",
            expression="urban_1987 = if(ndvi_1987 <= 0.1 && isnull(lakes), 1, if(isnull(lakes), 0, null()))",
        )
        cls.runModule(
            "r.mapcalc",
            expression="urban_2002 = if(ndvi_2002 <= 0.1 && isnull(lakes), 1, if(isnull(lakes), 0, null()))",
        )
        cls.runModule("r.slope.aspect", elevation="elevation", slope="slope")
        cls.runModule("r.grow.distance", input="lakes", distance="lakes_dist")
        cls.runModule("r.mapcalc", expression="lakes_dist_km = lakes_dist/1000.")
        cls.runModule("v.to.rast", input="streets_wake", output="streets", use="val")
        cls.runModule("r.grow.distance", input="streets", distance="streets_dist")
        cls.runModule("r.mapcalc", expression="streets_dist_km = streets_dist/1000.")
        cls.runModule(
            "r.futures.devpressure",
            input="urban_2002",
            output="devpressure",
            method="gravity",
            size=15,
            flags="n",
        )
        cls.runModule(
            "r.futures.pga",
            developed="urban_2002",
            development_pressure="devpressure",
            compactness_mean=0.4,
            compactness_range=0.05,
            discount_factor=0.1,
            patch_sizes="data/patches.txt",
            predictors=["slope", "lakes_dist_km", "streets_dist_km"],
            n_dev_neighbourhood=15,
            devpot_params="data/potential.csv",
            random_seed=2,
            num_neighbors=4,
            seed_search="random",
            development_pressure_approach="gravity",
            gamma=1.5,
            scaling_factor=1,
            subregions="zipcodes",
            demand="data/demand.csv",
            output="simulated",
        )

    @classmethod
    def tearDownClass(cls):
        cls.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=[
                "slope",
                "lakes_dist",
                "lakes_dist_km",
                "streets",
                "streets_dist",
                "streets_dist_km",
                "devpressure",
                "ndvi_2002",
                "ndvi_1987",
                "urban_1987",
                "urban_2002",
                "simulated",
                cls.actual,
            ],
        )

    def tearDown(self):
        self.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=[self.kappasim, self.quantity, self.allocation],
        )

    def test_kappasim_run(self):
        """Test if results is in expected limits"""
        self.runModule("g.region", flags="a", res=100)
        self.assertModule(
            "r.futures.gridvalidation",
            original="urban_2002",
            simulated="simulated",
            reference=self.actual,
            allocation_disagreement=self.allocation,
            quantity_disagreement=self.quantity,
            kappasimulation=self.kappasim,
        )


if __name__ == "__main__":
    test()
