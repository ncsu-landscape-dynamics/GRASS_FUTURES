#!/usr/bin/env python3
import pandas as pd

from grass.gunittest.case import TestCase
from grass.gunittest.main import test
import grass.script as gs


def compare_base_stats(reference, output_map):
    def stats(input_map):
        data = gs.read_command(
            "r.stats", flags="cn", input=["zipcodes", input_map], separator="comma"
        )
        cell_list = []
        for line in data.strip().splitlines():
            zipcode, step, cells = line.split(",")
            if step == "-1":
                cell_list.append(int(cells))
        return cell_list

    reference_stats = stats(reference)
    output_stats = stats(output_map)
    return [
        100 * abs(r - o) / (r if r > o else o)
        for r, o in zip(reference_stats, output_stats)
    ]


class TestSteering(TestCase):

    output = "pga_output"
    output_1 = "pga_output_1"
    output_2 = "pga_output_2"
    result = "result"
    result_steering = "result_steering"
    result_steering_distributed = "result_steering_distributed"
    output_demand = "output_demand_file.csv"
    output_demand_1 = "output_demand_file_1.csv"
    output_demand_2 = "output_demand_file_2.csv"
    output_development_pressure = "output_development_pressure"
    output_development_pressure_1 = "output_development_pressure_1"
    output_development_pressure_2 = "output_development_pressure_2"
    redistribution_1 = "redistribution_1"
    redistribution_2 = "redistribution_2"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule("g.region", raster="lsat7_2002_30@PERMANENT")
        cls.runModule("r.unpack", input="data/result.pack", output=cls.result)
        cls.runModule(
            "r.unpack", input="data/result_steering.pack", output=cls.result_steering
        )
        cls.runModule(
            "r.unpack",
            input="data/result_steering_distributed.pack",
            output=cls.result_steering_distributed,
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
            "r.watershed",
            elevation="elevation",
            drainage="drainage",
            stream="streams",
            threshold=1000,
        )
        cls.runModule(
            "r.watershed", elevation="elevation", basin="basin", threshold=5000
        )
        cls.runModule("r.null", map="basin", null=1000)
        cls.runModule(
            "r.stream.distance",
            stream_rast="streams",
            direction="drainage",
            elevation="elevation",
            method="downstream",
            difference="HAND",
        )
        cls.runModule("r.grow.distance", input="HAND", value="HAND_filled")
        cls.runModule(
            "r.lake",
            elevation="HAND_filled",
            water_level=8,
            lake="flood_1",
            seed="streams",
        )
        cls.runModule(
            "r.lake",
            elevation="HAND_filled",
            water_level=4,
            lake="flood_2",
            seed="streams",
        )
        cls.runModule(
            "r.lake",
            elevation="HAND_filled",
            water_level=2,
            lake="flood_3",
            seed="streams",
        )
        cls.runModule(
            "r.mapcalc",
            expression="flood_1 = if (flood_1, 0.05, null())",
            overwrite=True,
        )
        cls.runModule(
            "r.mapcalc",
            expression="flood_2 = if (flood_2, 0.1, null())",
            overwrite=True,
        )
        cls.runModule(
            "r.mapcalc",
            expression="flood_3 = if (flood_3, 0.2, null())",
            overwrite=True,
        )
        cls.runModule(
            "r.patch", input="flood_3,flood_2,flood_1", output="flood_probability"
        )
        cls.runModule("r.null", map="flood_probability", null=0)
        cls.runModule(
            "r.mapcalc", expression="acapacity = rand(-100, 100) * 0.01", seed=1
        )
        cls.runModule(
            "r.mapcalc", expression="urban_2002_steered = if (urban_2002 == 0, -1, 0)"
        )
        cls.runModule(
            "r.mapcalc", expression="urban_2002_flooded_steered = if (urban_2002 == 0, -1000, 0)"
        )
        rules_1 = """27511 = 27511
                    27513 = 27513
                    27518 = 27518
                    27539 = 27539
                    27606 = 27606
                    27607 = 27607"""
        gs.write_command(
            "r.reclass",
            input="zipcodes",
            output="subregions_1",
            rules="-",
            stdin=rules_1,
        )
        rules_2 = """27529 = 27529
                    27601 =27601
                    27603 = 27603
                    27604 = 27604
                    27605 = 27605
                    27608 = 27608
                    27610 = 27610"""
        gs.write_command(
            "r.reclass",
            input="zipcodes",
            output="subregions_2",
            rules="-",
            stdin=rules_2,
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
                "drainage",
                "HAND",
                "HAND_filled",
                "streams",
                "flood_1",
                "flood_2",
                "flood_3",
                "flood_probability",
                "acapacity",
                "basin",
                "urban_2002_steered",
                "urban_2002_flooded_steered",
                "subregions_1",
                "subregions_2",
                cls.result,
                cls.result_steering,
                cls.result_steering_distributed,
            ],
        )
        cls.del_temp_region()
        gs.try_remove(cls.output_demand)
        gs.try_remove(cls.output_demand_1)
        gs.try_remove(cls.output_demand_2)

    def tearDown(self):
        self.runModule(
            "g.remove",
            flags="f",
            type="raster",
            name=[
                self.output,
                self.output_1,
                self.output_2,
                self.output_development_pressure,
                self.output_development_pressure_1,
                self.output_development_pressure_2,
            ],
        )


    def test_pga_run_simple_steering(self):
        """Test if results is in expected limits"""
        for i in range(0, 7):
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_steered" if i == 0 else self.output,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="zipcodes",
                demand="data/demand.csv" if i == 0 else self.output_demand,
                output=self.output,
                flags="r",
                output_development_pressure=self.output_development_pressure,
                output_demand=self.output_demand,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )
        differences = compare_base_stats(self.result, self.output)
        for difference in differences:
            self.assertLess(difference, 5, "Cell count differs too much.")
        self.assertRastersNoDifference(
            actual=self.output, reference=self.result_steering, precision=1e-6
        )

    def test_pga_run_distributed_steering(self):
        """Test if results is in expected limits"""
        for i in range(0, 7):
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_steered" if i == 0 else self.output_1,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure_1,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="subregions_1",
                demand="data/demand.csv" if i == 0 else self.output_demand_1,
                output=self.output_1,
                flags="r",
                output_development_pressure=self.output_development_pressure_1,
                output_demand=self.output_demand_1,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_steered" if i == 0 else self.output_2,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure_2,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="subregions_2",
                demand="data/demand.csv" if i == 0 else self.output_demand_2,
                output=self.output_2,
                flags="r",
                output_development_pressure=self.output_development_pressure_2,
                output_demand=self.output_demand_2,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )
        self.runModule(
            "r.patch", input=[self.output_1, self.output_2], output=self.output
        )
        differences = compare_base_stats(self.result, self.output)
        for difference in differences:
            self.assertLess(difference, 6, "Cell count differs too much.")
        self.assertRastersNoDifference(
            actual=self.output,
            reference=self.result_steering_distributed,
            precision=1e-6,
        )

    def test_pga_run_flooding_steering(self):
        """Test if results is in expected limits"""
        for i in range(0, 7):
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_flooded_steered" if i == 0 else self.output,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="zipcodes",
                demand="data/demand.csv" if i == 0 else self.output_demand,
                output=self.output,
                flags="r",
                output_development_pressure=self.output_development_pressure,
                output_demand=self.output_demand,
                redistribution_matrix="data/matrix.csv",
                flood_maps_file="data/flood_depth_input_steered.csv",
                adaptive_capacity="acapacity",
                huc="basin",
                depth_damage_functions="data/damage_curves.csv",
                ddf_subregions="zipcodes",
                population_demand="data/population_demand.csv",
                response_func=[0.25, 0.5, 0.25, 0.5],
                response_stddev=0.1,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )

    def test_pga_run_flooding_distributed_steering(self):
        """Test if results is in expected limits"""
        for i in range(0, 7):
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_flooded_steered" if i == 0 else self.output_1,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure_1,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="subregions_1",
                demand="data/demand.csv" if i == 0 else self.output_demand_1,
                output=self.output_1,
                flags="r",
                output_development_pressure=self.output_development_pressure_1,
                output_demand=self.output_demand_1,
                redistribution_matrix="data/matrix.csv",
                flood_maps_file="data/flood_depth_input_steered.csv",
                adaptive_capacity="acapacity",
                huc="basin",
                depth_damage_functions="data/damage_curves.csv",
                ddf_subregions="zipcodes",
                population_demand="data/population_demand.csv",
                response_func=[0.25, 0.5, 0.25, 0.5],
                response_stddev=0.1,
                redistribution_external_output=self.redistribution_1,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )
            self.assertModule(
                "r.futures.simulation",
                developed="urban_2002_flooded_steered" if i == 0 else self.output_2,
                development_pressure="devpressure"
                if i == 0
                else self.output_development_pressure_2,
                compactness_mean=0.4,
                compactness_range=0.05,
                discount_factor=0.1,
                patch_sizes="data/patches.csv",
                predictors=["slope", "lakes_dist_km", "streets_dist_km"],
                n_dev_neighbourhood=15,
                devpot_params="data/potential.csv",
                random_seed=1,
                num_neighbors=4,
                seed_search="random",
                development_pressure_approach="gravity",
                gamma=1.5,
                scaling_factor=1,
                subregions="subregions_2",
                demand="data/demand.csv" if i == 0 else self.output_demand_2,
                output=self.output_2,
                flags="r",
                output_development_pressure=self.output_development_pressure_2,
                output_demand=self.output_demand_2,
                redistribution_matrix="data/matrix.csv",
                flood_maps_file="data/flood_depth_input_steered.csv",
                adaptive_capacity="acapacity",
                huc="basin",
                depth_damage_functions="data/damage_curves.csv",
                ddf_subregions="zipcodes",
                population_demand="data/population_demand.csv",
                response_func=[0.25, 0.5, 0.25, 0.5],
                response_stddev=0.1,
                redistribution_external_output=self.redistribution_2,
                num_steps=1,
                steering_step=i,
                overwrite=True,
            )
            if i < 6:
                part_1_df = pd.read_csv(f"{self.redistribution_1}_{i + 1}.csv", index_col=0)
                part_2_df = pd.read_csv(f"{self.redistribution_2}_{i + 1}.csv", index_col=0)
                migrants_row = (part_1_df + part_2_df).sum()
                demand_1_df = pd.read_csv(self.output_demand_1, index_col=0)
                demand_2_df = pd.read_csv(self.output_demand_2, index_col=0)
                # add all migrants to demand for next year
                demand_1_df.iloc[i + 1] += migrants_row
                demand_2_df.iloc[i + 1] += migrants_row
                demand_1_df.round(0).astype(int).to_csv(self.output_demand_1)
                demand_2_df.round(0).astype(int).to_csv(self.output_demand_2)

        self.runModule(
            "r.patch", input=[self.output_1, self.output_2], output=self.output
        )


if __name__ == "__main__":
    test()
