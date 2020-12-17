#!/usr/bin/env python3

from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestPGA(TestCase):

    output = 'pga_output'
    result = 'result'

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', raster="lsat7_2002_30@PERMANENT")
        cls.runModule('r.unpack', input='data/result.pack', output=cls.result)
        cls.runModule('r.mapcalc',
            expression="ndvi_2002 = double(lsat7_2002_40@PERMANENT - lsat7_2002_30@PERMANENT) / double(lsat7_2002_40@PERMANENT + lsat7_2002_30@PERMANENT)")
        cls.runModule('r.mapcalc',
            expression="ndvi_1987 = double(lsat5_1987_40@landsat - lsat5_1987_30@landsat) / double(lsat5_1987_40@landsat + lsat5_1987_30@landsat)")
        cls.runModule('r.mapcalc', expression="urban_1987 = if(ndvi_1987 <= 0.1 && isnull(lakes), 1, if(isnull(lakes), 0, null()))")
        cls.runModule('r.mapcalc', expression="urban_2002 = if(ndvi_2002 <= 0.1 && isnull(lakes), 1, if(isnull(lakes), 0, null()))")
        cls.runModule('r.slope.aspect', elevation='elevation', slope='slope')
        cls.runModule('r.grow.distance', input='lakes', distance='lakes_dist')
        cls.runModule('r.mapcalc', expression="lakes_dist_km = lakes_dist/1000.")
        cls.runModule('v.to.rast', input='streets_wake', output='streets', use='val')
        cls.runModule('r.grow.distance', input='streets', distance='streets_dist')
        cls.runModule('r.mapcalc', expression="streets_dist_km = streets_dist/1000.")
        cls.runModule('r.futures.devpressure', input='urban_2002', output='devpressure', method='gravity', size=15, flags='n')
        cls.runModule('r.watershed', elevation='elevation', drainage='drainage', stream='streams', threshold=1000)
        cls.runModule('r.watershed', elevation='elevation', basin='basin', threshold=5000)
        cls.runModule('r.null', map='basin', null=1000)
        cls.runModule('r.stream.distance', stream_rast='streams', direction='drainage',
                      elevation='elevation', method='downstream', difference='HAND')
        cls.runModule('r.grow.distance', input='HAND', value='HAND_filled')
        cls.runModule('r.lake', elevation='HAND_filled', water_level=8, lake='flood_1', seed='streams')
        cls.runModule('r.lake', elevation='HAND_filled', water_level=4, lake='flood_2', seed='streams')
        cls.runModule('r.lake', elevation='HAND_filled', water_level=2, lake='flood_3', seed='streams')
        cls.runModule('r.mapcalc', expression="flood_1 = if (flood_1, 0.05, null())", overwrite=True)
        cls.runModule('r.mapcalc', expression="flood_2 = if (flood_2, 0.1, null())", overwrite=True)
        cls.runModule('r.mapcalc', expression="flood_3 = if (flood_3, 0.2, null())", overwrite=True)
        cls.runModule('r.patch', input='flood_3,flood_2,flood_1', output='flood_probability')
        cls.runModule('r.null', map='flood_probability', null=0)
        cls.runModule('r.mapcalc', expression="acapacity = rand(-100, 100) * 0.01", seed=1)

    @classmethod
    def tearDownClass(cls):
        cls.runModule('g.remove', flags='f', type='raster',
                      name=['slope', 'lakes_dist', 'lakes_dist_km', 'streets',
                            'streets_dist', 'streets_dist_km', 'devpressure',
                            'ndvi_2002', 'ndvi_1987', 'urban_1987', 'urban_2002',
                            'drainage', 'HAND', 'HAND_filled', 'streams',
                            'flood_1', 'flood_2', 'flood_3',
                            'flood_probability', 'acapacity', 'basin', cls.result])
        cls.del_temp_region()

    def tearDown(self):
        self.runModule('g.remove', flags='f', type='raster', name=self.output)

    def test_pga_run(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2002', development_pressure='devpressure',
                          compactness_mean=0.4, compactness_range=0.05, discount_factor=0.1,
                          patch_sizes='data/patches.txt',
                          predictors=['slope', 'lakes_dist_km', 'streets_dist_km'],
                          n_dev_neighbourhood=15, devpot_params='data/potential.csv',
                          random_seed=1,
                          num_neighbors=4, seed_search='random', development_pressure_approach='gravity',
                          gamma=1.5, scaling_factor=1, subregions='zipcodes',
                          demand='data/demand.csv', output=self.output)
        self.assertRastersNoDifference(actual=self.output, reference=self.result, precision=1e-6)

    def test_pga_run_patch_library_multiple_columns(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2002', development_pressure='devpressure',
                          compactness_mean=0.4, compactness_range=0.05, discount_factor=0.1,
                          patch_sizes='data/patches.csv',
                          predictors=['slope', 'lakes_dist_km', 'streets_dist_km'],
                          n_dev_neighbourhood=15, devpot_params='data/potential.csv',
                          random_seed=1,
                          num_neighbors=4, seed_search='random', development_pressure_approach='gravity',
                          gamma=1.5, scaling_factor=1, subregions='zipcodes',
                          demand='data/demand.csv', output=self.output)
        self.assertRastersNoDifference(actual=self.output, reference=self.result, precision=1e-6)

    def test_pga_flooding(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2002', development_pressure='devpressure',
                          compactness_mean=0.4, compactness_range=0.05, discount_factor=0.1,
                          patch_sizes='data/patches.txt',
                          predictors=['slope', 'lakes_dist_km', 'streets_dist_km'],
                          n_dev_neighbourhood=15, devpot_params='data/potential.csv',
                          random_seed=1,
                          num_neighbors=4, seed_search='random', development_pressure_approach='gravity',
                          gamma=1.5, scaling_factor=1, subregions='zipcodes',
                          demand='data/demand.csv',
                          hand='HAND_filled', redistribution_matrix='data/matrix.csv',
                          flood_probability='flood_probability', adaptive_capacity='acapacity',
                          huc='basin', depth_damage_functions='data/damage_curves_single.csv', ddf_subregions='zipcodes',
                          population_demand='data/population_demand.csv',
                          output=self.output)

    def test_pga_flooding_DDF(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2002', development_pressure='devpressure',
                          compactness_mean=0.4, compactness_range=0.05, discount_factor=0.1,
                          patch_sizes='data/patches.txt',
                          predictors=['slope', 'lakes_dist_km', 'streets_dist_km'],
                          n_dev_neighbourhood=15, devpot_params='data/potential.csv',
                          random_seed=1,
                          num_neighbors=4, seed_search='random', development_pressure_approach='gravity',
                          gamma=1.5, scaling_factor=1, subregions='zipcodes',
                          demand='data/demand.csv',
                          hand='HAND_filled', redistribution_matrix='data/matrix.csv',
                          flood_probability='flood_probability', adaptive_capacity='acapacity',
                          huc='basin', depth_damage_functions='data/damage_curves.csv', ddf_subregions='zipcodes',
                          population_demand='data/population_demand.csv',
                          output=self.output)

if __name__ == '__main__':
    test()
