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

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()

    def tearDown(cls):
        # TODO: eventually, removing maps should be handled through testing framework fucntions
        cls.runModule('g.remove', flags='f', type='raster',
                      name=[cls.output, cls.result])

    def test_pga_run(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2002', development_pressure='devpressure',
                          compactness_mean=0.4, compactness_range=0.05, discount_factor=0.1,
                          patch_sizes='data/patches.txt',
                          predictors=['slope', 'lakes_dist_km', 'streets_dist_km'],
                          n_dev_neighbourhood=15, devpot_params='data/potential.csv',
                          random_seed=1,
                          num_neighbors=4, seed_search=2, development_pressure_approach='gravity',
                          gamma=1.5, scaling_factor=1, subregions='zipcodes',
                          demand='data/demand.csv', output=self.output)
        self.assertRastersNoDifference(actual=self.output, reference=self.result, precision=1e-6)


if __name__ == '__main__':
    test()
