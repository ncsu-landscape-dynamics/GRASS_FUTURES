from grass.gunittest.case import TestCase
from grass.gunittest.main import test


class TestPGA(TestCase):

    output = 'pga_output'

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule('g.region', n=1595760, s=1494780, e=1608600, w=1509630, res=30)

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()

    def tearDown(cls):
        """Remove viewshed map after each test method"""
        # TODO: eventually, removing maps should be handled through testing framework fucntions
        cls.runModule('g.remove', flags='f', type='raster',
                      name=cls.output)

    def test_pga_run(self):
        """Test if results is in expected limits"""
        self.assertModule('r.futures.pga', developed='urban_2010', development_pressure='devpressure_0_5',
                          compactness_mean=0.4, compactness_range=0.1, discount_factor=0.6,
                          patch_sizes='data/patches.txt',
                          predictors=['euc_interchanges_100', 'euc_protect', 'forest_1km_100', 'rdens_1km_100', 'water_dist'],
                          n_dev_neighbourhood=10, devpot_params='data/potential_FIPS.csv',
                          incentive_table='data/incentive_table.txt', random_seed=1,
                          num_neighbors=4, seed_search=2, development_pressure_approach='gravity',
                          gamma=0.5, scaling_factor=1, subregions='selected_counties_id',
                          demand='data/demand_noheader.csv', output=self.output)
        self.assertRastersNoDifference(actual=self.output, reference='urban_2030', precision=1e-6)


if __name__ == '__main__':
    test()
