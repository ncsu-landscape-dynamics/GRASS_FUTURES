## Required inputs by sub-model


##### All sub-models

1. **Study extent.** Specify the number of rows and columns for the entire study extent.

2.	**Number of sub-regions.** FUTURES is designed to capture variation across specified sub-regions within the full study extent. DEMAND and POTENTIAL can both be specified according to sub-regions. At this time, the same sub-regions must be used for all sub-models.

3.	**Sub-region map.** This input is a raster map that contains the sub-region index for each cell. For example, if you have 5 counties, each cell would be assigned a value between 1 and 5 where each value represents a separate county. If you do not wish to model by sub-region, all values in this map should be 1. See example file: 03_input\data\index.txt

##### DEMAND

1. **Demand table.** The demand table (.txt) tells the program how many cells to convert for each region at each time step. The first column should contain an identifier for the time step (e.g. the year) and subsequent columns should contain the number of cells to convert for each sub-region (region 1, region 2 … region n). See example file: 03_input\demand.txt

##### POTENTIAL

1.	**Parameter table.** The parameter table (.cfg) contains the coefficients for a statistical model that are used to calculate the value of development potential. These coefficients can be unique for each sub-region within the study extent. The first column should contain the sub-region id, the second column should contain the slope of the regression equation, and subsequent columns should contain the coefficients for each variable in the regression equation. The column order must match the order that the mapped predictors are listed in the configuration file. The sixth column must contain the coefficient for the development pressure variable. See example file: devpotParams.cfg

2.	**Mapped predictors (n).** For each variable that will be used to calculate development potential, provide a raster map of that variable. See example files: 03_input\data\d2urbkm.txt, 03_input\data\d2intkm.txt, 03_input\data\d2rdskm.txt, 03_input\data\gdp.txt, 03_input\data\slope.txt

3. **Constraint map.** A raster map containing the value of the constraint parameter used in scenario analysis. You can use this map to reduce the potential for development in specified locations. The value of this parameter must be between 0 and 1. This value will be multiplied by the development potential; as the value decreases from 1 to 0, the potential for development decreases. This map must be present for all scenarios. If the constraint parameter is not in use, set the value of all cells in this map to “1.” See example file: 03_input\data\weight_1.txt

4.	**Incentive table.** A table (.csv) that modifies the value of development potential, primarily used in scenario analysis. This table includes a list of possible values for POTENTIAL (column 1) followed by the appropriate transformation value (column 2). It was designed make high potential areas higher and low potential areas lower, or to make all areas more equally suitable (transformation = raise potential to the power of x). This table must be present for all scenarios, if no transformation is desired, these columns should contain the exact same values. See example files: 03_input\probLookup.csv (no transformation), 03_input\probLookup_x2.csv (transformation)

##### PGA
1.	**Land cover map T1.** A raster map representing the starting state of the landscape at the beginning of the simulation (Developed = 0; Available for Development = 1). Water, protected areas, and water should be excluded (NoData or -9999). Note: FUTURES is not designed to forecast road expansion, so existing roads are typically converted to gridlines and removed from land available for new development. See example file: 03_input\data\lc96.txt

2.	**Seed search method.** The location of a new development seed can be determined based on a uniform distribution (1) or based on the development potential (2).

3.	**Patch sizes.** A text file (.txt) that includes a list of numbers representing the distribution of patch sizes. The numbers indicate the number of cells belonging to a single patch. This can be any size distribution, or it can be calibrated according to historical patterns. See xxx for more information. See example file: 03_input\patch_sizes.txt

4.	**Discount factor.** The discount factor can be used to alter the size of simulated patches and must have a value between 0 and 1. This factor is multiplied by the patch sizes listed in the patch size file. (Note. Is it possible for this value be greater than 1 in order to make patches larger than those listed in the patch size file? I don’t see anything in the code that would force it to be between 0 and 1)

5.	**Neighbor rule.** The number of neighbors to consider in patch formation (4 or 8).

6.	**Patch compactness parameter.** The shape of development patches are controlled in the model by a “patch compactness parameter” which must have a value between 0 and 1. The mean and range of this parameter needs to be defined. As the value of the parameter increases, patches become more compact. Additional details and recommendations for calibrating the parameter can be found at xxx.

7.	**Development pressure.** The development pressure variable is recalculated at each time step. The correct approach, alpha (referred to as gamma in Meentemeyer et al. 2014), and scaling factor should be specified.  Additional details regarding development pressure calculation available at xxx.





