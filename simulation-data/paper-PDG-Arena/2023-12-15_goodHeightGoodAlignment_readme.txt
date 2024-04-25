Here are described some of the objects contained in 2023-12-15_goodHeightGoodAlignment.RData

# Vocabulary

PDG-Light : the previous name of PDG-Arena

Inventory types :
E0 = RegdemoMonosp = RM inventories
E1A = RegdemoPlurisp = R inventories
E1B = IrregdemoMonosp = irregular and monospecific inventories (not used in the publication)
E2 = IrregdemoPlurisp = O inventories


# Simulation tables

Tables beginning by 'standYearTable' have an entry for each site-year
Tables beginning by 'standPeriodTable' have an entry for each site-period
Tables beginning by 'treeYearTable' have an entry for each tree-year
Tables beginning by 'treePeriodTable' have an entry for each tree-period 

Defined periods are 1996-2002, 2003-2006, 2007-2013 and 1996-2013 (full range of simulated year)

Suffixes of the tables indicates the inventory type used in the simulation (see vocabulary).

For these four types of table, each entry contains variables corresponding to the simulations results, in particular:
- vegAbsorbance : the proportion of radiation absorbed by vegetation (from 0 to 1)
- GPPy_abs_sim, the yearly gross primary production, in g/year
- WVIy_sim, the yearly wood volume increment in m3 / year
- BAIy_sim, the yearly basal area increment in cm2/year
- transpiration, the yearly transpiration rate, in mm/year
- REWmin, the minimal reserve of extractable water reached in year (during the average of a period, if it is a period table)

Additionnaly, the variables WVIy_mes and BAIy_mes are given and do note represent simulated variables. They are respectively, the yearly wood volume increment and basal area increment measured at stand scale.

# Simulation lists

These following objects may be hard to manipulate; we then suggest to contact the first author of the publication if you want to go further.
Lists beginning by 'simuList' contains all the simulation results at a yearly level. The output are described on the CASTANEA documentation, which is not public yet and can be shared on request.

First suffixes of simulation lists indicates the inventory type used in the simulation (see vocabulary). Suffix '_st' indicates that results were converted from individual scale to stand scale.
