Here are described some of the objects contained in the .RData file

# Vocabulary

PDG-Light : the previous name of PDG-Arena

Inventory types :
E0 = RegdemoMonosp = RM inventories
E1A = RegdemoPlurisp = R inventories
E2 = IrregdemoPlurisp = O inventories


# Simulation tables

These tables are formatted simulation results. They are used to plot tree and stand scale plots of the publication.

Tables beginning by 'standYearTable' have an entry for each site-year
Tables beginning by 'standPeriodTable' have an entry for each site-period
Tables beginning by 'treeYearTable' have an entry for each tree-year
Tables beginning by 'treePeriodTable' have an entry for each tree-period 

Defined periods are 1996-2002, 2003-2006, 2007-2013 and 1996-2013 (full range of simulated year)
Suffixes of the tables indicates the inventory type used in the simulation (see vocabulary).

Example: standPeriodTable_E2 contains the formatted results for the original inventories (irregular and plurispecific) at the stand scale and averaged on periods.

For these four types of table, each entry contains variables corresponding to the simulations results, in particular:
- vegAbsorbance : the proportion of radiation absorbed by vegetation (from 0 to 1)
- GPPy_abs_sim, the yearly gross primary production, in g/year
- WVIy_sim, the yearly wood volume increment in m3 / year
- BAIy_sim, the yearly basal area increment in cm2/year
- transpiration, the yearly transpiration rate, in mm/year
- REWmin, the minimal reserve of extractable water reached in year (during the average of a period, if it is a period table)

Additionnaly, the variables WVIy_mes and BAIy_mes are given and do note represent simulated variables. They are respectively, the yearly wood volume increment and basal area increment measured at stand scale.


# Simulation lists

These lists are raw results of the simulation and are not available in the repository. We suggest to contact the first author of the publication if you want to go further.

Lists beginning by 'simuList' contains all the simulation results at a yearly level. The output are described on the CASTANEA documentation, which is not public yet and can be shared on request. First suffixes of simulation lists indicates the inventory type used in the simulation (see vocabulary). Suffix '_st' indicates that results were converted from individual scale to stand scale.

Example: simuListRegdemoMonosp_st contains raw simulation results for regular and monospecific inventories, converted at the stand scale. 