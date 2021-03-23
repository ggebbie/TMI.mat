# TMI
Total Matrix Intercomparison MATLAB diagnostic routines 

Following Gebbie, G., and P. Huybers:  "Total matrix intercomparison: A method for resolving the geometry of water mass pathways", J. Phys. Oceanogr., 40(8), doi:10.1175/2010JPO4272.1, 1710-1728, 2010. 

Gebbie, G., and P. Huybers. "How is the ocean filled?", Geophys. Res. Lett., 38, L06604, doi:10.1029/2011GL046769, 2011 

Gebbie, G., and P. Huybers, "The mean age of ocean waters inferred from radiocarbon observations", 2012, JPO.

Gebbie, G., "How much did Glacial North Atlantic Water shoal?", 2014, Paleoceanography.

# Codes by Geoffrey (Jake) Gebbie, WHOI, ggebbie@whoi.edu, 2009-2021.

Version 1, 07 May 2009.\
Version 2, 06 Aug 2010.\
Version 3, 21 Apr 2011 -- minor changes.\
Version 4, 13 July 2011, makes names consistent with papers.\
Version 5, 28 July 2011, add TMI transient tracer simulation model.\
Version 6, Nov 2012, bug fixes, use one LU decomp for both fwd and
                        adjoint, added global inversion example,
                        SynTraCE-21 workshop update \
Version 6.1, Jan 2013, added biogeochemical example, add
                       vector_to_field back into tarball.\
Version 6.2, July 2015, added sq.m function,
                        fixed d_all to properly divide Atlantic/Pacific and put White Sea into Arctic.\
Version 7, Sept. 2016, major improvements to transient run: 2 types of initial conditions and boundary conditions.\
Version 8, Jan. 2021, bug fixes, especially those found by Elaine McDonagh

# MAIN DIAGNOSTIC ROUTINES:

steadystate_diagnostics.m  : examples of analysis for the TMI pathways matrix.\
transient_driver.m : run a TMI transient tracer simulation model.

# DATA FILES

A_2deg_2010.mat : TMI pathways matrix with 2x2 degree horizontal
                  resolution and 33 levels  G & H 2010) \
A_4deg_2014.mat : TMI pathways matrix with 4x4 degree horizontal
                  resolution and 33 levels (reference case, Gebbie 2014 (Paleoceanography))  \
A_4deg_2012.mat : updated from 2010 TMI pathways matrix with 4x4 degree horizontal
                  resolution and 33 levels (unpublished)  \
A_4deg_2010.mat : TMI pathways matrix with 4x4 degree horizontal resolution and 33 levels \
L_4deg_2012.mat : TMI time tendency matrix (G & H 2012)\
d_all_4deg.mat:  predefined oceanographically-relevant surface dye patches.\
c_all_4deg.mat:  predefined oceanographically-relevant initial dye concentrations for transient simulation.\
tracerobs_4deg_33lev_woce.mat  : WOCE global hydrographic climatology + GISS O18, box-averaged at 4x4 resolution.\
tracerobs_4deg_33lev_TMI.mat  : TMI climatology (unpublised, produced from 2012 pathways matrix)\
tracerobs_2deg_33lev_woce.mat  : WOCE global hydrographic climatology + GISS O18, box-averaged at 2x2 resolution.


# SECONDARY MATLAB FUNCTIONS AND SCRIPTS: 

field_to_vector.m : transfer from a 3D field to a vector\
vector_to_field.m : and vice versa\
sw_dist.m  : calculate distance using the seawater toolbox.\
make_initial_conditions: for the TMI transient simulation.\
make_boundary_conditions: for the TMI transient simulation -- allow time varying surface concentration\
get_target: target surface concentration for transient run\
get_tendency: get time rate of change in transient simulation\
get_restoring_flux: calculate a Newtonian relaxation forcing at the sea surface (for transient run)\
mixit.m: used in make_initial_conditions\
objfun.m: used in fmincon (example 4)\
sq.m: replaces the squeeze function


