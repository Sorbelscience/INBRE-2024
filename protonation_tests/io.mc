##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Thu Jun 20 16:09:15 2024

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path k1_fixed_ph7.0.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 5, dime = (65, 65, 65)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 2793 atoms
Valist_getStatistics:  Max atom coordinate:  (27.596, 27.562, 25.492)
Valist_getStatistics:  Min atom coordinate:  (-19.59, -33.365, -25.864)
Valist_getStatistics:  Molecule center:  (4.003, -2.9015, -0.186)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 1.25, 1.25, 1.25
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.9375, 0.9375, 0.9375
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.75, 0.75, 0.75 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 44.003 37.0985 39.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -35.997 -42.9015 -40.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 34.003 27.0985 29.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -25.997 -32.9015 -30.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 34.003 27.0985 29.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -25.997 -32.9015 -30.186
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.3113
Vpbe_ctor2:  solute dimensions = 50.012 x 64.2494 x 54.4762
Vpbe_ctor2:  solute charge = -11
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127282
Vpbe_ctor2:  Debye length = 7.8566
Vpbe_ctor2:  zkappa2 = 1.27239
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.3 max radius
Vclist_setupGrid:  Grid lengths = (59.398, 73.139, 63.568)
Vclist_setupGrid:  Grid lower corner = (-25.696, -39.471, -31.97)
Vclist_assignAtoms:  Have 3562702 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 232.352
Vacc_storeParms:  Using 2352-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.374505
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 6.235030e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (065, 065, 065)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.089200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.536600e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 7.111220e-01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.870173e-02
Vprtstp: contraction number = 4.870173e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 5.939488e-03
Vprtstp: contraction number = 1.219564e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 8.439091e-04
Vprtstp: contraction number = 1.420845e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.329145e-04
Vprtstp: contraction number = 1.574986e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.150055e-05
Vprtstp: contraction number = 1.617622e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.723303e-06
Vprtstp: contraction number = 1.731724e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.959171e-07
Vprtstp: contraction number = 1.869085e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.157430e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 5.884290e-01
Vpmg_setPart:  lower corner = (-35.997, -42.9015, -40.186)
Vpmg_setPart:  upper corner = (44.003, 37.0985, 39.814)
Vpmg_setPart:  actual minima = (-35.997, -42.9015, -40.186)
Vpmg_setPart:  actual maxima = (44.003, 37.0985, 39.814)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.500000e-05
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.3113
Vpbe_ctor2:  solute dimensions = 50.012 x 64.2494 x 54.4762
Vpbe_ctor2:  solute charge = -11
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127282
Vpbe_ctor2:  Debye length = 7.8566
Vpbe_ctor2:  zkappa2 = 1.27239
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.3 max radius
Vclist_setupGrid:  Grid lengths = (59.398, 73.139, 63.568)
Vclist_setupGrid:  Grid lower corner = (-25.696, -39.471, -31.97)
Vclist_assignAtoms:  Have 3562702 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 232.352
Vacc_storeParms:  Using 2352-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -25.997, -32.9015, -30.186
VPMG::focusFillBound -- New mesh maxs = 34.003, 27.0985, 29.814
VPMG::focusFillBound -- Old mesh mins = -35.997, -42.9015, -40.186
VPMG::focusFillBound -- Old mesh maxs = 44.003, 37.0985, 39.814
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.376237
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 8.032840e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (065, 065, 065)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 8.779000e-03
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.393200e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.115663e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.596737e-02
Vprtstp: contraction number = 5.596737e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.761434e-03
Vprtstp: contraction number = 1.208103e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.780883e-04
Vprtstp: contraction number = 1.446569e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.546915e-04
Vprtstp: contraction number = 1.581570e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.461344e-05
Vprtstp: contraction number = 1.591130e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.887943e-06
Vprtstp: contraction number = 1.579602e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.954902e-07
Vprtstp: contraction number = 1.788838e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.502610e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 6.339150e-01
Vpmg_setPart:  lower corner = (-25.997, -32.9015, -30.186)
Vpmg_setPart:  upper corner = (34.003, 27.0985, 29.814)
Vpmg_setPart:  actual minima = (-25.997, -32.9015, -30.186)
Vpmg_setPart:  actual maxima = (34.003, 27.0985, 29.814)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.113000e-03
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.877909e+00
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Thu Jun 20 16:10:09 2024

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path k1_fixed_ph4.6.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 5, dime = (65, 65, 65)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 2793 atoms
Valist_getStatistics:  Max atom coordinate:  (27.596, 27.562, 25.492)
Valist_getStatistics:  Min atom coordinate:  (-19.59, -33.365, -25.864)
Valist_getStatistics:  Molecule center:  (4.003, -2.9015, -0.186)
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1855):  coarse grid center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO(./apbs/src/generic/nosh.c, 1860):  fine grid center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 1.25, 1.25, 1.25
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.9375, 0.9375, 0.9375
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.75, 0.75, 0.75 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1972):  coarse mesh center = 4.003 -2.9015 -0.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 44.003 37.0985 39.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -35.997 -42.9015 -40.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 34.003 27.0985 29.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -25.997 -32.9015 -30.186
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 34.003 27.0985 29.814
NOsh_setupCalcMGAUTO (./apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -25.997 -32.9015 -30.186
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.3113
Vpbe_ctor2:  solute dimensions = 50.012 x 64.2494 x 54.4762
Vpbe_ctor2:  solute charge = -11
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127282
Vpbe_ctor2:  Debye length = 7.8566
Vpbe_ctor2:  zkappa2 = 1.27239
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.3 max radius
Vclist_setupGrid:  Grid lengths = (59.398, 73.139, 63.568)
Vclist_setupGrid:  Grid lower corner = (-25.696, -39.471, -31.97)
Vclist_assignAtoms:  Have 3562702 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 232.352
Vacc_storeParms:  Using 2352-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.370109
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 6.262550e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (065, 065, 065)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 7.761000e-03
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 3.900800e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 6.955510e-01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 4.870173e-02
Vprtstp: contraction number = 4.870173e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 5.939488e-03
Vprtstp: contraction number = 1.219564e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 8.439091e-04
Vprtstp: contraction number = 1.420845e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.329145e-04
Vprtstp: contraction number = 1.574986e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.150055e-05
Vprtstp: contraction number = 1.617622e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.723303e-06
Vprtstp: contraction number = 1.731724e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.959171e-07
Vprtstp: contraction number = 1.869085e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.530650e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 6.133180e-01
Vpmg_setPart:  lower corner = (-35.997, -42.9015, -40.186)
Vpmg_setPart:  upper corner = (44.003, 37.0985, 39.814)
Vpmg_setPart:  actual minima = (-35.997, -42.9015, -40.186)
Vpmg_setPart:  actual maxima = (44.003, 37.0985, 39.814)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.200000e-05
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 37.3113
Vpbe_ctor2:  solute dimensions = 50.012 x 64.2494 x 54.4762
Vpbe_ctor2:  solute charge = -11
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127282
Vpbe_ctor2:  Debye length = 7.8566
Vpbe_ctor2:  zkappa2 = 1.27239
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.3 max radius
Vclist_setupGrid:  Grid lengths = (59.398, 73.139, 63.568)
Vclist_setupGrid:  Grid lower corner = (-25.696, -39.471, -31.97)
Vclist_assignAtoms:  Have 3562702 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 232.352
Vacc_storeParms:  Using 2352-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -25.997, -32.9015, -30.186
VPMG::focusFillBound -- New mesh maxs = 34.003, 27.0985, 29.814
VPMG::focusFillBound -- Old mesh mins = -35.997, -42.9015, -40.186
VPMG::focusFillBound -- Old mesh maxs = 44.003, 37.0985, 39.814
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.375424
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 7.665010e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (065, 065, 065)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.421500e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.166400e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.117769e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 5.596737e-02
Vprtstp: contraction number = 5.596737e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 6.761434e-03
Vprtstp: contraction number = 1.208103e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 9.780883e-04
Vprtstp: contraction number = 1.446569e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.546915e-04
Vprtstp: contraction number = 1.581570e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.461344e-05
Vprtstp: contraction number = 1.591130e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.887943e-06
Vprtstp: contraction number = 1.579602e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.954902e-07
Vprtstp: contraction number = 1.788838e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.472860e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 6.300500e-01
Vpmg_setPart:  lower corner = (-25.997, -32.9015, -30.186)
Vpmg_setPart:  upper corner = (34.003, 27.0985, 29.814)
Vpmg_setPart:  actual minima = (-25.997, -32.9015, -30.186)
Vpmg_setPart:  actual maxima = (34.003, 27.0985, 29.814)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 3.700000e-05
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.864954e+00
