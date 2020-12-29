FUNCTIONALITIES OF DIFFERENT FILES USED IN PIHM (in order of how they are used within the code)
pihm.c: This is the main file from which all other functions are called.
pihm.h: This declares and defines the variable type and data structure
read_alloc.c: This function is first called in pihm.c. It reads all the PIHM input files
inititalize.c: This function is next called in pihm.c. It initializes the physiographic and mteorogological properties on each grid.
f.c: This function forms the ODE that is solved in pihm.c
f_functions.h: declares the functions used in f.c
flux_cal.c:Calculate lateral and vertical fluxes. 
f_functions.c: Contains a bunch of functions that is called directly or indirectly called in flux_cal.c. These functions evaluate gradients, effective conductivities and overland fluxes, vertical fluxes (using the function fluxCalc_Ele()), and river fluxes (fluxCalc_Riv())
is_sm_et.c: Calculates interception storage, loss and melt in a weakly coupled manner
f_decoupled.c: If the subsurface vertical flux has to be decoupled
var.c: Evaluates some of the functions used in f_decoupled.c
print.c: Prints all relevant variables
update.c: update time counter for focing time series	
