
#---------------------------------------------------------
# Import the ParFlow TCL package
#---------------------------------------------------------
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

#source pftest.tcl
#set passed 1

pfset SILO.Filetype HDF5

set tmp0 [pfload ./manaus.out.overlandsum.00001.silo]
pfsave $tmp0 -pfb perm_z2.pfb
