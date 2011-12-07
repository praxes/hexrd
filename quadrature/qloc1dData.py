# DO-NOT-DELETE revisionify.begin() 
#
#   Copyright (c) 2007-2009 Lawrence Livermore National Security,
#   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
#   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
#   
#   Please also read the file NOTICES.
#   
#   This file is part of the mdef package (version 0.2) and is
#   free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   
#   A copy of the GNU Lesser General Public License may be found in the
#   file NOTICES. If this file is missing, see
#   <http://www.gnu.org/licenses/>.
#
# DO-NOT-DELETE revisionify.end() 
import numpy as num
from math import sqrt

'data'
ndim = 1

xi_a2  = sqrt(1.0e0/3.0e0) # 0.577350269189626e0
xi_a3  = sqrt(0.6e0)

xi_8a= 0.960289856497536e0
xi_8b= 0.796666477413627e0
xi_8c= 0.525532409916329e0
xi_8d= 0.183434642495650e0
w_8a = 0.101228536290376e0
w_8b = 0.222381034453374e0
w_8c = 0.313706645877887e0
w_8d = 0.362683783378362e0
