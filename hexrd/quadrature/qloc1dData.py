# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details on dowloading the source,
# see the file COPYING.
# 
# Please also see the file LICENSE.
# 
# This program is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this program (see file LICENSE); if not, write to
# the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA or visit <http://www.gnu.org/licenses/>.
# ============================================================
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
