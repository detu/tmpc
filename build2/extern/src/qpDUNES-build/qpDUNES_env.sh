##
##	This file is part of qpDUNES.
##
##	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
##	Copyright (C) 2012-2014 by Janick Frasch, Hans Joachim Ferreau et al. 
##	All rights reserved.
##
##	qpDUNES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpDUNES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpDUNES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##


################################################################################
#
# Description:
#	qpDUNES configuration file
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2013 - 2014.
#
# NOTE:
#	- Linux/Unix only.
#
# Usage:
#	- Linux - Ubuntu:
#		* Users are supposed to source this file into ~/.bashrc file.
#
################################################################################

################################################################################
#
# Definitions for both users and developers.
#
################################################################################

# 
# Tell the user project where to find our headers, libraries and external
# packages, etc.
#
export qpDUNES_ENV_INCLUDE_DIRS="/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES;/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES/include;/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES/include/qp;/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES/interfaces;/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES/interfaces/mpc"
export qpDUNES_ENV_LIBRARY_DIRS="/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-build/lib"

#
# List of qpDUNES libraries
#
export qpDUNES_ENV_STATIC_LIBRARIES="qpdunes"
