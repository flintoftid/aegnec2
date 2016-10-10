#
# aegnec2 - Dynamically Allocated Numerical Electromagnetics Code Version 2 
# Copyright (C) 1998-2016 Ian D. Flintoft <ian.flintoft@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

add_library ( nec2d STATIC arc.f atgn2.f blckin.f blckot.f cabc.f cang.f cmngf.f 
                           cmset.f cmss.f cmsw.f cmws.f cmww.f conect.f couple.f datagn.f 
                           db10.f db20.f efld.f eksc.f ekscx.f etmns.f facgf.f facio.f 
                           factr.f factrs.F fbar.f fblock.f fbngf.f ffld.f fflds.f 
                           filsuf.f gf.f gfil.f gfld.f gfout.f gh.f gwave.f gx.f gxx.f 
                           helix.f hfk.f hintg.f hsfld.f hsflx.f intrp.f intx.f isegno.f lfactr.f 
                           load.f ltsolv.f lunscr.f move.f nec2d.F nefld.f netwk.F nfpat.f nhfld.f 
                           parsit.f patch.f pcint.f prnt.f qdsrc.f rdpat.f readgm.f readmn.f 
                           reblk.f reflc.f rom2.f sbf.f sflds.f solgf.f solve.f solves.F subph.f 
                           tbf.f test.f timer.F trio.f unere.f upcase.f wire.f zint.f 
                           dataeq.inc dataj.inc ggrid.inc gnd.inc gwav.inc incom.inc
                           matpar.inc nec2d.inc smat.inc tmh.inc tmi.inc )

add_executable( fnec2 fdrive.f limits.inc nec2d.inc )
target_link_libraries( fnec2 nec2d ${NEC_LIBS} )
set_target_properties( fnec2 PROPERTIES LINKER_LANGUAGE Fortran )
