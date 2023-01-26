# This file is part of xtb.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

add_subdirectory("api")
add_subdirectory("coulomb")
add_subdirectory("disp")
add_subdirectory("docking")
add_subdirectory("extern")
add_subdirectory("freq")
add_subdirectory("gfnff")
add_subdirectory("io")
add_subdirectory("iff")
add_subdirectory("main")
add_subdirectory("mctc")
add_subdirectory("param")
add_subdirectory("prog")
add_subdirectory("solv")
add_subdirectory("tblite")
add_subdirectory("type")
add_subdirectory("xtb")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/aespot.f90"
  "${dir}/anharmlib.f90"
  "${dir}/aoparam.f90"
  "${dir}/approxrab.f90"
  "${dir}/axis_trafo.f90"
  "${dir}/basic_geo.f90"
  "${dir}/bfgs.f90"
  "${dir}/bias_path.f90"
  "${dir}/blowsy.f90"
  "${dir}/broyden.f90"
  "${dir}/charge_model.f90"
  "${dir}/coffee.f90"
  "${dir}/constrain_param.f90"
  "${dir}/constrain_pot.f90"
  "${dir}/constr.f90"
  "${dir}/cqpath.f90"
  "${dir}/cube.f90"
  "${dir}/david2.f90"
  "${dir}/david.f"
  "${dir}/define.f90"
  "${dir}/detrotra.f90"
  "${dir}/dipole.f90"
  "${dir}/drsp.f"
  "${dir}/dtrafo.f90"
  "${dir}/dynamic.f90"
  "${dir}/eeq_model.f90"
  "${dir}/elem.f90"
  "${dir}/embedding.f90"
  "${dir}/enso_printout.f90"
  "${dir}/eqrot.f90"
  "${dir}/esp.f"
  "${dir}/expire.f90"
  "${dir}/features.F90"
  "${dir}/filetools.F90"
  "${dir}/fixparam.f90"
  "${dir}/foden.f90"
  "${dir}/fragment.f90"
  "${dir}/generate_wsc.f90"
  "${dir}/geoopt_driver.f90"
  "${dir}/geosum.f90"
  "${dir}/getkopt.f90"
  "${dir}/getname.f90"
  "${dir}/getsymnum.f90"
  "${dir}/gfn0param.f90"
  "${dir}/gfn_paramset.f90"
  "${dir}/gfn_prparam.f90"
  "${dir}/grad_core.f90"
  "${dir}/header.f90"
  "${dir}/hessian.F90"
  "${dir}/hlex.f90"
  "${dir}/ifind.f90"
  "${dir}/iniq.f90"
  "${dir}/intgrad.f90"
  "${dir}/intmodes.f90"
  "${dir}/intpack.f90"
  "${dir}/lidep.f90"
  "${dir}/lindh.f90"
  "${dir}/lin.f90"
  "${dir}/lin_mod.f90"
  "${dir}/local.f90"
  "${dir}/locmode.f90"
  "${dir}/lopt.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/makel.f90"
  "${dir}/matinv.f90"
  "${dir}/mdoptim.f90"
  "${dir}/metadynamic.f90"
  "${dir}/modef.f90"
  "${dir}/model_hessian.f90"
  "${dir}/neighbor.f90"
  "${dir}/onetri.f90"
  "${dir}/oniom.f90"
  "${dir}/optimizer.f90"
  "${dir}/pbc.f90"
  "${dir}/pbc_tools.f90"
  "${dir}/peeq_module.f90"
  "${dir}/pocketscan.f90"
  "${dir}/pqn.f90"
  "${dir}/printmold.f90"
  "${dir}/printout.f90"
  "${dir}/prmat.f90"
  "${dir}/pseudodiag.f90"
  "${dir}/qpot.f90"
  "${dir}/qsort.f90"
  "${dir}/rdcoord2.f90"
  "${dir}/read_gfn_param.f90"
  "${dir}/readin.f90"
  "${dir}/readl2.f90"
  "${dir}/readl.f90"
  "${dir}/relaxation_engine.f90"
  "${dir}/restart.f90"
  "${dir}/rmrottr.f90"
  "${dir}/rmsd.f90"
  "${dir}/scan_driver.f90"
  "${dir}/scanparam.f90"
  "${dir}/scc_core.f90"
  "${dir}/scf_module.F90"
  "${dir}/screening.f90"
  "${dir}/set_module.f90"
  "${dir}/setparam.f90"
  "${dir}/shake_module.f90"
  "${dir}/shiftlp.f90"
  "${dir}/single.f90"
  "${dir}/slater.f90"
  "${dir}/sphereparam.f90"
  "${dir}/spline2.f90"
  "${dir}/spline3.f90"
  "${dir}/splitparam.f90"
  "${dir}/stm.f"
  "${dir}/surfac.f"
  "${dir}/thermo.f90"
  "${dir}/timing.f90"
  "${dir}/topology.f90"
  "${dir}/wrbas.f90"
  "${dir}/wrgbw.f90"
  "${dir}/wrmodef.f90"
  "${dir}/wrmo.f90"
  "${dir}/xhelp.f90"
  "${dir}/zmatpr.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
set(prog ${prog} PARENT_SCOPE)
