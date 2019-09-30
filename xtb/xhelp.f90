! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine citation(iunit)
integer,intent(in) :: iunit
write(iunit,'(3x,a)') &
   "Cite this work as:", &
   "* S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput., 2017,",&
   "  13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118", &
   "* C. Bannwarth, S. Ehlert and S. Grimme., J. Chem. Theory Comput., 2019,",&
   "  15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176",&
   "* P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, ChemRxiv, 2019, preprint.",&
   "  DOI: 10.26434/chemrxiv.8326202.v1",&
   "",&
   "for DFT-D4:",&
   "* E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017,",&
   "  147, 034112. DOI: 10.1063/1.4993215", &
   "* E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,", &
   "  C. Bannwarth and S. Grimme, J. Chem. Phys., 2019, 150, 154122.", &
   "  DOI: 10.1063/1.5090222", &
   "",&
   "for sTDA-xTB:",&
   "* S. Grimme and C. Bannwarth, J. Chem. Phys., 2016, 145, 054103.",&
   "  DOI: 10.1063/1.4959605",&
   "",&
   "in the mass-spec context:",&
   "* V. Asgeirsson, C. Bauer and S. Grimme, Chem. Sci., 2017, 8, 4879.",&
   "  DOI: 10.1039/c7sc00601b",&
   "",&
   "for metadynamics refer to:",&
   "* S. Grimme, J. Chem. Theory Comput., 2019, 155, 2847-2862", &
   "  DOI: 10.1021/acs.jctc.9b00143", &
   "",&
   "with help from (in alphabetical order)",&
   "C. Bannwarth, F. Bohle, G. Brandenburg, E. Caldeweyher, M. Checinski,",&
   "S. Dohm, S. Ehlert, S. Ehrlich, F. März, H. Neugebauer, J. Pisarek,",&
   "P. Pracht, P. Shushkov, and S. Spicher.",&
   ""
end subroutine citation

subroutine help
   use iso_fortran_env, only : istdout => output_unit
   write(istdout,'(a)') &
   "Usage: xtb [options] <geometry> [options]", &
   "",&
   "<geometry> may be provided as valid TM coordinate file (*coord in Bohr) or",&
   "in xmol format (*xyz in Ångström).",&
   "",&
   "Options:",&
   "",&
   "   -c,--chrg INT     specify molecular charge as INT,",&
   "                     overrides .CHRG file and xcontrol option",&
   "   -u,--uhf INT      specify Nalpha-Nbeta as INT,",&
   "                     overrides .UHF file and xcontrol option",&
   "",&
   "   -a,--acc REAL     accuracy for SCC calculation, lower is better (default = 1.0)",&
   "",&
   "      --iteration INT number of iterations in SCC (default = 250)",&
   "",&
   "      --cycles       number of cycles in ANCopt (default = automatic)",&
   "",&
   "      --gfn INT      specify parametrisation of GFN-xTB (default = 2)",&
   "",&
   "      --qmdff        use QMDFF for single point (needs solvent-file)",&
   "",&
   "      --tm           use TURBOMOLE for single point (needs control-file)",&
   "",&
   "      --orca         use ORCA for single point (writes ORCA input)",&
   "",&
   "      --mopac        use MOPAC for single point (writes MOPAC input)",&
   "",&
   "      --periodic     uses periodic boundary conditions (in developement)",&
   "",&
   "      --etemp REAL   electronic temperature (default = 300K)",&
   "",&
   "      --vparam FILE  Parameter file for vTB calculation",&
   "",&
   "      --xparam FILE  Parameter file for xTB calculation (not used)",&
   "",&
   "   -g,--gbsa [SOLVENT [STATE]] generalized born (GB) model with",&
   "                     solvent accessable surface area (SASA) model",&
   "",&
   "      --pop          requests printout of Mulliken population analysis",&
   "",&
   "      --molden       requests printout of molden file",&
   "",&
   "      --dipole      requests dipole printout",&
   "",&
   "      --wbo          requests Wiberg bond order printout",&
   "",&
   "      --lmo          requests localization of orbitals",&
   "",&
   "      --fod          requests FOD calculation, adjusts",&
   "                     electronic temperature to 12500 K if possible",&
   "",&
   "      --scc, --sp    performs a single point calculation",&
   "",&
   "      --vip          performs calculation of ionisation potential",&
   "",&
   "      --vea          performs calculation of electron affinity",&
   "",&
   "      --vipea        performs calculation of IP and EA",&
   "",&
   "      --vomega       performs calculation of electrophilicity index",&
   "",&
   "      --vfukui       calculate Fukui indicies using GFN-xTB",&
   "",&
   "      --esp          calculate electrostatic potential on VdW-grid",&
   "",&
   "      --stm          calculate STM image",&
   "",&
   "      --grad         performs a gradient calculation",&
   "",&
   "   -o,--opt [LEVEL]  call ancopt(3) to perform a geometry optimization,",&
   "                     levels from crude, sloppy, loose, normal (default),",&
   "                     tight, verytight to extreme can be chosen",&
   "",&
   "      --optts [LEVEL] [ROOT] call ancopt(3) to perform a transition state",&
   "                     optimization, may need to perform a hessian calculation first",&
   "",&
   "      --hess         perform a numerical hessian calculation on input geometry",&
   "",&
   "      --ohess [LEVEL] perform a numerical hessian calculation on",&
   "                     an ancopt(3) optimized geometry",&
   "",&
   "      --md           molecular dynamics simulation on start geometry",&
   "",&
   "      --omd          molecular dynamics simulation on ancopt(3) optimized",&
   "                     geometry, a loose optimization level will be chosen.",&
   "",&
   "      --metadyn [INT] meta dynamics simulation on start geometry",&
   "                     saving INT snapshots to bias the simulation",&
   "",&
   "      --siman        conformational search by simulated annealing based on",&
   "                     molecular dynamics. Conformers are optimized with ancopt.",&
   "",&
   "      --modef INT    modefollowing algorithm.  INT specifies the mode",&
   "                     that should be used for the modefollowing.",&
   "",&
   "   -I,--input FILE   use FILE as input source for xcontrol(7) instructions",&
   "",&
   "      --namespace STRING give this xtb(1) run a namespace.",&
   "                     All files, even temporary ones, will be named accordingly",&
   "                     (might not work everywhere).",&
   "",&
   "      --[no]copy     copies the xcontrol file at startup (default = true)",&
   "",&
   "      --[no]restart  restarts calculation from xtbrestart (default = true)",&
   "",&
!$ "   -P,--parallel INT number of parallel processes",&
!$ "",&
   "      --define       performs automatic check of input and terminate",&
   "",&
   "      --version      print version and terminate",&
   "",&
   "      --citation     print citation and terminate",&
   "",&
   "      --license      print license and terminate",&
   "",&
   "   -v,--verbose      be more verbose (not supported in every unit)",&
   "",&
   "   -h,--help         show this message",&
   "",&
   "Useful Maschine Settings:",&
   "",&
   "export MKL_NUM_THREADS=<NCORE>",&
   "export OMP_THREAD_LIMIT=<NCORE>",&
   "export OMP_NUM_THREADS=${OMP_THREAD_LIMIT},1",&
   "export OMP_STACKSIZE=500m",&
   "ulimit -s unlimited",&
   "",&
   "Output Conventions:",&
   "",&
   "total energies are given in atomic units (Eh)",&
   "gaps/HL energies are given in eV",&
!  this was a nice joke in the devel version, but I have to take it out for the main release
!  "Please read the manual carefully (just kidding... there is no manual)",&
   "",&
   "More information can be obtained by `man 1 xtb` and `man 7 xcontrol`",&
   "or at https://xtb-docs.readthedocs.io/en/latest/contents.html",&
   ""
end subroutine help

subroutine help_legacy
   use iso_fortran_env, only : id => output_unit
   write(id,'(''Usage: xtb <geometry> [options]'',/)')

   write(id,'(''<geometry> may be provided as'','//&
   &    'x,''valid TM coordinate file (*coord in Bohr) or'','//&
   &    '/,''in xmol format (*xyz in Ångström).'',/)')

   write(id,'(''Options:'',/)')

   write(id,'(3x,''-c, --chrg <int>  '','// &
   &          'x,''Set charge to molecule'','//&
   &      '/,22x,''overrides charge in .CHRG file'')')

   write(id,'(3x,''-u, --uhf <int>   '','// &
   &          'x,''Number pf unpaired electrons or 2S'','//&
   &      '/,22x,''overrides setting in .UHF file'')')

   write(id,'(3x,''    --nox/--nodiff'','// &
   &          'x,''skip second, extended part in sTDA-xTB'')')

   write(id,'(3x,''    --pop         '','// &
   &          'x,''do population analysis'')')

   write(id,'(3x,''    --mowr <real> '','// &
   &          'x,''cut MO write above (default 3 Eh)'')')

   write(id,'(3x,''    --mold[en]    '','// &
   &          'x,''write formatted molden file'')')

   write(id,'(3x,''    --dip[ole]    '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --wbo         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --lmo         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --rdchrg      '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --md          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --hessl       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --optts       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --ohessl      '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --mdopt       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --path        '','// &
   &          'x,''deactivated feature'')')

   write(id,'(3x,''    --vip         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --vea         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --tm          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --orca        '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --LS          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --mopac       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --qmdff       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --parx <file> '','// &
   &          'x,''read parameters for sTDA-xTB calculation'',' // &
   &      '/,22x,''(default: $XTBHOME/.param_stda2.xtb)'')')

   write(id,'(3x,''    --parv <file> '','// &
   &          'x,''read parameters for vTB part in sTDA-xTB'',' // &
   &      '/,22x,''(default: $XTBHOME/.param_stda1.xtb)'')')

   write(id,'(3x,''    --xtemp <real>'','// &
   &          'x,''electronic temperature for xTB (default 0 K)'')')

   write(id,'(3x,''    --etemp <real>'','// &
   &          'x,''electronic temperature for GFN (default 300 K)'')')

   write(id,'(3x,''    --fod         '','// &
   &          'x,''calculate the FOD and write molden.input file'','// &
   &      '/,22x,''with appropriate occupations for plotting'','// &
   &      '/,22x,''(defualt T = 12500 K)'')')

   write(id,'(3x,''    --gfn1        '','// &
   &          'x,''use GFN1-xTB parametrisation'')')

   write(id,'(3x,''    --gfn2        '','// &
   &          'x,''use GFN2-xTB parametrisation'')')

   write(id,'(3x,''    --gfn2d3      '','// &
   &          'x,''GFN2-xTB with D3 instead of D4 dispersion'')')

   write(id,'(3x,''    --scc         '','// &
   &          'x,''do single point only'')')

   write(id,'(3x,''    --grad        '','// &
   &          'x,''calculate GFNn-xTB gradient'')')

   write(id,'(3x,''    --acc <real>  '','// &
   &          'x,''GFNn-xTB accuracy (default 1.0)'')')

   write(id,'(3x,''    --opt [level] '','// &
   &          'x,''optimize at GFNn-xTB level, level can be one of'',' // &
   &      '/,22x,''crude, vloose, loose, tight, vtight, extreme'')')

   write(id,'(3x,''    --hess        '','// &
   &          'x,''compute Hessian at GFNn-xTB level'')')

   write(id,'(3x,''    --ohess       '','// &
   &          'x,''optimize and compute Hessian'')')

   write(id,'(3x,''    --omd         '','// &
   &          'x,''optimize and do MD'')')

   write(id,'(3x,''    --siman       '','// &
   &          'x,''conformational search'')')

   write(id,'(3x,''    --screen      '','// &
   &          'x,''opt. loop over ensemble'')')

   write(id,'(3x,''    --gmd         '','// &
   &          'x,''annealed MD for GMD procedure'')')

   write(id,'(3x,''    --modef       '','// &
   &          'x,''vibrational mode following'')')

   write(id,'(3x,''-g, --gbsa [solv] [state]'','// &
   &      '/,22x,''use GBSA implicit solvent for solvent [solv]'',' // &
   &      '/,22x,''and solvation state [state] = reference, bar1M'','// &
   &      '/,22x,''(default: bar1M)'')')

   write(id,'(3x,''    --help        '','// &
   &          'x,''show this message'')')

   write(id,'(/,''Useful Maschine Settings:'',/)')

   write(id,'(3x,''setenv MKL_NUM_THREADS <NCORE>'')')
   write(id,'(3x,''setenv OMP_STACKSIZE 500m'')')
   write(id,'(3x,''limit stacksize unlimited'')')

   write(id,'(/,''Output Conventions:'',/)')

   write(id,'(3x,''total energies are given in atomic units (Eh)'')')
   write(id,'(3x,''gaps/HL energies are given in eV'')')

   write(id,'(/,''Please read the manual carefully'',' // &
   &         'x,''(just kidding... there is no manual)'',/)')
end subroutine help_legacy

subroutine definebanner
   use iso_fortran_env, only : id => output_unit
   write(id,"(""           _"")")
   write(id,"(""          | |"")")
   write(id,"(""  _ __ ___| |_ _   _ _ __ _ __        __"")")
   write(id,"("" | '__/ _ \ __| | | | '__| '_ \      / _|"")")
   write(id,"("" | | |  __/ |_| |_| | |  | | | |___ | |_"")")
   write(id,"("" |_|  \___|\__|\__,_|_|  |_| |_/ _ \|  _|"")")
   write(id,"(""         _       __ _         | (_) | |"")")
   write(id,"(""        | |     / _(_)         \___/|_|"")")
   write(id,"(""      __| | ___| |_ _ _ __   ___"")")
   write(id,"(""     / _` |/ _ \  _| | '_ \ / _ \"")")
   write(id,"(""    | (_| |  __/ | | | | | |  __/"")")
   write(id,"(""     \__,_|\___|_| |_|_| |_|\___|"")")
end subroutine definebanner
