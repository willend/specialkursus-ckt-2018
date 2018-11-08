function [status, link] = sqw_phonons_requirements(options)
% sqw_phonons_requirements: check for availability of ASE and MD codes
%
% MPI, EMT, GPAW, Abinit, Elk, QE, VASP
%
% returns a structure with a field for each MD software being 1 when available.
status = [];

% set-up a temporary directory for the tests
d = tempname;
p = pwd;
mkdir(d);

link.python = 'http://www.python.org';
link.mpirun = 'http://www.openmpi.org';
link.ase    = 'https://wiki.fysik.dtu.dk/ase';
link.emt    = link.ase;
link.phonopy= 'https://atztogo.github.io';
link.gpaw   = 'http://wiki.fysik.dtu.dk/gpaw';
link.elk    = 'http://elk.sourceforge.net';
link.abinit = 'http://www.abinit.org/';
link.quantumespresso = 'http://www.quantum-espresso.org/';
link.qeutil = 'https://jochym.github.io/qe-doc/';
link.vasp   = 'http://www.vasp.at/';
link.octopus= 'http://octopus-code.org/';
link.cp2k   = 'http://www.cp2k.org/';
link.siesta = 'https://departments.icmab.es/leem/siesta/';
link.hdf5storage='https://pythonhosted.org/hdf5storage/';

cd(d)
status = sqw_phonons_requirements_safe(link, options);
cd(p)
if isempty(status)
  return
end

for f=fieldnames(options)'
  if isfield(status, f{1})
    status.(f{1}) = options.(f{1});
    disp([ 'Using ' f{1} ' = ' status.(f{1}) ' (User definition)' ])
  end
end

cd(p);
rmdir(d, 's');

% ==============================================================================
function status = sqw_phonons_requirements_safe(link, options)

  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
  else           precmd=''; end

  disp('Available packages (system):');
  status = [];
  % test for python
  if isfield(options, 'python')
    status.python = options.python;
  else
    status.python = '';
    for calc={'python'}
      % now test executable
      [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
      if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
          status.python=calc{1};
          st = 0;
          disp([ '  Python          (' link.python ') as "' status.python '"' ]);
          break;
      end
    end
  end
  if isempty(status.python)
    disp([ mfilename ': ERROR: Python not installed. This is required.' ]);
    status = [];
    return
  end

  % test for ASE in Python
  [status.ase, result] = system([ precmd status.python ' -c "import ase"' ]);
  if status.ase ~= 0
    disp([ mfilename ': ERROR: requires ASE to be installed.' ])
    disp([ '  Get it at <' link.ase '>.' ]);
    disp('  Packages exist for Debian/Mint/Ubuntu, RedHat/Fedora/SuSE, MacOSX and Windows.');
    disp([ mfilename ': ASE not installed. This is required.' ]);
    status = [];
    return
  end
    
  disp([ mfilename ': using ASE ' result ]);
  status.emt='ase-run';
  status.ase=sscanf(result,'%d.%d');
  
  % test for hdf5storage
  [st, result] = system([ precmd status.python ' -c "import hdf5storage"' ]);
  if any(st == 0)
    status.hdf5storage = 'hdf5storage';
    disp([ '  hdf5storage     (' link.hdf5storage ') as "' status.hdf5storage '"' ]);
  else
    disp([ mfilename ': ERROR: requires python-hdf5storage to be installed.' ])
    disp('  Packages exist for Debian/Mint/Ubuntu, RedHat/Fedora/SuSE, Python PIP and conda.');
    status = [];
    return
  end
  
  disp('  EMT             only for Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O');
  
  % test for mpirun
  status.mpirun = '';
  for calc={'mpirun','mpiexec'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        status.mpirun=calc{1};
        st = 0;
        disp([ '  MPI             (' link.mpirun ') as "' status.mpirun '"' ]);
        break;
    end
  end
  
  % test for PhonoPy
  [st, result] = system([ precmd status.python ' -c "from phonopy import Phonopy"' ]);
  if any(st == 0)
    status.phonopy = 'phonopy';
    disp([ '  PhonoPy         (' link.phonopy ') as "' status.phonopy '"' ]);
  else
    status.phonopy = '';
  end
  
  % test for GPAW
  [st, result] = system([ precmd status.python ' -c "from gpaw import GPAW"' ]);
  status.gpaw='';
  if any(st == 0:2)
    % now test executable
    for calc={'gpaw','gpaw-python'}
      [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
      if (st == 0 || st == 2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
          status.gpaw=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.gpaw)
    disp([ '  GPAW            (' link.gpaw ') as "' status.gpaw '"' ]);
  end
  
  % test for Elk
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.elk import ELK"' ]);
  status.elk='';
  if any(st == 0:2)
    % now test executable
    for calc={'elk','elk-lapw'}
      [st,result]=system([ precmd calc{1} ]);
      if (st == 0 || st == 2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
          status.elk=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.elk)
    disp([ '  Elk             (' link.elk ') as "' status.elk '"' ]);
  end
  
  % test for ABINIT
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.abinit import Abinit"' ]);
  status.abinit='';
  if any(st == 0:2)
    for calc={'abinit','abinis','abinip'}
      % now test executable
      [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
      if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
          status.abinit=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.abinit)
    disp([ '  ABINIT          (' link.abinit ') as "' status.abinit '"' ]);
  end
  
  % test for QuantumEspresso
  status.quantumespresso = '';
  for calc={'pw.x','pw.exe','pw','pwscf'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        status.quantumespresso=calc{1};
        st = 0;
        disp([ '  QuantumEspresso (' link.quantumespresso ') as "' status.quantumespresso '"' ]);
        break;
    end
  end
  
  % test for QE/ASE from QEutil
  status.qeutil = '';
  [st, result] = system([ precmd status.python ' -c "from qeutil import QuantumEspresso"' ]);
  if any(st == 0:2)
    status.qeutil='qeutil';
    disp([ '  QEutil          (' link.qeutil ') as "' status.qeutil '"' ]);
  else
    status.qeutil='';
  end
  
  % test for QE/ASE native (ASE >= 3.15)
  status.qease = '';
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.espresso import Espresso"' ]);
  if any(st == 0)
    status.qease='ase.calculators.espresso';
  else
    status.qease='';
  end
  
  % test for VASP
  [st, result] = system([ precmd 'vasp' ]);
  if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
    status.vasp = 'vasp';
  else
    status.vasp = '';
  end
  if ~isempty(status.vasp)
    disp([ '  VASP            (' link.vasp ') as "' status.vasp '"' ]);
  end
  % must set:
  % VASP_COMMAND=vasp
  % VASP_PP_PATH=/usr/share/vasp/pseudo
  
  % test for Octopus
  status.octopus = '';
  for calc={'octopus','octopus_mpi'}  % octopus_mpi is obsolete, 2nd choice.
    % now test executable
    [st,result]=system([ precmd calc{1} ' -v' ]);
    if any(st == 0:2) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        status.octopus=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.octopus)
    disp([ '  Octopus         (' link.octopus ') as "' status.octopus '"' ]);
  end
  
  % test for CP2K
  status.cp2k = '';
  for calc={'cp2k_shell','cp2k_shell.popt'}
    % now test executable
    [st,result]=system([ precmd 'echo EXIT | ' calc{1} ]);
    if any(st == 0) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        status.cp2k=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.cp2k)
    disp([ '  CP2K            (' link.cp2k ') as "' status.cp2k '"' ]);
  end
  
  % test for SIESTA
  status.siesta = '';
  for calc={'siesta'}
    % now test executable
    [st,result]=system([ precmd 'echo | ' calc{1} ]);
    if any(st == 0:1) && (~ispc || isempty(strfind(result, [ '''' calc{1} '''' ])))
        status.siesta=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.siesta)
    disp([ '  SIESTA          (' link.siesta ') as "' status.siesta '"' ]);
  end

  % disp('Calculator executables can be specified as ''options.command=exe'' when building a model.');

  %  lj (lenard-jones)
  %  morse
  %  eam

