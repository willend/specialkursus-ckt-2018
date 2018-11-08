function filename = iData_private_saveas_m(a, filename, name, options)
  NL = sprintf('\n');
  if ~isdeployed
    str = [ 'function this=' name NL ];
  else
    str = '';
  end
  str = [ str ...
          '% Original data: ' NL ...
          '%   class:    ' class(a) NL ...
          '%   variable: ' inputname(1) NL ...
          '%   tag:      ' a.Tag NL ...
          '%   label:    ' a.Label NL ...
          '%   source:   ' a.Source NL ... 
          '%' NL ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) ' with iData/saveas' NL ...
          '% To use/import data, type ''' name ''' at the matlab prompt.' NL ...
          '% You will obtain an iData object (if you have iData installed) or a structure.' NL ...
          '%' NL ...
          class2str('this', a, options) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
    disp(message)
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
  if isdeployed
    disp([ 'Warning: The standalone/deployed version of iFit does not allow to read back' NL ...
           '  function definitions. This m-file has been converted to a script that you can' NL ...
           '  import as "this" by typing: run ' filename ]);
  end
