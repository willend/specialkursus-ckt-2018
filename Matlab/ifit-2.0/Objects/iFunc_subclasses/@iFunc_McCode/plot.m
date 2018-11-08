function [comps, fig, model]=plot(model, p, options, match)
% iFunc_McCode/plot: runs the model in --trace mode and capture the output
% grab all MCDISPLAY lines, and render the TRACE information into a figure.
%
%  plot(model)
%     displays the given McCode model with its defaults/current parameters.
%  plot(model, parameters)
%     same as above, but uses the given parameters (vector, structure, cell, string)
%  plot(model, parameters, options)
%     same as above, and specifies an output file format to generate
%  plot(model, parameters, options, match)
%     same as above, and only plot component names that match 'match'
%  plot(iFunc_McCode)
%     displays the default instrument geometry (templateDIFF)
%  [comps, fig] = plot(...)
%     returns the list of components, figure used for display and updated model.
%     
% Example:
%   model = mccode('templateDIFF')
%   plot(model)
%     a figure is generated, which shows the instrument geometry.
%
% syntax:
%   [comps, fig, model]=plot(model, p, options, match)
%
% input:
%   model:      instrument simulation as obtained with 'mccode' (iFunc)
%               or path to an instrument definition to compile.
%   parameters: parameters to use, as a vector, cell, structure... or empty for default
%   options:    rendering/export options. Can contain 
%     'html','x3d','png','pdf','fig','tif','jpg','eps', 'svg' or empty for no output
%   match:      token to only render components that match it, or empty for all.
%   
% output:
%   comps: a list of component specifications (structure array)
%   fig:   figure handle
%   model: model, updated or created (iFunc)
%
% See also: mccode, iFunc, iFunc/feval, http://www.mcstas.org

  if nargin == 0 || isempty(model)
    model = mccode('defaults');
  end
  if nargin < 2, p=[]; end
  if nargin < 3, options=''; end
  if nargin < 4, match  =''; end
  comps = []; fig = [];

  % check if we have an iFunc McCode object
  if ischar(model) && ~isempty(dir(model))
    model = mccode(model);
  end
  
  if ~isa(model, 'iFunc')
    error([ mfilename ': ERROR: Usage: plot(model) with model=path_to_mccode_instr or mccode(path_to_mccode_instr)' ]);
  end
  if ~isfield(model.UserData, 'options') || ~isfield(model.UserData,'instrument_executable')
    error([ mfilename ': ERROR: Usage: plot(model) with model=mccode(path_to_mccode_instr)' ]);
  end
  
  % read monitors as iData
  if model.UserData.options.ncount > 1e3
    monitors       = model.UserData.monitors;
  else
    monitors       = [];
  end
  val = [];
  
  if isfield(model.UserData, 'display_cache')
    % to use cache, the same parameters must have been used
    flag = true;
    if ~isfield(model.UserData.display_cache, 'ParameterValues') ...
      || ~isequal(model.UserData.display_cache.ParameterValues, model.ParameterValues)
      flag = false;
    end
    if isfield(model.UserData, 'Parameters_Constant')
      if ~isfield(model.UserData.display_cache, 'Parameters_Constant') ...
      || ~isequal(model.UserData.display_cache.Parameters_Constant, model.UserData.Parameters_Constant)
        flag = false;
      end
    end
    
    if flag && isfield(model.UserData.display_cache, 'comps')
      disp([ mfilename ': re-using instrument ' strtok(model.Name) ' Trace from cache...' ])
      comps = model.UserData.display_cache.comps;
    end
  end
  
  if isempty(comps)
    % switch to trace mode when not in Cache
    ncount = model.UserData.options.ncount;
    raw    = model.UserData.options.raw;
    model.UserData.options.trace  = 1;
    model.UserData.options.ncount = 1e3;  % create a new set of output files with reduced statistics
    model.UserData.options.raw    = '--no-output-files';  % much faster without data import
    % execute and capture output (TRACE)
    disp([ mfilename ': running instrument ' strtok(model.Name) ' in Trace mode...' ])
    output = evalc('[val,model]=feval(model,[],nan);');
    model.UserData.options.trace = 0;
    model.UserData.options.ncount= ncount;
    model.UserData.options.raw   =raw;
  end
  
  if isempty(monitors)
    monitors       = model.UserData.monitors;
  end
  try
    monitors_names = get(monitors,'Component');
  catch
    monitors_names = [];
  end

  % first extract the portion 'start':'end'
  disp([ mfilename ': rendering geometry...' ])
  
  if isempty(comps)
    index_start = strfind(output, 'MCDISPLAY: start');
    index_end   = strfind(output, 'MCDISPLAY: end');
    if numel(index_start) ~= 1 || numel(index_end) ~= 1
      disp([ mfilename ': The MCDISPLAY section is invalid (incomplete or multiple). Aborting.' ]);
      comps=output; fig=[];
      return
    end

    output_mcdisplay_init = output(1:index_start);
    output_mcdisplay_init = textscan(output_mcdisplay_init, '%s','Delimiter',sprintf('\n\r'));
    output_mcdisplay_init = output_mcdisplay_init{1};

    % initiate component structures
    % we build a struct array, with one entry per component, and fields:
    %   name: name
    %   pos: pos(3)
    %   rot: rot(3,3)
    %   x,y,z: set of points
    comps = mcdisplay_get_components(output_mcdisplay_init);
    clear output_mcdisplay_init
    
    % restrict output to the MCDISPLAY so that all searches are faster
    output_mcdisplay_section = output(index_start:index_end);
    clear output
    output_mcdisplay_section = textscan(output_mcdisplay_section, '%s','Delimiter','\n\r');
    output_mcdisplay_section = output_mcdisplay_section{1};

    % get the components in order, and identify the output section/lines
    % which are separated by e.g. MCDISPLAY: component <blah>
    index_mcdisplay_comp = find(~cellfun(@isempty, strfind(output_mcdisplay_section, 'MCDISPLAY: component ')));
    if numel(index_mcdisplay_comp) ~= numel(comps)
      disp([ mfilename ...
        ': WARNING: not the same number of declared components (' num2str(numel(comps)) ...
        ') and MCDISPLAY sections ' num2str(numel(index_mcdisplay_comp)) ])
    end

    % extract the multiline and circle stuff in each component mcdisplay section
    for index=1:numel(index_mcdisplay_comp)
      if index < numel(index_mcdisplay_comp), 
        next = index_mcdisplay_comp(index+1);
      else 
        next = numel(output_mcdisplay_section); end
      % get the MCDISPLAY section for a single component
      section = output_mcdisplay_section(index_mcdisplay_comp(index):next);
      % then we get the multiline and circle commands in this section
      for token = {'multiline' ,'circle'}
        [x,y,z] = mcdisplay_get_token(section, token{1});
        comps(index).x = [ comps(index).x nan x ];
        comps(index).y = [ comps(index).y nan y ];
        comps(index).z = [ comps(index).z nan z ];
      end
    end
    clear output_mcdisplay_section
  end

  % PLOTTING: transform the points and plot them
  fig = gcf; 
  if ischar(match) && ~isempty(match), match = cellstr(match); end
  if ~iscellstr(match), match = {}; end
  if ~isempty(match)
    set(fig, 'Name',[ 'Instrument: ' model.Name ' ' sprintf('%s ', match{:}) ]);
  else
    set(fig, 'Name',[ 'Instrument: ' model.Name ]);
  end
  colors='bgrcmk';
  for index=1:numel(comps)
    comp = comps(index);

    if ~isempty(match)
      found = false;
      for m=1:numel(match)
        match{m} = strtrim(match{m});
        if isempty(match{m}) || ~isempty(strfind(comp.name, match{m}))
          found = true;
          break
        end
      end
    else found = true;
    end
    if ~found, continue; end % not found match -> skip this comp
    
    r = [ comp.x ; comp.y ; comp.z ];
    if all(isnan(r)), continue; end
    R = comp.rot*r;
    x = R(1,:)+comp.pos(1);
    y = R(2,:)+comp.pos(2);
    z = R(3,:)+comp.pos(3);
    disp([' Component: ' comp.name ' [' num2str(index) '] AT (' sprintf('%g,%g,%g', comp.pos) ')'  ])
    c = mod(comp.index, numel(colors)); c=colors(c+1);
    h = plot3(z,x,y, [ c '-' ], 'DisplayName', comp.name);
    popup=uicontextmenu;
    uimenu(popup,'label', comp.name,'ForeGroundColor',c);
    uimenu(popup,'label', [ 'AT: ' mat2str(comp.pos) ]);
    % attach an option to view a monitor data set when available
    index_mon = find(strcmp(monitors_names, comp.name));
    if ~isempty(index_mon)
      % add 'values' field
      if numel(monitors) > 1, mon=monitors(index_mon); else mon=monitors; end
      uimenu(popup,'label', [ '[I err N]: ' mat2str(get(mon(1),'values')) ]);
      mon = set(mon, 'Title', comp.name);
      % add uicontextmenu element to plot data set(s) with subplot(UserData)
      % store UserData so that we can plot it
      uimenu(popup,'label', 'Plot monitor...', ...
        'Callback','figure; subplot(get(gcbo,''UserData''),''tight'');', ...
        'UserData', mon);
      
      % outline the monitor components
      set(h,'LineWidth',4); 
    end
    set(h,'uicontextmenu',popup);
    hold on
  end
  % plot a red dashed line that follows the centre location of components
  centers = [ 0 0 0 ];
  for index=1:numel(comps)
    comp = comps(index);
    centers = [ centers ; comp.pos(:)' ];
  end
  plot3(centers(:,3), centers(:,1), centers(:,2), 'r:', 'DisplayName', model.Name);
  
  xlabel('Z [m]');
  ylabel('X [m]');
  zlabel('Y [m]');
  daspect([1 1 1]);
  box on;
  a0 = gca;
  
  mp    = model.ParameterValues;
  names = model.Parameters;
  t = [];
  for index=1:numel(mp)
    if numel(mp) < index, valp = []; else valp = mp(index); end
    t = [ t names{index} '=' num2str(valp) ' ' ];
  end
  % add static parameters
  if isfield(model.UserData, 'Parameters_Constant') && ...
     ~isempty(model.UserData.Parameters_Constant) && ...
     ~isempty(fieldnames(model.UserData.Parameters_Constant))
    for f=fieldnames(model.UserData.Parameters_Constant)
      name = f{1}; valp=model.UserData.Parameters_Constant.(name);
      t = [ t name '=' num2str(valp) ' ' ];
    end
  end
  
  t = { sprintf([ 'Instrument: ' model.Name ' ' sprintf('%s ', match{:}) '\n']) ; t };
  t = textwrap(t, 80);
  t = sprintf('%s\n', t{:});
  title(t ,'Interpreter','None');
  model.UserData.title = t;
  
  % store Display info in cache
  model.UserData.display_cache.ParameterValues = mp;
  if isfield(model.UserData, 'Parameters_Constant')
    model.UserData.display_cache.Parameters_Constant = model.UserData.Parameters_Constant;
  end
  model.UserData.display_cache.comps = comps;

  plot_contextmenu(gca, model.Name, t, monitors);
  [~,filename] = fileparts(strtok(model.Name));
  if isempty(filename), filename = 'instrument'; end
  if isdir(model.UserData.options.dir)
    filename = fullfile(model.UserData.options.dir, filename); 
  end
  if ~isempty(match)
    filename = [ filename sprintf('%s', match{:}) ];
  end
  
  % export options to x3d/xhtml
  if exist('figure2xhtml')
    if ~isempty(strfind(options, 'html'))
      t(t=='<')='[';
      t(t=='>')=']';
      figure2xhtml(filename, fig, ...
        struct('title', model.Name, 'Description',t,'interactive',true));
      plot_exportmessage([ filename '.xhtml' ])
      plot_exportmessage([ filename '.x3d' ])
    end
    if ~isempty(strfind(options, 'x3d'))
      t(t=='<')='[';
      t(t=='>')=']';
      figure2xhtml(filename, fig, struct('interactive',true, ...
        'output', 'x3d','title',model.Name,'Description',t));
      plot_exportmessage([ filename '.x3d' ])
    end
  end
  if exist('plot2svg') && ~isempty(strfind(options, 'svg'))
    plot2svg([ filename '.svg' ], fig);
    plot_exportmessage([ filename '.svg' ])
  end
  
  % add the model value in a small insert (statistics will be limited as ncount=1e3)
  % we get the monitor values (which may be with higher statistics)
  if ~isempty(model.UserData.options.monitor) && isnumeric(model.UserData.options.monitor) && numel(monitors) > 1
    monitors=monitors(model.UserData.options.monitor);
  end
  if isa(monitors, 'iData')
    if numel(monitors) > 1, monitors = monitors(end); end
    val=double(monitors);
  end
  if ~isempty(val)
    if isvector(val)
      a=axes('Parent',fig,'position',[.8 .05 .15 .15]);
      h=plot(a, val); axis tight
      set(a,'XTick',[]);
      if ~isempty(model.UserData.options.monitor)
        title(model.UserData.options.monitor,'Interpreter','none');
      end
    elseif ndims(val) == 2
      a=axes('Parent',fig,'position',[.8 .05 .15 .15]);
      h=surf(a, val); view(2); axis tight; set(a,'XTick',[], 'YTick',[]);
      set(h,'EdgeColor','none')
      if ~isempty(model.UserData.options.monitor)
        title(model.UserData.options.monitor,'Interpreter','none');
      end
    end
    clear val
  end
  
  % export to static images
  for f={'png','pdf','fig','tif','jpg','eps'}
    if ~isempty(strfind(options, f{1}))
      try
        saveas(fig, [ filename '.' f{1} ], f{1});
        plot_exportmessage([ filename '.' f{1} ])
      catch ME
        disp(getReport(ME))
      end
    end
  end
  axes(a0);
  
  if ~isempty(inputname(1)) && nargout < 3
    assignin('caller',inputname(1),model); % update in original object
  end

end % plot





% ------------------------------------------------------------------------------
% initialize the component structures by searching name and pos/rot, e.g.
%   COMPONENT: "collimador_radial"
%   POS: 1.82089, 0, 19.6314, 0.822317, 0, -0.569029, -0, 1, 0, 0.569029, -0, 0.822317
function comps = mcdisplay_get_components(output)
  token = 'COMPONENT: "';
  index_token = find(~cellfun(@isempty, strfind(output, token)));
  comps = [];
  for index = index_token'
    this_line = output{index};
    compname  = strtok(this_line(numel(token):end),'"');
    pos = [];
    if index < numel(output)
      next_line = output{index+1};
      if strncmp(next_line, 'POS: ', 5), pos = str2num(next_line(6:end)); end
    end

    if ~isempty(compname) && numel(pos) == 12
      comp.name = compname;
      comp.pos  = pos(1:3);                 % absolute position
      comp.rot  = reshape(pos(4:end),3,3);  % absolute rotation matrix
      comp.index= numel(comps)+1;
      comp.x=[]; comp.y=[]; comp.z=[];
      comps = [ comps comp ];
    end
  end
end

% search for a plot command: multiline or circle. Get the local coords back.
function [X,Y,Z] = mcdisplay_get_token(output, token)
  index_token = find(~cellfun(@isempty, strfind(output, [ 'MCDISPLAY: ' token ])));
  X=[]; Y=[]; Z=[];
  if isempty(index_token), return; end
  for index = index_token(:)'
    % in each multiline, we replace the search token and execute the remaining part in eval
    this_line = output{index};
    this_line = strrep(this_line, 'MCDISPLAY: ','[x,y,z]=');
    eval([ this_line ';' ]);  % executes multiline(..), return a set of points (local coords)
    X = [ X x ];
    Y = [ Y y ];
    Z = [ Z z ];
  end
end

% get the multiline and circle arguments back
function [X,Y,Z]=multiline(npoints, varargin)
  X=[]; Y=[]; Z=[];
  for index=1:3:numel(varargin)
    X=[ X varargin{index} ];
    Y=[ Y varargin{index+1} ];
    Z=[ Z varargin{index+2} ];
  end
end

function [X,Y,Z]=circle(plane, x0,y0,z0, radius)
  % we create a set of points
  if radius ~=0
    phi=linspace(0,2*pi, 36); % 36 points along the circle
  else phi=0; end
  x=radius*sin(phi);
  y=radius*cos(phi);
  zero = 0*x;
  switch plane
  case {'xy','yx'}
    X=x; Y=y; Z=zero;
  case {'xz','zx'}
    X=x; Y=zero; Z=y;
  case {'zy','yz'}
    X=zero; Y=x; Z=y;
  otherwise
    X=[]; Y=[]; Z=[]; 
    disp([ mfilename ': unknown plane: ' plane ]);
    return
  end
  X=X+x0; Y=Y+y0; Z=Z+z0;
end

function plot_contextmenu(a, name, pars, monitors)
  % install a context menu on the axis object of the figure
  
  % build the contextual menu
  if iscell(pars), pars = sprintf('%s ', pars{:}); end
  uicm = uicontextmenu('Tag','plot_contextmenu_gca');
  uimenu(uicm, 'Label', [ 'About ' name '...' ], ...
           'Callback', [ 'helpdlg(''' strrep(pars, sprintf('\n'),' ') ''',''' name ''')' ]);
  if ~isempty(monitors)
    uimenu(uicm, 'Label', 'Plot monitors...', ...
      'UserData',monitors, ...
      'Callback','figure; subplot(get(gcbo,''UserData''),''tight'');');
  end
    
  uimenu(uicm, 'Separator','on', 'Label', 'Duplicate View...', 'Callback', ...
         [ 'tmp_cb.g=gca;' ...
           'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
           'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
           'set(gcf,''Name'',''Copy of McCode display: ' name '''); ' ...
           'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
           'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
           'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');']);
           
  uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
  uimenu(uicm, 'Label','Toggle aspect ratio','Callback','if all(daspect == 1) daspect(''auto''); else daspect([ 1 1 1 ]); end');
  uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
  uimenu(uicm, 'Label','Toggle legend','Callback','tmp_h=legend(''toggle''); set(tmp_h,''Interpreter'',''None''); if strcmp(get(tmp_h,''Visible''),''off''), legend(gca,''off''); end; clear tmp_h;');
  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
      '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
      'clear tmp_a tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;legend off;' ]);
  
  uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', 'Callback', ...
    [ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
    
  % attach the contextual menu
  set(a, 'UIContextMenu', uicm);
  
  % ==============================================================================
  % contextual menu for the figure
  % add rotate/pan/zoom tools to the figure in case java machine is not started
  if ~usejava('jvm')
    if isempty(get(gcf, 'UIContextMenu'))
      uicmf = uicontextmenu('Tag','iData_plot_contextmenu_fig');
      uimenu(uicmf, 'Label','Zoom on/off', 'Callback','zoom');
      uimenu(uicmf, 'Label','Pan on/off',  'Callback','pan');
      if ndims(a) >= 2
        uimenu(uicmf, 'Label', 'Rotate on/off', 'Callback','rotate3d');
      end
      uimenu(uicmf, 'Label','Legend on/off', 'Callback','legend(gca, ''toggle'',''Location'',''Best'');');
      uimenu(uicmf, 'Label','Print...', 'Callback','printpreview');

      set(gcf, 'UIContextMenu', uicmf);
      set(gcf, 'KeyPressFcn', @(src,evnt) eval('if lower(evnt.Character)==''r'', lighting none;alpha(1);shading flat;axis tight;rotate3d off; zoom off; pan off; end') );
    end
  end
    
end

function plot_exportmessage(filename)
  if ~isdeployed
    disp([ mfilename ': exported instrument view as: <a href="' filename '">' filename '</a>'])
  else
    disp([ mfilename ': exported instrument view as: ' filename ])
  end
end
