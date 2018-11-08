function mifit_Optimizers_Criteria(varargin)
% Optimizers/Criteria: open dialogue to change the fit criteria
  
  crit = {'least_square = (|Signal-Model|/Error)^2 [non-robust]', ...
    'least_absolute = |Signal-Model|/Error [robust]',...
    'least_median = median(|Signal-Model|/Error) [robust]',...
    'least_mean = mean(|Signal-Model|/Error) [robust]',...
    'least_max = max(|Signal-Model|/Error) [non-robust]',...
    'least_rfactor = (|Signal-Model|/Error)^2/(Signal/Error)^2 [non-robust]',...
    'max_likelihood =  (|Signal-Model|/σ)^2 + 1/2 log(2πσ) [robust]', ...
    'max_corrcoef = 1-corrcoef(Signal, Model)'};
  selection = [];
  if isappdata(mifit_fig,'CurrentOptimizerCriteria')
    selection = find(strcmp(getappdata(mifit_fig,'CurrentOptimizerCriteria'), strtok(crit)),1);
  end
  if isempty(selection), selection = 1; end
  % show operator selection dialogue
  [selection, ok] = listdlg('ListString', crit, 'ListSize',[400 160], ...
    'Name',[ 'miFit: Select the Fit Criteria ti use' ], ...
      'PromptString', 'The Fit criteria is the quantity which is to be minimized. It is a function of the Data set Signal, and the Model value. It is then divided by the number of degrees of freedom.', ...
      'InitialValue',selection, 'SelectionMode','single');
  if isempty(selection) || ok ~= 1, return; end
  crit = strtok(crit{selection});
  setappdata(mifit_fig,'CurrentOptimizerCriteria', crit);
