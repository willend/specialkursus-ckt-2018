function [s,lambda,distance,chwidth,energy,wavevector] = Sqw_search_lambda(s)
  % search for the wavelength etc in the object parameters
  
  lambda = []; distance = []; chwidth=[]; energy=[]; wavevector=[];
  if ~isfield(s, 'parameters')
    parameters = [];
  else
    parameters = get(s, 'parameters');
  end
  
  % get lambda, energy, wavevector and check for distance and channel width ----
  % first solution:  search in 'parameters'
  % second solution: findfield in the object
  if isfield(parameters, 'Wavelength') && parameters.Wavelength   
    lambda     = parameters.Wavelength; 
  else
    lambda     = Sqw_getT(s, {'IncidentWavelength', 'wavelength' 'lambda'});
  end
  if isfield(parameters, 'IncidentEnergy') && parameters.IncidentEnergy
    energy     = parameters.IncidentEnergy; 
  else
    energy     = Sqw_getT(s, {'IncidentEnergy' 'fixed_energy' 'energy' 'ei'});
  end
  if isfield(parameters, 'IncidentWavevector') && parameters.IncidentWavevector
    wavevector = parameters.IncidentWavevector; 
  else
    wavevector = Sqw_getT(s, {'IncidentWavevector','wavevector' 'ki'});
  end
  if isfield(parameters,'Distance') && parameters.Distance
    distance = parameters.Distance;
  else
    distance = Sqw_getT(s, {'Distance_Det_Sample','detector_distance', 'distance'});
  end
  if isfield(parameters, 'ChannelWidth') && parameters.ChannelWidth
    chwidth = parameters.ChannelWidth;
  else
    chwidth = Sqw_getT(s, {'ChannelWidth', 'Channel_width'});
  end
  
  % now we test if the retrieved value are OK
  if isempty(lambda)
    if     ~isempty(energy)
      if ischar(energy), energy = get(s, energy); end
      if ~isempty(energy) && energy > 0
        lambda = sqrt(81.805./energy);
      end
    elseif ~isempty(wavevector)
      if ischar(wavevector), wavevector = get(s, wavevector); end
      if ~isempty(wavevector) && wavevector > 0
        lambda = 2*pi./wavevector; 
      end
    end
  end
  
  if isempty(lambda)
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' undefined incident neutron wavelength/wavevector/energy.' ]);
    disp('    Define e.g. s.Wavelength=<lambda in Angs> or s.IncidentEnergy=<energy in meV>');
  else
    lambda = mean(lambda(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <wavelength>               =' num2str(lambda) ' [Angs]']);
    if ~isfield(s, 'IncidentWavelength')
      setalias(s, 'IncidentWavelength', lambda, 'Incident neutron Wavelength [Angs-1]');
    end
  end
  
  if ~isempty(lambda) && lambda > 0 && (isempty(energy) || ischar(energy))
    energy      = 81.805./lambda^2; 
    if ~isfield(s, 'IncidentEnergy')
      setalias(s, 'IncidentEnergy', energy, 'Incident neutron Energy [meV]');
    end
  end
  if ~isempty(lambda) && lambda > 0 && (isempty(wavevector)  || ischar(wavevector))
    wavevector  = 2*pi./lambda;
    if ~isfield(s, 'IncidentWavevector')
      setalias(s, 'IncidentWavevector', wavevector, 'Incident neutron Wavevector [Angs-1]');
    end
  end
  
  % search for a sample-to-detector distance
  if ~isempty(distance)
    distance = mean(distance(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <sample-detector distance> =' num2str(distance) ' [m]' ]);
    if ~isfield(s, 'Distance')
      setalias(s, 'Distance', distance, '[m] Sample-Detector distance'); 
    end
  end
  
  % search for the Channel Width
  if ~isempty(chwidth)
    chwidth = mean(chwidth);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <channel width>            =' num2str(chwidth) ' [time unit: s, ms or us]']);
    if ~isfield(s, 'ChannelWidth')
      setalias(s, 'ChannelWidth', chwidth, '[time unit] ToF Channel Width');
    end
  end
 
