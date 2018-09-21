function [] = CreateAnalyserSpec_2(filename)
  
    if nargin==0,
        display('Filenumber not given')
        return
    else
        fid=fopen(char(filename),'r');
    end

    while 1
           fileline = fgetl(fid);
           if ~isstr(fileline),break, end
           if findstr('%Column names',fileline),break,end
    end
    
      fileline = fgetl(fid)
      count = 1;
      Area = 0;
      Angs(1) = 0;
      VertCov = [];
      Lsa = [];
      Lad = [];
      BladeL = [];
      DetLen = [];
      Energy = [];
       while 1
        fileline = fgetl(fid);
        if findstr('%NoMoreAnalyzers',fileline),break,end
        %datarray = regexp(fileline, ';', 'split');
    
        C = textscan(fileline,'%f')
        En = C{1}(1);
        Width = C{1}(2);
        L1 = C{1}(3);
        L2 = C{1}(4);
        delta = C{1}(5);
        Nbla = C{1}(6);           
        [A, ang, Acc, BladeL(count), DetLen(count)] = AsymmetricRowland_NewGeometry_JustCalcs_2(En, Width, L1, L2, delta, Nbla, 1);
        Area = Area + A;
        Angs(count) = ang;
        VertCov(count) = Acc;
        Lsa(count) = L1;
        Lad(count) = L2;
        Energy(count) = En;
        count = count + 1;

       end
       
    fid = fopen(['Bifrost_Analyzer_Summary.txt'], 'wt');
	h = length(Lad);
    fprintf(fid, 'Lsa [mm]      Lad [mm]        Blade length [mm]       Detector length [mm]\n');	
	for n = 1:h    
	fprintf(fid, '%12.3f \t %12.3f \t %12.3f \t %12.3f \n', [Lsa(n) Lad(n)  BladeL(n) DetLen(n)]);	
    end
    fclose(fid)
    edit 