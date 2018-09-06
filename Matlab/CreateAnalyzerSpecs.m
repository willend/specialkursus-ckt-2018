function [] = CreateAnalyserSpec(filename)
  
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
        [A, ang, Acc] = AsymmetricRowland_NewGeometry_JustCalcs(En, Width, L1, L2, delta, Nbla, 1);
        Area = Area + A;
        count = count + 1;
        Angs(count) = ang;
        VertCov(count) = Acc;

       end
    display(['Total analyzer area is: ' num2str(((3*Area)), 4) 'm^2'])
    display(['Maximum angular coverage is: ' num2str(max(Angs), 4) 'degrees'])
    Angs
    
    
    VertCov
    
    
    