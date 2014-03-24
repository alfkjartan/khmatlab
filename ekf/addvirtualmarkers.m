function nkm = addvirtualmarkers(km)
%  nkm = addvirtualmarkers(km)
% 
% Adds virtual markers to the km. Actually, the proximal and distal
% joint centers are added to the track markers.
%
 
% Kjartan Halvorsen
% 2004-03-17
  
nkm=km;

nkm.p0 = combinep0jc(nkm.p0, nkm.jcs, '');

  
  
  
function np0 = combinep0jc(p0,jc,prefix)
  np0=p0;
  if ~isempty(jc)
    if ~isempty(jc{1})
      np0{1} = putvalue(np0{1},[prefix,'prox_jc'], jc{1});
    end
    for brjc = 2:length(jc)
      if (~isempty(jc{brjc}) & ~isempty(jc{brjc}{1}))
	np0{1} = putvalue(np0{1},[prefix,'dist_jc', ...
		    int2str(brjc-1)], jc{2}{1});
      end
    end
  end
  
  for br=2:length(np0);
    np0{br} = combinep0jc(np0{br},jc{br},[prefix,'branch']);
  end

  
  
  