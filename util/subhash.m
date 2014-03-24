function sh=subhash(h,keys)
%  sh=subhash(h,keys)
% This function returns a hash (either a matlab 'struct' object, or
% a 2dimensional cell array) that consists of a subset of the 
% key value pairs of h. The subset is specified by the cell array of
% strings 'keys'.
% 
% See also keys, getvalue, putvalue
% Kjartan Halvorsen
% 2000-11-04

if (isempty(keys))
   sh={};
   return
end

if (ischar(keys)) % convert to cell
   keys={keys};
end

if (isa(h,'struct'))
   hkeys=fieldnames(h);
   [check,ind,hind]=intersect(keys,hkeys);
   try
      [ind,ind2]=sort(ind);
      check=keys(ind);
      hind=hind(ind2);
   catch
      sh={};
      return
   end

   if (~(length(check)==length(keys)))
      warning('Not all keys found in hash h.');
   end

   structarg=cell(2,length(check));
   for k=1:length(check)
      structarg{1,k}=check{k};
      structarg{2,k}=getfield(h,hind(k));
   end

   sh=struct(structarg{:});

elseif (isa(h,'cell'))
   hkeys=h(:,1);

   check=intersect(keys,hkeys);
   if (~(length(check)==length(keys)))
      warning('Not all keys found in hash h.');
   end

   [slask,keyi,hkeyi]=intersect(keys,hkeys);
   try
      [keyi,hind]=sort(keyi);
      check=keys(keyi);
      hkeyi=hkeyi(hind);
   catch
      sh={};
      return
   end

   sh=cat(2,check(:),h(hkeyi,2));
end
