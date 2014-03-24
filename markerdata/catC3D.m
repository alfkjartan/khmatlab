function catC3D(infile1, infile2, outfile) 
%  combineC3D(infile1, infile2, outfile) 
% Combines the data in infile1 and infile2 and writes a new c3d
% file in outfile.

% Kjartan Halvorsen
% 2007-03-14

[mrks1, vidFrRt1, analogs1, analogFrRt1, events1, parameterGr1,...
 caminfo1,residerr1] = readC3D(infile1);
[mrks2, vidFrRt, analogs2, analogFrRt, events2, parameterGr2,...
 caminfo2,residerr2] = readC3D(infile2);

% Change only number of frames in Point parameter group.

pg = getc3dparam(parameterGr1, 'POINT', 'FRAME');
pg.data = size(mrks1,1) + size(mrks2,1);

parameterGr = setc3dparam(parameterGr1, 'POINT', pg);

% Check if analogs2 is empty. create zeros as needed
if isempty(analogs2)
  analogs2 = zeros(size(mrks2, 1)*analogFrRt1/vidFrRt1, ...
		   size(analogs1, 2));
end


% Write the file


writeC3D(cat(1, mrks1, mrks2), vidFrRt1, cat(1, analogs1, analogs2), ...
	 analogFrRt1, events1, parameterGr,...
	 cat(1, caminfo1, caminfo2), ...
	 cat(1, residerr1, residerr2),outfile);

