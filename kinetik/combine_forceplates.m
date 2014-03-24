function [nanalog, npg] = combine_forceplates(analog, pg)
%  [nanalog, npg] = combine_forceplates(analog, pg)
%  Combines all force plate data in a c3d file to a single force
%  plate. Works so far only for type 1 platforms
%
%  Input
%     analog         ->  analog data from c3d file (nfrs x
%                        nchannels)
%     pg             ->  parameter group from c3d file
%  Output
%     nanalog        <-  new set of analog data
%                        ( nfrs x (nchannels + 6) )    
%     npg             <-  new parameter group

% Kjartan Halvorsen
% 2005-03-10



