function vw=getCoordinates(tw)
% function vw=getCoordinates(tw)
% 

% Kjartan Halvorsen
% 1999-09-21


vw=cat(1, tw(1:3,4), vect(tw(1:3, 1:3)));