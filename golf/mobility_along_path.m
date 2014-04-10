function [Wpath, Wnormal, evmagn] = mobility_along_path(Wep, endpoint)
%% [Wpath, Wnorm] = mobility_along_path(Wep, endpoint)
%% Will return the mobility in the direction of movement of the end point, and the mobility
%% in the normal plane along the path of the end point.
%%
%% Input
%%    Wep            ->  Mobility tensor (3x3xNfrs) of the end point in spatial coordinates
%%    endpoint       ->  End point trajectory (Nfrsx3)

%% Kjartan Halvorsen
%% 2014-01-31



%% Find vectors in direction of end point velocity, using central difference
ev = endpoint(3:end,:) - endpoint(1:end-2,:);
evmagn = sqrt(sum(ev.^2, 2));

for i=1:3
    ev(:,i) = ev(:,i) ./ evmagn;
end

%% Repeat first and last element to get a complete record
ev = cat(1, ev(1,:), ev, ev(end,:));

%% Loop for each time frame (better, vectorized ways?)

Nfrs = size(Wep, 3);

Wpath = zeros(Nfrs,1);
Wnormal = zeros(Nfrs, 1);

for i = 1:Nfrs
    evel = ev(i,:)';
    Wpath(i) =evel'*Wep(:,:,i)*evel;
    Pr = eye(3) - evel*evel'; % Projection onto plane normal to velocity vector

    Wepproj = Pr'*Wep(:,:,i)*Pr;
    Wnormal(i) = norm(Wepproj);
end


