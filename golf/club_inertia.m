function I = club_inertia(z)
  %%  I = club_inertia(z)
  %% Returns the inertia of the (rest of the) club with respect to a
  %% point that is a distance z meters from the upper end of the shaft.
  %% The moment of inertia is defined with respect to axes defined with
  %%   z-axis along the shaft pointing down
  %%   y-axis along the upper edge of the blade, but normal to z
  %%   x-axis normal to z and y.

  %% Kjartan Halvorsen
  %% 2011-06-28

  shaft_mass = 0.3; % kg
  shaft_length = 1; % m
  shaft_radius_max = 0.005; % m
  shaft_radius_min = 0.003; % m
  shaft_com = [0; 0; 0.35*shaft_length]; % wrt top of shaft

  head_mass = 0.5; % kg
  head_radius = 0.03; % m
  head_com_0 = [0; 0.02; shaft_length+0.15]; % wrt top of shaft
  head_com = head_com - [0;0;z];

  Ihead_com = head_mass * 2/5 * head_radius*2 * eye(3);
  Ihead = Ihead_com + head_mass*( (head_com'*head_com)*eye(3) - \
                                 head_com*head_com');
