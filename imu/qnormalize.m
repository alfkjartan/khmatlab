function qn = qnormalize(q)
  %% qn = qnormalize(q)
  %% Normalizes a set of quaternions

  %% Kjartan Halvorsen
  %% 2011-03-20

  qn = q;

  for i=1:size(q,2)
    qn(:,i) = q(:,i) / norm(q(:,i));
  end

