function [tid] = get_time(packet_number_time, packet)
  %% Returns the time corresponding to the given packet, using the
  %% matrix packet_number_time as a look up table. 

  %% Kjartan Halvorsen
  %% 2012-09-10

  [tf, indx] = ismember(packet, packet_number_time(:,1));
  if tf
    tid = packet_number_time(indx, 2);
    return
  end

  before = find(packet_number_time(:,1) < packet);
  after = find(packet_number_time(:,1) > packet);
  
  %% Interpolate
  pb = packet_number_time(before(end),:);
  pa = packet_number_time(after(1),:);
  
  s = (packet-pb)/(pa-pb);
  
  tid = pb(2) + s*(pa(2) - pb(2));
