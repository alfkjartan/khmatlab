% hand example
clear all;
load hand_data;

%model structure
% parent joint sensor axis category
segment = [...
   0 0 0 0 0;     %1  root
   1 2 0 0 1;     %2
   2 4 0 0 1;     %3  hand
   
   3 3 0 3 1;     %4
   4 3 0 1 1;     %5  index 1
   5 1 0 2 2;     %6
   6 3 0 1 1;     %7  index 2
   7 1 0 2 2;     %8
   8 3 0 1 1;     %9  index 3
   
   3 1 0 1 2;     %10
  10 3 0 3 1;     %11
  11 3 0 1 1;     %12 middle 1
  12 1 0 2 2;     %13
  13 3 0 1 1;     %14 middle 2
  14 1 0 2 2;     %15
  15 3 0 1 1;     %16 middle 3
  
   9 2 2 0 2;     %17 pos 1
   7 2 2 0 2;     %18 pos 2
   7 2 2 0 2;     %19 pos 3
   5 2 2 0 2;     %20 pos 4
   5 2 2 0 2;     %21 pos 5
   
  16 2 2 0 2;     %22 pos 6
  14 2 2 0 2;     %23 pos 7
  14 2 2 0 2;     %24 pos 8
  12 2 2 0 2;     %25 pos 9
  12 2 2 0 2;     %26 pos 10
   
   3 2 2 0 2;     %27 pos 11
   3 2 2 0 2;     %28 pos 12
   3 2 2 0 2;     %29 pos 13
   3 2 2 0 2;     %30 pos 14
];

% noise and dynamics model (compressed)
dt = 1/66;
datS = [10, 50];
datR = [0, 50];
datV = 1;
ratio = 1;

% prepare and estimate
[map, info, S, R, V] = prepare(segment,3,dt,datS,datR,datV,ratio);
X = estimate_ekf(data, zeros(map.nX,1),S,segment,map,info,R,V);

% plot estimated joint angles
hinges = [7 8 10 12 14 15 17 19];
figure(1);
clf;
plot(X(hinges,:)');
ylabel('estimated joint angles (rad)');
xlabel('samples (66Hz)');
