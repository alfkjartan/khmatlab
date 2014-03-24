% example of 3-link model corresponds to fig. 5 in the paper
clear all;

% parent joint sensor axis category
segment = [...
   0 0 0 0 0;  %1  root
   1 2 0 0 1;  %2
   2 4 0 0 1;  %3  link 1
   3 3 0 3 1;  %4  link 2
   4 1 0 1 2;  %5
   5 3 0 2 1;  %6  link 3
   3 1 2 1 2;  %7  pos 1
   7 4 4 0 2;  %8  quat 1
   6 1 2 1 2;  %9  pos 2
   9 4 4 0 2;  %10 quat 2
];

% noise and dynamics model (compressed)
ratio = 10;
dt = 1/20;
datR = [0, 20];
datS = [3, 20];
datV = 0.2;

% simulation parameters
% set to nD=400, nRep=50 to reproduce fig. 5 from the paper
nD = 200;                               % length of data sequence
nRep = 10;                              % number of data sequences


% prepare information structures
[map, info, S, R, V] = prepare(segment, 3, dt,datS,datR,datV,ratio);
map1 = prepare(segment, 1);

% define constant elements of true state vector
x_w = [20 10 .1 .2 .3 30 0 0.5 1]';

% assign true state
xTrue = zeros(map.nX,nD);
xTrue(info.type==2,:) = repmat(x_w,[1,nD]);

nrm = [];
r2 = [];
r21 = [];
SaveX1 = [];
SaveX2 = [];
SaveSte1 = [];
SaveSte2 = [];

figure(1);
clf;
select = [8 10 14; 15 16 17];

fprintf('Running %d sequences\n', nRep);
for it = 1:nRep
   % generate random smooth sequence
   pos = randn(map.nX,nD+40)*20;
   [b,a] = butter(4,2*dt);
   pos = filtfilt(b,a,pos')';
   pos = pos(:,21:nD+20);
   pos(info.spatial==2,:) = pos(info.spatial==2,:)/ratio;
   pos = pos(info.type==1,:);
   xTrue(info.type==1,:) = pos;

   x_w1 = x_w + randn(9,1).*[3 3 .3 .3 .3 3 .3 .3 .3]';
   x0 = zeros(map.nX,1);
   x0(info.type==2) = x_w1;

   % generate noisy observations
   Y = zeros(map.nY,nD);
   for d=1:nD
      [dummy, Y(:,d)] = residual(Y(:,d),xTrue(:,d),segment,map);
   end
   for n=1:size(segment,1)
      if segment(n,3),
         if segment(n,3)<3,        % translation
            ii = (map.aY(n):map.aY(n)+2);
            Y(ii,:) = Y(ii,:) + datV*randn(3,nD);
         else                       % rotation
            ii = (map.aY(n):map.aY(n)+3);
            Y(ii,:) = Y(ii,:) + datV/ratio*randn(4,nD);
            Y(ii,:) = Y(ii,:) ./ (ones(4,1)*sqrt(sum(Y(ii,:).^2)));
         end
      end
   end

   % get estimated sequence
   [X, Ste, S1] = estimate(Y, x0, S, segment, map, info, R, V);

   % get estimated sequence with fixed parameters
   x1 = x0;
   x1(info.type==2) = median(X(info.type==2,nD/2:end),2);
   flg = (info.type~=2);
   S1(flg,flg) = S(flg,flg);   
   X1 = estimate(Y, x1, S1, segment, map1, info, R, V);

   % compute R2
   good = 1;
   for n=1:map.nX
      if info.type(n)==1,
         cc = corrcoef(xTrue(n,:)',X1(n,:)');
         r21 = [r21; cc(1,2)^2];
         cc = corrcoef(xTrue(n,:)',X(n,:)');
         r2 = [r2; cc(1,2)^2];
         if cc(1,2)^2<0.95,
            good = 0;
         end
      end
   end

   if good,
      % normalized errors
      tmp = X(info.type==1,nD/2:end)-xTrue(info.type==1,nD/2:end);
      tmp = tmp ./ Ste(info.type==1,nD/2:end);
      nrm = [nrm; tmp(:)];

      xx = X(info.type==2,:);
      ss = Ste(info.type==2,:);
      SaveX1 = [SaveX1; xx([1 2 6],:)-x_w([1 2 6])*ones(1,nD)];
      SaveX2 = [SaveX2; xx([3 4 5 7 8 9],:)-x_w([3 4 5 7 8 9])*ones(1,nD)];
      SaveSte1 = [SaveSte1; ss([1 2 6],:)];
      SaveSte2 = [SaveSte2; ss([3 4 5 7 8 9],:)];
      
      % plot results
      for k = 1:size(select,1)
         subplot(1,2,k);
         set(plot(X(select(k,:),:)','k'),'linewidth',0.5);
         hold on;
      end
   end

   fprintf('\n');
end

subplot(1,2,1);
axis([0 nD 5 35]);
set(gca,'xtick',0:100:nD,'xticklabel',0:5:20,'ytick',[10 20 30]);
box off;

subplot(1,2,2);
axis([0 nD -0.2 1.2]);
set(gca,'xtick',0:100:nD,'xticklabel',0:5:nD/20,'ytick',[0 1]);
box off;

figure(2);
clf;
set(plot(xTrue([7 9],:)'),'linewidth',4,'color',[.75 .75 .75]);
hold on;
set(plot(X([7 9],:)','k'),'linewidth',1);
axis([0 nD -2 2]);
set(gca,'xtick',0:100:nD,'xticklabel',0:5:nD/20,'ytick',-2:2);
box off;

figure(3);
clf;
[nn,xx] = hist(nrm,-4.1:0.1:4.1);
plot(xx(2:end),nn(2:end)/max(nn(2:end)),'.-');
hold on;
xx = -4:0.1:4;
plot(xx,exp(-xx.^2/2),'r');
axis([-4 4 0 1.05]);
set(gca,'ytick',[]);
box off;

figure(4);
clf;
subplot(2,1,1);
set(plot(mean(SaveSte1)),'linewidth',4,'color',[.75 .75 .75]);
hold on;
set(plot(std(SaveX1),'k'),'linewidth',1);
set(gca,'xtick',0:100:nD,'xticklabel',0:5:nD/20,'ytick',0:3);
axis([0 nD 0 3]);
box off;

subplot(2,1,2);
set(plot(mean(SaveSte2)),'linewidth',4,'color',[.75 .75 .75]);
hold on;
set(plot(std(SaveX2),'k'),'linewidth',1);
set(gca,'xtick',0:100:nD,'xticklabel',0:5:nD/20,'ytick',0:0.1:0.3);
axis([0 nD 0 0.3]);
box off;
