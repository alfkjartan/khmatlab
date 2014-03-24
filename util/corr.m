function R2 = corr(x,y)
% Returns the square of the correlation coefficient

R = corrcoef(x,y);
R2 = R(2,1)*R(2,1);