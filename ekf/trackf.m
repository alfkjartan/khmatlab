function [pp,F]=trackf(p,dt,n)
% [pp,F]=trackf(p,dt)
% Tracking filter for n dimensional state vector
%
% Input
%    p      ->   current state vector ([x, xdot]') 
%    dt     ->   sampling time. x(t+1) = x(t) + dt * xdot(t)
%    n      ->   degrees of freedom
% Output
%    pp     <-   new state vector
%    F      <-   linearized system matrix.

% Kjartan Halvorsen
% 2002-05-15

A=[eye(n) dt*eye(n) ; zeros(n,n) eye(n)];
pp=A*p;

F=A;