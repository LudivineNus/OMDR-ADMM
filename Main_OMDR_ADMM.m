%---------------------------------------------------------------
% Main OMDR-ADMM
%---------------------------------------------------------------

clear; clc; close all;

% load data cube X of dimensions L x P x K 
% A: true abundances 
% S: true endmembers 

% L: number of wavelengths
% P: number of pixels 
% K: number of slices

L = size(X,1);
P = size(X,2);
K = size(X,3);

% Set the number of endmembers R 

% For X, R=3

% Matrices initialization 

S = rand(L,R); 
V = zeros(R,P); 
U = zeros(L,R); 
PI = zeros(R,P); 
Lamb = zeros(L,R); 
M = zeros(R,R); 
N = zeros(L,R); 


% Parameter estimates for each new slice k of the hyperspectral image X

for k = 1 : K

Xnew = squeeze(X(:,:,k)); % Xnew: each new slice of the hyperspectral image X

[S,A,V,U,PI,Lamb,M,N] = OMDR_ADMM(Xnew,S,V,U,PI,Lamb,M,N,R,'rho',0.01,'alpha',0.99,'mu',0.001,'N1',100,'N2',10);

A_t(:,:,k) = A; %concatenation of the abundances for each slice k
S_t(:,:,k) = S; %concatenation of the endmembers for each slice k


end

