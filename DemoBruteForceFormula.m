clear all
close all
clc
%% define necessary params
target = 506.995745; % this is the monoisotopic mass of ATP (C10H16N5O13P3)
ppm = 5; % distance in parts-per-million to search
load El_Mass_Abund
xx=seqmatch({'C','H','N','O','S','P'},ele);
libmass=mass(xx,1);
libminmax=[1 100; 1 400;0 10;0 20;0 10; 0 10]
maxCH=200;
f=BruteForceFormula(target,ppm,libmass,libminmax,maxCH);
%% see the distribution of masses
hist(f*libmass,100)
%% see how well the list is reduced using known C & N elemental count
xx=find(f(:,1)==10 & f(:,3)==5)
f(xx,:)
%% see how well it works with allcomb
load El_Mass_Abund
target=500;
tol=0.01
xx=seqmatch({'C','H','O','N','S'},ele);
libmass=mass(xx,1); %need to replace 1 with most likely isotopes
[libmass sx]=sort(libmass); %sort it now so we can deal with heavies elements first
bnds=zeros(2,length(libmass));
bnds(2,:)=ceil((target+tol)./libmass') %compute upper and lower bounds for atom counts
tic
fmat=allcomb([bnds(1,1):bnds(2,1)],[bnds(1,2):bnds(2,2)],[bnds(1,3):bnds(2,3)],...
    [bnds(1,4):bnds(2,4)],...
    [bnds(1,5):bnds(2,5)]);
m=fmat*libmass;
idx=abs(m-target)>tol;
m(idx)=[];
fmat(idx,:)=[];
size(fmat)
toc
size(fmat)
