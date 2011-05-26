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