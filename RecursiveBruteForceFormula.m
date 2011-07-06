% given a list of possible elements, a monoisotopic mass, and a tolerance,
% return the list of chemical formulae that are within this tolerance.
%
% definitions:
% target mass = t;
% tolerance = dt;
% elemental monoisotopic masses = [a1 a2 a3 a4 a5]
% elemental atom counts = [n1 n2 n3 n4 n5]
%
% the idea of recursion is to calculate the chemical formulae that would be
% subset against a particular heteoatom composition.
function guess = RecursiveBruteForceFormula(x,tol)
clc
set(0,'RecursionLimit',200)
load El_Mass_Abund
xx=seqmatch({'C','H','O'},ele);
libmass=mass(xx,1); %need to replace 1 with most likely isotopes
[libmass sx]=sort(libmass); %sort it now so we can deal with heavies elements first
bnds=zeros(2,length(libmass));
bnds(2,:)=ceil((x+tol+1)./libmass')+1 %compute upper and lower bounds for atom counts
bnds(2,bnds(2,:)>10)=10;
bnds(2,1)=100;
t=0;
tic
guess=bnds(2,:); % start with the low bounds
element1=size(bnds,2);
element2=size(bnds,2);
eshift=0;
while ~t
    
    [t,guess]=test(x,tol,libmass,bnds,element1,element2,eshift,guess);
    % x is the mass,
    % tol is the tolerance,
    % libmass is the mass of all elements,
    % bnds is the max/min for the elements,
    % element1 is the current element being investigated
    % element2 is the element in the loop
    % guess is the chemical formula being considered
end
toc

function [t,guess]=test(x,tol,libmass,bnds,element1,element2,eshift,guess)
if closeEnough(x,tol,libmass,guess)
    t=1;
    guess=guess;
else
    [guess,element1,element2,eshift,bnds]=betterGuess(bnds,element1,element2,eshift,guess);
    [t,guess]=test(x,tol,libmass,bnds,element1,element2,eshift,guess);
end

function F=closeEnough(x,tol,libmass,guess)
abs(guess*libmass-x)
F=abs(guess*libmass-x)<tol

function [guess,element1,element2,eshift,bnds]=betterGuess(bnds,element1,element2,eshift,guess)
if guess(element2)==0
    eshift=eshift+1;
end

if (element2-eshift)==1
    eshift=0;
    element1=size(bnds,2);
    element2=size(bnds,2);
end

if guess(element1)==0
    eshift=eshift+1;
    bnds(2,element2-eshift)=bnds(2,element2-eshift)-1;
    guess=bnds(2,:);
end
guess(element1)=guess(element1)-1
