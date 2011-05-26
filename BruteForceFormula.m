function fmat=BruteForceFormula(target,ppm,elmass,elminmax,maxCH)%% Brute Force Chemical Formula Generator
mmu=target*ppm/1e6;
a=nchoosek(0:maxCH,2);
a=cat(1,a,[a(:,2) a(:,1)]);
a=cat(1,a,[0:maxCH; 0:maxCH]');
mass=a*elmass(1:2);
a(mass>(target+mmu) | (mass+elmass(3:end)'*elminmax(3:end,2))<(target-mmu),:)=[]; %remove any that are too big already
c1=1;
fmat=zeros(1e6,size(elmass,1));

        for m = elminmax(6,1):elminmax(6,2)
            for k = elminmax(5,1):elminmax(5,2)
                for j = elminmax(4,1):elminmax(4,2)
                    for i = elminmax(3,1):elminmax(3,2)
                        newa=zeros(size(a,1),length(elmass));
                        newa(:,1:2)=a;
                        newa(:,3)=i;
                        newa(:,4)=j;
                        newa(:,5)=k;
                        newa(:,6)=m;
                        mass=newa*elmass;
                        newa((mass>=(target+mmu) | mass<=(target-mmu)),:)=[]; %remove any that are too big or small already
                        L=size(newa,1);
                        c2=c1+L-1;
                        fmat(c1:c2,:)=newa;
                        c1=c1+L;
                    end
                end
            end
        end
fmat(c1:end,:)=[];
% str=sprintf('%d formulas in %6.2f seconds',size(fmat,1),toc)
