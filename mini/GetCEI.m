clear

path='E:\NTU\[R1.2]\Clouds and Environment\';

R=1;

rain=ncread([path,'3B42.2015.3hr.nc'],'pcp',[1121-R 201-R 1],[120+2*R 120+2*R Inf]);

tn=size(rain,3);
rainYES=zeros(122,122,tn);
rainYES(rain>1)=1;     % has rain = 1, no rain = 0

CEI=nan(120,120,tn);
for ti=1:tn
    for xi=1:120 
        for yi=1:120
            if rainYES(xi+1,yi+1,ti)==1     % calculate CEI when the pixel has rain
                CEI(xi,yi,ti)=nansum(nansum(rainYES(xi:xi+2,yi:yi+2,ti)))-1;
            end
        end
    end
end

rain=rain(1+R:120+R,1+R:120+R,:);

save('CEI2015.mat','rain','CEI')



