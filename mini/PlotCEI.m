clear
clf

load('CEI2006.mat')     % 120*120*2920
cloud=ncread('cloud2006.nc','cc');     % data from ECMWF, 121*121*2920
w1=ncread('w2006.nc','w');     % data from ECMWF, 121*121*2920
RH=ncread('RH2006.nc','r');     % data from ECMWF, 121*121*2920
tn=2920;

cloudL1=reshape((cloud(:,:,1,:)-cloud(:,:,2,:))./cloud(:,:,1,:),121,121,tn);
cloudH1=reshape((cloud(:,:,3,:)-cloud(:,:,2,:))./cloud(:,:,3,:),121,121,tn);
cloudL1(cloudL1<0)=0;
cloudH1(cloudH1<0)=0;
RHm1=reshape(RH(:,:,2,:),121,121,tn);

cloudL=zeros(120,120,tn);
cloudH=zeros(120,120,tn);
RHm=zeros(120,120,tn);
w=zeros(120,120,tn);
for xi=1:120
    for yi=1:120
        % change from 121*121 to 120*120
        cloudL(xi,yi,:)=mean(mean(cloudL1(xi:xi+1,yi:yi+1,:),1),2);
        cloudH(xi,yi,:)=mean(mean(cloudH1(xi:xi+1,yi:yi+1,:),1),2);
        RHm(xi,yi,:)=mean(mean(RHm1(xi:xi+1,yi:yi+1,:),1),2);
        w(xi,yi,:)=mean(mean(w1(xi:xi+1,yi:yi+1,:),1),2);
    end
end
clear cloudH1 cloudL1 RHm1 w1



%%

cloudL(w>0)=NaN;    % only consider Subsiding area
cloudH(w>0)=NaN;    % only consider Subsiding area
cloudLmean=reshape(nanmean(nanmean(cloudL,1),2),tn,1);
cloudHmean=reshape(nanmean(nanmean(cloudH,1),2),tn,1);

RHmstd=nan(tn,1);
rainstd=nan(tn,1);
sub=nan(tn,1);
for ti=1:tn
    RHmstd(ti)=std(reshape(RHm(:,:,ti),120*120,1),'omitnan');
    sub(ti)=nansum(nansum(w(:,:,ti)<0))/120/120;
end




scatter(MDI(1:2920),sub,'k.')
title('MDI vs Subsiding (w<0) fraction [at 2006]','FontSize',12)
xlabel('mass distance index')
ylabel('w<0 area fraction')


scatter(MDI(1:2920),RHmstd,'k.')
title('MDI vs standard deviation of RH [at 2006]','FontSize',12)
xlabel('mass distance index')
ylabel('standard deviation of RH at 500mb')

scatter(MDI(1:2920),cloudLmean,'k.')
title('MDI vs low cloud [at 2006]','FontSize',12)
xlabel('mass distance index')
ylabel('low cloud fraction in w<0 areas')


%%


ti=1449;

worldmap([0 30],[100 130])
load coastlines
pcolorm(0.125:0.25:29.875,100.125:0.25:129.875,rain(:,:,ti)')
plotm(coastlat,coastlon,'w')

colorbar
title('rain (mm/hr) [2006.7.1 00Z]','FontSize',12)
colormap('pink')





% 
% rainYES=zeros(120,120,tn);
% rainYES(rain>1)=1;
% RI=reshape(nansum(nansum(rain,1),2)./nansum(nansum(rainYES,1),2),tn,1);
% 
% MDI=nan(tn,1);
% for ti=1:tn
%     [yi,xi]=find(rain(:,:,ti)>1);
%     di=((xi-mean(xi)).^2+(yi-mean(yi)).^2).^0.5;
%     MDI(ti)=mean(di);
% end
% 



