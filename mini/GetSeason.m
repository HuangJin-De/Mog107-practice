clear

% devide 2006-2015 tim series into DJA and JJA
cloudall=ncread('cloud2006-2015.nc','cc');
tn=size(cloudall,4);
cloudL=reshape((cloudall(:,:,1,:)-cloudall(:,:,2,:))./cloudall(:,:,1,:),121,121,tn);
cloudL(cloudL<0)=0;
clear cloudall

wall=ncread('w2006-2015.nc','w');

dayn=repmat([31,28,31,30,31,30,31,31,30,31,30,31],10,1);
dayn(3,2)=29;
dayn(6,2)=29;

tDJF=((31+31+28)*10+2)*8;
tJJA=(30+31+31)*10*8;
cloudLDJF=nan(121,121,tDJF);
cloudLJJA=nan(121,121,tJJA);
wDJF=nan(121,121,tDJF);
wJJA=nan(121,121,tJJA);

tiDJF=0;
tiJJA=0;
tiALL=0;
for mi=1:10
    tiDJF=tiDJF+sum(dayn(mi,1:2))*8;
    tiALL=tiALL+sum(dayn(mi,1:2))*8;
    cloudLDJF(:,:,tiDJF-sum(dayn(mi,1:2))*8+1:tiDJF)=cloudL(:,:,tiALL-sum(dayn(mi,1:2))*8+1:tiALL);
    wDJF(:,:,tiDJF-sum(dayn(mi,1:2))*8+1:tiDJF)=wall(:,:,tiALL-sum(dayn(mi,1:2))*8+1:tiALL);
    tiJJA=tiJJA+sum(dayn(mi,6:8))*8;
    tiALL=tiALL+sum(dayn(mi,3:8))*8;
    cloudLJJA(:,:,tiJJA-sum(dayn(mi,6:8))*8+1:tiJJA)=cloudL(:,:,tiALL-sum(dayn(mi,6:8))*8+1:tiALL);
    wJJA(:,:,tiJJA-sum(dayn(mi,6:8))*8+1:tiJJA)=wall(:,:,tiALL-sum(dayn(mi,6:8))*8+1:tiALL);
    tiDJF=tiDJF+dayn(mi,12)*8;
    tiALL=tiALL+sum(dayn(mi,9:12))*8;
    cloudLDJF(:,:,tiDJF-dayn(mi,12)*8+1:tiDJF)=cloudL(:,:,tiALL-dayn(mi,12)*8+1:tiALL);
    wDJF(:,:,tiDJF-dayn(mi,12)*8+1:tiDJF)=wall(:,:,tiALL-dayn(mi,12)*8+1:tiALL);
end

cloudLDJF(wDJF>0)=NaN;
cloudLDJFmean=reshape(nanmean(nanmean(cloudLDJF,1),2),tDJF,1);
cloudLJJA(wJJA>0)=NaN;
cloudLJJAmean=reshape(nanmean(nanmean(cloudLJJA,1),2),tJJA,1);


save('cloudLseason','cloudLDJFmean','cloudLJJAmean')


[soundY,soundFs] = audioread('D:\Dropbox\[Tools]\ended.mp3');
sound(soundY,soundFs);
clear soundY soundFs

