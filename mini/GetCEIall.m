clear

% merge CEI date into 10 year time series
year=2006:2015;

tn=29216;
rainall=nan(120,120,tn);
CEIall=nan(120,120,tn);
ti=0;
for yi=1:10
    load(['CEI',num2str(year(yi)),'.mat'])
    ti=ti+size(rain,3);
    rainall(:,:,ti-size(rain,3)+1:ti)=rain;
    CEIall(:,:,ti-size(rain,3)+1:ti)=CEI;
end

CEImean=reshape(nanmean(nanmean(CEIall,1),2),tn,1);

MDI=nan(tn,1);
rainstd=nan(tn,1);
for ti=1:tn
    [yi,xi]=find(rainall(:,:,ti)>1);
    di=((xi-mean(xi)).^2+(yi-mean(yi)).^2).^0.5;
    MDI(ti)=mean(di);
    
    rainstd(ti)=std(reshape(rainall(:,:,ti),120*120,1),'omitnan');
end



save('CEIall.mat','CEImean','MDI','rainstd')

[soundY,soundFs] = audioread('D:\Dropbox\[Tools]\ended.mp3');
sound(soundY,soundFs);
clear soundY soundFs

