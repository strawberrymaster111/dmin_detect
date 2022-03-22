clc;
BG_idx = 2;
kz=6;
mz=5;
z=13;
thresh=30;
opt=1; %0-Berrou,1 default;
inum=3;%fastest
rate=kz/(kz+mz-2);
Param1=[' -code ','D:/LDPC/dmin_detect/data/ldpc_5GNR_BG' num2str(BG_idx) '_mz' num2str(mz) '_kz' num2str(kz) '_z' num2str(z) 'H.dat'];
Param2=[' -outFileCW ','ldpc_5GNR_BG' num2str(BG_idx) '_mz' num2str(mz) '_kz' num2str(kz) '_z' num2str(z) 'H.dat.txt'];
Param3=[' -thresh ',num2str(thresh)];
Param4=[' -speedup ',num2str(inum)];
Param5=[' -switch ',num2str(opt)];
Exefilename='Main.exe';
Exepath=fullfile('.\',[Exefilename,Param1,Param2,Param3,Param4,Param5]);
system(Exepath);

