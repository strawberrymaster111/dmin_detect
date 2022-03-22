%The program will transfer any parity check matrix into compact format for
%Approximate Nearest Words calculation;
%该程序计算由LDPC BG2->ANC程序的H矩阵；
BG_idx = 2;
z_idx = 2; %1,2,...,8 for {2,3,4,7,9,11,13,15}*2^n
z =13;
A0 = LDPC_5GNR(BG_idx,z_idx); %full base graph with lifting
if BG_idx==2, kz_max=10; mz_max=42; vz=2; gz=4; end
kz = 6; mz = 14;
A0 = (A0>0) .* (mod(A0-1,z)+1);
A = A0(1:mz, [1:kz kz_max+1:kz_max+mz]);
rate=kz/(kz+mz-2);
%rate=1/2,1/3;
%PMG 6*13
% A = [
%  1   1   1   1   1   1   1   1   0   0   0
%  9   7  13   8   3   6   2   1   1   0   0
% 11   0   5   2   4  12  10   0   1   0   0
%  6   1   0   0   0   0   0   0   0   1   0
%  1   0  12   0   0   3  10   0   0   0   1
%    ];
%num2str(A)
H = A2H(A,z);
 fid = fopen(['D:/LDPC/LDPC_Analysis/results/''ldpc_5GNR_BG' num2str(BG_idx) '_mz' num2str(mz) '_kz' num2str(kz) '_z' num2str(z) 'H.dat'],'w');
% fid = fopen('PMG6_13Hmac.dat','w');
 %以下三组为test H
%H=[1 0 0 1 0 0 1 0 0;0 1 0 0 1 0 0 1 0;0 0 1 0 0 1 0 0 1;1 0 0 1 0 0 1 0 0;0 1 0 0 1 0 0 1 0;0 0 1 0 0 1 0 0 1];
% H=[1 0 0 1 0 1 1;0 1 0 1 1 1 0;0 0 1 0 1 1 1];
%  H=[0 0 1 0 0 0 0 1 0 1 0 0 1 0 0 0
%  0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0
%  0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0
%  0 0 0 1 0 1 0 0 0 0 1 0 0 1 0 0
%  0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1
%  1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0
%  0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0
%  0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1
%  1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1
%  0 0 1 0 1 0 0 0 0 0 0 1 0 1 0 0
%  0 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0
%  1 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0];


[m,n]=size(H);
[r,c]=find(H);
index=sortrows([r,c],1);
rowmax=max(sum(H,2));
rowmax=full(rowmax);r=full(r);c=full(c);
fprintf(fid,'%d %d\n',n,m);
fprintf(fid,'%d ',rowmax);
 num=length(r);
dat=zeros(m,rowmax);j=1;
for i=1:num
    count=sum(r(:)==index(i,1));
    if j>count j=1; end
    dat(index(i,1),j)=index(i,2);
    j=j+1;
end
dat=sort(dat,2);
fprintf(fid,'\n');
for i=1:m
    count=0;
    for j=1:rowmax
        if(dat(i,j)) fprintf(fid,'%d ',dat(i,j));count=count+1; end;
    end
    while(count~=rowmax) fprintf(fid,'%d ',0); count=count+1;end
    fprintf(fid,'\n');
end