

fid=fopen('ldpc_5GNR_BG2_mz5_kz6_z13H.dat.txt');
z = 13;
v = 2;
puncture = 1:1:z*v;

dspect = zeros(1, 30);

while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);     
    
    % 以 ',' 作为分割数据的字符,结果为cell数组
    s = regexp(str,' ','split');    
    s = (s > z*v);
    
    idx = nnz(str);
    dsprct(idx) = dpsect(idx) + 1; 
end
fclose(fid);
