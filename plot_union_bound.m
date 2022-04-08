dmin_spect = 5;
figure;
snr = -1:0.1:10;
union_bound = zeros(1, length(snr));
idx = 1;
n = 11*13;
for isnr_db = -1:0.1:10
    isnr = 10^(isnr_db/10);
    union_bound(idx) = 3 * 13 * qfunc(sqrt(5*isnr));
    idx = idx + 1;
end

semilogy(snr, union_bound);
hold on;
plot(6.81,2.84*10^(-5),'o','color','r', 'MarkerSize', 5)
