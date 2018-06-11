Nseq = [10,8,4,2,1];
for i = Nseq
    fejer(i, 1)
    hold on
end

Legend=cell(5,1);
Legend{1}='F_{10}(f)';
Legend{2}='F_8(f)';
Legend{3}='F_4(f)';
Legend{4}='F_2(f)';
Legend{5}='F_1(f)';
legend(Legend);
saveas(gcf, 'simfejer.jpg')

function [kernel] = fejer(N, delta)
fmax = 1/(2*delta);
frange = -fmax:0.001:fmax;
kernelnum = delta * (sin(N*pi*frange*delta)).^2;
kerneldenom = N * (sin(pi*frange*delta)).^2;
kernel = kernelnum./kerneldenom;
plot(frange, kernel)
legend('N')
end





