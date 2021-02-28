
%    C, vector to be transformed
%    dt, time sampling

[Chat,freq]=centeredFFT(C,dt)
ipos = find(freq>0);
Ghat2 = abs(Ghat).^2;
figure
loglog(freq(ipos),Ghat2(ipos),'bo','MarkerFaceColor','b')
xlabel('frequency [cycles/yr]')
% The frequency range and x-tick labels may have to be changed.
set(gca,'XTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
set(gca,'XTickLabel',{'1/100000','1/10000','1/1000','1/100','1/10'})
ylabel('Spectral Power Density [put proper units]')
axis([5e-5 1e-2 1e-8 1e0])
set(gca,'YTick',[1e-8 1e-6 1e-4 1e-2 1e0])
% Consider putting your own title.
%title('TR163-31B, 4{^\circ}S, 84{^\circ}W, 3210m')
grid on





































