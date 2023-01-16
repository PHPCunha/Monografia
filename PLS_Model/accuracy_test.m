function [pvalue,dist_tt,meandiff,Resul] = accuracy_test(y,yhatA,yhatB,teste,niter,alpha)
% Teste randomino para comparacao da acuracias de dois modelos de
% calibracao multivariada.
% input
%   y   : valores y de referencia (do vetor de teste);
% yhatA : valores de y previsto por um dos modelos;
% yhatB : valores de y previsto pelo outro modelo;
% teste : tipo de teste:
%    randbi  : teste randomico bicaudal.
%    randuni : teste randomico unicaudal. yhatA > yhatB
%  tpareado  : teste-t pareado bicaudal.
% niter : numero de permutacoes; 500000
% alpha : nivel de significancia adotado; 0.05
%
% output
% pvalue : estatistica de teste.
% Se o pvalue < alpha : a diferenca na acuracia dos modelos e diferente.
%
% [pvalue,dist_tt,meandiff] = accuracy_test(y,ypA,ypB,randuni,500000,0.05)
%
% Paulo R. Filgueiras  - 26/06/2014
%
% Referencia
% Van der Voet, H. Chemom. Intel. Lab. Syst. 25, 1994, 313-323.
%
if nargin==3
    teste='randbi';
    niter=100000;
    alpha=0.05;
elseif nargin==4
    niter=100000;
    alpha=0.05;
elseif nargin==5
    alpha=0.05;
end

eA=y-yhatA;
eB=y-yhatB;
diff=eA.^2-eB.^2;
meandiff=mean(diff);
n=length(diff);
%niter=199;
sum=0;
dist_tt=[];
if strcmp(teste,'randbi');   % teste randomico bicaudal

for k=1:niter
    randomsign=2*round(rand(n,1))-1;
    signeddiff=randomsign.*diff;
    meansigneddiff=mean(signeddiff);
    dist_tt=[dist_tt;meansigneddiff];
    sum=sum+(abs(meansigneddiff)>=abs(meandiff));
end
pvalue=(sum+1)/(niter+1);

figure(1);
axes('FontSize',16,'FontName','Arial');
hist(dist_tt,50), vline(meandiff,'r');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k');
ylabel('Frequencia','FontSize',24,'FontName','arial');
xlabel('Distribuicao aleatoria','FontSize',24,'FontName','arial');


if pvalue < alpha
    s = sprintf('tcal = %g < alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com DIFERENCAS na acuracia.')
    Resul = "Diferente";
else
    s = sprintf('tcal = %g > alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com acuracias IGUAIS.')
    Resul = "Iguais";
end

elseif strcmp(teste,'randuni');    % teste randomico unicaudal
for k=1:niter
    randomsign=2*round(rand(n,1))-1;
    signeddiff=randomsign.*diff;
    meansigneddiff=mean(signeddiff);
    dist_tt=[dist_tt;meansigneddiff];
    sum=sum+(meansigneddiff>=meandiff);
end
pvalue=(sum+1)/(niter+1);


figure(1);
axes('FontSize',16,'FontName','Arial');
hist(dist_tt,50), vline(meandiff,'r');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','k');
ylabel('Frequencia','FontSize',24,'FontName','arial');
xlabel('Distribuicao aleatoria','FontSize',24,'FontName','arial');

if pvalue < alpha
    s = sprintf('p-cal = %g < alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com DIFERENCAS na acuracia.')
    Resul = "Diferente";
else
    s = sprintf('p-cal = %g > alpha = %g',pvalue,alpha); disp(s)
    disp('Modelos com acuracias IGUAIS.')
    Resul = "Iguais";
end

elseif strcmp(teste,'tpareado');
% teste-t em para medias (teste bicaudal)
diff1 = eA-eB;
meandiff1=mean(diff1);
stddiff1=std(diff1);
tcalc = (meandiff1*sqrt(length(diff1)))/(stddiff1);
tstat = tinv((1-alpha/2),length(diff1)-1);
pvalue = 1-tcdf(tcalc,length(diff1)-1);

disp(' Teste-t em para para medias ');

if tcalc < tstat
    s = sprintf('tcal = %g < t-tabelado = %g',tcalc,tstat); disp(s)
    disp('Modelos com acuracias IGUAIS')
    s = sprintf('pvalor = %g',pvalue); disp(s)
    Resul = "Iguais";
else
    s = sprintf('tcal = %g > t-tabelado = %g',tcalc,tstat); disp(s)
    disp('Modelos com DIFERENCAS na acuracia')
    s = sprintf('pvalor = %g',pvalue); disp(s);
    Resul = "Diferente";
end

end

function h = vline(x,lc)
[m,n] = size(x);
if m>1&n>1
  error('Error - input must be a scaler or vector')
elseif n>1
  x   = x';
  m   = n;
end

v     = axis;
if ishold
  for ii=1:m
    h(ii) = plot([1 1]*x(ii,1),v(3:4),lc);
  end
else
  hold on
  for ii=1:m
    h(ii) = plot([1 1]*x(ii,1),v(3:4),lc);
  end
  hold off
end

if nargout == 0;
  clear h
end
