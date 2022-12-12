function modelo=eparamenter(modelo,Xcal,Xtest,ycal,ytest)
%------------------------------------------------------------------------------------------
% Funcao criada para gerar parametros de avaliacao (evaluation paramenter) 
% de forma simples e facil, acompanhado do pldmodel do artigo da
% IFES-Ciencia. 
% Compativel em Octave e Matlab.
%
%
% Pedro H. Pereida da Cunha 31/08/2022
%
% Entradas:
% modelo = modelo; Deriavo de plsmodel (IFES-Ciencia)
% Xcal   = matriz X de calibracao;
% Xtest  = matriz X de test; 
% ycal   = vetor y de calibração;
% ytest  = vetor y de teste;      
%     
% Saidas:
%      seletividade = vetor que contem a seletividade de cada amostra. 
%      sensibilidade = sensibilidade do modelo.
%      sen_analitica = sensibilidade analítica do modelos;
%      ruido = ruido espectral
%      sinal_ruido = razão sinal ruído.
%      LD = limite de detecção do modelo.
%      LQ = limite de quantificação.
%      
%---------------------------------------------------------------------------------------------

[Xcal,para1,para2]=pretreat(Xcal,'center');   

x  = Xcal;
y  = ycal;
b  = modelo.coef.B;
T  = modelo.coef.T;
P  = modelo.coef.P;
vl = modelo.options.vl;
yp = modelo.Ytest(:,2);
yc = modelo.Ycal(:,2);

[L1,C1]=size(x);
[Rc,mx]=center(x,1);
[cc,my]=center(y,1);
cest=Rc*b;

nnas_c=cest/norm(b);
Z(1:L1,1)=1;
Z(1:L1,2)=nnas_c;
ccc=cc(:,1)+ones(L1,1)*my;
cte=pinv(Z)*ccc;
constante=-cte(1,1)/cte(2,1);
nnas_cm(:,1)=nnas_c(:,1)-constante;

% Sensibilidade(escalar)
sen=1/norm(b);

% seletividade
for i=1:L1
    sel=nnas_cm./norm(x(i,:));
end

xp=T*P'; 
xpm=xp+ones(size(x,1),1)*mx;
er=x-xpm;
va=var(er);
vva=mean(va); %média da variancia do ruido
dx=sqrt(vva); % ruido espectral

sr=nnas_cm./dx;

%++++++++++ Sensibilidade analitica ++++++++++
sena=sen/dx; %  (SENSIBILIDADE ANALÍTICA) = 410.5

% Seletividade e Sensibilidade
EParamenter.seletividade=sel;
EParamenter.sensibilidade=sen;
EParamenter.sen_analitica =sen/dx;
EParamenter.inv_sen_analitiva = 1/EParamenter.sen_analitica;

% Sinal Ruido
EParamenter.ruido_espectral = dx;
EParamenter.sinal_ruido = sr;
if min(sr)<0; 
EParamenter.sinal_ruido_min =0; 
else
EParamenter.sinal_ruido_min = min(sr);
end
EParamenter.sinal_ruido_max = max(sr);

% LD & LQ
EParamenter.LD = 3*dx*(1/sen);;
EParamenter.LQ = 10*dx*(1/sen);

modelo.EParamenter = EParamenter;


%% SUB-ROTINA

% CENTER;
function [cdata,me,cnewdata]=center(data,opt,newdata)

%#										
%#  function [cdata,me,ctest]=center(data,opt,newdata);			
%#										
%#  AIM: 	Centering along columns, rows or double centering		
%#										
%#  PRINCIPLE:  Removal of the column-, row- or overall mean from 		
%#              each column, row or both, respectively 		 		
%# 				 If a test data set is available it can ONLY be 
%#              column centered using the mean from the calibration
%#              data set.
%#
%#
%#  INPUT:	data: (m x n) matrix with m rows and n variables		
%#				opt: optional							
%#		     1 = column centering					
%#		     2 = row centering						
%#		     3 = double centering					
%#          newdata: (mt x n) test matrix with mt rows and n variables
%#					
%#			 							
%#  OUTPUT:	cdata: (m x n) matrix containing centered data			
%#				me: mean vector, overall mean (scalar)
%#              newdata: (mt*n) test matrix centered with the mean of data
%#	
%#										
%#  AUTHOR: 	Andrea Candolfi				 			
%#	    			Copyright(c) 1997 for ChemoAc					
%#          	FABI, Vrije Universiteit Brussel            			
%#          	Laarbeeklaan 103 1090 Jette					
%#   										
%# VERSION: 1.2 (25/02/2002)							
%#										
%#  TEST:   	I. Stanimirova	& S. Gourvénec & M. Zhang
%#										
	

[m,n]=size(data);

if nargin==1;
  opt=[4];
  while opt>3 || opt<=0 
    opt=input('column centering(1), row centering(2), double centering(3):');
  end
end


if opt==1			% column centering 
   me=mean(data);
   cdata=data-ones(m,1)*me;
end

if opt==2			% row centering
   me=mean(data')';
   cdata=data-me*ones(1,n);
end

if opt==3 	% double centering
   me=mean(mean(data));
   mej=mean(data');
   mei=mean(data);
   cdata=data-(ones(m,1)*mei)-(ones(n,1)*mej)'+(ones(m,n)*me);
end

if exist('newdata')==1			% center new data
    [mt,n]=size(newdata);
    
    if opt==1				% column centering 
        me=mean(data);
        cnewdata=newdata-ones(mt,1)*me;
    else
        error('Row centering and double centering are impossible to perform on a test set');
    end
    
end