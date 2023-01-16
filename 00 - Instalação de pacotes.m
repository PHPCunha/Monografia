%% 00 - Instalação de pacotes

% Como estamos falando de praticas que envolvem calculos estatisticos e
% matematicos, e interessante que certos pacotes sejam instalados no Octave e
% puxados toda vez que iniciamos uma rotina no Octave, lembrando que isso não é
% necessario no matlab.


% Localização dos Pacotes
cd('C:\Users\Exemplo\...\Monografia\Pacotes');

% Pacote voltado para estatistica.
pkg install statistics-1.4.3.tar.gz

% Caso a instalacao ocorra corretamente o pacote irá aparecer na lista.
pkg list

% Pacote voltado para uso do excel.
pkg install io-2.6.4.tar.gz

% Esses sao os pacotes necessarios para aplicacao do nosso tutorial. Para
% puxa-los em uma rotina você só precisa fazer os seguintes comandos.
pkg load statistics
pkg load io

% Esse e a versao prologo do nosso tutorial, espero que aproveitem.
