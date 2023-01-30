global Tsampling t
Tsampling = 1e-7 ;
fsampling = 1/Tsampling ;
% alphabet = [0, 1] ;
% alphabet = [-1, 1] ;
 alphabet = [-3, -1, 1, 3] ;
Tb = 1e-6 ;
Ts = Tb * log2(length(alphabet)) ;
nbits = 1000 ; nsymbols = nbits / log2(length(alphabet)) ; 
% N0 = 0.02 ;
% N0 = 0.04 ;
 N0 = 0.1 ;
t = (0:Tsampling:(nsymbols*Ts-Tsampling)) ;

% Transmitted bits
tx_bits = round(rand(1,nbits)) ;

% Transmitted symbols
tx_symbols = map(tx_bits,alphabet) ;

% Transmitted signal (channel input)
tx_signal = modulate(tx_symbols,Ts,nsymbols) ;

% Received signal (channel output)
 rx_signal = tx_signal + sqrt(fsampling * N0/2) * randn(1,length(t)) ;
% rx_signal = tx_signal;

% Received symbols
rx_symbols = demodulate(rx_signal,Ts,nsymbols) ;

% Received bits
rx_bits = demap(rx_symbols,alphabet) ;


plotfig = scatterplot(rx_symbols);
hold on;
scatterplot(tx_symbols, 1,0,'r*',plotfig);
hold off;

 figure
  plot(rx_symbols);
  plot(tx_bits);







% Task 1.1 - Mapper
%function tx_symbols = map(tx_bits,alphabet)
    % Initialise the tx_symbols with zeros
%    tx_symbols = zeros(1, length(tx_bits)/log2(length(alphabet)));
    % Create a matrix aux 
%    aux = reshape(tx_bits, length(tx_symbols), log2(length(alphabet)));
    % map the bits to alphabets
%    for k = 1:length(tx_symbols)
%        tx_symbols(k) = alphabet(bi2de(aux(k,:))+1);
%    end
%end


function tx_symbols = map(tx_bits,alphabet)
    nbits_tx_sy = log2(length(alphabet));
    length_tx_symbols = length(tx_bits)/nbits_tx_sy;
    aux_0 = reshape(tx_bits,[nbits_tx_sy,length_tx_symbols]);
    aux=aux_0';
    for k = 1:length_tx_symbols
        tx_symbols(k) = alphabet(bi2de(aux(k,:))+1);
    end
end 

% Task 1.2 Baseband Modulator
function tx_signal = modulate(tx_symbols,Ts,nsymbols) 
    global t
    tx_signal = zeros(1,length(t)) ;
    for k = 1:nsymbols
        tx_signal = tx_signal + tx_symbols(k) * (1/sqrt(Ts) * rectpuls(t-Ts/2-(k-1)*Ts,Ts)) ;
    end
end

% Task 2.1 Demodulator
function rx_symbols = demodulate(rx_signal,Ts,nsymbols)
    global t
    global Tsampling
    rx_symbols = zeros(1, length(nsymbols));
    for k = 1:nsymbols
        rx_symbols(k) = Tsampling * sum(rx_signal .* (1/sqrt(Ts) * rectpuls(t- Ts/2-(k-1)*Ts,Ts))) ;
    end
end

% Task 2.2 Demapper
function rx_bits = demap(rx_symbols,alphabet)
    rx_bits = zeros(log2(length(alphabet)),length(rx_symbols)) ; 
    for k1 = 1:length(rx_symbols)
        aux = Inf ;
        for k2 = 1:length(alphabet)
            if (rx_symbols(k1) - alphabet(k2))^2 < aux
                aux = (rx_symbols(k1) - alphabet(k2))^2 ;
                rx_bits(:,k1) = de2bi(k2-1,log2(length(alphabet)))' ;
            end 
        end
    end 
    rx_bits = reshape(rx_bits,length(rx_bits)*log2(length(alphabet)),1)';
end

















