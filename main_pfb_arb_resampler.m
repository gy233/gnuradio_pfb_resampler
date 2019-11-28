% book: Multirate Signal Processing for Communications Systems
% chapter 7.1, 7.2
% in GNU Radio, sampsPerSym = 32 (spb = sampling_freq/symbol_rate = 32)
%               Nsym = 11 * int(self._samples_per_symbol)

% use with /home/guyu/my_gnuradio_projects/cma/verify_block.grc
clc;clear;close all;

samples_per_symbol=2;

%% verify block data
data_length=10000;
filename='/home/guyu/my_gnuradio_projects/cma/pfb_resampler_in.bin';
[fid]=fopen(filename,'rb');
input_raw=fread(fid,data_length,'float32');
input_raw=input_raw(1:2:end)+1i*input_raw(2:2:end);
fclose(fid);

filename='/home/guyu/my_gnuradio_projects/cma/pfb_resampler_out.bin';
[fid]=fopen(filename,'rb');
output_verify=fread(fid,samples_per_symbol*data_length,'float32');
output_verify=output_verify(1:2:end)+1i*output_verify(2:2:end);
fclose(fid);

global d_taps_per_filter

%% PART1: polyphase filter bank %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RRC
% firdes.root_raised_cosine
nfilts=32;

gain=nfilts;
sampling_freq=nfilts;
symbol_rate=1;
alpha=0.35;
ntaps = nfilts * 11 * samples_per_symbol;    % make nfilts filters of ntaps each
[taps,t,num,den] = root_raised_cosine(gain,sampling_freq,symbol_rate,alpha,ntaps);

% RaisedCosineTransmitFilter
Nsym=11*samples_per_symbol;
sampsPerSym=nfilts;
rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', alpha, ...
    'FilterSpanInSymbols', Nsym, ...
    'OutputSamplesPerSymbol', sampsPerSym);
b1=coeffs(rctFilt);
rctFilt.Gain = 1/max(b1.Numerator);


figure(1)
plot(b1.Numerator/max(b1.Numerator),'b')
hold on
plot(taps,'r')
legend('RaisedCosineTransmitFilter','firdes.root raised cosine')

%% pfb_arb_resampler

% init pfb_arb_resampler
d_acc = 0; 

d_int_rate = 32;
[d_dec_rate,d_flt_rate]=set_rate(samples_per_symbol,d_int_rate);

d_last_filter = mod(floor(length(taps)/2),32)+1;

dtaps = [taps(2:end)-taps(1:end-1),0];
d_taps=create_taps(taps);
d_dtaps=create_taps(dtaps);

delay = samples_per_symbol * (d_taps_per_filter - 1.0) / 2.0;
d_delay = round(delay);

accum = d_delay * d_flt_rate;
accum_int = floor(accum);
accum_frac = accum - accum_int;
end_filter = round(mod(d_last_filter + d_delay * d_dec_rate + accum_int, d_int_rate));

d_est_phase_change = d_last_filter - (end_filter + accum_frac);

input=[zeros(d_taps_per_filter-1,1);input_raw];

% pfb_arb_resampler_ccf::filter
output=zeros(length(input),1);
i_out = 1;
i_in = 1;
j = d_last_filter;
j_rec=[];
d_acc_rec=[];
while(i_in <= length(input)-d_taps_per_filter)
    while(j <= d_int_rate)
        j_rec=[j_rec;j];
        d_acc_rec=[d_acc_rec;d_acc];
        o0=fliplr(d_taps(j,:))*input(i_in:i_in+d_taps_per_filter-1);
        o1=fliplr(d_dtaps(j,:))*input(i_in:i_in+d_taps_per_filter-1);
        output(i_out) = o0 + o1*d_acc;     % linearly interpolate between samples
        i_out=i_out+1;
        
        d_acc = d_acc + d_flt_rate;
        j = j + d_dec_rate + floor(d_acc);
        d_acc = mod(d_acc, 1);
    end
    i_in = i_in + floor(j / d_int_rate);
    j = mod(j,d_int_rate);
end


%% figures
figure
scatter(real(input),imag(input));
title('input (constellation)')

figure
subplot(121)
scatter(real(output_verify),imag(output_verify))
title('GNU radio')
subplot(122)
scatter(real(output),imag(output))
title('mine')
suptitle('output of pfb arb resampler (constellation)')

figure
plot(real(output),'b','LineWidth',2)
hold on
plot(real(output_verify),'r')
legend('GNU radio','mine')

%% PART2: equivalent filter %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% taps
nfilts=2;
gain=nfilts;
sampling_freq=nfilts;
symbol_rate=1;
alpha=0.35;
ntaps = nfilts * 11 * samples_per_symbol;    % make nfilts filters of ntaps each
[taps1,t1,num1,den1] = root_raised_cosine(gain,sampling_freq,symbol_rate,alpha,ntaps);

figure
scatter(t,taps,'b','LineWidth',2)
hold on
scatter(t1,taps1,'r')
legend('polyphase filter bank','equivalent filter')

%% filtering
input=reshape([input_raw,zeros(size(input_raw))].',2*length(input_raw),1);
output1=conv(input,taps1);

figure
plot(real(output),'b','LineWidth',2)
hold on
plot(real(output1),'r')
legend('polyphase filter bank','equivalent filter')

if isempty(find(abs(output-output1(1:length(output)))>0.01,1))
    disp('here is the equivalent filter')
end

%% PART3: matched filter %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RRC
Nsym=11*samples_per_symbol;
sampsPerSym=2;
rctFilt2 = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', alpha, ...
    'FilterSpanInSymbols', Nsym, ...
    'OutputSamplesPerSymbol', sampsPerSym);
b1=coeffs(rctFilt2);

figure
stem(t1,b1.Numerator/max(b1.Numerator),'b','LineWidth',2)
hold on
stem(t1,taps1,'r')
legend('RaisedCosineTransmitFilter','firdes.root raised cosine')

tx = 0: length(input_raw) - 1;
to = (0: length(input_raw)*sampsPerSym - 1) / sampsPerSym;
fltDelay = Nsym / 2;
yc = rctFilt2([input_raw; zeros(Nsym/2,1)]);
yc = yc(fltDelay*sampsPerSym+1:end);

figure
t_plot=100;
subplot(121)
stem(tx(1:t_plot), real(input_raw(1:t_plot)), 'kx'); 
hold on;
plot(to(1:t_plot*sampsPerSym), real(yc(1:t_plot*sampsPerSym)), 'b-'); 
hold off;

%% matched filter taps
rcrFilt = comm.RaisedCosineReceiveFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', alpha, ...
    'FilterSpanInSymbols', Nsym, ...
    'InputSamplesPerSymbol', sampsPerSym, ...
    'DecimationFactor',1);
yr = rcrFilt([yc; zeros(Nsym*sampsPerSym/2, 1)]); % equivalent to conv(yc,taps1,'same'); 
yr = yr(fltDelay*sampsPerSym+1:end);

subplot(122)
stem(tx(1:t_plot), real(input_raw(1:t_plot)), 'kx'); 
hold on;
plot(to(1:t_plot*sampsPerSym), real(yr(1:t_plot*sampsPerSym)), 'b-'); 
hold off;
