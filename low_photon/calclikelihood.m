function lik = calclikelihood(time,data,time_irf,irf,param_mat,...
    nexpo,fit_start,fit_end)

%   function lik = calclikelihood(time,data,irf,time_irf,param_mat,...
%    nexpo,fit_start,fitend)
%
%   Calculates likelihood function for fluorescence decay model based on
%   IRF and data given.
%
%   time: time vector in data (row vector)
%   data: decay data  (row vector, should be in the same size as time vector)
%   irf: irf vector  (row vector)
%   time_irf: irf time vector (row vector, should be in the same size as
%   irf vector)
%   param_mat: (Number of parameter)x(Number of grids in likelihood
%   function) matrix, contains trial parameters in each grid.
%   nexpo:  Number of exponential
%   fit_start, fit_end : time channels defining the time channels of your interest


dt_irf = time_irf(2)-time_irf(1);
irf = double(irf/(sum(irf)*dt_irf));

y_dat = data(fit_start:fit_end);

adc_ratio = round(length(irf)/length(time));

T = 12.58;
t = (0:dt_irf:T)';

% lengths
Nt = length(t);
Nirf = length(irf);
Ndata = length(data);

Nsubs = size(param_mat,2);

nfft = 2^nextpow2(2*length(t)+length(irf)-1);

%Total counts between fit_start and fit_end
totcount = sum(y_dat);

llik_vec = zeros(Nsubs,1);

tau1 = -999;
shift = -999;
if length(unique(param_mat(1,:))) == 1
    % shift IRF and Fourier transform
    shift = round(param_mat(1,1));
    Firf = fft(mycircshift(irf,shift),nfft);
    
    if nexpo == 1
        for i = 1:Nsubs
            
            param = param_mat(:,i);
            
            if tau1 ~= param(3)
                tau1 = param(3);
                decay1 = exp(-t/tau1);
            end
            
            decay_model = param(2)*decay1+(1-param(2));
            decay_model(end+1:end+Nt) = decay_model(1:Nt);
            
            %Fourier Transform of model
            Fmodel = fft(decay_model,nfft);
            
            Fmodel = Fmodel.*Firf;
            
            convmodel = real(ifft(Fmodel));
            
            convmodel = convmodel(Nt+1:Nt+Nirf);
            
            convmodel = reshape(convmodel,round(adc_ratio),Ndata);
            convmodel = sum(convmodel,1);
            p_model = convmodel(fit_start:fit_end)';
            
            %model normalization (max to total counts)
            p_model = p_model/sum(p_model)*totcount;
            
            %likelihood function calculation
            %vectorized N-param-dimensional likelihood function
            llik_vec(i) = sum((y_dat).*log(p_model));
            
            switch i
                case 1
                    disp('progress: 0%')
                case round(Nsubs/5)
                    disp('progress: 20%')
                case round(2*Nsubs/5)
                    disp('progress: 40%')
                case round(3*Nsubs/5)
                    disp('progress: 60%')
                case round(4*Nsubs/5)
                    disp('progress: 80%')
                case round(5*Nsubs/5)
                    disp('progress: 100%, fine search complete')
            end
        end
    elseif nexpo == 2
        
        for i = 1:Nsubs
            
            param = param_mat(:,i);
            
            if tau1 ~= param(3)
                tau1 = param(3);
                decay1 = exp(-t/tau1);
            end
            
            decay_model = param(2)*param(4)*decay1...
                +param(2)*(1-param(4))*exp(-t/(param(5)*param(3)))+(1-param(2));
            decay_model(end+1:end+Nt) = decay_model(1:Nt);
            
            %Fourier Transform of model
            Fmodel = fft(decay_model,nfft);
            
            Fmodel = Fmodel.*Firf;
            
            convmodel = real(ifft(Fmodel))*dt_irf;
            
            convmodel = convmodel(Nt+1:Nt+Nirf);
            
            convmodel = reshape(convmodel,round(adc_ratio),Ndata);
            convmodel = sum(convmodel,1);
            p_model = convmodel(fit_start:fit_end)';
            
            %model normalization (max to total counts)
            p_model = p_model/sum(p_model)*totcount;
            
            %likelihood function calculation
            %vectorized N-param-dimensional likelihood function
            llik_vec(i) = sum((y_dat).*log(p_model));
            
            switch i
                case 1
                    disp('progress: 0%')
                case round(Nsubs/5)
                    disp('progress: 20%')
                case round(2*Nsubs/5)
                    disp('progress: 40%')
                case round(3*Nsubs/5)
                    disp('progress: 60%')
                case round(4*Nsubs/5)
                    disp('progress: 80%')
                case round(5*Nsubs/5)
                    disp('progress: 100%, fine search complete')
            end
        end
        
    end
else
    if nexpo == 1
        for i = 1:Nsubs
            
            param = param_mat(:,i);
            
            if shift ~= round(param(1))
                shift = round(param(1));
                Firf = fft(mycircshift(irf,shift),nfft);
            end
            
            if tau1 ~= param(3)
                tau1 = param(3);
                decay1 = exp(-t/tau1);
            end
            
            decay_model = param(2)*decay1+(1-param(2));
            decay_model(end+1:end+Nt) = decay_model(1:Nt);
            
            %Fourier Transform of model
            Fmodel = fft(decay_model,nfft);
            
            Fmodel = Fmodel.*Firf;
            
            convmodel = real(ifft(Fmodel));
            
            convmodel = convmodel(Nt+1:Nt+Nirf);
            
            convmodel = reshape(convmodel,round(adc_ratio),Ndata);
            convmodel = sum(convmodel,1);
            p_model = convmodel(fit_start:fit_end)';
            
            %model normalization (max to total counts)
            p_model = p_model/sum(p_model)*totcount;
            
            %likelihood function calculation
            %vectorized N-param-dimensional likelihood function
            llik_vec(i) = sum((y_dat).*log(p_model));
            
            switch i
                case 1
                    disp('progress: 0%')
                case round(Nsubs/5)
                    disp('progress: 20%')
                case round(2*Nsubs/5)
                    disp('progress: 40%')
                case round(3*Nsubs/5)
                    disp('progress: 60%')
                case round(4*Nsubs/5)
                    disp('progress: 80%')
                case round(5*Nsubs/5)
                    disp('progress: 100%, fine search complete')
            end
        end
    elseif nexpo == 2
        
        for i = 1:Nsubs
            
            param = param_mat(:,i);
            
            if shift ~= round(param(1))
                shift = round(param(1));
                Firf = fft(mycircshift(irf,shift),nfft);
            end
            
            if tau1 ~= param(3)
                tau1 = param(3);
                decay1 = exp(-t/tau1);
            end
            
            decay_model = param(2)*param(4)*decay1...
                +param(2)*(1-param(4))*exp(-t/(param(5)*param(3)))+(1-param(2));
            decay_model(end+1:end+Nt) = decay_model(1:Nt);
            
            %Fourier Transform of model
            Fmodel = fft(decay_model,nfft);
            
            Fmodel = Fmodel.*Firf;
            
            convmodel = real(ifft(Fmodel))*dt_irf;
            
            convmodel = convmodel(Nt+1:Nt+Nirf);
            
            convmodel = reshape(convmodel,round(adc_ratio),Ndata);
            convmodel = sum(convmodel,1);
            p_model = convmodel(fit_start:fit_end)';
            
            %model normalization (max to total counts)
            p_model = p_model/sum(p_model)*totcount;
            
            %likelihood function calculation
            %vectorized N-param-dimensional likelihood function
            llik_vec(i) = sum((y_dat).*log(p_model));
            
            switch i
                case 1
                    disp('progress: 0%')
                case round(Nsubs/5)
                    disp('progress: 20%')
                case round(2*Nsubs/5)
                    disp('progress: 40%')
                case round(3*Nsubs/5)
                    disp('progress: 60%')
                case round(4*Nsubs/5)
                    disp('progress: 80%')
                case round(5*Nsubs/5)
                    disp('progress: 100%, fine search complete')
            end
        end
        
    end
end


lik = exp(llik_vec-max(llik_vec));

