classdef DTMF_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        SNRSlider                     matlab.ui.control.Slider
        SNRSliderLabel                matlab.ui.control.Label
        Label_4                       matlab.ui.control.Label
        Label_3                       matlab.ui.control.Label
        Label_2                       matlab.ui.control.Label
        Label                         matlab.ui.control.Label
        OutputEditField               matlab.ui.control.EditField
        OutputLabel                   matlab.ui.control.Label
        DetectButton                  matlab.ui.control.Button
        DetectionMethodsButtonGroup   matlab.ui.container.ButtonGroup
        AllpoleLPCButton              matlab.ui.control.RadioButton
        GoertzelAlgorithmButton       matlab.ui.control.RadioButton
        FFTButton                     matlab.ui.control.RadioButton
        Hz697                         matlab.ui.control.Label
        Hz770                         matlab.ui.control.Label
        Hz852                         matlab.ui.control.Label
        Hz941                         matlab.ui.control.Label
        Hz1209                        matlab.ui.control.Label
        Hz1336                        matlab.ui.control.Label
        Hz1477                        matlab.ui.control.Label
        Hz1633                        matlab.ui.control.Label
        ButtonD                       matlab.ui.control.StateButton
        ButtonPound                   matlab.ui.control.StateButton
        Button0                       matlab.ui.control.StateButton
        ButtonStar                    matlab.ui.control.StateButton
        ButtonC                       matlab.ui.control.StateButton
        Button9                       matlab.ui.control.StateButton
        Button8                       matlab.ui.control.StateButton
        Button7                       matlab.ui.control.StateButton
        ButtonB                       matlab.ui.control.StateButton
        Button6                       matlab.ui.control.StateButton
        Button5                       matlab.ui.control.StateButton
        Button4                       matlab.ui.control.StateButton
        ButtonA                       matlab.ui.control.StateButton
        Button3                       matlab.ui.control.StateButton
        Button2                       matlab.ui.control.StateButton
        InputLabel                    matlab.ui.control.Label
        Button1                       matlab.ui.control.StateButton
        InEditField                   matlab.ui.control.EditField
        GenerateButton                matlab.ui.control.Button
        GenerationMethodsButtonGroup  matlab.ui.container.ButtonGroup
        DirectDigitalFrequencySynthesisButton  matlab.ui.control.RadioButton
        DigitalOscillatorButton       matlab.ui.control.RadioButton
        BuiltInSinusoidalButton       matlab.ui.control.RadioButton
        ClrButton                     matlab.ui.control.Button
        Button_D                      matlab.ui.control.Button
        Button_C                      matlab.ui.control.Button
        Button_B                      matlab.ui.control.Button
        Button_A                      matlab.ui.control.Button
        Button_pound                  matlab.ui.control.Button
        Button_0                      matlab.ui.control.Button
        Button_star                   matlab.ui.control.Button
        Button_9                      matlab.ui.control.Button
        Button_8                      matlab.ui.control.Button
        Button_7                      matlab.ui.control.Button
        Button_6                      matlab.ui.control.Button
        Button_5                      matlab.ui.control.Button
        Button_4                      matlab.ui.control.Button
        Button_3                      matlab.ui.control.Button
        Button_2                      matlab.ui.control.Button
        Button_1                      matlab.ui.control.Button
        ReferenceUIAxes               matlab.ui.control.UIAxes
        DecodeUIAxes                  matlab.ui.control.UIAxes
        SpectrumUIAxes                matlab.ui.control.UIAxes
        GenerationUIAxes              matlab.ui.control.UIAxes
    end

    properties (Constant)
         F_ROW = [697, 770, 852, 941]; % low frequency in DTMF
         F_COL = [1209, 1336, 1477, 1633]; % high frequency in DTMF
         F = [697, 779, 852, 941, 1209, 1336, 1477, 1633]; % all the involved frequency
         F_sample = 8000; % Sampling Frequency    Default = 8000 cycles/sec
         Goe_N = 205; % the number of samples in Goertzel Algorithm
         N = 1600; % Tones exists 0.2s
         PAUSE_TIME = 800; % Gap between two keys is 0.1s
         MAP = ['1', '2', '3', 'A'; 
                '4', '5', '6', 'B'; 
                '7', '8', '9', 'C'; 
                '*', '0', '#', 'D']; % MAP stores the keypad
         sintablen = 1024; % the length of look-up table in DDFS
         lpc_rank = 12; % the order of lpc model
    end
    
    properties (Access = public)
         METHOD_GENERATION = '';
         METHOD_DETECTION = '';
         input_cached = ''; 
         input_current = '';
         x = [];
         t = [];
         total_t = [];
         pit = [];
         tone = [];
         index = [];
         SINTAB = [];
         K = [];
    end
    
    methods (Access = private)

        function plot_reference(app)
            a697=[1 -2*cos(2*pi*18/205) 1];
            a770=[1 -2*cos(2*pi*20/205) 1];
            a852=[1 -2*cos(2*pi*22/205) 1];
            a941=[1 -2*cos(2*pi*24/205) 1];
            a1209=[1 -2*cos(2*pi*31/205) 1];
            a1336=[1 -2*cos(2*pi*34/205) 1];
            a1477=[1 -2*cos(2*pi*38/205) 1];
            a1633=[1 -2*cos(2*pi*42/205) 1];
            
            [w1, ~]=freqz([1 -exp(-2*pi*18/205)],a697,512,app.F_sample);
            [w2, ~]=freqz([1 -exp(-2*pi*20/205)],a770,512,app.F_sample);
            [w3, ~]=freqz([1 -exp(-2*pi*22/205)],a852,512,app.F_sample);
            [w4, ~]=freqz([1 -exp(-2*pi*24/205)],a941,512,app.F_sample);
            [w5, ~]=freqz([1 -exp(-2*pi*31/205)],a1209,512,app.F_sample);
            [w6, ~]=freqz([1 -exp(-2*pi*34/205)],a1336,512,app.F_sample);
            [w7, ~]=freqz([1 -exp(-2*pi*38/205)],a1477,512,app.F_sample);
            [w8, f]=freqz([1 -exp(-2*pi*42/205)],a1633,512,app.F_sample);
            plot(app.ReferenceUIAxes,f,abs(w1)/1000,f,abs(w2)/1000,f,abs(w3)/1000,f,abs(w4)/1000,f,abs(w5)/1000,f,abs(w6)/1000,f,abs(w7)/1000,f,abs(w8)/1000);
        end

        function dtmf = super_position(app, f1, f2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:                             %
            %   f1 - the low frequency in DTMF   %
            %   f2 - the high frequency in DTMF  %
            % Output:                            %
            %   dtmf - the DTMF signal in domain %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dtmf = sin(f1 * app.pit) + sin(f2 * app.pit);
        end
        
        function dtmf = digital_oscillator(app, f1, f2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:                             %
            %   f1 - the low frequency in DTMF   %
            %   f2 - the high frequency in DTMF  %
            % Output:                            %
            %   dtmf - the DTMF signal in domain %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            dtmf = filter([0 sin(2 * pi * f1 / app.F_sample)], [1 -2 * cos(2 * pi * f1 / app.F_sample) 1], app.x)...
                 + filter([0 sin(2 * pi * f2 / app.F_sample)], [1 -2 * cos(2 * pi * f2 / app.F_sample) 1], app.x);
            % x is the unit impulse signal
        end

        function dtmf = direct_digital_frequency_synthesis(app, f1, f2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input:                             %
            %   f1 - the low frequency in DTMF   %
            %   f2 - the high frequency in DTMF  %
            % Output:                            %
            %   dtmf - the DTMF signal in domain %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            step1 = (f1/app.F_sample) * app.sintablen; % calculate the tuning word for f1
            step2 = (f2/app.F_sample) * app.sintablen; % calculate the tuning word for f2
            index1 = mod(round(app.index .* step1), app.sintablen) + 1; % calculate the address for f1
            index2 = mod(round(app.index .* step2), app.sintablen) + 1; % calculate the address for f2
            dtmf = app.SINTAB(index1) + app.SINTAB(index2);
        end
        
        function dtmf = generate(app, f1, f2)
            
            switch app.METHOD_GENERATION
                case 'BuiltInFunction'
                    dtmf = super_position(app, f1, f2);
                case 'DigitalOscillator'
                    dtmf = digital_oscillator(app, f1, f2);
                case 'DirectDigitalFrequencySynthesis'
                    dtmf = direct_digital_frequency_synthesis(app, f1, f2);
            end
        end

        function dtmf_noise = GaussNoise(app, dtmf)
            % Input:
            %   dtmf - dtmf signal in time domain
            % Output: 
            %   dtmf_noise - dtmf signal plus a noise
            noise = randn(size(dtmf));           
            len = length(dtmf);                     
            signal_power = 1 / len * sum(dtmf.*dtmf);   
            noise_power = 1 / len * sum(noise.*noise);
            % app.SNRSlider.Value is the SNR that we can adjust
            noise_variance = signal_power / (10^(app.SNRSlider.Value/10)); 
            noise = sqrt(noise_variance / noise_power) * noise;     
            dtmf_noise = dtmf + noise;  
        end
        
        function char = find_char(app, f1, f2)
            frow = abs(app.F_ROW - f1);
            fcol = abs(app.F_COL - f2); 
            row = find(frow == min(frow));
            col = find(fcol == min(fcol));
            char = app.MAP(row, col);
        end

        function [char, f1, f2] = FourierTransform(app, dtmf)
            % generate unilateral spectrum of DTMF signal
            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            
            % plot the spectrum
            f = (0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.DecodeUIAxes, f, dtmf_spectrum);
            app.DecodeUIAxes.XTick = [0 500 1000 1500 2000];
            
            % find the frequencies corresponding to two peak values
            [~, ind] = findpeaks(dtmf_spectrum,'SortStr','descend', 'MinPeakHeight', 0.2);
            f1 = min(f(ind(1)),f(ind(2)));
            f2 = max(f(ind(1)),f(ind(2)));
            
            % find_char is a self-defined function return corresponding 
            % character according to given f1 and f2
            char = find_char(app, f1, f2);
        end
        
        function [char, f1, f2] = Goertzel(app, dtmf)
            % sample the first N signal points
            dtmf = dtmf(1 : app.Goe_N); 

            % apply goertzel algorithm
            goe_dtmf = goertzel(dtmf, app.K); 
            
            % plot the spectrum
            stem(app.DecodeUIAxes, app.F, abs(goe_dtmf)); 
            app.DecodeUIAxes.XTick = app.F;
            app.DecodeUIAxes.XTickLabelRotation = 80;

            % find f1 and f2
            [~, f1] = max(goe_dtmf(1 : 4));
            [~, f2] = max(goe_dtmf(5 : 8));
            f1 = app.F_ROW(f1);
            f2 = app.F_COL(f2);

            % find the corresponding character
            char = find_char(app, f1, f2);
        end
        
        function [char, f1, f2] = LPC(app, dtmf)
            % LPC coefficients
            [a, ~] = lpc(dtmf, app.lpc_rank); 

            % generate frequency vector
            f = (0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf); 

            % generate frequency response of the LPC filter
            dtmf_spectrum = abs(freqz(1, a, f, app.F_sample));

            % plot the spectrum of LPC
            plot(app.DecodeUIAxes, f, dtmf_spectrum);
            app.DecodeUIAxes.XTick = [0 500 1000 1500 2000];
            
            % find the frequencies corresponding to two peak values
            [~, ind] = findpeaks(dtmf_spectrum,'SortStr','descend', 'MinPeakHeight', 2);
            f1 = min(f(ind(1)),f(ind(2)));
            f2 = max(f(ind(1)),f(ind(2)));

            % find the corresponding character
            char = find_char(app, f1, f2);
        end

        function [char, f1, f2] = detect(app, dtmf)
            
            switch app.METHOD_DETECTION
                case 'FourierTransform'
                    [char, f1, f2] = FourierTransform(app, dtmf);
                case 'GoertzelAlgorithm'
                    [char, f1, f2] = Goertzel(app, dtmf);
                case 'AllpoleLPC'
                    [char, f1, f2] = LPC(app, dtmf);
            end
        end

        function Extinguish(app, char, f1, f2)
            frow = abs(app.F_ROW - f1);
            switch find(frow == min(frow))
                case 1
                    app.Hz697.FontColor = 'Black';
                case 2
                    app.Hz770.FontColor = 'Black';
                case 3
                    app.Hz852.FontColor = 'Black';
                case 4
                    app.Hz941.FontColor = 'Black';
            end
            
            fcol = abs(app.F_COL - f2); 
            switch find(fcol == min(fcol))
                case 1
                    app.Hz1209.FontColor = 'Black';
                case 2
                    app.Hz1336.FontColor = 'Black';
                case 3
                    app.Hz1477.FontColor = 'Black';
                case 4
                    app.Hz1633.FontColor = 'Black';
            end
            
            switch char
                case '1'
                    app.Button1.BackgroundColor = [0.96 0.96 0.96];
                case '2'
                    app.Button2.BackgroundColor = [0.96 0.96 0.96];
                case '3'
                    app.Button3.BackgroundColor = [0.96 0.96 0.96];
                case '4'
                    app.Button4.BackgroundColor = [0.96 0.96 0.96];
                case '5'
                    app.Button5.BackgroundColor = [0.96 0.96 0.96];
                case '6'
                    app.Button6.BackgroundColor = [0.96 0.96 0.96];
                case '7'
                    app.Button7.BackgroundColor = [0.96 0.96 0.96];
                case '8'
                    app.Button8.BackgroundColor = [0.96 0.96 0.96];
                case '9'
                    app.Button9.BackgroundColor = [0.96 0.96 0.96];
                case '0'
                    app.Button0.BackgroundColor = [0.96 0.96 0.96];
                case '#'
                    app.ButtonPound.BackgroundColor = [0.96 0.96 0.96];
                case '*'
                    app.ButtonStar.BackgroundColor = [0.96 0.96 0.96];
                case 'A'
                    app.ButtonA.BackgroundColor = [0.96 0.96 0.96];
                case 'B'
                    app.ButtonB.BackgroundColor = [0.96 0.96 0.96];
                case 'C'
                    app.ButtonC.BackgroundColor = [0.96 0.96 0.96];
                case 'D'
                    app.ButtonD.BackgroundColor = [0.96 0.96 0.96];
            end
        end

        function Display(app, char, f1, f2)
            frow = abs(app.F_ROW - f1);
            switch find(frow == min(frow))
                case 1
                    app.Hz697.FontColor = 'Red';
                case 2
                    app.Hz770.FontColor = 'Red';
                case 3
                    app.Hz852.FontColor = 'Red';
                case 4
                    app.Hz941.FontColor = 'Red';
            end
            
            fcol = abs(app.F_COL - f2); 
            switch find(fcol == min(fcol))
                case 1
                    app.Hz1209.FontColor = 'Red';
                case 2
                    app.Hz1336.FontColor = 'Red';
                case 3
                    app.Hz1477.FontColor = 'Red';
                case 4
                    app.Hz1633.FontColor = 'Red';
            end
            
            switch char
                case '1'
                    app.Button1.BackgroundColor = 'Red';
                case '2'
                    app.Button2.BackgroundColor = 'Red';
                case '3'
                    app.Button3.BackgroundColor = 'Red';
                case '4'
                    app.Button4.BackgroundColor = 'Red';
                case '5'
                    app.Button5.BackgroundColor = 'Red';
                case '6'
                    app.Button6.BackgroundColor = 'Red';
                case '7'
                    app.Button7.BackgroundColor = 'Red';
                case '8'
                    app.Button8.BackgroundColor = 'Red';
                case '9'
                    app.Button9.BackgroundColor = 'Red';
                case '0'
                    app.Button0.BackgroundColor = 'Red';
                case '#'
                    app.ButtonPound.BackgroundColor = 'Red';
                case '*'
                    app.ButtonStar.BackgroundColor = 'Red';
                case 'A'
                    app.ButtonA.BackgroundColor = 'Red';
                case 'B'
                    app.ButtonB.BackgroundColor = 'Red';
                case 'C'
                    app.ButtonC.BackgroundColor = 'Red';
                case 'D'
                    app.ButtonD.BackgroundColor = 'Red';
            end
        end

    end
    


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.t = (0:app.N-1) / app.F_sample;
            app.pit = 2 * pi * app.t;
            app.x = zeros(1, length(app.t));
            app.x(1) = 1;
            app.METHOD_GENERATION = 'BuiltInFunction';
            app.METHOD_DETECTION = 'FourierTransform';
            app.SINTAB = sin(2*pi*(0:app.sintablen-1)./app.sintablen);
            app.index = 0:app.N - 1;
            app.K = round(app.F / app.F_sample * app.Goe_N) + 1;
            plot_reference(app);
        end

        % Button pushed function: Button_1
        function Button_1Pushed(app, event)
            app.input_current = '1';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_2
        function Button_2Pushed(app, event)
            app.input_current = '2';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_3
        function Button_3Pushed(app, event)
            app.input_current = '3';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_4
        function Button_4Pushed(app, event)
            app.input_current = '4';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_5
        function Button_5Pushed(app, event)
            app.input_current = '5';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_6
        function Button_6Pushed(app, event)
            app.input_current = '6';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_7
        function Button_7Pushed(app, event)
            app.input_current = '7';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_8
        function Button_8Pushed(app, event)
            app.input_current = '8';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_9
        function Button_9Pushed(app, event)
            app.input_current = '9';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_0
        function Button_0Pushed(app, event)
            app.input_current = '0';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_A
        function Button_APushed(app, event)
            app.input_current = 'A';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_B
        function Button_BPushed(app, event)
            app.input_current = 'B';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_C
        function Button_CPushed(app, event)
            app.input_current = 'C';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_D
        function Button_DPushed(app, event)
            app.input_current = 'D';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Selection changed function: GenerationMethodsButtonGroup
        function GenerationMethodsButtonGroupSelectionChanged(app, event)
            selectedButton = app.GenerationMethodsButtonGroup.SelectedObject;
            switch selectedButton
                case app.BuiltInSinusoidalButton
                    app.METHOD_GENERATION = 'BuiltInFunction';
                case app.DigitalOscillatorButton
                    app.METHOD_GENERATION = 'DigitalOscillator';
                case app.DirectDigitalFrequencySynthesisButton
                    app.METHOD_GENERATION = 'DirectDigitalFrequencySynthesis';
            end
        end

        % Selection changed function: DetectionMethodsButtonGroup
        function DetectionMethodsButtonGroupSelectionChanged(app, event)
            selectedButton = app.DetectionMethodsButtonGroup.SelectedObject;
            switch selectedButton
                case app.FFTButton
                    app.METHOD_DETECTION = 'FourierTransform';
                case app.GoertzelAlgorithmButton
                    app.METHOD_DETECTION = 'GoertzelAlgorithm';
                case app.AllpoleLPCButton
                    app.METHOD_DETECTION = 'AllpoleLPC';
            end
        end

        % Button pushed function: GenerateButton
        function GenerateButtonPushed(app, event)
            time = (0 : app.N * 3 / 2 * length(app.input_cached) - 1) ./ app.F_sample;
            cla(app.GenerationUIAxes);
            time_maximum = 0.3 * length(app.input_cached);
            plot(app.GenerationUIAxes, time, app.tone);
            ttick = (0:0.03:0.30) .* length(app.input_cached);
            app.GenerationUIAxes.XTick = ttick;
            app.GenerationUIAxes.XTickLabel = ttick;
            axis(app.GenerationUIAxes, [0 time_maximum -1 1]);
        end

        % Button pushed function: DetectButton
        function DetectButtonPushed(app, event)
            len = app.N + app.PAUSE_TIME;
            
            for i = 0 : ((length(app.tone) / len) - 1)               
                st = i * len;
                ed = (i + 1) * len - len / 3;
                
                current = app.tone(st + 1 : ed);
                [char, f1, f2] = detect(app, current);
                app.OutputEditField.Value = [app.OutputEditField.Value, char];
                Display(app, char, f1, f2);
                pause(0.5);
                Extinguish(app, char, f1, f2);
                pause(0.5);
            end
        end

        % Button pushed function: ClrButton
        function ClrButtonPushed(app, event)
            cla(app.GenerationUIAxes);
            cla(app.DecodeUIAxes);
            cla(app.SpectrumUIAxes);
            app.input_current = '';
            app.input_cached = '';
            app.InEditField.Value = '';
            app.OutputEditField.Value = '';
            app.tone = [];
            app.GenerationUIAxes.XTick = (0:0.02:0.20);
            app.GenerationUIAxes.XTickLabel = (0:0.02:0.20);
            axis(app.GenerationUIAxes, [0 0.2 -1 1]);
        end

        % Button pushed function: Button_star
        function Button_starPushed(app, event)
            app.input_current = '*';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end

        % Button pushed function: Button_pound
        function Button_poundPushed(app, event)
            app.input_current = '#';
            app.input_cached = [app.input_cached, app.input_current];
            app.InEditField.Value = app.input_cached;
            [row, col] = find(app.MAP == app.input_current);
            [f1, f2] = deal(app.F_ROW(row), app.F_COL(col));
            
            dtmf = generate(app, f1, f2);
            dtmf = dtmf / max(abs(dtmf));
            dtmf = GaussNoise(app, dtmf);
            sound(dtmf, app.F_sample);

            plot(app.GenerationUIAxes, app.t, dtmf);
            pause_signal = zeros(1, app.PAUSE_TIME) + 0.01;
            pause_signal = GaussNoise(app, pause_signal);
            app.tone = [app.tone, dtmf, pause_signal];

            dtmf_spectrum = abs(fft(dtmf)) / length(dtmf) * 2;
            dtmf_spectrum(1) = dtmf_spectrum(1) / 2;
            dtmf_spectrum = dtmf_spectrum(1:(length(dtmf)+1)/2);
            f=(0:1:((length(dtmf)-1)/2))*app.F_sample/length(dtmf);
            plot(app.SpectrumUIAxes, f, dtmf_spectrum);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 842 667];
            app.UIFigure.Name = 'MATLAB App';

            % Create GenerationUIAxes
            app.GenerationUIAxes = uiaxes(app.UIFigure);
            title(app.GenerationUIAxes, 'Generated DTMF')
            xlabel(app.GenerationUIAxes, 'Time(seconds)')
            ylabel(app.GenerationUIAxes, 'Amplitude')
            zlabel(app.GenerationUIAxes, 'Z')
            app.GenerationUIAxes.XLim = [0 0.2];
            app.GenerationUIAxes.YLim = [-1 1];
            app.GenerationUIAxes.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2];
            app.GenerationUIAxes.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'; '0.20'};
            app.GenerationUIAxes.Position = [241 444 557 195];

            % Create SpectrumUIAxes
            app.SpectrumUIAxes = uiaxes(app.UIFigure);
            title(app.SpectrumUIAxes, 'Generated Spectrum')
            xlabel(app.SpectrumUIAxes, 'Frequency(Hz)')
            ylabel(app.SpectrumUIAxes, 'Amplitude')
            zlabel(app.SpectrumUIAxes, 'Z')
            app.SpectrumUIAxes.XLim = [0 2000];
            app.SpectrumUIAxes.Position = [249 233 279 198];

            % Create DecodeUIAxes
            app.DecodeUIAxes = uiaxes(app.UIFigure);
            title(app.DecodeUIAxes, 'Decoded Spectrum')
            xlabel(app.DecodeUIAxes, 'Frequency(Hz)')
            ylabel(app.DecodeUIAxes, 'Amplitude')
            zlabel(app.DecodeUIAxes, 'Z')
            app.DecodeUIAxes.XLim = [0 2000];
            app.DecodeUIAxes.Position = [519 233 279 198];

            % Create ReferenceUIAxes
            app.ReferenceUIAxes = uiaxes(app.UIFigure);
            title(app.ReferenceUIAxes, 'Reference Spectrum')
            xlabel(app.ReferenceUIAxes, 'Frequency(Hz)')
            ylabel(app.ReferenceUIAxes, 'Amplitude')
            zlabel(app.ReferenceUIAxes, 'Z')
            app.ReferenceUIAxes.XLim = [0 2000];
            app.ReferenceUIAxes.Position = [519 9 279 198];

            % Create Button_1
            app.Button_1 = uibutton(app.UIFigure, 'push');
            app.Button_1.ButtonPushedFcn = createCallbackFcn(app, @Button_1Pushed, true);
            app.Button_1.Position = [49 555 30 30];
            app.Button_1.Text = '1';

            % Create Button_2
            app.Button_2 = uibutton(app.UIFigure, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @Button_2Pushed, true);
            app.Button_2.Position = [78 555 30 30];
            app.Button_2.Text = '2';

            % Create Button_3
            app.Button_3 = uibutton(app.UIFigure, 'push');
            app.Button_3.ButtonPushedFcn = createCallbackFcn(app, @Button_3Pushed, true);
            app.Button_3.Position = [107 555 30 30];
            app.Button_3.Text = '3';

            % Create Button_4
            app.Button_4 = uibutton(app.UIFigure, 'push');
            app.Button_4.ButtonPushedFcn = createCallbackFcn(app, @Button_4Pushed, true);
            app.Button_4.Position = [49 526 30 30];
            app.Button_4.Text = '4';

            % Create Button_5
            app.Button_5 = uibutton(app.UIFigure, 'push');
            app.Button_5.ButtonPushedFcn = createCallbackFcn(app, @Button_5Pushed, true);
            app.Button_5.Position = [78 526 30 30];
            app.Button_5.Text = '5';

            % Create Button_6
            app.Button_6 = uibutton(app.UIFigure, 'push');
            app.Button_6.ButtonPushedFcn = createCallbackFcn(app, @Button_6Pushed, true);
            app.Button_6.Position = [107 526 30 30];
            app.Button_6.Text = '6';

            % Create Button_7
            app.Button_7 = uibutton(app.UIFigure, 'push');
            app.Button_7.ButtonPushedFcn = createCallbackFcn(app, @Button_7Pushed, true);
            app.Button_7.Position = [49 497 30 30];
            app.Button_7.Text = '7';

            % Create Button_8
            app.Button_8 = uibutton(app.UIFigure, 'push');
            app.Button_8.ButtonPushedFcn = createCallbackFcn(app, @Button_8Pushed, true);
            app.Button_8.Position = [78 497 30 30];
            app.Button_8.Text = '8';

            % Create Button_9
            app.Button_9 = uibutton(app.UIFigure, 'push');
            app.Button_9.ButtonPushedFcn = createCallbackFcn(app, @Button_9Pushed, true);
            app.Button_9.Position = [107 497 30 30];
            app.Button_9.Text = '9';

            % Create Button_star
            app.Button_star = uibutton(app.UIFigure, 'push');
            app.Button_star.ButtonPushedFcn = createCallbackFcn(app, @Button_starPushed, true);
            app.Button_star.Position = [49 468 30 30];
            app.Button_star.Text = '*';

            % Create Button_0
            app.Button_0 = uibutton(app.UIFigure, 'push');
            app.Button_0.ButtonPushedFcn = createCallbackFcn(app, @Button_0Pushed, true);
            app.Button_0.Position = [78 468 30 30];
            app.Button_0.Text = '0';

            % Create Button_pound
            app.Button_pound = uibutton(app.UIFigure, 'push');
            app.Button_pound.ButtonPushedFcn = createCallbackFcn(app, @Button_poundPushed, true);
            app.Button_pound.Position = [107 468 30 30];
            app.Button_pound.Text = '#';

            % Create Button_A
            app.Button_A = uibutton(app.UIFigure, 'push');
            app.Button_A.ButtonPushedFcn = createCallbackFcn(app, @Button_APushed, true);
            app.Button_A.Position = [136 555 30 30];
            app.Button_A.Text = 'A';

            % Create Button_B
            app.Button_B = uibutton(app.UIFigure, 'push');
            app.Button_B.ButtonPushedFcn = createCallbackFcn(app, @Button_BPushed, true);
            app.Button_B.Position = [136 526 30 30];
            app.Button_B.Text = 'B';

            % Create Button_C
            app.Button_C = uibutton(app.UIFigure, 'push');
            app.Button_C.ButtonPushedFcn = createCallbackFcn(app, @Button_CPushed, true);
            app.Button_C.Position = [136 497 30 30];
            app.Button_C.Text = 'C';

            % Create Button_D
            app.Button_D = uibutton(app.UIFigure, 'push');
            app.Button_D.ButtonPushedFcn = createCallbackFcn(app, @Button_DPushed, true);
            app.Button_D.Position = [136 468 30 30];
            app.Button_D.Text = 'D';

            % Create ClrButton
            app.ClrButton = uibutton(app.UIFigure, 'push');
            app.ClrButton.ButtonPushedFcn = createCallbackFcn(app, @ClrButtonPushed, true);
            app.ClrButton.Position = [49 584 116 29];
            app.ClrButton.Text = 'Clr';

            % Create GenerationMethodsButtonGroup
            app.GenerationMethodsButtonGroup = uibuttongroup(app.UIFigure);
            app.GenerationMethodsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @GenerationMethodsButtonGroupSelectionChanged, true);
            app.GenerationMethodsButtonGroup.TitlePosition = 'centertop';
            app.GenerationMethodsButtonGroup.Title = 'Generation Methods';
            app.GenerationMethodsButtonGroup.Position = [16 300 212 93];

            % Create BuiltInSinusoidalButton
            app.BuiltInSinusoidalButton = uiradiobutton(app.GenerationMethodsButtonGroup);
            app.BuiltInSinusoidalButton.Text = 'Built-In Sinusoidal';
            app.BuiltInSinusoidalButton.Position = [7 46 118 22];
            app.BuiltInSinusoidalButton.Value = true;

            % Create DigitalOscillatorButton
            app.DigitalOscillatorButton = uiradiobutton(app.GenerationMethodsButtonGroup);
            app.DigitalOscillatorButton.Text = 'Digital Oscillator';
            app.DigitalOscillatorButton.Position = [8 25 109 22];

            % Create DirectDigitalFrequencySynthesisButton
            app.DirectDigitalFrequencySynthesisButton = uiradiobutton(app.GenerationMethodsButtonGroup);
            app.DirectDigitalFrequencySynthesisButton.Text = 'Direct Digital Frequency Synthesis';
            app.DirectDigitalFrequencySynthesisButton.Position = [9 4 207 22];

            % Create GenerateButton
            app.GenerateButton = uibutton(app.UIFigure, 'push');
            app.GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateButtonPushed, true);
            app.GenerateButton.Position = [50 430 117 30];
            app.GenerateButton.Text = 'Generate';

            % Create InEditField
            app.InEditField = uieditfield(app.UIFigure, 'text');
            app.InEditField.Position = [73 617 92 22];

            % Create Button1
            app.Button1 = uibutton(app.UIFigure, 'state');
            app.Button1.Text = '1';
            app.Button1.Position = [308 152 30 30];

            % Create InputLabel
            app.InputLabel = uilabel(app.UIFigure);
            app.InputLabel.HorizontalAlignment = 'center';
            app.InputLabel.Position = [48 617 16 22];
            app.InputLabel.Text = 'In:';

            % Create Button2
            app.Button2 = uibutton(app.UIFigure, 'state');
            app.Button2.Text = '2';
            app.Button2.Position = [361 152 30 30];

            % Create Button3
            app.Button3 = uibutton(app.UIFigure, 'state');
            app.Button3.Text = '3';
            app.Button3.Position = [414 152 30 30];

            % Create ButtonA
            app.ButtonA = uibutton(app.UIFigure, 'state');
            app.ButtonA.Text = 'A';
            app.ButtonA.Position = [467 152 30 30];

            % Create Button4
            app.Button4 = uibutton(app.UIFigure, 'state');
            app.Button4.Text = '4';
            app.Button4.Position = [308 105 30 30];

            % Create Button5
            app.Button5 = uibutton(app.UIFigure, 'state');
            app.Button5.Text = '5';
            app.Button5.Position = [361 105 30 30];

            % Create Button6
            app.Button6 = uibutton(app.UIFigure, 'state');
            app.Button6.Text = '6';
            app.Button6.Position = [414 105 30 30];

            % Create ButtonB
            app.ButtonB = uibutton(app.UIFigure, 'state');
            app.ButtonB.Text = 'B';
            app.ButtonB.Position = [467 105 30 30];

            % Create Button7
            app.Button7 = uibutton(app.UIFigure, 'state');
            app.Button7.Text = '7';
            app.Button7.Position = [308 58 30 30];

            % Create Button8
            app.Button8 = uibutton(app.UIFigure, 'state');
            app.Button8.Text = '8';
            app.Button8.Position = [361 58 30 30];

            % Create Button9
            app.Button9 = uibutton(app.UIFigure, 'state');
            app.Button9.Text = '9';
            app.Button9.Position = [414 58 30 30];

            % Create ButtonC
            app.ButtonC = uibutton(app.UIFigure, 'state');
            app.ButtonC.Text = 'C';
            app.ButtonC.Position = [467 58 30 30];

            % Create ButtonStar
            app.ButtonStar = uibutton(app.UIFigure, 'state');
            app.ButtonStar.Text = '*';
            app.ButtonStar.Position = [308 11 30 30];

            % Create Button0
            app.Button0 = uibutton(app.UIFigure, 'state');
            app.Button0.Text = '0';
            app.Button0.Position = [361 11 30 30];

            % Create ButtonPound
            app.ButtonPound = uibutton(app.UIFigure, 'state');
            app.ButtonPound.Text = '#';
            app.ButtonPound.Position = [414 11 30 30];

            % Create ButtonD
            app.ButtonD = uibutton(app.UIFigure, 'state');
            app.ButtonD.Text = 'D';
            app.ButtonD.Position = [467 11 30 30];

            % Create Hz1633
            app.Hz1633 = uilabel(app.UIFigure);
            app.Hz1633.Position = [458 190 46 22];
            app.Hz1633.Text = '1633Hz';

            % Create Hz1477
            app.Hz1477 = uilabel(app.UIFigure);
            app.Hz1477.Position = [405 190 46 22];
            app.Hz1477.Text = '1477Hz';

            % Create Hz1336
            app.Hz1336 = uilabel(app.UIFigure);
            app.Hz1336.Position = [352 190 46 22];
            app.Hz1336.Text = '1336Hz';

            % Create Hz1209
            app.Hz1209 = uilabel(app.UIFigure);
            app.Hz1209.Position = [299 190 46 22];
            app.Hz1209.Text = '1209Hz';

            % Create Hz941
            app.Hz941 = uilabel(app.UIFigure);
            app.Hz941.Position = [258 15 40 22];
            app.Hz941.Text = '941Hz';

            % Create Hz852
            app.Hz852 = uilabel(app.UIFigure);
            app.Hz852.Position = [258 62 40 22];
            app.Hz852.Text = '852Hz';

            % Create Hz770
            app.Hz770 = uilabel(app.UIFigure);
            app.Hz770.Position = [258 109 40 22];
            app.Hz770.Text = '770Hz';

            % Create Hz697
            app.Hz697 = uilabel(app.UIFigure);
            app.Hz697.Position = [258 156 40 22];
            app.Hz697.Text = '697Hz';

            % Create DetectionMethodsButtonGroup
            app.DetectionMethodsButtonGroup = uibuttongroup(app.UIFigure);
            app.DetectionMethodsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @DetectionMethodsButtonGroupSelectionChanged, true);
            app.DetectionMethodsButtonGroup.TitlePosition = 'centertop';
            app.DetectionMethodsButtonGroup.Title = 'Detection Methods';
            app.DetectionMethodsButtonGroup.Position = [16 156 212 93];

            % Create FFTButton
            app.FFTButton = uiradiobutton(app.DetectionMethodsButtonGroup);
            app.FFTButton.Text = 'FFT';
            app.FFTButton.Position = [7 46 44 22];
            app.FFTButton.Value = true;

            % Create GoertzelAlgorithmButton
            app.GoertzelAlgorithmButton = uiradiobutton(app.DetectionMethodsButtonGroup);
            app.GoertzelAlgorithmButton.Text = 'Goertzel Algorithm';
            app.GoertzelAlgorithmButton.Position = [8 25 121 22];

            % Create AllpoleLPCButton
            app.AllpoleLPCButton = uiradiobutton(app.DetectionMethodsButtonGroup);
            app.AllpoleLPCButton.Text = 'All-pole LPC';
            app.AllpoleLPCButton.Position = [9 4 207 22];

            % Create DetectButton
            app.DetectButton = uibutton(app.UIFigure, 'push');
            app.DetectButton.ButtonPushedFcn = createCallbackFcn(app, @DetectButtonPushed, true);
            app.DetectButton.Position = [50 401 117 30];
            app.DetectButton.Text = 'Detect';

            % Create OutputLabel
            app.OutputLabel = uilabel(app.UIFigure);
            app.OutputLabel.HorizontalAlignment = 'center';
            app.OutputLabel.Position = [45 131 44 22];
            app.OutputLabel.Text = 'Output:';

            % Create OutputEditField
            app.OutputEditField = uieditfield(app.UIFigure, 'text');
            app.OutputEditField.Position = [96 131 102 22];

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [55 84 131 22];
            app.Label.Text = '2020190502030 刘代宸';

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.Position = [55 63 131 22];
            app.Label_2.Text = '2020190502017 李敬元';

            % Create Label_3
            app.Label_3 = uilabel(app.UIFigure);
            app.Label_3.Position = [55 42 131 22];
            app.Label_3.Text = '2020190501035 丁麒晔';

            % Create Label_4
            app.Label_4 = uilabel(app.UIFigure);
            app.Label_4.Position = [55 21 131 22];
            app.Label_4.Text = '2020190501019 王博禹';

            % Create SNRSliderLabel
            app.SNRSliderLabel = uilabel(app.UIFigure);
            app.SNRSliderLabel.HorizontalAlignment = 'right';
            app.SNRSliderLabel.Position = [19 270 31 22];
            app.SNRSliderLabel.Text = 'SNR';

            % Create SNRSlider
            app.SNRSlider = uislider(app.UIFigure);
            app.SNRSlider.Limits = [0 50];
            app.SNRSlider.Position = [71 279 150 3];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DTMF_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end