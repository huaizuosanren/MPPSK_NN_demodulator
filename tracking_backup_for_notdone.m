function [trackResults, channel]= tracking(data, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%                    - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $
% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, 10);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, 10);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, 10);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, 10);
trackResults.I_E            = zeros(1, 10);
trackResults.I_L            = zeros(1, 10);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, 10);
trackResults.Q_P            = zeros(1, 10);
trackResults.Q_L            = zeros(1, 10);

% Loop discriminators
trackResults.dllDiscr       = inf(1, 10);
trackResults.dllDiscrFilt   = inf(1, 10);
trackResults.pllDiscr       = inf(1, 10);
trackResults.pllDiscrFilt   = inf(1, 10);
trackResults.fllDiscr       = inf(1, 10);
min3000Address              = 0;
%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);%对不同数量的信道进行排序

%///////////////////The following variables was used to initialized SNR////////////////////////
KKNWPR                         =   0;
WBPI                           =   0;
WBPQ                           =   0;
NBPI                           =   0;
NBPQ                           =   0;
%timeSamples                    =  floor(settings.msToProcess/20);
timeSamplesCount               =   1;
estimation                     =   20;
%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);%对不同数量的信道进行排序



hwb = waitbar(0,'Synchornizing...');
%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels
    
    channel(channelNr).acquiredFreq=settings.IF+80;
    % Only process if PRN is non zero (acquisition was successful)
        % Save additional information - each channel's tracked PRN
        
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample) 
       % fseek(fid, ...
       %       settings.skipNumberOfBytes + channel(channelNr).codePhase-1, ...
       %       'bof');


        % Get a vector with the C/A code sampled 1x/chip
        %caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        %caCode = [caCode(1023) caCode caCode(1)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        %carrier/Costas loop parameters
%        oldCarrNco   = 0.0;
%        oldCarrError = 0.0;
        
        old_d_pll_x = 0;
        old_d_pll_w = 0;
        
        old_d_dll_w = 0;
        %% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;        %单位是秒

% Calculate filter coefficient values都没使用呀，估计以后调参的时候用吧
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
%[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
%                                    settings.pllDampingRatio, ...
%                                    0.25);

% Parameters for 2nd-order FLL aided 3rd-order PLL
d_pll_b3 = 2.400; 
d_pll_a3 = 1.100; 
d_pll_a2 = 1.414; 
d_pll_w0p = settings.pllNoiseBandwidth / 0.7845; 
d_pll_w0p2 = d_pll_w0p * d_pll_w0p; 
d_pll_w0p3 = d_pll_w0p2 * d_pll_w0p; 

d_pll_w0f = settings.fllNoiseBandwidth / 0.53; 
d_pll_w0f2 = d_pll_w0f * d_pll_w0f;
% Parameters for 1st-order FLL aided 2nd-order PLL
% d_pll_a2 = 1.414;
% d_pll_w0p = settings.pllNoiseBandwidth / 0.53;
% d_pll_w0p2 = d_pll_w0p * d_pll_w0p;
% d_pll_w0f = settings.fllNoiseBandwidth / 0.25;

% Parameters for 2nd-order DLL
d_dll_a2 = 1.414;
d_dll_w0p = settings.dllNoiseBandwidth / 0.53;
d_dll_w0p2 = d_dll_w0p * d_dll_w0p;

old_I_P = 1;
old_Q_P = 1;


presentAddress                 =   0;
basicAddressPositive           =   0;
basicAddressNegative           =   0;
FF                             =   0;
loopCount                      =  20;

%get the location of data
remblkloc=0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1 : codePeriods
            
%% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 50) == 0)
                try
                    waitbar(loopCnt/codePeriods, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; Completed ',int2str(loopCnt), ...
                            ' of ', int2str(codePeriods), ' msec']);                       
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

%% Read next block of data ------------------------------------------------            
            % Find the size of a "block" or code period in whole samples
            
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = settings.codeFreqBasis / settings.samplingFreq;
            
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this
            % interation 
            %[rawSignal, samplesRead] = fread(fid, ...
            %                                 blksize, settings.dataType);
            %rawSignal = rawSignal';  %transpose vector
            
            rawSignal=zeros(1,blksize);
            rawSignal(1:blksize)= data(remblkloc+1:remblkloc+blksize);

            %remember which data was read 
            remblkloc=remblkloc+blksize;
            
            
            % If did not read in enough samples, then could be out of 
            % data - better exit 
           % if (samplesRead ~= blksize)
           %     disp('Not able to read the specified number of samples  for tracking, exiting!')
           %     fclose(fid);
           %    return
           % end

%% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            %tcode       = (remCodePhase-earlyLateSpc) : ...
            %              codePhaseStep : ...
            %              ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            %tcode2      = ceil(tcode) + 1;
            tcode2_length= blksize;
            %earlyCode   = caCode(tcode2);
            
            % Define index into late code vector
            %tcode       = (remCodePhase+earlyLateSpc) : ...
            %              codePhaseStep : ...
            %              ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            %tcode2      = ceil(tcode) + 1;
            %lateCode    = caCode(tcode2);
            
            % Define index into prompt code vector
            %tcode       = remCodePhase : ...
             %             codePhaseStep : ...
              %            ((blksize-1)*codePhaseStep+remCodePhase);
            %tcode2      = ceil(tcode) + 1;
            %promptCode  = caCode(tcode2);
            promptCode= ones(1,tcode2_length); 
           % remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = carrCos .* rawSignal;
            iBasebandSignal = carrSin .* rawSignal;

            % Now get early, late, and prompt values for each
            %I_E = sum(earlyCode  .* iBasebandSignal);
            %Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Q_P = sum(promptCode .* qBasebandSignal);
            %I_L = sum(lateCode   .* iBasebandSignal);
            %Q_L = sum(lateCode   .* qBasebandSignal);
            trackResults.I_PMatrix(loopCnt)       =    I_P;   
            trackResults.Q_PMatrix(loopCnt)       =    Q_P;   
            
            
%% Try to make the I_P bit synchornized ----------------------------------
            if (I_P>0)&&(loopCount~=0)
                basicAddressPositive = basicAddressPositive + 1;
                presentAddress       = presentAddress       + 1;
                if (mod(basicAddressNegative,20)~=0)&&(basicAddressNegative~=0)
                    basicAddressNegative   =   0;
                    loopCount              =   20;
                elseif((basicAddressNegative~=0)&&mod(basicAddressPositive,20)==0)
                    if(mod(basicAddressPositive,20)==0)
                        basicAddressPositive=0;
                    end
                    loopCount              =   loopCount    -   1;
                elseif((basicAddressPositive~=0)&&mod(basicAddressPositive,20)==0)
                    if(mod(basicAddressPositive,20)==0)
                        basicAddressPositive=0;
                    end
                    loopCount              =   loopCount    -   1;
                end
            elseif (I_P<0)&&(loopCount~=0)
                basicAddressNegative = basicAddressNegative + 1;
                presentAddress       = presentAddress       + 1;
                if(mod(basicAddressPositive,20)~=0)&&(basicAddressPositive~=0)
                    basicAddressPositive   =   0;
                    loopCount              =   20;
                elseif((basicAddressPositive~=0)&&mod(basicAddressNegative,20)==0)
                    if(mod(basicAddressNegative,20)==0)
                        basicAddressNegative=0;
                    end
                    loopCount              =   loopCount     -  1;
                elseif((basicAddressNegative~=0)&&mod(basicAddressNegative,20)==0)
                    if(mod(basicAddressNegative,20)==0)
                        basicAddressNegative=0;
                    end
                    loopCount              =   loopCount     -  1;
                end
            end
            if(loopCount==0)&&(FF==0)
                presentAddress    =  presentAddress     -    399;
                FF=1;
            end
               %trackResults(channelNr).basicAddressNegative(loopCnt) = basicAddressNegative;
               %trackResults(channelNr).basicAddressPositive(loopCnt) = basicAddressPositive;
               
               
%% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            PLL_discriminator = atan(Q_P / I_P) / (2.0 * pi);
            
            dot   = old_I_P * I_P + old_Q_P * Q_P;
            cross = Q_P * old_I_P - I_P * old_Q_P; 
            
            %FLL_discriminator = atan2(cross , dot) / (2.0 * pi * PDIcarr);
            FLL_discriminator = atan(cross / dot) / (2.0 * pi * PDIcarr);
%             if(loopCnt > 2000)
%                 if(FLL_discriminator < -30)
%                     FLL_discriminator = -30;
%                 elseif(FLL_discriminator > 30)
%                     FLL_discriminator = 30;
%                 end
%             end
            old_I_P = I_P;
            old_Q_P = Q_P;
%             
%             % Implement carrier loop filter and generate NCO command
%             carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
%                 (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
%             oldCarrNco   = carrNco;
%             oldCarrError = carrError;
% 
%             % Modify carrier freq based on NCO command
%             carrFreq = carrFreqBasis + carrNco;
% 
%             trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

%%%%&&&&& 3rd-order PLL %%%%%%%%%%
%%%%             d_pll_w = d_pll_w + PDIcode * (d_pll_w0p3 * carrError + d_pll_w0f2 * FLL_discriminator); 
%%%%             d_pll_x = d_pll_x + PDIcode * (0.5*d_pll_w + d_pll_a2 * d_pll_w0f * FLL_discriminator + d_pll_a3 * d_pll_w0p2 * carrError); 
            d_pll_w = old_d_pll_w + PDIcarr * (d_pll_w0p3 * PLL_discriminator + d_pll_w0f2 * FLL_discriminator); 
            d_pll_x = old_d_pll_x + PDIcarr * (0.5 * (d_pll_w + old_d_pll_w) + d_pll_a2 * d_pll_w0f * FLL_discriminator + d_pll_a3 * d_pll_w0p2 * PLL_discriminator); 
            carrNco  = 0.5 * (d_pll_x + old_d_pll_x) + d_pll_b3 * d_pll_w0p * PLL_discriminator; 
            carrFreq = carrFreqBasis + carrNco;
            old_d_pll_w = d_pll_w;
            old_d_pll_x = d_pll_x;
            %% chongji filter to get the baseband signals
            
            if loopcnt>=400
                
            end
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%&&&&& 2nd-order PLL %%%%%%%%%%
%             d_pll_w = old_d_pll_w + PLL_discriminator * d_pll_w0p2 * PDIcarr + FLL_discriminator * d_pll_w0f * PDIcarr;
%             carrNco = 0.5 * (d_pll_w + old_d_pll_w) + d_pll_a2 * d_pll_w0p * PLL_discriminator;
%             old_d_pll_w = d_pll_w;
%             carrFreq = carrFreqBasis + carrNco;
            

            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;


%% Find DLL error and update code NCO -------------------------------------
%%no DLL loop -----------------------------------------           
            
            if  loopCnt >1000
                if trackResults.I_PMatrix(loopCnt) * trackResults.I_PMatrix(loopCnt-1) < 0
                    min3000Address   =   loopCnt;
                end
            end
            
%% Record various measures to show in postprocessing ----------------------
            % Record sample number (based on 8bit samples)
            %trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid);

           % trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
           % trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = PLL_discriminator;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            trackResults(channelNr).fllDiscr(loopCnt)       = FLL_discriminator;
            trackResults(channelNr).fll_cross(loopCnt)      = cross;
            trackResults(channelNr).fll_dot(loopCnt)        = dot;

            %trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
           % trackResults(channelNr).I_L(loopCnt) = I_L;
            %trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            %trackResults(channelNr).Q_L(loopCnt) = Q_L;
            if  I_P>30
                trackResults(channelNr).I_Ptoone(loopCnt)=1;
            elseif I_P<-30
                trackResults(channelNr).I_Ptoone(loopCnt)=-1;
            end
       
        end % for loopCnt

end
        close(hwb)
        

hwb = waitbar(0,'Tracking...');
%% Initialize tracking variables ==========================================



PDIcode = 0.100;%0.875; %码环

PDIcarr = 0.02; %载波环

originalAddress    =             1;
timeCoherent       =             20;
I_E1 = 0;
Q_E1 = 0;
I_P1 = 0;
Q_P1 = 0;
I_L1 = 0;
Q_L1 = 0;

settings.pllNoiseBandwidth       = 0;%12;     %[Hz]
settings.fllNoiseBandwidth       = 1.56402;%1.5;%7;      %[Hz]
% Parameters for 2nd-order FLL aided 3rd-order PLL
d_pll_w0p = settings.pllNoiseBandwidth / 0.7845; 
d_pll_w0f = settings.fllNoiseBandwidth / 0.53; 






d_pll_w0p2 = d_pll_w0p * d_pll_w0p; 
d_pll_w0p3 = d_pll_w0p2 * d_pll_w0p; 

d_pll_w0f2 = d_pll_w0f * d_pll_w0f;
% Parameters for 1st-order FLL aided 2nd-order PLL
% d_pll_a2 = 1.414;
% d_pll_w0p = settings.pllNoiseBandwidth / 0.53;
% d_pll_w0p2 = d_pll_w0p * d_pll_w0p;
% d_pll_w0f = settings.fllNoiseBandwidth / 0.25;

% Parameters for 2nd-order DLL
d_dll_a2 = 1.414;
d_dll_w0p = settings.dllNoiseBandwidth / 0.53;
d_dll_w0p2 = d_dll_w0p * d_dll_w0p;





controlVariable    =   timeCoherent/estimation;

       
a                                     =  size (trackResults.I_E);
trackResults(channelNr).numbers       =            a(2);

        %% NWPR C/NO Estimator ----------------------------------------------------
presentAddressEnd   =   floor((settings.msToProcess-presentAddress)/20);
nwpr_cnr  =  zeros(1,presentAddressEnd);
for countNWPR        =    presentAddress  :  25000
    if timeSamplesCount<presentAddressEnd
        KKNWPR        =        KKNWPR   +    1;
        WBPI          =        WBPI     +    trackResults.I_P(countNWPR)^2;
        WBPQ          =        WBPQ     +    trackResults.Q_P(countNWPR)^2;
        NBPI          =        NBPI     +    trackResults.I_P(countNWPR);
        NBPQ          =        NBPQ     +    trackResults.Q_P(countNWPR);
        while KKNWPR == 20
            NBP       =        NBPI^2   +    NBPQ^2;
            WBP       =        WBPI     +    WBPQ;
            NP        =        NBP      /    WBP;
            nwpr_snr  =        (NP - 1) / (20 - NP);
            nwpr_cnr(timeSamplesCount) = abs(10*log10(nwpr_snr) + 10*log10(settings.samplingFreq/2)-10*log10(10230));
            if NP >= 20 || NP <= 1
                countNWPR           =    countNWPR        -   1;
                presentAddressEnd   =  presentAddressEnd  -   1;
            end
            KKNWPR           =        0;
            WBPI             =        0;
            WBPQ             =        0;
            WBP              =        0;
            NBPI             =        0;
            NBPQ             =        0;
            timeSamplesCount = timeSamplesCount + 1;
        end      
    end
end
figure (301);
plot(nwpr_cnr)
end

% Close the waitbar
%close(hwb)
