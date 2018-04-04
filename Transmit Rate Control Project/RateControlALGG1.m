
%% Waveform Configuration

cfgVHT = wlanVHTConfig;         
cfgVHT.ChannelBandwidth = 'CBW40'; % 40 MHz channel bandwidth
cfgVHT.MCS = 1;                    % QPSK rate-1/2
cfgVHT.APEPLength = 4096;          % APEP length in bytes

% Set random stream for repeatability of results
s = rng(21);

%% Channel Configuration

tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-D';
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth;
tgacChannel.NumTransmitAntennas = 1;
tgacChannel.NumReceiveAntennas = 1;
tgacChannel.TransmitReceiveDistance = 3; % Distance in meters for NLOS
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 1794;

% Set the sampling rate for the channel
sr = wlanSampleRate(cfgVHT);
tgacChannel.SampleRate = sr;

%% Rate Control Algorithm Parameters
%Setting the maximum supported MCS, depending on the chosen bandwidth
BW = sscanf(cfgVHT.ChannelBandwidth,'CBW%d');
if BW < 40
    maxMCS = 8;
else
    maxMCS = 9;
end
tempMCS = cfgVHT.MCS;
minMcsLookup = [0 10 15 17 21 25 27 29 32 35];
mcsLookup = [0 10 15 17 21 25 27 29 32 35];
changeLookupCounter = 0;
packetsCounter = 0;
faultyBound = 1;
sparseFaultyBound = 8;
changeLookupCounterBound = 2;
packetsCounterBound = 5;
snrDecreasingBound = 2;
sparseFaultyBoundCheck = 2;
faultyBoundCheck = 1;
deviation = 3 * 0.01;

%% Simulation Parameters

numPackets = 100; % Number of packets transmitted during the simulation 
walkSNR = true; 

% Select SNR for the simulation
if walkSNR
    meanSNR = 26;   % Mean SNR
    amplitude = 10; % Variation in SNR around the average mean SNR value
    % Generate varying SNR values for each transmitted packet
    baseSNR = sin(linspace(1,10,numPackets))*amplitude+meanSNR;
    snrWalk = baseSNR(1); % Set the initial SNR value
    % The maxJump controls the maximum SNR difference between one
    % packet and the next 
    maxJump = 1;
else
    % Fixed mean SNR value for each transmitted packet. All the variability
    % in SNR comes from a time varying radio channel
    snrWalk = 22; %#ok<UNRCH>
end

% To plot the equalized constellation for each spatial stream set
% displayConstellation to true
displayConstellation = false;
if displayConstellation
    ConstellationDiagram = comm.ConstellationDiagram; %#ok<UNRCH>
    ConstellationDiagram.ShowGrid = true;
    ConstellationDiagram.Name = 'Equalized data symbols';
end

% Define simulation variables
snrMeasured = zeros(1,numPackets);
MCS = zeros(1,numPackets);
ber = zeros(1,numPackets);
packetLength = zeros(1,numPackets);


%% Processing Chain

for numPkt = 1:numPackets 
    if walkSNR
        % Generate SNR value per packet using random walk algorithm biased
        % towards the mean SNR
        snrWalk = 0.9*snrWalk+0.1*baseSNR(numPkt)+rand(1)*maxJump*2-maxJump;
    end
    
    % Generate a single packet waveform
    txPSDU = randi([0,1],8*cfgVHT.PSDULength,1,'int8');
    txWave = wlanWaveformGenerator(txPSDU,cfgVHT,'IdleTime',5e-4);
    
    % Receive processing, including SNR estimation
    y = processPacket(txWave,snrWalk,tgacChannel,cfgVHT);
    
    % Plot equalized symbols of data carrying subcarriers
    if displayConstellation && ~isempty(y.EstimatedSNR)
        release(ConstellationDiagram);
        ConstellationDiagram.ReferenceConstellation = helperReferenceSymbols(cfgVHT);
        ConstellationDiagram.Title = ['Packet ' int2str(numPkt)];
        ConstellationDiagram(y.EqDataSym(:));
        drawnow 
    end
    
    % Store estimated SNR value for each packet
    if isempty(y.EstimatedSNR) 
        snrMeasured(1,numPkt) = 0;
    else
        snrMeasured(1,numPkt) = y.EstimatedSNR;
    end
    
    % Determine the length of the packet in seconds including idle time
    packetLength(numPkt) = y.RxWaveformLength/sr;
    
    % Calculate packet error rate (PER)
    if isempty(y.RxPSDU)
        % Set the PER of an undetected packet to NaN
        ber(numPkt) = NaN;
    else
        [~,ber(numPkt)] = biterr(y.RxPSDU,txPSDU);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compare the estimated SNR to the threshold, and adjust the MCS value
    % used for the next packet
    faulty = 0;
    sparseFaulty = 0;
    snrDecreasingCounter = 0;
    
    if numPkt > 1
        ind = MCS(numPkt - 1) + 1;
    else
        ind = cfgVHT.MCS + 1;
    end

    %Implementation of the lookup table. Depending on the estimated SNR
    %the corresponding MCS is chosen
    if y.EstimatedSNR > mcsLookup(ind)
        for i=ind:1:maxMCS
            if y.EstimatedSNR < mcsLookup(i + 1)
                tempMCS = i - 1;
                break;
            elseif i == maxMCS
                tempMCS = maxMCS;
                break;
            end
        end
    else
        for i=ind:-1:2
            if y.EstimatedSNR > mcsLookup(i - 1)
                tempMCS = i - 2;
                break;
            end
        end
    end
    
    if numPkt > 1      
        %Counting how many packets in a row are received faulty
        for i = numPkt:-1:numPkt-faultyBound
            if i < 1
                break;
            end
            if ber(1, i) > 0
                faulty = faulty + 1;
            end
        end
        
        %Counting how many packets in a specified window, just before the
        %received one, are faulty
        for i = numPkt:-1:numPkt-sparseFaultyBound
            if i < 1
                break;
            end
            if ber(1, i) > 0
                sparseFaulty = sparseFaulty + 1;
            end
        end       
        
        %Counting how many times in a row the SNR is keep on dropping
        for i = numPkt:-1:numPkt-snrDecreasingBound
            if i < 2
                break;
            end
            if snrMeasured(1, i) < snrMeasured(1, i - 1)
                snrDecreasingCounter = snrDecreasingCounter + 1;
            else
                snrDecreasingCounter = 0;
                break;
            end
        end
        
        %If the look-up table thresholds were increased or they are not
        %equal to their initial values, we use this part of the code either
        %to set the thresholds available to be increased again or to
        %decrease them step by step back to their initial values
        if changeLookupCounter > 0 || ~isequal(mcsLookup, minMcsLookup)
            packetsCounter = packetsCounter + 1;
            if changeLookupCounter > changeLookupCounterBound
                changeLookupCounter = 0;
                packetsCounter = 0;
            elseif packetsCounter == packetsCounterBound
                changeLookupCounter = 0;
                packetsCounter = 0;
                if ~isequal(mcsLookup, minMcsLookup)
                    for i=2:1:10
                        mcsLookup(1, i) = mcsLookup(1, i) - 1;
                    end
                end
            end
        end
        
        %We decrease the MCS by one, or we decrease the MCS by one and
        %increase the thresholds, depending on the condition
        if sparseFaulty > sparseFaultyBoundCheck && ber(1, numPkt) > 0
            if tempMCS >= cfgVHT.MCS
                if cfgVHT.MCS > 0
                    tempMCS = cfgVHT.MCS - 1;
                else
                    tempMCS = 0;
                end          
            end
            if changeLookupCounter == 0
                for i=2:1:10
                    mcsLookup(1, i) = mcsLookup(1, i) + 2;
                end
            end            
            changeLookupCounter = changeLookupCounter + 1;
        elseif faulty > faultyBoundCheck
            if tempMCS >= cfgVHT.MCS
                if cfgVHT.MCS > 0
                    tempMCS = cfgVHT.MCS - 1;
                else
                    tempMCS = 0;
                end          
            end
        end
        
        %If the estimated SNR is dropping three times in a row and it is
        %very close to the lowest threshold of the current MCS, then we
        %decrease the MCS by one
        if snrDecreasingCounter == 3 && tempMCS >= cfgVHT.MCS && ...
                snrMeasured(1, numPkt) < mcsLookup(1,cfgVHT.MCS + 1) * (1 + deviation)
            if cfgVHT.MCS > 0
                tempMCS = cfgVHT.MCS - 1;
            else
                tempMCS = 0;
            end            
        end
        
        %If we increased the MCS and we received a faulty packet, then we
        %decrease the MCS by one back to its previous value
        if cfgVHT.MCS > MCS(numPkt - 1) && ber(1, numPkt) > 0 && tempMCS >= cfgVHT.MCS
            if cfgVHT.MCS > 0
                tempMCS = cfgVHT.MCS - 1;
            else
                tempMCS = 0;
            end
        end
    end
 
    MCS(numPkt) = cfgVHT.MCS;
    %Setting the MCS for the next transmission
    cfgVHT.MCS = tempMCS;
end

%% Display and Plot Simulation Results
overallDataRate = 8*cfgVHT.APEPLength*(numPackets-numel(find(ber)))/sum(packetLength)/1e6;
overallPacketErrorRate = numel(find(ber))/numPackets;
throughput = movsum(8*cfgVHT.APEPLength.*(ber==0),3)./movsum(packetLength,3)/1e6;
overallThroughput = mean(throughput);
%Display and plot simulation results
disp(['Overall data rate: ' num2str(overallDataRate) ' Mbps']);
disp(['Overall packet error rate: ' num2str(overallPacketErrorRate)]);
disp(['Overall throughput: ' num2str(overallThroughput)]);
plotResults(ber,packetLength,snrMeasured,MCS,cfgVHT);
disp(['min SNR: ' num2str(min(snrMeasured))]);
disp(['max SNR: ' num2str(max(snrMeasured))]);
disp(['Average SNR: ' num2str(mean(snrMeasured))]);
%Restore default stream
rng(s);

% * |processPacket|: Add channel impairments and decode receive packet
% * |plotResults|: Plot the simulation results

function Y = processPacket(txWave,snrWalk,tgacChannel,cfgVHT)
    % Pass the transmitted waveform through the channel, perform
    % receiver processing, and SNR estimation.
    
    chanBW = cfgVHT.ChannelBandwidth; % Channel bandwidth
    % Set the following parameters to empty for an undetected packet
    estimatedSNR = [];
    eqDataSym = [];
    noiseVarVHT = [];
    rxPSDU = [];
    
    % Get the number of occupied subcarriers in VHT fields
    [vhtData,vhtPilots] = helperSubcarrierIndices(cfgVHT,'VHT');
    Nst_vht = numel(vhtData)+numel(vhtPilots);
    Nfft = helperFFTLength(cfgVHT); % FFT length
    
    % Pass the waveform through the fading channel model
    rxWave = tgacChannel(txWave);
    
    % Create an instance of the AWGN channel for each transmitted packet
    awgnChannel = comm.AWGNChannel;
    awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
    % Normalization
    awgnChannel.SignalPower = 1/tgacChannel.NumReceiveAntennas;
    % Account for energy in nulls
    awgnChannel.SNR = snrWalk-10*log10(Nfft/Nst_vht);
    
    % Add noise
    rxWave = awgnChannel(rxWave);
    rxWaveformLength = size(rxWave,1); % Length of the received waveform
    
    % Recover packet
    ind = wlanFieldIndices(cfgVHT); % Get field indices
    pktOffset = wlanPacketDetect(rxWave,chanBW); % Detect packet
    
    if ~isempty(pktOffset) % If packet detected
        % Extract the L-LTF field for fine timing synchronization
        LLTFSearchBuffer = rxWave(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
    
        % Start index of L-LTF field
        finePktOffset = wlanSymbolTimingEstimate(LLTFSearchBuffer,chanBW);
     
        % Determine final packet offset
        pktOffset = pktOffset+finePktOffset;
        
        if pktOffset<15 % If synchronization successful
            % Extract L-LTF samples from the waveform, demodulate and
            % perform noise estimation
            LLTF = rxWave(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            demodLLTF = wlanLLTFDemodulate(LLTF,chanBW);

            % Estimate noise power in non-HT fields
            noiseVarVHT = helperNoiseEstimate(demodLLTF,chanBW,cfgVHT.NumSpaceTimeStreams,'Per Antenna');

            % Extract VHT-LTF samples from the waveform, demodulate and
            % perform channel estimation
            VHTLTF = rxWave(pktOffset+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
            demodVHTLTF = wlanVHTLTFDemodulate(VHTLTF,cfgVHT);
            chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);

            % Recover equalized symbols at data carrying subcarriers using
            % channel estimates from VHT-LTF
            [rxPSDU,~,eqDataSym] = wlanVHTDataRecover( ...
                rxWave(pktOffset + (ind.VHTData(1):ind.VHTData(2)),:), ...
                chanEstVHTLTF,mean(noiseVarVHT),cfgVHT);
            
            % SNR estimation per receive antenna
            powVHTLTF = mean(VHTLTF.*conj(VHTLTF));
            estSigPower = powVHTLTF-noiseVarVHT;
            estimatedSNR = 10*log10(mean(estSigPower./noiseVarVHT));
        end
    end
    
    % Set output
    Y = struct( ...
        'RxPSDU',           rxPSDU, ...
        'EqDataSym',        eqDataSym, ...
        'RxWaveformLength', rxWaveformLength, ...
        'NoiseVar',         noiseVarVHT, ...
        'EstimatedSNR',     estimatedSNR);
    
end

function plotResults(ber,packetLength,snrMeasured,MCS,cfgVHT)
    % Visualize simulation results

    figure('Outerposition',[60 15 1000 850])
    subplot(4,1,1);
    plot(MCS);
    xlabel('Packet Number')
    ylabel('MCS')
    title('MCS selected for transmission')

    subplot(4,1,2);
    plot(snrMeasured);
    xlabel('Packet Number')
    ylabel('SNR')
    title('Estimated SNR')

    subplot(4,1,3);
    plot(find(ber==0),ber(ber==0),'x') 
    hold on; stem(find(ber>0),ber(ber>0),'or') 
    if any(ber)
        legend('Successful decode','Unsuccessful decode') 
    else
        legend('Successful decode') 
    end
    xlabel('Packet Number')
    ylabel('BER')
    title('Instantaneous bit error rate per packet')

    subplot(4,1,4);
    windowLength = 3; % Length of the averaging window
    movDataRate = movsum(8*cfgVHT.APEPLength.*(ber==0),windowLength)./movsum(packetLength,windowLength)/1e6;
    plot(movDataRate)
    xlabel('Packet Number')
    ylabel('Mbps')
    title(sprintf('Throughput over last %d packets',windowLength))
end