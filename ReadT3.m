function ReadT3(par,trunc,skip) % Read PicoQuant Unified TTTR Files (without using uifileopen)
% This is demo code. Use at your own risk. No warranties.ead
% Marcus Sackrow, PicoQuant GmbH, December 2013
% Peter Kapusta, PicoQuant GmbH, November 2016
% Edited script: text output formatting changed by KAP.

% Note that marker events have a lower time resolution and may therefore appear
% in the file slightly out of order with respect to regular (photon) event records.
% This is by design. Markers are designed only for relatively coarse
% synchronization requirements such as image scanning.

% T Mode data are written to an output file [filename].out
% We do not keep it in memory because of the huge amout of memory
% this would take in case of large files. Of course you can change this,
% e.g. if your files are not too big.
% Otherwise it is best process the data on the fly and keep only the results.

% All HeaderData are introduced as Variable to Matlab and can directly be
% used for further analysis
    arguments
        par.ONthreshold = 0;
        par.OFFp = 0.2;
        par.bin = 10;
        par.roll = 8;
        par.DetEff = 0.09;
        par.output = 'Test';
        par.path = 'no imput';
        par.file = '.ptu';
        trunc.range = [0,3600];
        trunc.gate = [-65535,65535];
        skip.read = '.ptu';
        skip.blinking = 1;
        skip.g2 = 1;
        skip.lifetime = 1;
        skip.duration = 1;
        skip.plot = 1;
    end

    clc
    close all

    % some constants
    tyEmpty8      = hex2dec('FFFF0008');
    tyBool8       = hex2dec('00000008');
    tyInt8        = hex2dec('10000008');
    tyBitSet64    = hex2dec('11000008');
    tyColor8      = hex2dec('12000008');
    tyFloat8      = hex2dec('20000008');
    tyTDateTime   = hex2dec('21000008');
    tyFloat8Array = hex2dec('2001FFFF');
    tyAnsiString  = hex2dec('4001FFFF');
    tyWideString  = hex2dec('4002FFFF');
    tyBinaryBlob  = hex2dec('FFFFFFFF');
    % RecordTypes
    rtHydraHarp2T3   = hex2dec('01010304');% (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
    TTResultFormat_TTTRRecType = 0;
    TTResult_NumberOfRecords = 0;
    MeasDesc_Resolution = 0;
    MeasDesc_GlobalResolution = 0;
    TTResult_StopAfter = 0;

    if strcmp(par.file,'.ptu')
        [file, path]=uigetfile('*.ptu','Select data files:','Multiselect','on');
        if iscell(file)
            numberofFiles=length(file);
        else
            numberofFiles=1;
        end
    else
        file=par.file;
        if strcmp(par.path,'no imput')
            path='C:\Data Temporary\Test Runs\';
        else
            path=par.path;
        end
        numberofFiles=1;
    end

    for i=1:numberofFiles
        if numberofFiles==1
            filename=file(1:end-4);
        else
            filename=file{1,i};
            filename=filename(1:end-4);
        end
        
        outdir=[path,filename,'\'];

        if strcmp(skip.read,'.ptu')
            fid=fopen([path,filename,'.ptu']); % read .ptu files
            Magic = fread(fid, 8, '*char');
            if not(strcmp(Magic(Magic~=0)','PQTTTR'))
                error('Magic invalid, this is not an PTU file.');
            end
            Version = fread(fid, 8, '*char');
            fprintf(1,'Tag Version: %s\n', Version);

            % read headers
            % there is no repeat.. until (or do..while) construct in matlab so we use
            % while 1 ... if (expr) break; end; end;
            while 1
                % read Tag Head
                TagIdent = fread(fid, 32, '*char'); % TagHead.Ident
                TagIdent = (TagIdent(TagIdent ~= 0))'; % remove #0 and more more readable
                TagIdx = fread(fid, 1, 'int32');    % TagHead.Idx
                TagTyp = fread(fid, 1, 'uint32');   % TagHead.Typ
                % TagHead.Value will be read in the
                % right type function
                TagIdent = genvarname(TagIdent);    % remove all illegal characters
                if TagIdx > -1
                    EvalName = [TagIdent '(' int2str(TagIdx + 1) ')'];
                else
                    EvalName = TagIdent;
                end
                fprintf(1,'\n   %-40s', EvalName);
                % check Typ of Header
                switch TagTyp
                    case tyEmpty8
                        fread(fid, 1, 'int64');
                        fprintf(1,'<Empty>');
                    case tyBool8
                        TagInt = fread(fid, 1, 'int64');
                        if TagInt==0
                            fprintf(1,'FALSE');
                            eval([EvalName '=false;']);
                        else
                            fprintf(1,'TRUE');
                            eval([EvalName '=true;']);
                        end
                    case tyInt8
                        TagInt = fread(fid, 1, 'int64');
                        fprintf(1,'%d', TagInt);
                        eval([EvalName '=TagInt;']);
                    case tyBitSet64
                        TagInt = fread(fid, 1, 'int64');
                        fprintf(1,'%X', TagInt);
                        eval([EvalName '=TagInt;']);
                    case tyColor8
                        TagInt = fread(fid, 1, 'int64');
                        fprintf(1,'%X', TagInt);
                        eval([EvalName '=TagInt;']);
                    case tyFloat8
                        TagFloat = fread(fid, 1, 'double');
                        fprintf(1, '%e', TagFloat);
                        eval([EvalName '=TagFloat;']);
                    case tyFloat8Array
                        TagInt = fread(fid, 1, 'int64');
                        fprintf(1,'<Float array with %d Entries>', TagInt / 8);
                        fseek(fid, TagInt, 'cof');
                    case tyTDateTime
                        TagFloat = fread(fid, 1, 'double');
                        fprintf(1, '%s', datestr(datenum(1899,12,30)+TagFloat)); % display as Matlab Date String
                        eval([EvalName '=datenum(1899,12,30)+TagFloat;']); % but keep in memory as Matlab Date Number
                        File_CreatingDateTime=datestr(File_CreatingTime);
                    case tyAnsiString
                        TagInt = fread(fid, 1, 'int64');
                        TagString = fread(fid, TagInt, '*char');
                        TagString = (TagString(TagString ~= 0))';
                        fprintf(1, '%s', TagString);
                        if TagIdx > -1
                            EvalName = [TagIdent '{' int2str(TagIdx + 1) '}'];
                        end
                        eval([EvalName '=[TagString];']);
                    case tyWideString
                        % Matlab does not support Widestrings at all, just read and
                        % remove the 0's (up to current (2012))
                        TagInt = fread(fid, 1, 'int64');
                        TagString = fread(fid, TagInt, '*char');
                        TagString = (TagString(TagString ~= 0))';
                        fprintf(1, '%s', TagString);
                        if TagIdx > -1
                            EvalName = [TagIdent '{' int2str(TagIdx + 1) '}'];
                        end
                        eval([EvalName '=[TagString];']);
                    case tyBinaryBlob
                        TagInt = fread(fid, 1, 'int64');
                        fprintf(1,'<Binary Blob with %d Bytes>', TagInt);
                        fseek(fid, TagInt, 'cof');
                    otherwise
                        error('Illegal Type identifier found! Broken file?');
                end
                if strcmp(TagIdent, 'Header_End')
                    break
                end
            end
            fprintf(1, '\n----------------------\n');

            % Check recordtype
            if TTResultFormat_TTTRRecType == rtHydraHarp2T3
                fprintf(1,'HydraHarp V2 T3 data\n');
            else
                error('Illegal RecordType!');
            end

            % Initialization
            T3WRAPAROUND = 1024;
            RepTime=round(MeasDesc_GlobalResolution*1e9,3);
            Res=round(MeasDesc_Resolution*1e9,3);
            reliab=1;
            reason='_';
            nbins=ceil(RepTime/Res); % number of bins for TCSPC
            tdini=-0.12*RepTime;
            tdfin=0.88*RepTime;
            [~,~] = mkdir(outdir);

            % read time tags, obtain 3 waves CHN, Sync, PhotonTime, and 
            fprintf(1,'Reading T3 data from %s...\n',filename);
            T3Record = fread(fid, TTResult_NumberOfRecords, 'ubit32');     % read all rows, all 32 bits in each row:
            %   +-------------------------------+  +-------------------------------+
            %   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
            %   +-------------------------------+  +-------------------------------+
            Sync = bitand(T3Record,1023);       % the lowest 10 bits:
            %   +-------------------------------+  +-------------------------------+
            %   | | | | | | | | | | | | | | | | |  | | | | | | |x|x|x|x|x|x|x|x|x|x|
            %   +-------------------------------+  +-------------------------------+
            PhotonTime = bitand(bitshift(T3Record,-10),32767);   % the next 15 bits:
            %   the dtime unit depends on "Resolution" that can be obtained from header
            %   +-------------------------------+  +-------------------------------+
            %   | | | | | | | |x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x| | | | | | | | | | |
            %   +-------------------------------+  +-------------------------------+
            CHN = bitand(bitshift(T3Record,-25),63);   % the next 6 bits:
            %   +-------------------------------+  +-------------------------------+
            %   | |x|x|x|x|x|x| | | | | | | | | |  | | | | | | | | | | | | | | | | |
            %   +-------------------------------+  +-------------------------------+
            special = bitand(bitshift(T3Record,-31),1);   % the last bit:
            %   +-------------------------------+  +-------------------------------+
            %   |x| | | | | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
            %   +-------------------------------+  +-------------------------------+
            fclose(fid);
            OverflowCorrection=Sync;
            OverflowCorrection(special==0)=0;
            OverflowCorrection=cumsum(OverflowCorrection);
            cnt_ov=OverflowCorrection(end);
            Sync=Sync+OverflowCorrection*T3WRAPAROUND;
            Sync=Sync(special==0);
            CHN=CHN(special==0);
            PhotonTime=PhotonTime(special==0)*Res;

            % correct delay time between two channels
            [c0,t0]=histcounts(PhotonTime(CHN==0),nbins);
            [c1,t1]=histcounts(PhotonTime(CHN==1),nbins);
            [~,d0]=rollmax(c0,par.roll);
            [~,d1]=rollmax(c1,par.roll);
            PhotonTime(CHN==0)=PhotonTime(CHN==0)-t0(d0+1);
            PhotonTime(CHN==1)=PhotonTime(CHN==1)-t1(d1+1);
            Sync(PhotonTime<tdini)=Sync(PhotonTime<tdini)+1;
            Sync(PhotonTime>=tdfin)=Sync(PhotonTime>=tdfin)-1;
            PhotonTime(PhotonTime<tdini)=PhotonTime(PhotonTime<tdini)+RepTime;
            PhotonTime(PhotonTime>=tdfin)=PhotonTime(PhotonTime>=tdfin)-RepTime;
            TrueTime=round(Sync*RepTime+PhotonTime,3);
            PhotonTime=round(PhotonTime,3);
            i0=length(CHN(CHN==0));
            i1=length(CHN(CHN==1));
            cnt_ph=i0+i1;

            % write translated data to file
            outfile = [path,filename,'.dat'];
            fpout = fopen(outfile,'W');
            fprintf(1,'Delay corrected (d0=%g ns,d1=%g ns), writing data to %s...\n',t0(d0+1),t1(d1+1),outfile);
            fprintf(fpout,'%s      %g\n','RepTime(ns)',RepTime);
            fprintf(fpout,'%s   %g\n','Resolution(ns)',Res);
            fprintf(fpout,'%s      %g\n','Duration(s)',TTResult_StopAfter/1000);
            fprintf(fpout,'%s %u\n','Channel 0 Counts',i0);
            fprintf(fpout,'%s %u\n','Channel 1 Counts',i1);
            fprintf(fpout,'%s\n','------Header End------');
            fprintf(fpout,'CHN,tTrue_ns,tPh_ns');
            fclose(fpout);
            writematrix([CHN,TrueTime,PhotonTime],outfile,'WriteMode','append');
        else
            error('function coming soon...')
        end
        
        % truncate data for statistical analysis
        if trunc.range(1)<0
            trunc.range(1)=0;
        end
        if trunc.range(2)>TTResult_StopAfter/1000
            trunc.range(2)=TTResult_StopAfter/1000;
        end
        duration=trunc.range(2)-trunc.range(1);
        if duration<=0
            error('Zero duration!')
        end
        recordfile=[filename,'_',num2str(trunc.range(1)),'-',num2str(trunc.range(2)),'s'];
        if duration~=TTResult_StopAfter/1000
            CHN=CHN((TrueTime>=trunc.range(1)*1e9)&(TrueTime<=trunc.range(2)*1e9));
            Sync=Sync((TrueTime>=trunc.range(1)*1e9)&(TrueTime<=trunc.range(2)*1e9));
            PhotonTime=PhotonTime((TrueTime>=trunc.range(1)*1e9)&(TrueTime<=trunc.range(2)*1e9));
            TrueTime=TrueTime((TrueTime>=trunc.range(1)*1e9)&(TrueTime<=trunc.range(2)*1e9));
        end
        fprintf(1,'Data saved, analyzing data within the range %s\n',recordfile);
        fprintf(1,'Statistics obtained from the data:\n%i photons (%i in CHN0,%i in CHN1), %i overflows.\n',cnt_ph,i0,i1,cnt_ov);

        % Blinking and ON/OFF statistics
        if skip.blinking==1
            tBlink=(trunc.range(1)*1e9:par.bin*1e6:ceil(trunc.range(2)*1e3/par.bin)*par.bin*1e6)';
            [cBlink,~,bindx]=histcounts(TrueTime,tBlink);
            cBlink=cBlink';
            tBlink=tBlink/1e9;
            cBlink_max=max(cBlink);
            counts=(0:1:cBlink_max+1)';
            occurrence=histcounts(cBlink,counts)';
            max_occurrence=max(occurrence);
            if par.ONthreshold==0
                stdev_ON=(sqrt(2*log(max_occurrence)+4*cBlink_max)-sqrt(2*log(max_occurrence)))/2;
                ON_level=round(stdev_ON^2,1);
                ON_threshold=round(ON_level-stdev_ON*3,1);
                if ON_threshold<1
                    ON_threshold=0.5;
                    reliab=reliab*0;
                    reason=[reason,'noON_']; %#ok<*AGROW,*NASGU> 
                    reducedwidth=0;
                    ON_level=1;
                    stdev_ON=0;
                    pw=[sum(occurrence),ON_level,1.2];
                else
                    pw=[sum(occurrence(fix(ON_threshold)+1:cBlink_max+1)),ON_level,1.2];
                    occurrence_fit=fitpoiss(counts,occurrence,pw,min(fix(ON_threshold)+1,fix(cBlink_max/2)+1),cBlink_max+1);
                    pw(1:3)=coeffvalues(occurrence_fit);
                    reducedwidth=round(pw(3),2);
                    ON_level=round(pw(2),1);
                    stdev_ON=sqrt(pw(2));
                    ON_threshold=round(ON_level-stdev_ON*3,1);
                end
            else
                ON_threshold=Par.ONthreshold;
            end
            OFF_level=round(ON_level*par.OFFp,1);
            OFF_threshold=round(OFF_level+3*sqrt(OFF_level),1);
            if OFF_threshold>ON_threshold
                OFF_threshold=ON_threshold;
                reliab=reliab*0;
                reason=[reason,'lowON_'];
            end
            pON=round(length(cBlink(cBlink>ON_threshold))/duration*par.bin/10,1);
            pON_poiss=round(pw(1)/duration*par.bin/10,1);
            if pON_poiss<pON*0.5
                reliab=reliab*0;
                reason=[reason,'nonPoiss_'];
            end
            pOFF=round(length(cBlink(cBlink<OFF_threshold))/duration*par.bin/10,1);
            pGrey=round(100-pON-pOFF,1);
            if pGrey<0
                pGrey=0;
            end
            if pON>=90
                BlinkingType='NonBlinking';
            elseif pOFF>=90
                BlinkingType='Dark';
                reliab=reliab*0;
                reason=[reason,'dark_'];
            elseif pGrey>=10&&pOFF>=10
                BlinkingType='Flickering';
            elseif pOFF>pGrey
                BlinkingType='BinaryBlinking';
            elseif pOFF<10
                BlinkingType='Fluctuating';
            else
                BlinkingType='Error!';
                reliab=reliab*0;
                reason=[reason,'error_'];
            end
            AvgN=round(ON_level*RepTime/par.bin/1e6/par.DetEff,3);
            fprintf(1,'Estimated <N> = %g\n',AvgN);
            % Write blinking trace to file
            outfile=[outdir,'blinking_',recordfile,'.txt'];
            fpout=fopen(outfile,'W');
            fprintf(fpout,'%s      %g\n','Duration(s)',duration);
            fprintf(fpout,'%s     %g\n','Bin time(ms)',par.bin);
            fprintf(fpout,'%s            %g\n','ONlvl',ON_level);
            fprintf(fpout,'%s          %g\n','ONthrld',ON_threshold);
            fprintf(fpout,'%s           %g\n','OFFlvl',OFF_level);
            fprintf(fpout,'%s         %g\n','OFFthrld',OFF_threshold);
            fprintf(fpout,'%s        %g\n','ONtime(%)',pON);
            fprintf(fpout,'%s       %g\n','OFFtime(%)',pOFF);
            fprintf(fpout,'%s\n','------Header End------');
            fprintf(fpout,'tBlink,cBlink');
            fclose(fpout);
            writematrix([tBlink(2:end),cBlink],outfile,'WriteMode','append');
        else
            error('Function coming soon...')
        end

        % Power-law duration analysis
        if skip.duration==1
            [tONavg,tnONavg,ParON,ParnON,PdurON,PdurnON]=DurationAnalysis(cBlink,(tBlink-trunc.range(1)),ON_threshold);
            fprintf(1,'<t_ON>=%gs,<t_~ON>=%gs,a_ON=%g,a_~ON=%g\n',tONavg,tnONavg,round(ParON(2),2),round(ParnON(2),2));
            [tnOFFavg,tOFFavg,ParnOFF,ParOFF,PdurnOFF,PdurOFF]=DurationAnalysis(cBlink,(tBlink-trunc.range(1)),OFF_threshold);
            fprintf(1,'<t_~OFF>=%gs,<t_OFF>=%gs,a_~OFF=%g,a_OFF=%g\n',tnOFFavg,tOFFavg,round(ParnOFF(2),2),round(ParOFF(2),2));
            % Write lifetime trace and fits to file
            outfile=[outdir 'durationanalysis_' recordfile '.txt'];
            fpout=fopen(outfile,'W');
            fprintf(fpout,'%s        %g\n','<t_ON> = ',tONavg);
            fprintf(fpout,'%s       %g\n','<t_~ON> = ',tnONavg);
            fprintf(fpout,'%s          %g\n','a_ON = ',round(ParON(2),2));
            fprintf(fpout,'%s         %g\n','a_~ON = ',round(ParnON(2),2));
            fprintf(fpout,'%s      %g\n','<t_~OFF> = ',tnOFFavg);
            fprintf(fpout,'%s       %g\n','<t_OFF> = ',tOFFavg);
            fprintf(fpout,'%s        %g\n','a_~OFF = ',round(ParnOFF(2),2));
            fprintf(fpout,'%s         %g\n','a_OFF = ',round(ParOFF(2),2));
            fprintf(fpout,'%s\n','------Header End------');
            fprintf(fpout,'tBlink,PdurON,PdurnON,PdurnOFF,PdurOFF');
            fclose(fpout);
            writematrix([tBlink(2:end),PdurON,PdurnON,PdurnOFF,PdurOFF],outfile,'WriteMode','append');
        else
            error('Function coming soon...');
        end

        % Lifetime statistics and fittings
        if skip.lifetime==1
            tTCSPC=(tdini:Res:tdfin+Res)';
            cTCSPC=histcounts(PhotonTime,tTCSPC)';
            tTCSPC=tTCSPC(1:end-1);
            cTCSPC_max=max(cTCSPC);
            ParDecay=zeros(1,7);
            ParDecay(1)=mean(cTCSPC(1:fix(0.1*nbins))); %Background
            ParDecay(2)=cTCSPC_max;
            ParDecay(3)=4;
            nexp=1;
            while nexp<=3
                [decay_fit,DecayType]=fitdecay(tTCSPC,cTCSPC,nexp,ParDecay);
                ParDecay(1:2*nexp+1)=round(coeffvalues(decay_fit),1);
                decay_fit=ParDecay(1)+ParDecay(2)*exp(-tTCSPC/ParDecay(3));
                if nexp==2
                    decay_fit=decay_fit+ParDecay(4)*exp(-tTCSPC/ParDecay(5));
                end
                if nexp==3
                    decay_fit=decay_fit+ParDecay(6)*exp(-tTCSPC/ParDecay(7));
                end
                decay_fit(tTCSPC<0)=ParDecay(1);
                residue=decay_fit-cTCSPC;
                residue(tTCSPC<0)=0;
                residue(tTCSPC>48)=0;
                chisq=sum((residue./decay_fit).^2)/nbins;
                if chisq<2
                    break
                elseif nexp<3
                    ParDecay(2*nexp)=ParDecay(2*nexp)/2;
                    ParDecay(2*nexp+3)=ParDecay(2*nexp+1)*2;
                    ParDecay(2*nexp+2)=ParDecay(2*nexp);
                    nexp=nexp+1;
                else
                    break
                end
            end
            Atot=round(ParDecay(2)+ParDecay(4)+ParDecay(6));
            % Write lifetime trace and fits to file
            outfile=[outdir,'lifetime_',recordfile,'.txt'];
            fpout=fopen(outfile,'W');
            fprintf(fpout,'%s          %g\n','Atot = ',Atot);
            fprintf(fpout,'%s            %g\n','BL = ',ParDecay(1));
            fprintf(fpout,'%s            %g\n','A1 = ',round(ParDecay(2)/Atot,3));
            fprintf(fpout,'%s     %g\n','tau1 (ns) = ',ParDecay(3));
            fprintf(fpout,'%s            %g\n','A2 = ',round(ParDecay(4)/Atot,3));
            fprintf(fpout,'%s     %g\n','tau2 (ns) = ',ParDecay(5));
            fprintf(fpout,'%s            %g\n','A3 = ',round(ParDecay(4)/Atot,3));
            fprintf(fpout,'%s     %g\n','tau3 (ns) = ',ParDecay(7));
            fprintf(fpout,'%s      %s\n','Decay Type:',DecayType);
            fprintf(fpout,'%s\n','------Header End------');
            fprintf(fpout,'tTCSPC,cTCSPC,decay_fit,residue');
            fclose(fpout);
            writematrix([tTCSPC,cTCSPC,decay_fit,residue],outfile,'WriteMode','append');
        else
            error('Function coming soon...');
        end

        % g2 statistics and write g2 to file
        if skip.g2==1
            if trunc.gate(1)<tdini
                trunc.gate(1)=tdini;
            end
            if trunc.gate(2)>tdfin
                trunc.gate(2)=tdfin;
            end
            if trunc.gate(1)>tdini||trunc.gate(2)<tdfin
                CHN=CHN((PhotonTime>=trunc.gate(1))&(PhotonTime<=trunc.gate(2)));
                Sync=Sync((PhotonTime>=trunc.gate(1))&(PhotonTime<=trunc.gate(2)));
                PhotonTime=PhotonTime((PhotonTime>=trunc.gate(1))&(PhotonTime<=trunc.gate(2)));
                outfile=[outdir,'g2_',recordfile,'_gated',num2str(trunc.gate(1)),'-',num2str(trunc.gate(2)),'ns.txt'];
            else
                outfile=[outdir,'g2_',recordfile,'.txt'];
            end
            [g2,gtime,g2zero,g2zero_corr,g2_err]=g2_pulseresolved(PhotonTime(CHN==0),PhotonTime(CHN==1),Sync(CHN==0),Sync(CHN==1),1,RepTime,'pBG',ParDecay(1)*nbins/(duration*ON_level*1000/par.bin));
            g2_max=max(g2);
            fprintf('Calculating ON state g2:\n')
            [~,ONPhotonTime,PhotonMrk]=RetrievePhotons(ON_threshold,cBlink_max,cBlink,bindx,PhotonTime);
            ONSync=Sync(PhotonMrk);
            ONCHN=CHN(PhotonMrk);
            [g2_ON,gtime_ON,g2zero_ON,g2zero_corr_ON,g2_err_ON]=g2_pulseresolved(ONPhotonTime(ONCHN==0),ONPhotonTime(ONCHN==1),ONSync(ONCHN==0),ONSync(ONCHN==1),1,RepTime,'pBG',ParDecay(1)*nbins/(duration*ON_level*1000/par.bin));
            % Write g2 and g2_ON to file
            fid=fopen(outfile,'w');
            fprintf(fid,'g2(0) =          %g\ng2(0)_corr =     %g\ng2(0)_ON_corr =  %g\nGate open (ns)   %g\nGate close (ns)  %g\n%s\n%s,%s,%s',g2zero,g2zero_corr,g2zero_corr_ON,trunc.gate(1),trunc.gate(2),'------Header End------','gtime_ns','g2','g2_ON');
            fclose(fid);
            writematrix([gtime;g2;g2_ON]',outfile,'WriteMode','append');
            if g2zero>=0.4
                reliab=reliab*0;
                reason=[reason,'hig2'];
            end
        else
            error('Function coming soon...');
        end

        % Plot blinking, blinking stats, g2 and lifetime
        if skip.plot==1
            figure('Position',[600 400 640 720]);
            subplot(3,2,1); % blinking
            histogram('BinEdges',tBlink','BinCounts',cBlink','DisplayStyle','stairs'); title([recordfile,' ',BlinkingType]);ylabel(['Count rate (ct/',num2str(par.bin),' ms)']);xlabel('Time (s)');xlim(trunc.range);ylim([0,cBlink_max*1.2+1]);
            text(trunc.range(1)+1,cBlink_max*1.1,['<N> = ',num2str(AvgN),' ON Lvl:',num2str(round(ON_level/10,1)),'kHz']);
            hold on;
            plot(trunc.range,[ON_threshold,ON_threshold],'color','r');
            plot(trunc.range,[OFF_threshold,OFF_threshold],'color','b');
            hold off;
            subplot(3,2,2); % blinking stats
            histogram('BinEdges',counts','BinCounts',occurrence','Orientation','horizontal'); title(['Measured on ',File_CreatingDateTime]);ylabel(['Count rate (ct/',num2str(par.bin),' ms)']);xlabel('Occurrence');xlim([0,max_occurrence*1.2+1]);ylim([0,cBlink_max*1.2+1]);text(max_occurrence,cBlink_max*1.1,[num2str(pON),'%ON'],'color','r');text(max_occurrence,0.1*cBlink_max,[num2str(pOFF),'%OFF'],'color','b');
            hold on;
            plot([0;max_occurrence*1.2],[ON_threshold;ON_threshold],'color','r');
            plot([0;max_occurrence*1.2],[OFF_threshold;OFF_threshold],'color','b');
            hold off;
            subplot(3,2,3); % g2
            plot(gtime,g2,gtime_ON,g2_ON);title(['g2 = ',num2str(g2zero)]);ylabel('g2 (counts)');xlabel('Delay (ns)');xlim([-2.5*RepTime,2.5*RepTime]);ylim([0,g2_max*1.2+1]);text(-2*RepTime,g2_max*1.1,['g2_c_o_r_r = ',num2str(g2zero_corr) ' g2_O_N = ',num2str(g2zero_corr_ON)]);
            subplot(3,2,4); % lifetime
            plot(tTCSPC,cTCSPC,'.'); set(gca,'Yscale','log'); title(DecayType);ylabel('PL intensity (a.u.)');xlabel('Time (ns)');xlim([-0.1*RepTime,0.9*RepTime]);ylim([1,Atot*5+2]);
            text(0.5*RepTime,Atot*2,['A=',num2str(Atot),' BL=' num2str(round(ParDecay(1),2))]);
            text(0.5*RepTime,Atot*0.8,['T_1=',num2str(ParDecay(3)),'ns(',num2str(round(ParDecay(2)/Atot,3)),')']);
            if ParDecay(4)>0
                text(0.5*RepTime,Atot*0.32,['T_2=',num2str(ParDecay(5)),'ns(' num2str(round(ParDecay(4)/Atot,3)),')']);
            end
            if ParDecay(6)>0
                text(0.5*RepTime,Atot*0.128,['T_3=',num2str(ParDecay(7)),'ns(' num2str(round(ParDecay(6)/Atot,3)),')']);
            end
            hold on;
            plot(tTCSPC,decay_fit,'LineStyle','-','Color','r');
            hold off;
            subplot(3,2,5); % power-law duration analysis using ON threshold
            plot(tBlink(2:end),PdurON,'o','Color','r')
            hold on;
            plot(tBlink(2:end),PdurnON,'o','Color','b')
            plot(tBlink(2:end),ParON(1).*(tBlink(2:end).^(-ParON(2))),'Color','r')
            plot(tBlink(2:end),ParnON(1).*(tBlink(2:end).^(-ParnON(2))),'Color','b')
            set(gca,'XScale','log','YScale','log','Xlim',[par.bin/1000,duration],'Ylim',[par.bin/(duration*1000),1]);ylabel('Probability');xlabel('ON/~ON Duration');title(['a_O_N=',num2str(round(ParON(2),2)),' a_~_O_N=',num2str(round(ParnON(2),2))]);
            text(0.02*duration,0.5,['<t_O_N>=',num2str(round(tONavg,3)),'s']);
            text(0.02*duration,0.2,['<t_~_O_N>=',num2str(round(tnONavg,3)),'s']);
            hold off;
            subplot(3,2,6); % power-law duration analysis using OFF threshold
            plot(tBlink(2:end),PdurnOFF,'o','Color','r')
            hold on;
            plot(tBlink(2:end),PdurOFF,'o','Color','b')
            plot(tBlink(2:end),ParnOFF(1).*(tBlink(2:end).^(-ParnOFF(2))),'Color','r')
            plot(tBlink(2:end),ParOFF(1).*(tBlink(2:end).^(-ParOFF(2))),'Color','b')
            set(gca,'XScale','log','YScale','log','Xlim',[par.bin/1000,duration],'Ylim',[par.bin/(duration*1000),1]);ylabel('Probability');xlabel('OFF/~OFF Duration');title(['a_~_O_F_F=',num2str(round(ParnOFF(2),2)),' a_O_F_F=',num2str(round(ParOFF(2),2))]);
            text(0.02*duration,0.5,['<t_~_O_F>=',num2str(round(tnOFFavg,3)),'s']);
            text(0.02*duration,0.2,['<t_O_F_F>=',num2str(round(tOFFavg,3)),'s']);
            hold off;
            saveas(gcf,fullfile([path,'Plot_',recordfile,'_N',num2str(fix(AvgN)),'p',num2str(fix(AvgN*10-fix(AvgN)*10)),num2str(fix(AvgN*100-fix(AvgN*10)*10)),num2str(fix(AvgN*1000-fix(AvgN*100)*10)),'_',BlinkingType,'.jpg']));
        end

        outmat={File_CreatingDateTime,path,filename,trunc.range,trunc.gate,AvgN,duration,cnt_ph,ON_level,ON_threshold,pON,pON_poiss,OFF_level,OFF_threshold,pOFF,pGrey,tONavg,tnONavg,ParON(2),ParnON(2),tnOFFavg,tOFFavg,ParnOFF(2),ParOFF(2),g2zero,g2zero_corr,g2_err,g2zero_ON,g2zero_corr_ON,g2_err_ON,reducedwidth,ParDecay,BlinkingType,DecayType,reliab,reason};
        writecell(outmat,'C:\Data Temporary\Data g2 Stat Paper\Single QD Statistics.xlsx','FileType','spreadsheet','Sheet',par.output,'WriteMode','append');
    end

function [ft,DecayType]=fitdecay(wx,wy,n,par)
if n==1
    fitfunc=fittype('y0+a*exp(-1/tau*x)','coefficients',{'y0','a','tau'},'independent','x');
    ft=fit(wx,wy,fitfunc,'exclude',wx<0,'Weights',wy,'startpoint',par(1:3),'lower',[par(1),sqrt(par(2))-par(1),0.032],'upper',[par(1),par(2)+sqrt(par(2))-par(1),128]);
    DecayType='monoexp';
elseif n==2
    fitfunc=fittype('y0+a1*exp(-1/tau1*x)+a2*exp(-1/tau2*x)','coefficients',{'y0','a1','tau1','a2','tau2'},'independent','x');
    ft=fit(wx,wy,fitfunc,'exclude',wx<0,'Weights',wy,'startpoint',par(1:5),'lower',[par(1),par(2)/16,0.032,par(2)/16,par(3)],'upper',[par(1),par(2)+sqrt(par(2))-par(1),par(3),par(2)+sqrt(par(2))-par(1),128]);
    DecayType='biexp';
else
    fitfunc=fittype('y0+a1*exp(-1/tau1*x)+a2*exp(-1/tau2*x)+a3*exp(-1/tau3*x)','coefficients',{'y0','a1','tau1','a2','tau2','a3','tau3'},'independent','x');
    ft=fit(wx,wy,fitfunc,'exclude',wx<0,'Weights',wy,'startpoint',par(1:7),'lower',[par(1),(par(2)+par(4))/16,0.032,(par(2)+par(4))/16,par(3),(par(2)+par(4))/16,par(5)],'upper',[par(1),par(2)+par(4)+sqrt(par(2)+par(4))-par(1),par(3),par(2)+par(4)+sqrt(par(2)+par(4))-par(1),par(5),par(2)+par(4)+sqrt(par(2)+par(4))-par(1),128]);
    DecayType='multiexp';
end

function ft=fitpoiss(xw,yw,pw,LwB,UpB)
    fitfunc=fittype('Area/(sqrt(2*pi*mean)*width_red)*exp(-(x-mean)^2/(2*mean*width_red^2))','coefficients',{'Area','mean','width_red'},'independent','x');
    ft=fit(xw(LwB:UpB),yw(LwB:UpB),fitfunc,'startpoint',pw,'lower',[pw(1)/3,LwB,1],'upper',[pw(1),UpB,1.2]);
