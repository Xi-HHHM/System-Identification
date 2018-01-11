%19.12.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

DateiMenge = 693;
Warning = 0.05; % Schwerwert
%%  Initialisieren
pathname = '';
filename = '';

% --- Merkmale Initialisieren ---
MaxAbweichungOn = zeros(DateiMenge,3);  % maximale Abeichung
MaxAbweichungOff = zeros(DateiMenge,3);

QuaErrOn = zeros(DateiMenge,3);     % quadratische Fehler
QuaErrOff = zeros(DateiMenge,3);

StdQErrOn = zeros(DateiMenge,1);    % Standardabweichung der Quadratischen Fehler
StdQErrOff = zeros(DateiMenge,1);

StdMAbwOff = zeros(DateiMenge,1);   % Std der max. Abweichung
StdMAbwOn = zeros(DateiMenge,1);

Ttlist = zeros(DateiMenge,1);       % Zeitkonstante
T1list = zeros(DateiMenge,1);       % Totzeit

Drehzahl = zeros(DateiMenge,1);     % Drehzahl

MittelDruck = zeros(DateiMenge,1);
FriLuft = zeros(DateiMenge,1);      % Frischluftfüllung, die den Arbeitspunkt beschreibt

FriLuftOn = zeros(DateiMenge,3);    % Frischluftfüllung kurz vor den Sprüngen
FriLuftOff = zeros(DateiMenge,3);

Diff_Boost_On = zeros(DateiMenge,3);    % Boost kurz vor den Sprüngen
Diff_Boost_Off = zeros(DateiMenge,3);

Mittel_bCrit = zeros(DateiMenge,3);     % Mittelwert der bCrit

Diff_Des_On = zeros(DateiMenge,1);      % Des kurz vor den Sprüngen
Diff_Des_Off = zeros(DateiMenge,1);

%% Daten Bearbeiten
FehlerListe = 46;
for i = 1 : DateiMenge
    if ~ismember(i,FehlerListe)      
    datafile = sprintf('.%03i.dat', i);
    datafile = strcat(filename,datafile);
    
    Daten = mdfload2([pathname '\' datafile],...
        'UEGO_rLamS1B1\ETKC:1',...         
        'GEMChCont_rAirCh_VW\ETKC:1',...   
        'GECGsl_rEgrDesDyn_VW\ES910/Simulation Controller:1',...  % Sollwert eAGR
        'GEM_EgrLp_rEgrCyl_VW\ES910/Simulation Controller:1',...  % Modellwert eAGR
        '?GEM_EgrLp_tiFilMdl_VW',...  % Zeitkonstante eAGR
        '?GEM_EgrLp_tiDlyMdl_VW',...  % Totzeit eAGR
        'Epm_nEng_VW\ETKC:1',...      % Drehzahl
        'rl_w_msg\ETKC:1',...         % Frischluftfüllung
        '?PI0',...
        '?Air_pInMnfBoostBas_VW',...
        '?GEMChCont_pInMnf_VW',...
        '?GECGsl_pInMnfDes_VW',...
        '?GECGsl_bCritAcvCtlTvh_pInMnf_VW'...
        ); 
    
    tIn = Daten.GECGsl_rEgrDesDyn_VW_ES910_SimulationController_1(:,1);
    dataIn = Daten.GECGsl_rEgrDesDyn_VW_ES910_SimulationController_1(:,2);
    
    tOut = Daten.UEGO_rLamS1B1_ETKC_1(:,1);
    dataOut = Daten.UEGO_rLamS1B1_ETKC_1(:,2);
    
    Tt = Daten.GEM_EgrLp_tiDlyMdl_VW_GEM_EgrLp_Ausgangsadapter_ES910_Simulatio(:,2);
    T1 = Daten.GEM_EgrLp_tiFilMdl_VW_GEM_EgrLp_Ausgangsadapter_ES910_Simulatio(:,2);
    
    Modell_eAGR = Daten.GEM_EgrLp_rEgrCyl_VW_ES910_SimulationController_1(:,2);
    
    Ttlist(i) = median(Tt);
    T1list(i) = median(T1);

    % --------------- Datenauswerten ------------------
    minDataIn = min(dataIn);
    
    TtIdx = floor(Tt/0.01);
    
    % --- Positionen der Sprünge(Soll-AGR-Füllung) erkennen ---
    posOn  = find(diff((dataIn- minDataIn)>(0.1*max(dataIn-minDataIn)))>0);
    posOff  = find(diff((dataIn- minDataIn)>(0.1*max(dataIn-minDataIn)))<0);
    
    posOnTt = posOn + TtIdx(posOn(:)',1);
    posOffTt = posOff + TtIdx(posOff(:)',1);
    
    % Die Arbeitpunkt wird mit der Messung t = 5s definiert.
    Drehzahl(i,:) = Daten.Epm_nEng_VW_ETKC_1(5,2);
    FriLuft(i,:) = Daten.rl_w_msg_ETKC_1(5,2);       
    MittelDruck(i,:) = Daten.PI0_IndiCom_1(15,2);
    
    % Punkte, an dennen die Messung von einander getrennt wird
    Schnittpunkt = floor((posOnTt + posOffTt)/2);        % Array
    Schnittpunkt2 = floor((posOffTt(1) + posOnTt(2))/2); % Skalar
    Fenster = Schnittpunkt2 - Schnittpunkt(1); % Zeitspanne einzeles Sprungs
    
    % Soll, Sprünge 
    dataInOn = zeros(3,Fenster+1);
    dataInOff = zeros(3,Fenster+1);
    
    % Gemischabweichung
    dataOutOn = zeros(3,Fenster+1);
    dataOutOff = zeros(3,Fenster+1);
    
    % Modellierte eAGR-Sprünge
    dataMeAGROn = zeros(3,Fenster+1);
    dataMeAGROff = zeros(3,Fenster+1);
    
    % Frischluftfüllung kurz vor den Sprüngen
    dataFriLuftOn = zeros(3,Fenster+1);
    dataFriLuftOff = zeros(3,Fenster+1);
    
    % Boost des Sprüngen Mittelwert
    dataDiff_Boost_On = zeros(3,Fenster+1);
    dataDiff_Boost_Off = zeros(3,Fenster+1);
    
    dataDiff_Des_On = zeros(3,Fenster+1);
    dataDiff_Des_Off = zeros(3,Fenster+1);
    
    % Boost_KF und Des_KF für die ganze Messung
    Diff_Boost = Daten.Air_pInMnfBoostBas_VW_ETKC_1(:,2) - Daten.GEMChCont_pInMnf_VW_ETKC_1(:,2);
    Diff_Des = Daten.GECGsl_pInMnfDes_VW_ETKC_1(:,2) - Daten.GEMChCont_pInMnf_VW_ETKC_1(:,2);
    

    % Messung in geeignete Länge schneiden
    for ip = 1 : 3
        if length(Schnittpunkt) >= ip
            if ip == 1 && (Schnittpunkt(ip)-Fenster < 0)
            % Nur für den Fall, dass der erste Sprung ganz am Anfang
            % der Messung sich befindet.
                dataInOn(ip,end-Schnittpunkt(ip)+1:end) = dataIn(1:Schnittpunkt(ip));
                dataOutOn(ip,end-Schnittpunkt(ip)+1:end) = dataOut(1:Schnittpunkt(ip));
                
                dataMeAGROn(ip,end-Schnittpunkt(ip)+1:end) = Modell_eAGR(1:Schnittpunkt(ip));         
                dataFriLuftOn(ip,end-Schnittpunkt(ip)+1:end) = Daten.GEMChCont_rAirCh_VW_ETKC_1(1:Schnittpunkt(ip),2);
                
                dataDiff_Boost_On(ip,end-Schnittpunkt(ip)+1:end) = Diff_Boost(1:Schnittpunkt(ip));
                dataDiff_Des_On(ip,end-Schnittpunkt(ip)+1:end) = Diff_Des(1:Schnittpunkt(ip));
                
                Mittel_bCrit(i,ip) = mean(Daten.GECGsl_bCritAcvCtlTvh_pInMnf_VW_ETKC_1(1:Schnittpunkt(ip),2));
            else
                dataMeAGROn(ip,:) = Modell_eAGR(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip));
                dataFriLuftOn(ip,:) = Daten.GEMChCont_rAirCh_VW_ETKC_1(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip),2);
                
                dataInOn(ip,:) = dataIn(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip));
                dataOutOn(ip,:) = dataOut(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip));
                
                dataDiff_Boost_On(ip,:) = Diff_Boost(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip));
                dataDiff_Des_On(ip,:) = Diff_Des(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip));
                
                Mittel_bCrit(i,ip) = mean(Daten.GECGsl_bCritAcvCtlTvh_pInMnf_VW_ETKC_1(Schnittpunkt(ip)-Fenster:Schnittpunkt(ip),2));
            end
        
        % Negativer Spruenge
        dataInOff(ip,:) = dataIn(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster);
        dataOutOff(ip,:) = dataOut(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster); 
        
        dataMeAGROff(ip,:) = Modell_eAGR(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster);
        dataFriLuftOff(ip,:) = Daten.GEMChCont_rAirCh_VW_ETKC_1(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster,2);
        
        dataDiff_Boost_Off(ip,:) = Diff_Boost(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster);
        dataDiff_Des_Off(ip,:) = Diff_Des(Schnittpunkt(ip):Schnittpunkt(ip)+Fenster);
        end
    end
    
    IdxSprungOn = zeros(size(dataInOn,1),1);
    IdxSprungOff = zeros(size(dataInOff,1),1);
    
    % Ergebnisse von DES
    Diff_Des_On(i) = mean(max(dataDiff_Des_On,[],2));
    Diff_Des_Off(i) = mean(min(dataDiff_Des_Off,[],2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------- Postiver Sprung ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AdjDauer = 55*round(3500/mean(Drehzahl(i,:))); % Zeitfenster variiert sich adaptiv 
    for ip = 1: length(IdxSprungOn)
        if length(posOnTt) >= ip
        IdxSprungOn(ip)  = posOnTt(ip)- (Schnittpunkt(ip)-Fenster)+1; 
        [MaxSprungOn, Imax] = max(dataOutOn(ip,IdxSprungOn(ip):IdxSprungOn(ip)+AdjDauer));
        [MinSprungOn, Imin] = min(dataOutOn(ip,IdxSprungOn(ip):IdxSprungOn(ip)+AdjDauer));
        if abs(MaxSprungOn - dataOutOn(ip,IdxSprungOn(ip))) > abs(MinSprungOn - dataOutOn(ip,IdxSprungOn(ip)))
            SpritzIdx = Imax;
        else
            SpritzIdx = Imin;
        end
        
        FriLuftOn(i,ip) = mean(dataFriLuftOn(ip,1:IdxSprungOn));
        Diff_Boost_On(i,ip) = mean(dataDiff_Boost_On(ip,IdxSprungOn-200:IdxSprungOn)); % ---> Speichern
        
        n = IdxSprungOn(ip) + SpritzIdx;
        MaxAbweichungOn(i,ip) = dataOutOn(ip,n) - dataOutOn(ip,IdxSprungOn(ip));
        
        while abs(dataOutOn(ip,n)- dataOutOn(ip, IdxSprungOn(ip))) >5e-4 &&...
              abs(dataOutOn(ip,n)- 1) >5e-4 && n < size(dataOutOn,2)
            n = n + 1;
        end                
        FensterOnEnd = n;
 
        % --- Quadratische Fehler ---
        QuaErrOn(i,ip) = sum((dataOutOn(ip,IdxSprungOn(ip):FensterOnEnd)-1).^2);
        end
    end
    
    if dataMeAGROn(:,n+10) < Warning
        warning('Die %03i. Messung hat eventuell keine Daten drin',i);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------------  Negativer Sprung ----------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ip = 1 : length(IdxSprungOff)
        if length(posOffTt) >= ip
        IdxSprungOff(ip)  = posOffTt(ip)- Schnittpunkt(ip)+1;
        [MaxSprungOff, Imax] = max(dataOutOff(ip,IdxSprungOff(ip):IdxSprungOff(ip)+AdjDauer));
        [MinSprungOff, Imin] = min(dataOutOff(ip,IdxSprungOff(ip):IdxSprungOff(ip)+AdjDauer));
        if abs(MaxSprungOff - dataOutOff(ip,IdxSprungOff(ip))) > abs(MinSprungOff - dataOutOff(ip,IdxSprungOff(ip)))
            SpritzIdx = Imax;
        else
            SpritzIdx = Imin;
        end
        
        FriLuftOff(i,ip) = mean(dataFriLuftOff(ip,1:IdxSprungOff));   % ---> Speichern
        Diff_Boost_Off(i,ip) = mean(dataDiff_Boost_Off(ip,IdxSprungOff-200:IdxSprungOff)); % ---> Speichern
        
        n = IdxSprungOff(ip) + SpritzIdx;
        MaxAbweichungOff(i,ip) = dataOutOff(ip,n) - dataOutOff(ip,IdxSprungOff(ip));
        
        while (n < size(dataOutOff,2) &&...
                abs(dataOutOff(ip,n)- dataOutOff(ip, IdxSprungOff(ip))) >1e-3 &&...
              abs(dataOutOff(ip,n)- 1) >1e-3)
            n = n + 1;
        end
        FensterOffEnd = n;
        % --- Quadratische Fehler ---
        QuaErrOff(i,ip) = sum((dataOutOff(ip,IdxSprungOff(ip):FensterOffEnd)-1).^2);
        end
    end

    
    % --- Standardabweichung ---
    StdQErrOn(i) = std(QuaErrOn(i,:)); % Standardabweichung der Quadratischen Fehler
    StdQErrOff(i) = std(QuaErrOff(i,:));
    StdMAbwOff(i) = std(MaxAbweichungOff(i,:)); %Std der max. Abweichung
    StdMAbwOn(i) = std(MaxAbweichungOn(i,:));
    end
end


%% Beste Parameterkombination auswerten
Idx = (1: length(Drehzahl))';

% --- Arbeispunkt genegieren
DrehzahlOnr = round(Drehzahl(:,1)/100)*100;

%
diffFriLuft = FriLuft(1:end-1) - FriLuft(2:end);

% Hier die Index von Arbeitspunktswechsel finden.
Wechsel = [0;Idx(abs(diffFriLuft)>=7)];
IdxDrehzahlWecksel = [0;Idx(abs(diff(DrehzahlOnr))>150)];

Wechsel(end+1:end+length(IdxDrehzahlWecksel)) = IdxDrehzahlWecksel;
Wechsel = unique(Wechsel); % zwei mal vorgekommene Punkte weg

% Index der besten Parameterkombination
IdxBestOn = zeros(length(Wechsel),3);
IdxBestOff = zeros(length(Wechsel),3);

% Verteilung der Parameterkombintaion in 5*3 Matrix
% 13 14 15
% 10 11 12
% 7  8  9
% 4  5  6
% 1  2  3 
Idx_in15_BestOn = zeros(length(Wechsel),3);
Idx_in15_BestOff = zeros(length(Wechsel),3);

% 
BestMaxAbwOn = zeros(length(Wechsel),3);
BestMaxAbwOff = zeros(length(Wechsel),3);

%
for k = 1 : length(Wechsel)
    if k < length(Wechsel)
        Bereich = Wechsel(k)+1:Wechsel(k+1);
    else
        Bereich = Wechsel(k)+1:length(FriLuft);
    end
    for j = 1 :3
        [~, idxBestOn] = min(abs(MaxAbweichungOn(Bereich,1)));
        [~, idxBestOff] = min(abs(MaxAbweichungOff(Bereich,1)));
        Idx_in15_BestOn(k,j) = idxBestOn;
        Idx_in15_BestOff(k,j) = idxBestOff;
        IdxBestOn(k,j) = idxBestOn + Bereich(1) -1;
        IdxBestOff(k,j) = idxBestOff + Bereich(1) -1;
        
        BestMaxAbwOn(k,j) = MaxAbweichungOn(IdxBestOn(k,j),j);
        BestMaxAbwOff(k,j) = MaxAbweichungOff(IdxBestOff(k,j),j);
    end
    
    % Datainame der besten Parameterkombination speichern
    IdxOn = sprintf('.%03i.dat',IdxBestOn(k,1));
    IdxOff = sprintf('.%03i.dat',IdxBestOff(k,1));
    DateiNameOn(k,:) = strcat(filename,IdxOn);
    DateiNameOff(k,:) = strcat(filename,IdxOn);
    
end

m_BestMaxAbwOn = median(BestMaxAbwOn,2);
m_BestMaxAbwOff = median(BestMaxAbwOff,2);

DrehzahlOnr = round(Drehzahl/100)*100;

% Drehzahl
Best_DrehzahlOn = DrehzahlOnr(IdxBestOn(:,1));
Best_DrehzahlOff = DrehzahlOnr(IdxBestOff(:,1));

% Frischfüllung kurz vor den Sprüngen von beste Parameterkombination
Best_FriLuftOn = FriLuftOn(IdxBestOn(:,1),:);
Best_FriLuftOff = FriLuftOff(IdxBestOff(:,1),:);

% Totzeit 
Best_TtOn = Ttlist(IdxBestOn);
Best_TtOff = Ttlist(IdxBestOff);

% Zeitkonstante
Best_T1On = T1list(IdxBestOn);
Best_T1Off = T1list(IdxBestOff);

% Boost kurz vor den Sprüngen von beste Parameterkombination
Best_Boost_On = mean(Diff_Boost_On(IdxBestOn),2);
Best_Boost_Off = mean(Diff_Boost_Off(IdxBestOn),2);

% Des von beste Parameterkombination
Best_Des_On = Diff_Des_On(IdxBestOn(:,1));
Best_Des_Off = Diff_Des_Off(IdxBestOff(:,1));

Best_Mittel_bCrit = Mittel_bCrit(IdxBestOn);
Best_MittelDruck = MittelDruck(IdxBestOn(:,1));

% Wichtige Merkmale Speichern
Datenfilename = strcat(pathname,'\Daten_',datafile,'.mat');
save(Datenfilename,...
    'Best_DrehzahlOn',...
    'Best_DrehzahlOff',...
    'Best_FriLuftOn',...
    'Best_FriLuftOff',...
    'Best_TtOn',...
    'Best_TtOff',...
    'Best_T1On',...
    'Best_T1Off',...
    'm_BestMaxAbwOn',...
    'm_BestMaxAbwOff',...
    'BestMaxAbwOn',...
    'BestMaxAbwOff',...
    'Idx_in15_BestOn',...
    'Idx_in15_BestOff',...
    'Best_Mittel_bCrit',...
    'Best_MittelDruck',...
    'Best_Boost_On',...
    'Best_Boost_Off',...
    'Best_Des_On',...
    'Best_Des_Off',...
    'DateiNameOn',...
    'DateiNameOff'...
    );


% Tabelle erzeugen
   Drehzahl =  Best_DrehzahlOn;
   FriLuftOn = mean(Best_FriLuftOn,2);
   FriLuftOff = mean(Best_FriLuftOff,2);
   TtOn =  Best_TtOn(:,1);
   TtOff = Best_TtOff(:,1);
   T1On = Best_T1On(:,1);
   T1Off = Best_T1Off(:,1);
   m_MaxAbwOn = m_BestMaxAbwOn;
   m_MaxAbwOff = m_BestMaxAbwOff;
   MaxAbwOn = BestMaxAbwOn;
   MaxAbwOff = BestMaxAbwOff;
   Mittel_bCrit = Best_Mittel_bCrit;
   MittelDruck = Best_MittelDruck;
   Boost_On = Best_Boost_On;
   Boost_Off = Best_Boost_Off;
   Des_On = Best_Des_On;
   Des_Off = Best_Des_Off;
   DateiNameOn = cellstr(DateiNameOn);
   DateiNameOff = cellstr(DateiNameOff);
    
   Header = {'Drehzahl',...
       'FriLuftOn',...
       'FriLuftOff',...
       'TtOn',...
       'TtOff',...
       'T1On',...
       'T1Off',...
       'MittelDruck',...
       'Boost_On',...
       'Boost_Off',...
       'Des_On',...
       'Des_Off',...
       'DateiNameOn',...
       'DateiNameOff'}; 
   
   
   Data2write = [Drehzahl,...
       FriLuftOn,...
       FriLuftOff,...
       TtOn,...
       TtOff,...
       T1On,...
       T1Off,...
       MittelDruck,...
       Boost_On,...
       Boost_Off,...
       Des_On,...
       Des_Off];
   
   
   Filename2write = [DateiNameOn,...
       DateiNameOff];
       
   system('taskkill /F /IM EXCEL.EXE');
   xlsname = strcat('D:',pathname,'\',filename,'.xlsx');
   xlswrite(xlsname,Header)
   xlswrite(xlsname,Data2write,'Tabelle1','A2');
   xlswrite(xlsname,Filename2write,'Tabelle1','M2');
   
   
