%%%%% From "Select_Mag_Phase_Enhanced_2FocalLength_New.m" 2016.10.17
%%%%% Use the manually selected parameters

clear all;
clc;
close all;

% filename1=['./data/SingleLayer_P120_075THz_4'];   %%%
% load(strcat(filename1,'.mat'))
% Mag_S1  =Mag_Phase1(:,1);
% Deg_S1  =Mag_Phase1(:,2);
%
% filename2=['./data/SingleLayer_P120_04THz_1'];  %%%
% load(strcat(filename2,'.mat'))
% Mag_S2  =Mag_Phase2(:,1);
% Deg_S2  =Mag_Phase2(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Begin: Calculate Required Phase for 1-D lens with total N cells: P=120um,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Total_N=151;
P=120;

% x1=0.5*P+P*[0:floor(Total_N/2)-1];   %%% um
x1=P*[0:floor(Total_N/2)];   %%% um

F1=5.25e3; %%%% um: focal length 1  7  Lbd  %%%% lambda=750 um @0.4T
F2=13.5e3;  %%% um: focal length 2: 18 Lbd

TwoFocal=0;  %%% =0: 1 focal point; =1: 2 focal point

%%%%  For frequency 1: high %%%%%%%%%%%%%
f1=0.75; %THz
Lam1=300/f1;   %%% um
% Phi_1=2*pi/Lam1*(sqrt(x1.^2+F^2)-F);

Phi1_1=2*pi/Lam1*(sqrt(x1.^2+F1^2)-F1);
Phi2_1=2*pi/Lam1*(sqrt(x1.^2+F2^2)-F2);

if(TwoFocal==1)
    AE1=0.5*(exp(1i*Phi1_1)+exp(1i*Phi2_1));
else
    AE1=exp(1i*Phi1_1);
end
%AE1=0.5*(exp(-1i*Phi1_1)+exp(-1i*Phi2_1));  %%% NOT WORK with "-" sign %%%%%%%[important]

Amp1=abs(AE1)';
The1=asind(Amp1)/2;  %%% rotation angle NEEDED
Phi_1=angle(AE1)';
Phi_deg1=wrapTo360(Phi_1*180/pi);

%%%%  For frequency 2: Low %%%%%%%%%%%%%
f2=0.4; % THz
Lam2=300/f2;   %%% um @0.4 THz

Phi1_2=2*pi/Lam2*(sqrt(x1.^2+F1^2)-F1);
Phi2_2=2*pi/Lam2*(sqrt(x1.^2+F2^2)-F2);

if(TwoFocal==1)
    AE2=0.5*(exp(1i*Phi1_2)+exp(1i*Phi2_2));  %%% For Two foc
else
    AE2=exp(1i*Phi1_2);
end

%AE1=0.5*(exp(-1i*Phi1_1)+exp(-1i*Phi2_1));  %%% NOT WORK with "-" sign %%%%%%
Amp2=abs(AE2)';
The=asind(Amp2)/2;
Phi_2=angle(AE2)';
Phi_deg2=wrapTo360(Phi_2*180/pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End: Calculate Required Phase for 1-D lens with total 151 cells: P=120um,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phase1_R=Phi_deg1;   %%% Required phase for HIGH frequency
Phase2_R=Phi_deg2;   %%% Required phase for LOW frequency
TotalN=length(Phase1_R);
N=3;
Mag_DigitArray=linspace(0,1,N);

indx=1:1:length(Amp1);
indx=indx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For HIGH frequency: 0.75 THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=0;

Phase1_DigitArray=k1+[0:45:360];

phase1=[];
idx1=[];
mag1=[];

for II=1:length(Phi_deg1)
    phase_t=Phi_deg1(II);
    [~,idx]=min(abs(phase_t-Phase1_DigitArray));
    phase_d=Phase1_DigitArray(idx);
    phase1=[phase1;phase_d];
    idx1=[idx1;idx];
end

idx1_4=find(idx1==9); %%% find out the "-": phase index >4
idx1(idx1_4)=idx1(idx1_4)-8;

idx1_4=find(idx1>4); %%% find out the "-": phase index >4
idx1(idx1_4)=idx1(idx1_4)-4;


for II=1:length(Amp1)
    mag_t=Amp1(II);
    [~,idx]=min(abs(mag_t-Mag_DigitArray));
    mag_d=Mag_DigitArray(idx);
    mag1=[mag1;mag_d];
end

theta1=asind(Amp1)/2;

theta1(idx1_4)=-1*theta1(idx1_4);

all_1=[indx Phi_deg1 phase1 idx1 Amp1 mag1 theta1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For LOW frequency: 0.4 THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2=0; %%% xxxxxxxx

Phase2_DigitArray=k2+[0:45:360];

phase2=[];
idx2=[];
mag2=[];

for II=1:length(Phi_deg2)
    phase_t=Phi_deg2(II);
    [~,idx]=min(abs(phase_t-Phase2_DigitArray));
    phase_d=Phase2_DigitArray(idx);
    phase2=[phase2;phase_d];
    idx2=[idx2;idx];
end

idx2_4=find(idx2==9); %%% find out the "-": phase index >4
idx2(idx2_4)=idx2(idx2_4)-8;

idx2_4=find(idx2>4); %%% find out the "-": phase index >4

idx2(idx2_4)=idx2(idx2_4)-4;


for II=1:length(Amp2)
    
    mag_t=Amp2(II);
    [~,idx]=min(abs(mag_t-Mag_DigitArray));
    mag_d=Mag_DigitArray(idx);
    mag2=[mag2;mag_d];
    
end

theta2=asind(Amp2)/2;

theta2(idx2_4)=-1*theta2(idx2_4);

all_2=[indx Phi_deg2 phase2 idx2 Amp2 mag2 theta2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileID_Para1 = fopen('Required_phase_1Freq_2Focal.txt','w');
% 
% for II=1:length(Phi_deg1)
%     fprintf(fileID_Para1,'The index number: %3d \r\n \r\n',II);
%     fprintf(fileID_Para1,'scshape_Angle = "%4.1f" \r\n',alpha);
%     fprintf(fileID_Para1,'scshape_Width = "%3.1f" \r\n', 5);
%     fprintf(fileID_Para1,'scshape_OuterRadius = "%3.1f" \r\n',r);
%     fprintf(fileID_Para1,'scshape_RotationAngle = "%4.1f" \r\n',the_sign*The1_Cal);
%     fprintf(fileID_Para1,'scslot_Angle = "%4.1f" \r\n', 60);
%     fprintf(fileID_Para1,'scslot_Width = "%3.1f" \r\n', 5);
%     fprintf(fileID_Para1,'scslot_OuterRadius = "%3.1f" \r\n',55);
%     fprintf(fileID_Para1,'scslot_RotationAngle = "%4.1f" \r\n \r\n',the1_sign*45);
% end
% 
% fclose(fileID_Para1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% fileID=fopen('ToCST_8phase_Nmag.txt','w');
% 
% for II=1:length(idx1)
%     fprintf(fileID,'%d,',idx1(II));
% end
% 
% fprintf(fileID,'\r\n');
% for II=1:length(idx1)
%     fprintf(fileID,'%.1f,',theta1(II));
% end
% fprintf(fileID,'\r\n');
% fclose(fileID);



% SWITHCH=1;   %%% 1 for HIGH frequency; 2 for LOW frequency
% HIGH=0;
% LOW=1;
% %%%========================================================
% %%%  DO NOT DELETE: FOR HIGH Frequency
% %%%========================================================
%
% if (HIGH==1)
%
%     fileID = fopen('Required_phase_High.txt','w');
%
%     fileID_Para1 = fopen('Required_phase_High_Para.txt','w');
%
%     fprintf(fileID,'Required_Phase_For_High_Frequency \r\n');
%     fprintf(fileID,'----------------DigitizedP----Phase---Mag--IDX-------- --- r ---- alpha --- the ---------------the1---------\r\n');
%
%     for k=16
%
%         display(['=================   new iteration 1: k=' num2str(k) '  ======================'])
%
%         Phase1_DigitArray=k+[0:45:315];
%
%         phase1_d=wrapTo360(k+Phase1_R);
%         % phase2_d=k+Phase2_R;
%
%         for II=1:length(phase1_d)
%
%             phase_d=phase1_d(II);
%
%             [~,idx]=min(abs(phase_d-Phase1_DigitArray));
%
%             phase_t=Phase1_DigitArray(idx);
%
%             dif=2;
%
%             X1=find(Deg_S1<phase_t+dif & Deg_S1>phase_t-dif);  %% index
%
%             Y1=Mag_S1(X1)';
%             Idx1=find( Y1==(max(Y1)) );
%             Angle1=Deg_S1(X1(Idx1));
%
%             fprintf('%3d reqP=%6.2f reqD=%6.2f\t',II,phase_d,phase_t);
%             %             fprintf('(%8.3f %8.3f  %d)\n',Angle1, max(Y1), X1(Idx1));
%             fprintf('(%6.2f %6.2f  %3d  reqM=%.3f)\n',Angle1, Y1(Idx1), X1(Idx1),Amp1(II));
%
%             if(II==1)
%                 A1=Y1(Idx1);
%             end
%
%             The1_Cal=real(asind(A1/Y1(Idx1)*Amp1(II))/2); %%% A1 is the first one of the Transmission amplitude, corresponding to Max Amp=1;
%
%             fprintf(fileID,'%3d reqP=%6.2f reqD=%6.2f\t',II,phase_d,phase_t);
%             % fprintf(fileID,'%8.3f %8.3f  %3d \t |==%3d==  req=%.2f\r\n',Angle1, Y1(Idx1), X1(Idx1), II, Amp(II));
%
%             [alpha,r,the_sign,the1_sign]=FindParameter075(X1(Idx1));
%
%             fprintf(fileID,'%6.2f %6.2f  %3d  reqM=%.2f r=%2d \t alp= %3d |the=%.1f=|\t%3d |the1=%.1f| \r\n',Angle1, Y1(Idx1), X1(Idx1),  Amp1(II),r,alpha,the_sign*The1_Cal, II,the1_sign*45);
%             fprintf(fileID,'-------------------------------------------------------------------------------------------------------------\r\n');
%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             fprintf(fileID_Para1,'The index number: %3d \r\n \r\n',II);
%             fprintf(fileID_Para1,'scshape_Angle = "%4.1f" \r\n',alpha);
%             fprintf(fileID_Para1,'scshape_Width = "%3.1f" \r\n', 5);
%             fprintf(fileID_Para1,'scshape_OuterRadius = "%3.1f" \r\n',r);
%             fprintf(fileID_Para1,'scshape_RotationAngle = "%4.1f" \r\n',the_sign*The1_Cal);
%             fprintf(fileID_Para1,'scslot_Angle = "%4.1f" \r\n', 60);
%             fprintf(fileID_Para1,'scslot_Width = "%3.1f" \r\n', 5);
%             fprintf(fileID_Para1,'scslot_OuterRadius = "%3.1f" \r\n',55);
%             fprintf(fileID_Para1,'scslot_RotationAngle = "%4.1f" \r\n \r\n',the1_sign*45);
%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         end
%     end
%
%     fclose(fileID);
%
%     fclose(fileID_Para1);
%
% end %%% For HIGH frequency
%
% %%%========================================================
% %%%  DO NOT DELETE: FOR LOW Frequency
% %%%========================================================
%
% if (LOW==1)
%
%     for k=5
%
%         display(['=================   new iteration 2: k=' num2str(k) '  ======================'])
%
%         fileID = fopen('Required_phase_LOW.txt','w');
%
%         fileID_Para1 = fopen('Required_phase_LOW_Para.txt','w');
%
%         fprintf(fileID,'Required_Phase_For_High_Frequency \r\n');
%         fprintf(fileID,'----------------DigitizedP----Phase---Mag--IDX-------- --- r ---- alpha --- the ---------------the1---------\r\n');
%
%         Phase2_DigitArray=k+[0:45:315];
%
%         phase2_d=wrapTo360(k+Phase2_R);  %%% 55 Cells: 0.40 THz
%
%         for II=1:length(phase2_d)
%
%             phase_d=phase2_d(II);
%
%             [~,idx]=min(abs(phase_d-Phase2_DigitArray));
%
%             phase_t=Phase2_DigitArray(idx);
%
%             dif=2;
%
%             X2=find(Deg_S2<phase_t+dif & Deg_S2>phase_t-dif);  %% index
%
%             Y2=Mag_S2(X2)';
%
%
%             if(II==1)
%                 Idx2=find( Y2==(max(Y2)) );
%                 A2=Y2(Idx2);
%             else
%                 [~,Idx2]=min(abs(A2-Y2));
%             end
%
%             Angle2=Deg_S2(X2(Idx2));
%
%             fprintf('Phase == %2d ==: %8.3f\t',II, phase_t);
%             fprintf('(%8.3f %8.3f  %3d  reqM=%.3f)\n',Angle2, Y2(Idx2), X2(Idx2),Amp(II)); %%% Last one is the required magnitude
%
%
%
%             The_Cal=real(asind(A2/Y2(Idx2)*Amp(II))/2); %%% A2 is the first one of the Transmission amplitude, corresponding to Max Amp=1;
%
%             fprintf(fileID,'==%3d== %8.3f\t',II, phase_t);
%
%
%             [alpha1,the_sign,the1_sign]=FindParameter04(X2(Idx2));
%
%             fprintf(fileID,'%6.2f %6.2f  %3d  reqM=%.2f r1=%2d \t alp= %3d |the=%.1f=|\t%3d |the1=%.1f| \r\n',Angle2, Y2(Idx2), X2(Idx2),  Amp(II),55,alpha1,the_sign*45, II,the1_sign*The_Cal); %%r1=55
%             fprintf(fileID,'-------------------------------------------------------------------------------------------------------------\r\n');
%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             fprintf(fileID_Para1,'The index number: %3d \r\n \r\n',II);
%             fprintf(fileID_Para1,'scshape_Angle = "%4.1f" \r\n',62);
%             fprintf(fileID_Para1,'scshape_Width = "%3.1f" \r\n', 5);
%             fprintf(fileID_Para1,'scshape_OuterRadius = "%3.1f" \r\n',30);
%             fprintf(fileID_Para1,'scshape_RotationAngle = "%4.1f" \r\n',the_sign*45);
%             fprintf(fileID_Para1,'scslot_Angle = "%4.1f" \r\n', alpha1);
%             fprintf(fileID_Para1,'scslot_Width = "%3.1f" \r\n', 5);
%             fprintf(fileID_Para1,'scslot_OuterRadius = "%3.1f" \r\n',55);
%             fprintf(fileID_Para1,'scslot_RotationAngle = "%4.1f" \r\n \r\n',the1_sign*The_Cal);
%
%             %
%         end
%         %     pause;
%     end
%
%     fclose(fileID);
%     fclose(fileID_Para1);
%
% end
%
%

