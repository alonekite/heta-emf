clear
clc
close all

tic
% Defining Variables & Statement
% Subscripts use underscores between symbols & Vectors use repeating letters (eg: HH, some are not)


% Defining the parameters of the antenna array (number of row elements, number of column elements, interval)
N_row = 1;   % The number of elements in a row of an antenna array
N_col = 4;   % The number of elements in a column of an antenna array
N_ant_ele = N_row * N_col;   % The number of antenna array elements
Tilt_deg = 0;   % The tilt angle of the antenna array (eg: 4 deg)
R_Ant_Ein = 0.3;    % Distance from the antenna element closest to the center to the center

method = 'MRT'; % MRT, ZF, CB
% n - Tx; k - Rx
% h_kn : channel coefficient
% u_r : voltage amplitude
% HH : normalizing channel vector
% WW: procoding matrix
% WW_CB : codebook steering matrix
% Dist_BStoUE : distance from the center of BS to UE
% dd : distance vector form the BS array center to the nth element
% cc : unit vector in the codebook direction assigned to the kth UE (cc_k)
% Supplementary information for cc: Unit_c is cc_k (not the cc_k in Fomula (6))
% cc_k_Max : the chosen beamforming direction (Fomula (6))
% tt : transmit vector
% ss : transmitted symbol

% define the numbers of header lines
N_BasicInfo_S = 15; % Where the basic information starts (row)
%N_BasicInfo_E = 22; % Where the basic information ends (row)
N_header = 51; % number of header lines (for the ray)
% define the numbers of header lines of 'Pattern'
N_header_Patt = 7;

file_name = 'Site  3 Antenna 1 Rays.str';
[info,point] = get_winprop_ray(file_name,N_BasicInfo_S,N_header);
% read elements' patterns
[Pattern_Eles] = get_element_pattern(N_header_Patt);
clear freq key n_path n_point nn sss
clear tline a val

Ant_Position = info.Ant_Position;
Ant_H = info.Ant_H;
Ant_V = info.Ant_V;
Ant_outPow = info.Ant_outPow;
freq_M = info.freq_M;
freq = info.freq ;
lambda = info.lambda;

% Calculate the sum of the channel coefficient h_kn, the normalizing ?H column-wise and direction vector cc_i
num_UEs = size(point);    % num_UEs(2) is  how many UEs (or points)
h_kn_noTrans = zeros(num_UEs(2),N_ant_ele);
HH = zeros(num_UEs(2),N_ant_ele);
Unit_c = zeros(num_UEs(2),3);
for ii = 1 : num_UEs(2)
    num_Path_on1Point = size(point(ii).path);   % num_Path_on1Point(2) is num of path for one UE(point)
    Mat_Path_on1Point = zeros(num_Path_on1Point(2),1);    % Given an empty column vector for quick sum
    for jj = 1 :  num_Path_on1Point(2)
        Mat_Path_on1Point(jj) = point(ii).path(jj).h_kn_on1path;
    end
    h_kn_point(ii) = sum(Mat_Path_on1Point);   % Formula(1)
    % Given the matrix of h_kn
    h_kn_noTrans(ii,:) = deal(h_kn_point(ii));   % I'm not sure that h_k1=h_k2=......=h_kn?????
    h_kn_Frob_point(ii) = norm((h_kn_noTrans(ii,:))','fro');
    HH(ii,:) = h_kn_noTrans(ii,:) / h_kn_Frob_point(ii);   % Formula(2)
    Dist_BStoUE(ii) = pdist2(Ant_Position,point(ii).position);
    Unit_c(ii,:) = Ant_Position - point(ii).position;   % Calculate the direction vector between the center of BS and point(ii) (not the unit vector)
    % Calculate the unit vector in rectangular coordinate system between the center of BS and point(ii)
    Norm_Unit_c(ii,:) = norm(Unit_c(ii,:));
    Unit_c(ii,:) = Unit_c(ii,:) / Norm_Unit_c(ii,:);   % Unit_c is cc_k (the unit vector)
end
h_kn = (h_kn_noTrans)';
clear ii jj num_Path_on1Point Mat_Path_on1Point
clear Norm_Unit_c Theta_point Phi_point h_kn_point h_kn_Frob_point
clear h_kn_noTrans

% Given the physical distance D_ele and the electrical distance Dr_ele
D_ele = 0.5; % physical distance (Unit:m)
% Electrical Distance
Dr_ele = 2*pi*D_ele/lambda;

% Calculate the procoding matrices ?W (Fomula (3) & (5))
alfa = norm(HH);   % alfa is the real-valued normalization coefficient
HH_H = HH';   % hermitian transpose: conj(A.') (also conj(A).')
% Given a transition GtoE matrix & synthesis matrix
[Row_Pattern_Eles, Col_Pattern_Eles] = size(Pattern_Eles);
Pattern_Eles_GtoE = zeros(Row_Pattern_Eles, Col_Pattern_Eles);
Pattern_syn = zeros(Row_Pattern_Eles, 3);
switch method
    case 'MRT'
        WW_MRF = alfa * HH_H;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First-Calculate the GtoE Transitional Pattern
        % Step1-transform G(gain) to E(field strength or amplitude of field component)
        % Step2-E*phase information (The phase information of each element is different, with the center as the origin)
        % Step3-E*WW_MRF
        % Step4-sum(E*WW_MRF) and to G finally
            % Select the XOZ-plane(phi=0) and YOZ-plane(phi=90)
            % Find the coordinate and values (Theta & Gain)
            % Note: the First column is THETA, the second column is PHI
            Pattern_XOZ = zeros(181, 2);   % 0-180 (181 numbers)
            Pattern_YOZ = zeros(181, 2);
        % Combine the 4 patterns
        for ii = 1:Row_Pattern_Eles
            Pattern_Eles_GtoE(ii,1) = Pattern_Eles(ii,1);
            Pattern_Eles_GtoE(ii,2) = Pattern_Eles(ii,2);
            Pattern_Eles_GtoE(ii,3) = 20^(Pattern_Eles(ii,3)/10);
            Pattern_Eles_GtoE(ii,3) = Pattern_Eles_GtoE(ii,3)*exp(-3/2*j*Dr_ele);
            Pattern_Eles_GtoE(ii,3) = Pattern_Eles_GtoE(ii,3)*WW_MRF(1);
            Pattern_Eles_GtoE(ii,4) = 20^(Pattern_Eles(ii,4)/10);
            Pattern_Eles_GtoE(ii,4) = Pattern_Eles_GtoE(ii,4)*exp(-1/2*j*Dr_ele);
            Pattern_Eles_GtoE(ii,4) = Pattern_Eles_GtoE(ii,4)*WW_MRF(2);
            Pattern_Eles_GtoE(ii,5) = 20^(Pattern_Eles(ii,5)/10);
            Pattern_Eles_GtoE(ii,5) = Pattern_Eles_GtoE(ii,5)*exp(1/2*j*Dr_ele);
            Pattern_Eles_GtoE(ii,5) = Pattern_Eles_GtoE(ii,5)*WW_MRF(3);
            Pattern_Eles_GtoE(ii,6) = 20^(Pattern_Eles(ii,6)/10);
            Pattern_Eles_GtoE(ii,6) = Pattern_Eles_GtoE(ii,6)*exp(3/2*j*Dr_ele);
            Pattern_Eles_GtoE(ii,6) = Pattern_Eles_GtoE(ii,6)*WW_MRF(4);
            % Second-Calculate the Synthetic Pattern
              Pattern_syn(ii,1) = Pattern_Eles(ii,1);
              Pattern_syn(ii,2) = Pattern_Eles(ii,2);
              Pattern_syn(ii,3) = sum(Pattern_Eles_GtoE(ii,3:6));
              Pattern_syn(ii,3) = abs(Pattern_syn(ii,3));
              Pattern_syn(ii,3) = 20*log10(Pattern_syn(ii,3));
              % Select the XOZ-plane(phi=0) and YOZ-plane(phi=90) and Gains
                 Pattern_XOZ(:,1) = Pattern_syn(1:360:Row_Pattern_Eles,1);
                 Pattern_XOZ(:,2) = Pattern_syn(1:360:Row_Pattern_Eles,3);
                 Pattern_YOZ(:,1) = Pattern_syn(91:360:Row_Pattern_Eles,1);
                 Pattern_YOZ(:,1) = Pattern_syn(91:360:Row_Pattern_Eles,3);
        end
        clear ii Col_Pattern_Eles
        clear Pattern_Eles_GtoE
        
    case 'ZF'
        WW_ZF = alfa * HH_H * inv(HH*HH_H);   % Warning: 矩阵接近奇异值，或者缩放错误。结果可能不准确。        
    case 'CB'
        % dd : distance vector from the BS array center to the nth element
        % Noted that the arrangement order of antenna elements should be consistent with the element order of matrix in MATLAB (Columnar distribution or Counter-clockwise direction)
        % The distance between the two antenna element arrays is: 2*R_Ant_Ein (VERY IMPORTANT!!!)
        % Function rot() is counter-clockwise direction
        dd = zeros(N_ant_ele,3);
        d_El_even1 = R_Ant_Ein*cosd(45);        % When the number of elements is even, the radius from the FIRST circle element to the center
        d_El_even2 = R_Ant_Ein*3*cosd(45);    % When the number of elements is even, the radius from the SECOND circle element to the center
        d_El_odd1 = R_Ant_Ein*2;                     % When the number of elements is odd, the radius from the FIRST/SECOND circle element to the center
        if Tilt_deg == 0
            if N_row == N_col
                switch N_row
                    case 2
                        dd(1,:) = [0, d_El_even1, -d_El_even1];
                        dd(2,:) = [0, -d_El_even1, -d_El_even1];
                        dd(3,:) = [0, d_El_even1, -d_El_even1];
                        dd(4,:) = [0, d_El_even1, d_El_even1];
                    case 3
                        dd(1,:) = [0, -d_El_odd1, d_El_odd1];
                        dd(2,:) = [0, -d_El_odd1, 0];
                        dd(3,:) = [0, -d_El_odd1, -d_El_odd1];
                        dd(4,:) = [0, 0, d_El_odd1];
                        dd(5,:) = [0, 0, 0];
                        dd(6,:) = [0, 0, -d_El_odd1];
                        dd(7,:) = [0, d_El_odd1, d_El_odd1];
                        dd(8,:) = [0, d_El_odd1, 0];
                        dd(9,:) = [0, d_El_odd1, -d_El_odd1];
                    otherwise
                        dd = 0;
                end
            elseif (N_row ~= N_col)
                dd(1,:) = [0, -R_Ant_Ein*3, 0];
                dd(2,:) = [0, -R_Ant_Ein, 0];
                dd(3,:) = [0, R_Ant_Ein, 0];
                dd(4,:) = [0, R_Ant_Ein*3, 0];
            else
                dd = 0;
            end
        else
            if N_row == N_col
                switch N_row
                    case 2
                        dd(1,:) = [0, d_El_even1, -d_El_even1];
                        dd(2,:) = [0, -d_El_even1, -d_El_even1];
                        dd(3,:) = [0, d_El_even1, -d_El_even1];
                        dd(4,:) = [0, d_El_even1, d_El_even1];
                        dd_T = (dd)';   % Function roty() act on column vectors
                        dd_T_rot = roty(Tilt_deg)*dd_T;
                        dd = (dd_T_rot)';
                    case 3
                        dd(1,:) = [0, -d_El_odd1, d_El_odd1];
                        dd(2,:) = [0, -d_El_odd1, 0];
                        dd(3,:) = [0, -d_El_odd1, -d_El_odd1];
                        dd(4,:) = [0, 0, d_El_odd1];
                        dd(5,:) = [0, 0, 0];
                        dd(6,:) = [0, 0, -d_El_odd1];
                        dd(7,:) = [0, d_El_odd1, d_El_odd1];
                        dd(8,:) = [0, d_El_odd1, 0];
                        dd(9,:) = [0, d_El_odd1, -d_El_odd1];
                        dd_T = (dd)';   % Function roty() act on column vectors
                        dd_T_rot = roty(Tilt_deg)*dd_T;
                        dd = (dd_T_rot)';
                    otherwise
                        dd = 0;
                end
            elseif (N_row ~= N_col)
                dd(1,:) = [0, -R_Ant_Ein*3, 0];
                dd(2,:) = [0, -R_Ant_Ein, 0];
                dd(3,:) = [0, R_Ant_Ein, 0];
                dd(4,:) = [0, R_Ant_Ein*3, 0];
            else
                dd = 0;
            end
        end
        clear d_El_even1 d_El_even2 d_El_odd1
        clear Ant_Elr Ant_Elc R_Ant_Ein
        clear d_El_even1 d_El_even2 d_El_odd1
        clear dd_T dd_T_rot
        
        % Calculate the steering column-vectors bb_k (Unit_c is cc_k)
        bb = zeros(num_UEs(2),N_ant_ele);   % Now is non-transposed
        for kk = 1:num_UEs(2)
            for nn = 1:N_ant_ele
                bb(kk,nn) = exp(2*pi*i*dot(dd(nn,:),Unit_c(kk,:)) / lambda);
            end
        end
        bb = (bb)';
        clear nn kk
        WW_CB = alfa * bb;
        % Calculate the chosen beamforming direction cc_k_Max (Fomula (6))
        % cc_k_Max_Dot : maximum product value
        % cc_k_Max_i : the position of the maximum product value
        [cc_k_Max_Dot,cc_k_Max_i] = max(dot(h_kn,bb));
        % Calculate the unit vector in polar and azimuthal angles
        Theta_cc_k_Max = acosd(Unit_c(cc_k_Max_i,3));
        Phi_cc_k_Max = asind(Unit_c(cc_k_Max_i,2) / sind(Theta_cc_k_Max));
        cc_k_Max_Angle = [Theta_cc_k_Max,Phi_cc_k_Max];
        cc_k_Max(1,:) = Unit_c(cc_k_Max_i,:);
        clear cc_k_Max_Dot Theta_cc_k_Max  Phi_cc_k_Max
end

% Show the 2D & 3D pattern
% Show the XOZ-plane(phi=0) and YOZ-plane(phi=90)
% Find the coordinate and values (Theta & Gain)
% Note: the First column is THETA, the second column is PHI
Pattern_XOZ = zeros(360, 2);
Pattern_YOZ = zeros(360, 2);
Pattern_XOZ = zeros(360, 2);
Pattern_YOZ = zeros(360, 2);
tt = 1;
if tt >= 359
   for pp = 1:Row_Pattern_Eles
      if Pattern_syn(pp,2) == 0
            Pattern_XOZ(tt,1) = Pattern_syn(pp,1);
            Pattern_XOZ(tt,2) = Pattern_syn(pp,3);
            tt = tt + 1;
            break
      elseif Pattern_syn(pp,2) == 90
          Pattern_YOZ(tt,:) = [Pattern_syn(pp,1) Pattern_syn(pp,3)];
          tt = tt + 1;
          break
      else
          break
      end
   end
end


clear alfa HH_H
clear N_BasicInfo_S N_col N_header N_header_Patt N_row
clear Ant_H Ant_outPow Ant_Position Ant_V
clear Dist_BStoUE
clear Dr_ele
clear freq_M R_Ant_Ein

% % Calculate the transmit vector tt
% ss = zeros(num_UEs(2),1);
% ss(:,1) = deal(sqrt(1/num_UEs(2)));
% tt_MRF= WW_MRF * ss;
% tt_ZF = WW_ZF * ss;
% 
% % Average the antenna patterns
% % read radiation pattern information
% N_Pattern_Header = 7; % number of header lines for the radiation pattern file
% fid = fopen('far field result_2500.000_MHz_Total(1).apa','rt');
% % read the columns of radiation pattern files
% lines_Pattern = 0;
% while ~feof(fid)
%     fgetl(fid);
%     lines_Pattern = lines_Pattern +1;
% end
% num_Points = lines_Pattern - N_Pattern_Header;   % the number of points
% fclose(fid);
% fid = fopen('far field result_2500.000_MHz_Total(1).apa','rt');
% FormatString = '%s %s %s';
% Pattern_S = textscan(fid,FormatString,'HeaderLines',N_Pattern_Header);
% Pattern_S = [Pattern_S{1,1} Pattern_S{1,2} Pattern_S{1,3}];   % Pattern_S is the pattern for each element
% fclose(fid);
% %%%%%% Copy this part when the num_Pattern changed (except changing the file name & function name)
% fid = fopen('far field result_2500.000_MHz_Total(2).apa','rt');
% FormatString = '%s %s %s';
% Pattern_Add1(1,:) = textscan(fid,FormatString,'HeaderLines',N_Pattern_Header);
% Pattern_Add1 = Pattern_Add1{1,3};
% fclose(fid);
% %%%%%%% Copy this part
% %%%%%%% Copy this part when the num_Pattern changed (except changing the file name & function name)
% fid = fopen('far field result_2500.000_MHz_Total(3).apa','rt');
% FormatString = '%s %s %s';
% Pattern_Add2(1,:) = textscan(fid,FormatString,'HeaderLines',N_Pattern_Header);
% Pattern_Add2 = Pattern_Add2{1,3};
% fclose(fid);
% %%%%%%% Copy this part
% %%%%%%% Copy this part when the num_Pattern changed (except changing the file name & function name)
% fid = fopen('far field result_2500.000_MHz_Total(4).apa','rt');
% FormatString = '%s %s %s';
% Pattern_Add3(1,:) = textscan(fid,FormatString,'HeaderLines',N_Pattern_Header);
% Pattern_Add3 = Pattern_Add3{1,3};
% fclose(fid);
% %%%%%%% Copy this part
% Pattern_S = [Pattern_S Pattern_Add1 Pattern_Add2 Pattern_Add3];
% Pattern_S = cellfun(@str2num,Pattern_S);   % Pattern_S is double
% Pattern_C = Pattern_S(:,3:6);   % Pattern_C is the pattern for calculation (only the ampitude, without the theta & phi)
% % Calculate the sum of the patterns (Fomula (10))
% Pattern_Final = zeros(num_Points,1);
% A = zeros(1,N_ant_ele);
% for p = 1 : num_Points
%     for n = 1 : N_ant_ele
%         A(n) = Pattern_C(p,n)*tt_MRF(n)*exp(-2*pi*i*dot(dd(n,:),cc_k_Max) / lambda);
%     end
%     Pattern_Final(p) = sum(A);
% end
% % why the complex number?
% 
% clear N_Pattern_Header lines_Pattern_read FormatString lines_Pattern
% clear Pattern_Add1 Pattern_Add2 Pattern_Add3
% clear A p n
% 
% 
% clear c0 Angle_Ant lamda N_header Num_Ant_El num_UEs
toc
fprintf('%f\n', toc)