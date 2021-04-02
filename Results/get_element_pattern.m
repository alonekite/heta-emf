function [Pattern_Eles] = get_element_pattern(N_header_Patt)
c0 = 3e8;

% read the 3D pattern of 'Element_1' 
Pattern_name_1 = '3D_Pattern_Port1_3100.000_MHz_Total.apa';
fid = fopen(Pattern_name_1,'r');
for nn = 1:N_header_Patt+1
    tline = fgetl(fid);
end
n_azimuth = 0; % Given the original number of azimuth
while tline~=-1   % tline~=-1?
    [key,val]=strtok(tline);
    Pattern_InfoA = str2num(key);
    Pattern_InfoB = str2num(val);
    n_azimuth = n_azimuth+1;
    Pattern1st(n_azimuth).Azimuth = [Pattern_InfoA(1) Pattern_InfoB(1)];   % The 3D pattern for the FIRST element
    Pattern1st(n_azimuth).RealizedGain = Pattern_InfoB(2);        
    % read next line
    tline = fgetl(fid);
end
fclose(fid);
[Row_pattern_1, Col_pattern_1] = size(Pattern1st);
clear fid key n_azimuth nn 
clear Pattern_InfoA Pattern_InfoB tline val
clear Row_pattern_1 Pattern_name_1

% Pre storage an azimuth matrix
    Azimuth = zeros(Col_pattern_1,2);

% read the 3D pattern of 'Element_2' 
Pattern_name_2 = '3D_Pattern_Port2_3100.000_MHz_Total.apa';
fid = fopen(Pattern_name_2,'r');
for nn = 1:N_header_Patt+1
    tline = fgetl(fid);
end
n_azimuth = 0; % Given the original number of azimuth
while tline~=-1   % tline~=-1?
    [key,val]=strtok(tline);
    Pattern_InfoA = str2num(key);
    Pattern_InfoB = str2num(val);
    n_azimuth = n_azimuth+1;
    Pattern2nd(n_azimuth).Azimuth = [Pattern_InfoA(1) Pattern_InfoB(1)];   % The 3D pattern for the FIRST element
    Pattern2nd(n_azimuth).RealizedGain = Pattern_InfoB(2);
            Azimuth(n_azimuth,1) = [Pattern_InfoA(1)];
            Azimuth(n_azimuth,2) = [Pattern_InfoB(1)];
    % read next line
    tline = fgetl(fid);
end
fclose(fid);
[Row_pattern_2, Col_pattern_2] = size(Pattern2nd);
clear fid key n_azimuth nn 
clear Pattern_InfoA Pattern_InfoB tline val
clear Row_pattern_2 Pattern_name_2

% read the 3D pattern of 'Element_3' 
Pattern_name_3 = '3D_Pattern_Port3_3100.000_MHz_Total.apa';
fid = fopen(Pattern_name_3,'r');
for nn = 1:N_header_Patt+1
    tline = fgetl(fid);
end
n_azimuth = 0; % Given the original number of azimuth
while tline~=-1   % tline~=-1?
    [key,val]=strtok(tline);
    Pattern_InfoA = str2num(key);
    Pattern_InfoB = str2num(val);
    n_azimuth = n_azimuth+1;
    Pattern3rd(n_azimuth).Azimuth = [Pattern_InfoA(1) Pattern_InfoB(1)];   % The 3D pattern for the FIRST element
    Pattern3rd(n_azimuth).RealizedGain = Pattern_InfoB(2);
    % read next line
    tline = fgetl(fid);
end
fclose(fid);
[Row_pattern_3, Col_pattern_3] = size(Pattern3rd);
clear fid key n_azimuth nn 
clear Pattern_InfoA Pattern_InfoB tline val
clear Row_pattern_3 Pattern_name_3

% read the 3D pattern of 'Element_4' 
Pattern_name_4 = '3D_Pattern_Port4_3100.000_MHz_Total.apa';
fid = fopen(Pattern_name_4,'r');
for nn = 1:N_header_Patt+1
    tline = fgetl(fid);
end
n_azimuth = 0; % Given the original number of azimuth
while tline~=-1   % tline~=-1?
    [key,val]=strtok(tline);
    Pattern_InfoA = str2num(key);
    Pattern_InfoB = str2num(val);
    n_azimuth = n_azimuth+1;
    Pattern4th(n_azimuth).Azimuth = [Pattern_InfoA(1) Pattern_InfoB(1)];   % The 3D pattern for the FIRST element
    Pattern4th(n_azimuth).RealizedGain = Pattern_InfoB(2);
    % read next line
    tline = fgetl(fid);
end
fclose(fid);
[Row_pattern_4, Col_pattern_4] = size(Pattern4th);
clear fid key n_azimuth nn 
clear tline val
clear Row_pattern_4 Pattern_name_4

% Determine whether the dimensions of each element's pattern are equal
if Col_pattern_1 == Col_pattern_2
   if Col_pattern_2 == Col_pattern_3
      if Col_pattern_3 == Col_pattern_4
          disp('Patterns Legal');
      else 
         error('Patterns Illegal');  
      end
   else 
      error('Patterns Illegal');  
   end
else
   error('Patterns Illegal');  
end

Gain1st = [Pattern1st.RealizedGain];
Gain1st = Gain1st(:);
Gain2nd = [Pattern2nd.RealizedGain];
Gain2nd = Gain2nd(:);
Gain3rd = [Pattern3rd.RealizedGain];
Gain3rd = Gain3rd(:);
Gain4th = [Pattern4th.RealizedGain];
Gain4th = Gain4th(:);
Pattern_Eles = [Azimuth(:,1), Azimuth(:,2), Gain1st, Gain2nd, Gain3rd, Gain4th];

clear Azimuth 
clear Gain1st Gain2nd Gain3rd Gain4th
clear Pattern_InfoA Pattern_InfoB

end