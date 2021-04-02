function [info,point] = get_winprop_ray(file_name,N_BasicInfo_S,N_header)
c0 = 3e8;
% read basic information
fid = fopen(file_name,'r');
C = textscan(fid, '%s %s %s %[^\n]',7,'Headerlines',N_BasicInfo_S);
B = {C{1, 2}{1, 1}, C{1, 3}{1, 1}, C{1, 4}{1, 1}, C{1, 2}{3, 1}, C{1, 2}{4, 1}, C{1, 2}{5, 1}, C{1, 2}{7, 1}};
b = cellfun(@str2num,B);       % b is double
info.Ant_Position = [b(1) b(2) b(3)];
info.Ant_H = b(4);
info.Ant_V = b(5);
info.Ant_outPow = b(6);
info.freq_M = b(7);
info.freq = info.freq_M*1e6;
info.lambda = c0 / info.freq;
fclose(fid);
clear C B b N_BasicInfo_S

% read ray information
fid = fopen(file_name,'r');
for nn = 1:N_header+1
    tline = fgetl(fid);
end

n_point = 0;

while tline~=-1   % tline~=-1?
    [key,val]=strtok(tline); 
    val=strrep(val,'D','0');   % 0 = D = diffraction
    val=strrep(val,'R','1');   % 1 = R = reflexion
    val=strrep(val,'T','2');   % 2 = T = transmission []
    val=val(2:end);  
    a = str2num(val);   
    sss = size(a);   
    if sss(2) <= 9
        switch key    
            case 'POINT'   
                n_point = n_point + 1;
                point(n_point).position = [a(1) a(2) a(3)]; 
                n_path = 0;
            case 'PATH' 
                n_path  = n_path+1;
                point(n_point).path(n_path).Delay = a(1)*1e-9;   % convert unit from ns to s
                point(n_point).path(n_path).FieldStrength = (10^(a(2)/20))*1e-6;   % convert unit from dBuV/m to V/m
                point(n_point).path(n_path).TypofPath = a(3); % 0 = deterministic, 5 = empirical
                point(n_point).path(n_path).NumofInteractions  = a(4); 
                point(n_point).path(n_path).Phase = a(5); % [rad]
                point(n_point).path(n_path).DoD = [a(6) a(7)]; % [azimuth elevation]
                point(n_point).path(n_path).DoA = [a(8) a(9)]; % [azimuth elevation]
                % Calculate the channel coefficient h_kn_on1path for point(n_point).path(n_path)
                point(n_point).path(n_path).h_kn_on1path = (10^(a(2)/20))*1e-6*exp(-2*pi*i*info.freq*a(1)*1e-9);
            otherwise
                error('key invalid')
        end
    elseif (sss(2) > 9) & (sss(2) <= 15)
        switch key    
            case 'POINT'   
                n_point = n_point + 1;
                point(n_point).position = [a(1) a(2) a(3)]; 
                n_path = 0;
            case 'PATH' 
                n_path  = n_path+1;
                point(n_point).path(n_path).Delay = a(1)*1e-9;   % convert unit from ns to s
                point(n_point).path(n_path).FieldStrength = (10^(a(2)/20))*1e-6;   % convert unit from dBuV/m to V/m
                point(n_point).path(n_path).TypofPath = a(3); % 0 = deterministic, 5 = empirical
                point(n_point).path(n_path).NumofInteractions  = a(4); 
                point(n_point).path(n_path).Phase = a(5); % [rad]
                point(n_point).path(n_path).DoD = [a(6) a(7)]; % [azimuth elevation]
                point(n_point).path(n_path).DoA = [a(8) a(9)]; % [azimuth elevation]
                point(n_point).path(n_path).InteractionPoint_0 = [a(10) a(11) a(12)]; % Coordinate of interaction point_0 (x,y,z)
                point(n_point).path(n_path).InteractionType_0 = a(13); % 0-diffraction,1-reflexion,2-transmission []
                point(n_point).path(n_path).IDofObject_0  = a(14);
                point(n_point).path(n_path).IDofMeterial_0 = a(15);
                % Calculate the channel coefficient h_kn_on1path for point(n_point).path(n_path)
                point(n_point).path(n_path).h_kn_on1path = (10^(a(2)/20))*1e-6*exp(-2*pi*i*info.freq*a(1)*1e-9);
            otherwise
                error('key invalid')
        end
    else
        switch key    
            case 'PATH' 
                n_path  = n_path+1;
                point(n_point).path(n_path).Delay = a(1)*1e-9;   % convert unit from ns to s
                point(n_point).path(n_path).FieldStrength = (10^(a(2)/20))*1e-6;   % convert unit from dBuV/m to V/m
                point(n_point).path(n_path).TypofPath = a(3); % 0 = deterministic, 5 = empirical
                point(n_point).path(n_path).NumofInteractions  = a(4); 
                point(n_point).path(n_path).Phase = a(5); % [rad]
                point(n_point).path(n_path).DoD = [a(6) a(7)]; % [azimuth elevation]
                point(n_point).path(n_path).DoA = [a(8) a(9)]; % [azimuth elevation]
                point(n_point).path(n_path).InteractionPoint_0 = [a(10) a(11) a(12)]; % Coordinate of interaction point_0 (x,y,z)
                point(n_point).path(n_path).InteractionType_0 = a(13); % 0-diffraction,1-reflexion,2-transmission []
                point(n_point).path(n_path).IDofObject_0  = a(14);
                point(n_point).path(n_path).IDofMeterial_0 = a(15);
                point(n_point).path(n_path).InteractionPoint_1 = [a(16) a(17) a(18)]; % Coordinate of interaction point_1 (x,y,z)
                point(n_point).path(n_path).InteractionType_1 = a(19); % 0-diffraction,1-reflexion,2-transmission []
                point(n_point).path(n_path).IDofObject_1  = a(20);
                point(n_point).path(n_path).IDofMeterial_1 = a(21);
                % Calculate the channel coefficient h_kn_on1path for point(n_point).path(n_path)
                point(n_point).path(n_path).h_kn_on1path = (10^(a(2)/20))*1e-6*exp(-2*pi*i*info.freq*a(1)*1e-9);
                
            otherwise
                error('key invalid')
        end
    end
    
    % read next line
    tline = fgetl(fid);
end
fclose(fid);
end