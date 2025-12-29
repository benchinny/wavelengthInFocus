% addpath(genpath(fullfile('toolboxes')));
% KbName('UnifyKeyNames');
function OPNT()
global name_map l_trombone_f r_trombone_f l_opto_f r_opto_f enable_optotunes enable_trombones zaber opto
addpath(genpath(fullfile('toolboxes')));
enable_optotunes = true;
enable_trombones = true;

%% Setup anonymous functions
% mag_to_d = @(x) 1/(x*0.1);


measured_optopower = [8, 9, 10, 11, 12];
measured_trombonepos_l = [1174521, 850677, 595023, 379309, 202719]; %updated 1/24 using shear plate and updated using sharp psf measurments after lens replacement
l_trombone_f = @(x) round(interp1(measured_optopower, measured_trombonepos_l, x, 'pchip')); % Note that this uses shape-preserving interpolation
l_opto_f = @(x) interp1(measured_trombonepos_l, measured_optopower, x, 'pchip'); % Note that this uses shape-preserving interpolation

% Setup interpolation function for right trombone
measured_trombonepos_r = [451223, 772527, 1018330, 1217662, 1381313];
r_trombone_f = @(x) round(interp1(measured_optopower, measured_trombonepos_r, x, 'pchip')); % Note that this uses shape-preserving interpolation
r_opto_f = @(x) interp1(measured_trombonepos_r, measured_optopower, x, 'pchip'); % Note that this uses shape-preserving interpolation

% Plot left trombone measurements
% x = linspace(8, 14, 1000);
% y = l_trombone_f(x);
% plot(x, y);
% hold on;
% plot(measured_optopower, measured_trombonepos_l, 'or');
% 
% % Plot right trombone measurements
% y = r_trombone_f(x);
% plot(x, y);
% plot(measured_optopower, measured_trombonepos_r, 'ob');
% legend('Left','Right');

name_map = containers.Map;

%% Open translation stage ports
if enable_trombones
    portName = 'COM9'; % Name of the serial port to use.
    baudRate = 115200; % Baud rate the Zaber device is configured to use.
    
    % Initialize port.
    port = serial(portName);
    set(port, ...
        'BaudRate', baudRate,   'DataBits', 8, ...
        'FlowControl', 'none',  'Parity', 'none', ...
        'StopBits', 1,          'Terminator','CR/LF');
    
    set(port, 'Timeout', 0.5)
    warning off MATLAB:serial:fgetl:unsuccessfulRead
    
    % Open the port.
    fopen(port);
    
    % In this example we know we're using the ASCII protocol, so just
    % instantiate it directly.
    protocol = Zaber.AsciiProtocol(port);
    
    % Rotation stage
    deviceAddress = 1; % Address the Zaber device is configured to use.
    ident = 'rotation';
    name_map(ident) = deviceAddress;
    zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
    zaber(name_map(ident)).unit_scale = 768000/360;
    zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
    zaber(name_map(ident)).move_deg = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
    fprintf('Device %d is a %s with firmware version %f\n', ...
        deviceAddress, zaber(deviceAddress).control.Name, ...
        zaber(deviceAddress).control.FirmwareVersion);
    
    % Right translation stage
    deviceAddress = 2; % Address the Zaber device is configured to use.
    ident = 'r_trombone';
    name_map(ident) = deviceAddress;
    zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
    zaber(name_map(ident)).unit_scale = 1584048/75.440286;
    zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
    zaber(name_map(ident)).move_mm = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
    fprintf('Device %d is a %s with firmware version %f\n', ...
        deviceAddress, zaber(deviceAddress).control.Name, ...
        zaber(deviceAddress).control.FirmwareVersion);
    
    % Left translation stage
    deviceAddress = 3; % Address the Zaber device is configured to use.
    ident = 'l_trombone';
    name_map(ident) = deviceAddress;
    zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
    zaber(name_map(ident)).unit_scale = 1584048/75.440286;
    zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
    zaber(name_map(ident)).move_mm = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
    fprintf('Device %d is a %s with firmware version %f\n', ...
        deviceAddress, zaber(deviceAddress).control.Name, ...
        zaber(deviceAddress).control.FirmwareVersion);
    
    zaber(name_map('rotation')).control.home();
    zaber(name_map('r_trombone')).control.home();
    zaber(name_map('l_trombone')).control.home();
    zaber(name_map('rotation')).control.waitforidle();
    zaber(name_map('r_trombone')).control.waitforidle();
    zaber(name_map('l_trombone')).control.waitforidle();
    
    % fclose(port);
    % delete(port);
end

%% Open optotune ports
try
    if enable_optotunes
        serial_ports = seriallist;
        opto = [];
        for p = [1 2 3 4 5 6]
            fprintf('Optotune serial port: %s\n', serial_ports(p))
            opto(p).port = serial_ports(p);
            opto(p).control = Optotune(opto(p).port);
            opto(p).control.Open();
            fprintf('Optotune serial number: %s\n', opto(p).control.serial_number)
            switch opto(p).control.serial_number
                case "AQAA2859"
                    fprintf('Optotune ID: Left display\n'); %
                    name_map('l_disp') = p;
                case "AQAA3868"
                    fprintf('Optotune ID: Left trombone, near\n'); %
                    name_map('l_t_near') = p;
                case "AQAA4262"
                    fprintf('Optotune ID: Left trombone, far\n'); %
                    name_map('l_t_far') = p;
                case "AQAA4190"
                    fprintf('Optotune ID: Right display\n'); %
                    name_map('r_disp') = p;
                case "AQAA3948"
                    fprintf('Optotune ID: Right trombone, near\n'); %
                    name_map('r_t_near') = p;
                case "AQAA3864"
                    fprintf('Optotune ID: Right trombone, far\n');
                    name_map('r_t_far') = p;
            end
            opto(p).control.modeFocalPower();
            fprintf('\n');
        end
    end
catch ERROR
    a = instrfind();
    if ~isempty(a) %isempty(Obj) returns logical 1 (true) if Obj is an empty ExptData object. Otherwise, it returns logical 0 (false). An empty ExptData object contains no data elements.
        fclose(a);
        delete(a)
        clear a
    end
    if enable_optotunes
        for p = 1:6
            opto(p).control.Close();
        end
    end
    if enable_trombones
        fclose(port);
        delete(port);
    end
    rethrow(ERROR)
end



% % addpath(genpath(fullfile('toolboxes')));
% % KbName('UnifyKeyNames');
% 
% enable_optotunes = true;
% enable_trombones = true;
% 
% %% Setup anonymous functions
% mag_to_d = @(x) 1/(x*0.1);
% measured_optopower = [8, 9, 10, 11, 12];
% measured_trombonepos_l =[1118902, 797630, 536757, 342365 ,176524];%updated 1/24 using shear plate but psf is not sharpest
% l_trombone_f = @(x) round(interp1(measured_optopower, measured_trombonepos_l, x, 'pchip')); % Note that this uses shape-preserving interpolation
% l_opto_f = @(x) interp1(measured_trombonepos_l, measured_optopower, x, 'pchip'); % Note that this uses shape-preserving interpolation
% 
% % Setup interpolation function for right trombone
% measured_trombonepos_r = [451223, 772527, 1018330, 1217662, 1381313];
% r_trombone_f = @(x) round(interp1(measured_optopower, measured_trombonepos_r, x, 'pchip')); % Note that this uses shape-preserving interpolation
% r_opto_f = @(x) interp1(measured_trombonepos_r, measured_optopower, x, 'pchip'); % Note that this uses shape-preserving interpolation
% 
% % Plot left trombone measurements
% x = linspace(8, 12, 100);
% y = l_trombone_f(x);
% % plot(x, y);
% % hold on
% % plot(measured_optopower, measured_trombonepos_l, 'or');
% 
% % Plot right trombone measurements
% y = r_trombone_f(x);
% % plot(x, y);
% % plot(measured_optopower, measured_trombonepos_r, 'or');
% 
% name_map = containers.Map;
% 
% %% Open translation stage ports
% if enable_trombones
%     portName = 'COM9'; % Name of the serial port to use.
%     baudRate = 115200; % Baud rate the Zaber device is configured to use.
%     
%     % Initialize port.
%     port = serial(portName);
%     set(port, ...
%         'BaudRate', baudRate,   'DataBits', 8, ...
%         'FlowControl', 'none',  'Parity', 'none', ...
%         'StopBits', 1,          'Terminator','CR/LF');
%     
%     set(port, 'Timeout', 0.5)
%     warning off MATLAB:serial:fgetl:unsuccessfulRead
%     
%     % Open the port.
%     fopen(port);
%     
%     % In this example we know we're using the ASCII protocol, so just
%     % instantiate it directly.
%     protocol = Zaber.AsciiProtocol(port);
%     
%     % Rotation stage
%     deviceAddress = 1; % Address the Zaber device is configured to use.
%     ident = 'rotation';
%     name_map(ident) = deviceAddress;
%     zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
%     zaber(name_map(ident)).unit_scale = 768000/360;
%     zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
%     zaber(name_map(ident)).move_deg = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
%     fprintf('Device %d is a %s with firmware version %f\n', ...
%         deviceAddress, zaber(deviceAddress).control.Name, ...
%         zaber(deviceAddress).control.FirmwareVersion);
%     
%     % Right translation stage
%     deviceAddress = 2; % Address the Zaber device is configured to use.
%     ident = 'r_trombone';
%     name_map(ident) = deviceAddress;
%     zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
%     zaber(name_map(ident)).unit_scale = 1584048/75.440286;
%     zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
%     zaber(name_map(ident)).move_mm = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
%     fprintf('Device %d is a %s with firmware version %f\n', ...
%         deviceAddress, zaber(deviceAddress).control.Name, ...
%         zaber(deviceAddress).control.FirmwareVersion);
%     
%     % Left translation stage
%     deviceAddress = 3; % Address the Zaber device is configured to use.
%     ident = 'l_trombone';
%     name_map(ident) = deviceAddress;
%     zaber(name_map(ident)).control = Zaber.AsciiDevice.initialize(protocol, deviceAddress);
%     zaber(name_map(ident)).unit_scale = 1584048/75.440286;
%     zaber(name_map(ident)).move = @(x) zaber(name_map(ident)).control.moveabsolute(x);
%     zaber(name_map(ident)).move_mm = @(x) zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale);
%     fprintf('Device %d is a %s with firmware version %f\n', ...
%         deviceAddress, zaber(deviceAddress).control.Name, ...
%         zaber(deviceAddress).control.FirmwareVersion);
%     
%     zaber(name_map('rotation')).control.home();
%     zaber(name_map('r_trombone')).control.home();
%     zaber(name_map('l_trombone')).control.home();
%     zaber(name_map('rotation')).control.waitforidle();
%     zaber(name_map('r_trombone')).control.waitforidle();
%     zaber(name_map('l_trombone')).control.waitforidle();
%     
%     % fclose(port);
%     % delete(port);
% end
% 
% %% Open optotune ports
% try
%     if enable_optotunes
%         serial_ports = seriallist;
%         opto = [];
%         for p = 1:6
%             fprintf('Optotune serial port: %s\n', serial_ports(p))
%             opto(p).port = serial_ports(p);
%             opto(p).control = Optotune(opto(p).port);
%             opto(p).control.Open();
%             fprintf('Optotune serial number: %s\n', opto(p).control.serial_number)
%             switch opto(p).control.serial_number
%                 case "AQAA2859"
%                     fprintf('Optotune ID: Left display\n');
%                     name_map('l_disp') = p;
%                 case "AQAA3868"
%                     fprintf('Optotune ID: Left trombone, near\n');
%                     name_map('l_t_near') = p;
%                 case "AQAA4262"
%                     fprintf('Optotune ID: Left trombone, far\n');
%                     name_map('l_t_far') = p;
%                 case "AQAA4190"
%                     fprintf('Optotune ID: Right display\n');
%                     name_map('r_disp') = p;
%                 case "AQAA3948"
%                     fprintf('Optotune ID: Right trombone, near\n');
%                     name_map('r_t_near') = p;
%                 case "AQAA3864"
%                     fprintf('Optotune ID: Right trombone, far\n');
%                     name_map('r_t_far') = p;
%             end
%             opto(p).control.modeFocalPower();
%             fprintf('\n');
%         end
%     end
% catch ERROR
%     a = instrfind();
%     if ~isempty(a) %isempty(Obj) returns logical 1 (true) if Obj is an empty ExptData object. Otherwise, it returns logical 0 (false). An empty ExptData object contains no data elements.
%         fclose(a);
%         delete(a)
%         clear a
%     end
%     if enable_optotunes
%         for p = 1:6
%             opto(p).control.Close();
%         end
%     end
%     if enable_trombones
%         fclose(port);
%         delete(port);
%     end
%     rethrow(ERROR)
% end
% 
