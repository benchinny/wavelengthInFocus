
    function send_tcp0fiat(tcp_socket,stage)
    global log
        %{
            6 byte serial output format:
            Byte 1: Trial
            Byte 2: Stage
            Byte 3: Condition
            Byte 4: Stimulus Diameter
            Byte 5: Stimulus Frequency
            Byte 6: Stimulus Diopters
        %}
        % tcp_socket = tcpserver("169.229.228.24",31000);
        scene.enable_tcp = 1;
        scene.trial_num = 121;
        if scene.enable_tcp
            if stage > 0
                trial_ = scene.trial_num;
                condition_ = 0; %scene.trials(trial_).condition;
                stimdiam_ = 0; %scene.trials(trial_).stimulus_diameter;
                stimfreq_ = 0; %find(stimulus_speeds == scene.trials(trial_).stimulus_speed);
                stimdiop_ = 0; %scene.stim_diopters_ind;
            else
                trial_ = 0; %0
                condition_ = 0; %0;
                stimdiam_ = 0;
                stimfreq_ = 0;
                stimdiop_ = 0;
            end
            assert((stage <= 255) && (stage >= 0))
            assert((condition_ <= 255) && (condition_ >= 0))
            assert((stimdiam_ <= 255) && (stimdiam_ >= 0))
            assert((stimfreq_ <= 255) && (stimfreq_ >= 0))
            assert((stimdiop_ <= 255) && (stimdiop_ >= 0))

            data = [trial_ stage condition_ stimdiam_ stimfreq_ stimdiop_];
            write(tcp_socket,data,"uint8")
        end
        % clear (tcp_socket);

    end