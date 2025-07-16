
    function send_tcp0(scene, stage)
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
            cmsg(['TCP sending ' num2str(data(1)) ':' num2str(data(2)) ...
                  ':' num2str(data(3)) ':' num2str(data(4)) ...
                  ':' num2str(data(5)) ':' num2str(data(6))], ...
                  log.DEBUG, log.LEVEL);
            fwrite(scene.tcp_socket, data);
        end
    end
