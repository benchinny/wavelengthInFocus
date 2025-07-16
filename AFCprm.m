%%210604 VRCMprm save varichrome parameters into struct

        if ey(1)=='R'
        AFCp.dspl_pwr=opto(name_map('r_disp')).control.getFocalPower.focal_power;
        AFCp.near_pwr=opto(name_map('r_t_near')).control.getFocalPower.focal_power;
        AFCp.far_pwr=opto(name_map('r_t_far')).control.getFocalPower.focal_power;
        AFCp.trmb_pos=zaber(name_map('r_trombone')).control.getposition;
        elseif ey(1)=='L'
        AFCp.dspl_pwr=opto(name_map('l_disp')).control.getFocalPower.focal_power;
        AFCp.near_pwr=opto(name_map('l_t_near')).control.getFocalPower.focal_power;
        AFCp.far_pwr=opto(name_map('l_t_far')).control.getFocalPower.focal_power;
        AFCp.trmb_pos=zaber(name_map('l_trombone')).control.getposition;
        else
        AFCp.dspl_pwr=[opto(name_map('l_disp')).control.getFocalPower.focal_power   opto(name_map('r_disp')).control.getFocalPower.focal_power];
        AFCp.near_pwr=[opto(name_map('l_t_near')).control.getFocalPower.focal_power opto(name_map('r_t_near')).control.getFocalPower.focal_power];
        AFCp.far_pwr=[opto(name_map('l_t_far')).control.getFocalPower.focal_power opto(name_map('r_t_far')).control.getFocalPower.focal_power];
        AFCp.trmb_pos=[zaber(name_map('l_trombone')).control.getposition zaber(name_map('l_trombone')).control.getposition];
        end
        