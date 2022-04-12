%% 3D Ray Tracing to calculate the attenuation and activity weighting

function [ACV, w] = RayTracing3D_ACV_w(phi,theta,xs,ys,zs,ActivityMap,AttenuationMap,attenuation_table,energy_index,R_FOV)
    
    % constants and start values
    [xdim,ydim,zdim] = size(ActivityMap);
    LineIntegralAtt = 0.0;
    LineIntegralWeight = 0.0;
    
    % assuming isotropic voxels for calculation of path length (in cm)
    voxel_dim = R_FOV / (.5 * xdim);

    % x,y,z of the start position
    x = xs;
    y = ys;
    z = zs;

    % direction of ray in x,y,z directions
    e_x = cosd(phi);
    e_y = sind(phi);
    e_z = sind(theta);

    % prevent dividing by zero
    if e_x == 0
        e_x = 1e-15;
    end
    if e_y == 0
        e_y = 1e-15;
    end
    if e_z == 0
        e_z = 1e-15;
    end

    % directional constants used in scan path (delta_i(x,y,z)) and
    % coordinates first voxel-boundaries for ray tracing (ix,iy,iz)
    if e_x >= 0
        delta_ix = 1;
        ix = xs + 1;
    elseif e_x < 0
        delta_ix = -1;
        ix = xs - 1;
    end
    if e_y >= 0
        delta_iy = 1;
        iy = ys + 1;
    elseif e_y < 0
        delta_iy = -1;
        iy = ys - 1;
    end
    if e_z >= 0
        delta_iz = 1;
        iz = zs + 1;
    elseif e_z < 0
        delta_iz = -1;
        iz = zs - 1;
    end

    % prevent divisions inside inner loops, so pre-calculate
    inv_e_x = 1/e_x;
    inv_e_y = 1/e_y;
    inv_e_z = 1/e_z;

    while round(x) >= 1 && round(x) <= xdim && round(y) >= 1 && round(y) <= ydim && round(z) >= 1 && round(z) <= zdim

        % calculate distance to each voxel boundary
        a_x = abs((ix-x)*inv_e_x);
        a_y = abs((iy-y)*inv_e_y);
        a_z = abs((iz-z)*inv_e_z);

        % set change in all directions to 0
        d_ix = 0; d_iy = 0; d_iz = 0;

        % check which voxel boundary is nearest and update a and d_x, d_y and d_z accordingly
        if a_x <= a_y && a_x <= a_z
            if a_x == a_y 
                if a_x == a_z
                    d_ix = delta_ix;
                    d_iy = delta_iy;
                    d_iz = delta_iz;

                    a = a_x;
                else
                    d_ix = delta_ix;
                    d_iy = delta_iy;

                    a = a_x;
                end
            else
                if a_x == a_z
                    d_ix = delta_ix;
                    d_iz = delta_iz;

                    a = a_x;
                else
                    d_ix = delta_ix;

                    a = a_x;
                end
            end
        elseif a_y <= a_x && a_y <= a_z
            if a_y == a_x
                d_ix = delta_ix;
                d_iy = delta_iy;

                a = a_y;
            else
                if a_y == a_z
                    d_iy = delta_iy;
                    d_iz = delta_iz;

                    a = a_y;
                else
                    d_iy = delta_iy;

                    a = a_y;
                end
            end
        elseif a_z <= a_x && a_z <= a_x
            if a_z == a_x
                d_ix = delta_ix;
                d_iz = delta_iz;

                a = a_z;
            else
                if a_z == a_y
                    d_iy = delta_iy;
                    d_iz = delta_iz;

                    a = a_z;
                else
                    d_iz = delta_iz;

                    a = a_z;
                end
            end
        end

        % calculate the length through current voxel
        l = (((1+(e_z^2))*(a^2))^.5)*voxel_dim;
        
        % get the attenuation coefficient
        if AttenuationMap(round(x),round(y),round(z)) > 0
            tissue = AttenuationMap(round(x),round(y),round(z));
            atn_coef = attenuation_table(energy_index,tissue);
        else
            atn_coef = 0;
        end
        
        % update Lineintegrals
        LineIntegralAtt = LineIntegralAtt + l*atn_coef;
        LineIntegralWeight = LineIntegralWeight + l*ActivityMap(round(x),round(y),round(z));

        % update x, y, and z
        x = x + a*e_x;
        y = y + a*e_y;
        z = z + a*e_z;
        
        % update ix, iy, and iz
        ix = ix + d_ix;
        iy = iy + d_iy;
        iz = iz + d_iz;

    end
    
    %% Outputs
    
    ACV = exp(-LineIntegralAtt);
    w = LineIntegralWeight;
    
end