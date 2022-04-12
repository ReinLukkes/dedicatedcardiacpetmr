%% 2D Ray Tracing to calculate the attenuation and activity weighting

function [ACV, w] = RayTracing2D(ThetaRay1, ActivityMap, AttenuationMap , i, j, R_FOV)               
    
    %% Global values throughout function
    sz = size(AttenuationMap);
    dim = sz(1);
    voxel_dim = R_FOV / (.5 * dim); % size voxel in cm (assume isotropic voxel)
    
    %% Initial values
    s_curr = 0.0;               % current value of s (current position on path)
    s_next = 0.0;               % next value of s (next position on path)
    LineIntegralAtt = 0.0;      % Line Integral value
    LineIntegralWeight = 0.0;   % Activity weight line integral

    %%_____________________________________________________________________
    %% Ray 1 - attenuation + activity weighting
    
    % compute unit directions of path
    e_x = cosd(ThetaRay1);
    e_y = sind(ThetaRay1);

    % prevent dividing by zero
    if (e_x == 0) 
        e_x = 1e-15; 
    end
    if (e_y == 0) 
        e_y = 1e-15; 
    end

    % directional constants used in scan path
    if (e_x>=0) 
        delta_ix = 1; 
        d_x = 1;
    elseif (e_x<0) 
        delta_ix = -1;
        d_x = -1;
    end
    
    if (e_y >= 0) 
        delta_iy = 1;
        d_y = 1;
    elseif (e_y < 0) 
        delta_iy = -1;
        d_y = -1;
    end

    % start voxel
    ix = i;
    iy = j;

    % distance from start voxel
    Dx = d_x - i;
    Dy = d_y - j;

    % prevent divisions inside inner loops, so pre-calculate
    inv_e_x = 1/e_x;
    inv_e_y = 1/e_y;

    % compute line integral;
    while ((ix > 1) && (iy > 1) && (ix < dim) && (iy < dim))
        
        % possible intersections of path with voxel boundary
        % (x, y boundary)
        % direction of the path is taken into account
        s_x = double(ix+Dx)*inv_e_x;      % s_x is total path to x-intersection
        s_y = double(iy+Dy)*inv_e_y;      % s_y is total path to y-intersection

        % only the closest intersection is really encountered
        % find this intersection and update voxel index for this direction
        % in some cases more boundaries are possible (45 degree paths)
        if(s_x <= s_y)              % intersection at x-boundary
            s_next = s_x;
            ix = ix + delta_ix;     % x-index next voxel
        end

        if(s_y <= s_x)              % intersection at y-boundary
            s_next = s_y;
            iy = iy + delta_iy;     % y-index next voxel
        end

        % length through the voxel
        l = (s_next-s_curr)*voxel_dim;

        % calculate line integral
        LineIntegralAtt = double(LineIntegralAtt) + AttenuationMap(ix,iy)*l;
        LineIntegralWeight = double(LineIntegralWeight) + ActivityMap(ix,iy)*l;

        % update voxelcount and current position
        s_curr = s_next;
    end

    %%_____________________________________________________________________
    %% Ray 2 - attenuation
    
%     if LineIntegralWeight ~= 0.0
%         s_curr = 0.0;           % current value of s (current position on path)
%         s_next = 0.0;           % next value of s (next position on path)
% 
%         % compute unit directions of path
%         e_x = cosd(ThetaRay2);
%         e_y = sind(ThetaRay2);
% 
%         % prevent dividing by zero
%         if (e_x == 0) 
%             e_x = 1e-15; 
%         end
%         if (e_y == 0) 
%             e_y = 1e-15; 
%         end
% 
%         % directional constants used in scan path
%         if (e_x>=0) 
%             delta_ix = 1; 
%             d_x = 1;
%         elseif (e_x<0) 
%             delta_ix = -1;
%             d_x = -1;
%         end
% 
%         if (e_y >= 0) 
%             delta_iy = 1;
%             d_y = 1;
%         elseif (e_y < 0) 
%             delta_iy = -1;
%             d_y = -1;
%         end
% 
%         % start voxel
%         ix = i;
%         iy = j;
% 
%         % distance from start voxel
%         Dx = d_x - i;
%         Dy = d_y - j;
% 
%         % prevent divisions inside inner loops, so pre-calculate
%         inv_e_x = 1/e_x;
%         inv_e_y = 1/e_y;
% 
%         % compute line integral;
%         while ((ix > 1) && (iy > 1) && (ix < dim) && (iy < dim))
% 
%             % possible intersections of path with voxel boundary
%             % (x, y boundary)
%             % direction of the path is taken into account
%             s_x = double(ix+Dx)*inv_e_x;      % s_x is total path to x-intersection
%             s_y = double(iy+Dy)*inv_e_y;      % s_y is total path to y-intersection
% 
%             % only the closest intersection is really encountered
%             % find this intersection and update voxel index for this direction
%             % in some cases more boundaries are possible (45 degree paths)
%             if(s_x <= s_y)              % intersection at x-boundary
%                 s_next = s_x;
%                 ix = ix + delta_ix;     % x-index next voxel
%             end
% 
%             if(s_y <= s_x)                % intersection at y-boundary
%                 s_next = s_y;
%                 iy = iy + delta_iy;         % y-index next voxel
%             end
% 
%             % length through the voxel
%             l = (s_next-s_curr)*voxel_dim;
% 
%             % calculate line integral
%             LineIntegralAtt = double(LineIntegralAtt) + AttenuationMap(ix,iy)*l;
% 
%             % update voxelcount and current position
%             s_curr = s_next;
%         end
%     end
    
    %%_____________________________________________________________________
    %% Outputs
    
    ACV = exp(-LineIntegralAtt);    % Calculate attenuation correction value
    w = LineIntegralWeight;
 
end