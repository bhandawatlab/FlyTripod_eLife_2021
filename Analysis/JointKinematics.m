classdef JointKinematics
    %JointKinematics Determine joint kinematics from leg points.
    %   Detailed explanation goes here
    
    properties
        ThC_U
        CTr_U
        FTi_U
        TiTa_U
        Ta_U
    end
    
    methods
        function obj = JointKinematics(ap_axis_points)
            %JointKinematics Construct an instance of this class
            
            % Update fields
            obj.ThC_U = ap_axis_points(1).ThC_U;
            obj.CTr_U = ap_axis_points(1).CTr_U;
            obj.FTi_U = ap_axis_points(1).FTi_U;
            obj.TiTa_U = ap_axis_points(1).TiTa_U;
            obj.Ta_U = ap_axis_points(1).Ta_U;
        end
        
        function joint_angles = run(obj)
            % run Determine joint kinematics
            joint_angles = struct();  % Output structure
            
            % The points are relative to the AP-axis:
            % u1=PA
            % u2=LR
            % u3=VD
            
            %% Get the normalized femur vector
            % Subtract the attachment point
            femur_vector = obj.FTi_U - obj.CTr_U;
            
            % Normalize
            femur_lengths = vecnorm(femur_vector,2,2);
            femur_norm = femur_vector ./ femur_lengths;
            
            %% Get the normalized tibia vector
            % Subtract the attachment point
            tibia_vector = obj.TiTa_U - obj.FTi_U;
            
            % Normalize
            tibia_lengths = vecnorm(tibia_vector,2,2);
            tibia_norm = tibia_vector ./ tibia_lengths;
            
            %% Get the normalized tarsus vector
            % Subtract the attachment point
            tarsus_vector = obj.Ta_U - obj.TiTa_U;
            
            % Normalize
            tarsus_lengths = vecnorm(tarsus_vector,2,2);
            tip_norm = tarsus_vector ./ tarsus_lengths;

            % Coxa levation-depression is the angle between the femur vector
            % and the vertical axis (u3) [0,pi]. The angle between vectors is
            % θ=arccos(dot(x,y)/|x||y|)
            point_count = size(obj.CTr_U, 1);
            lev_dep = acos(dot(femur_norm, repmat([0 0 1], point_count, 1),2));
            joint_angles(1).lev_dep = lev_dep;

            % We also want the coxa's horizontal angle (azimuth) which
            % represents retraction-protraction (RP).
            % First, use the dot product to get the femur vector's projection
            % onto the X and Y axes.
            xfemur_proj = dot(femur_norm, repmat([1 0 0], point_count, 1),2);
            yfemur_proj = dot(femur_norm, repmat([0 1 0], point_count, 1),2);

            % Next, obtain the femur vector's components only in the X
            % (cos(RP)) and Y (sin(RP)) axis directions by dividing the
            % projections by the Z component (sin(LD)).
            zfemur_comp = sin(lev_dep);
            xfemur_no_z = xfemur_proj ./ zfemur_comp;
            yfemur_no_z = yfemur_proj ./ zfemur_comp;

            % Then use the tangent function to solve for RP [-pi, pi].        
            ret_pro = atan2(yfemur_no_z, xfemur_no_z);
            joint_angles(1).ret_pro = ret_pro;

            % For the remaining angles, we establish a reference frame at the
            % femur-tibia (FTi) joint.
            z_fti_norm = femur_norm;  % Already normalized
            y_fti = cross(z_fti_norm, repmat([0 0 1], point_count, 1),2);
            y_fti_norm = y_fti ./ vecnorm(y_fti,2,2);  % Normalize
            x_fti = cross(y_fti, z_fti_norm,2);
            x_fti_norm = -(x_fti ./ vecnorm(x_fti,2,2));  % Normalize

            % The extension-flexion angle is the angle between the femur and
            % tibia vectors.
            ext_flex = acos(dot(femur_norm, tibia_norm, 2));
            joint_angles(1).ext_flex = ext_flex;

            % The pronation-supination angle (PS) is the angle between the
            % tibia vector and the FTi X-axis.
            % First, use the dot product to get the femur vector's projection
            % onto the X and Y axes.
            xtibia_proj = dot(tibia_norm, x_fti_norm,2);
            ytibia_proj = dot(tibia_norm, y_fti_norm,2);

            % Next, obtain the femur vector's components only in the X
            % (cos(PS)) and Y (sin(PS)) axis directions by dividing the
            % projections by the Z component (sin(EF)).
            ztibia_comp = sin(ext_flex);
            xtibia_no_z = xtibia_proj ./ ztibia_comp;
            ytibia_no_z = ytibia_proj ./ ztibia_comp;

            % Then use the tangent function to solve for PS.        
            pron_sup = atan2(ytibia_no_z, xtibia_no_z);
            joint_angles(1).pron_sup = pron_sup;
        end
    end
end
