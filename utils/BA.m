function F = BA(X)

load BA.mat

% Fix the first camera pose and Z position of second camera to get the unique solution 
X1_p = state.camState{1,1}.p_C_G;
X1_a = state.camState{1,1}.eul'; 
X2_Z = groundT.camState{1,2}.p_C_G(3,:);

X2_XY = X(1:2);
X2_p = [X(1:2); X2_Z];
X2_a = X(3:5); 

object = 0; 
r_o = []; 

for i = 1 : length(state.camState)
    for j = 1 : length(feat_ob)
        
        k = length(state.camState); 
        p_f_G = X(6+6*(k-2)+3*(j-1) :8+6*(k-2)+3*(j-1));
        
        if i == 1
            
            if isnan(measurement{1,feat_ob(:,j)}(:,i)) 
                break; 
            end 
            
            C_CG = eul2rotm(X1_a')';
            p_f_C = C_CG*(p_f_G - X1_p);
            
        elseif i == 2
            
            if isnan(measurement{1,feat_ob(:,j)}(:,i))
                break;
            end
            
            C_CG = eul2rotm(X2_a')'; 
            p_f_C = C_CG*(p_f_G - X2_p); 
            
        elseif i > 2
            
            if isnan(measurement{1,feat_ob(:,j)}(:,i))
                break;
            end
            C_CG = eul2rotm(X(9+6*(k-3):11+6*(k-3))')'; 
            p_f_C = C_CG*(p_f_G - X(6+6*(k-3):8+6*(k-3))');  
            
        end 
        
        
        z_hat = [p_f_C(1)/p_f_C(3); p_f_C(2)/p_f_C(3)];
       
        r = measurement{1,feat_ob(:,j)}(:,i) - z_hat;
        
        r_o = [r_o, r]; 
        
        object = object + (r'*r);
    end
end


% for i = 1 : length(state.camState)
%     for j = 1 : length(state.featState)
%         
%         k = length(state.camState); 
%         p_f_G = X(6+6*(k-2)+3*(j-1) :8+6*(k-2)+3*(j-1));
%         
%         if i == 1
%             C_CG = eul2rotm(X1_a')';
%             p_f_C = C_CG*(p_f_G - X1_p);
%         elseif i == 2
%             C_CG = eul2rotm(X2_a')'; 
%             p_f_C = C_CG*(p_f_G - X2_p);          
%         elseif i > 2
%             C_CG = eul2rotm(X(9+6*(k-3):11+6*(k-3))')'; 
%             p_f_C = C_CG*(p_f_G - X(6+6*(k-3):8+6*(k-3))');          
%         end 
%         
%         
%         z_hat = [p_f_C(1)/p_f_C(3); p_f_C(2)/p_f_C(3)];
%         
%         
%         r = measurement{1,j}(:,i) - z_hat;
%         
%         r_o = [r_o, r]; 
%         
%         object = object + (r'*r);
%     end
% end

r_o;
F = object; 








