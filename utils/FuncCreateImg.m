function FuncCreateImg(groundT, state, measurements,state_k,indices, input_value)


GT_P0 = [];
for i = 1 : length(groundT.camState)
    GT_P0(end+1, :) = [v_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];
end
image = [];
for i = 1:length(state_k)
    image = [];
    if i == 1
        image(1,:) = [v_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];
        c = 1;
    else
        if strcmp(input_value, 'Ground_truth') || strcmp(input_value, 'Initialization_4')
            image(1,:) = [v_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];
        end 
        if strcmp(input_value, 'Estimated') || strcmp(input_value, 'Initialization_5') 
            image(1,:) = [v_InvRotMatrixYPR22(quatToRotMat(state.camState{1,i}.q_C_G)), state.camState{1,i}.p_C_G'];
        end 
        c = 1;
    end
    for j = 1 : length(measurements{1,state_k(1,i)}.y)
        if ~isempty(indices{1,j})
            if ~isnan(measurements{1,state_k(1,i)}.y(1,j))
                image(c+1,:) = [j, measurements{1,state_k(1,i)}.y(:,j)', zeros(1,3)];
                c = c+1;
            end
        end
    end
    imageName = strcat('DataPrepareBA/Starry/Image',int2str(i),'.mat');
    save(imageName, 'image');
end

save('DataPrepareBA/Starry/GT_PO_PA.mat', 'GT_P0')