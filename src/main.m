function main()
orijin_image = imread('../resources/horse.jpg');
to_black_image = double(orijin_image);
% to_black_image = double(imread('../resources/horse_inpaint.jpg'));
global inpainting_image obj_similarity_vector obj_Neighbour_patch_vector
to_black_image(120:169, 255:264, 1:3) = 0;
inpainting_image = to_black_image;
PATCH_ORDER = 7;
N_ORDER = 35;
while 1
    [black_area_row_array, ~] = FIND_BLACK(inpainting_image, 1);    % row_index_lst, column_index_lst
    if numel(black_area_row_array) == 0
        break;
    end
    [edge_row, edge_column] = Extract_edge(PATCH_ORDER, 0);
    black_pixel_num = numel(black_area_row_array);
    fprintf('%d pixels are black same\n', black_pixel_num);
    max_priority = 0;
    for index = 1:numel(edge_row)     
        image_row =  edge_row(index); image_column =  edge_column(index);
        [Patch_priority, similarity_vector, Neighbour_patch_vector] = ...
        Patch_select(image_row, image_column, PATCH_ORDER, N_ORDER);
        if Patch_priority > max_priority             % find maximum priority
            max_priority = Patch_priority;
            obj_similarity_vector = similarity_vector;
            obj_Neighbour_patch_vector = Neighbour_patch_vector;
            priority_row = image_row;
            priority_column = image_column;
        end
    end
    fprintf('row = %d, column = %d \n', priority_row, priority_column);
    Patch_inpaint(priority_row, priority_column, PATCH_ORDER);
end
inpainting_image = uint8(inpainting_image);
figure;
subplot(131); imshow(uint8(orijin_image)), title('orijinal image');
subplot(132); imshow(uint8(to_black_image)), title('pre-restoration image'); 
subplot(133); imshow(inpainting_image), title('restored image'); 
end


% Function: To find black area pixel(index or row/column sub)
% Call format: [black_area_row, black_area_column] = FIND_BLACK(Patch, row_or_not)
% Input parameters: Patch, row_or_not, patch to calculate and output mode
% Output: sub_row, sub_column(row_or_not = 1), index(row_or_not = 0)
function[black_area_row, black_area_column] = FIND_BLACK(Patch, row_or_not)
if row_or_not == 1
    [black_area_row, black_area_column]  = find(Patch(:, :, 1)== 0 & Patch(:, :, 2)== 0 & Patch(:, :, 3) == 0);
else
    index  = find(Patch(:, :, 1)== 0 & Patch(:, :, 2)== 0 & Patch(:, :, 3) == 0);
    num = numel(Patch(:, :, 1));
    black_area_row = [index; index + num; index + 2 * num];
    black_area_column = 0;
end
end

% Function: Calculate distance between two patch(reshape :, 1)
% Call format: distance = Calculate_distance(Patch, Neighbor_patch) 
% Input para: Patch1, Patch2(Regardless of the order)
% Output para: distance
function[distance] = Calculate_distance(Patch, Neighbor_patch) 
distance = sum((Patch - Neighbor_patch).^2)/numel(Patch);
end

% Function: Extract edge of the black area
% Call format: [edge_row, edge_column] = Extract_edge(Black_row, Black_column)
% Input para: 
% Output para:
function [edge_row, edge_column] = Extract_edge(PATCH_ORDER, tiaohe_coefficient)
global inpainting_image
[inpainting_image_row, inpainting_image_column,~] = size(inpainting_image);
index_lst = find(inpainting_image(:, :, 1)== 0 & inpainting_image(:, :, 2)== 0 & inpainting_image(:, :, 3) == 0);
grey_image = ones(inpainting_image_row, inpainting_image_column); grey_image(index_lst) = 0;
grey_image = imerode(grey_image, strel('square', double(uint8(PATCH_ORDER / 2 - tiaohe_coefficient))));   % 腐蚀
BW = edge(grey_image, 'sobel');
[edge_row, edge_column] = find(BW ~= 0);
end

function[patch_priority, similarity_vector, Neighbour_patch_vector] = Patch_select(image_row, image_column, PATCH_ORDER, N_order)
global inpainting_image
patch_width = (PATCH_ORDER - 1) / 2;
patch = inpainting_image((image_row - patch_width): (image_row + patch_width), ...
(image_column - patch_width): (image_column + patch_width), 1:3);
[patch_black_area_index, ~] = FIND_BLACK(patch, 0);      % index
row = numel(patch_black_area_index);
if row ~= 0
    patch = reshape(patch, PATCH_ORDER * PATCH_ORDER * 3, 1);
    patch_confidence = 1 - row / PATCH_ORDER / PATCH_ORDER;
    similarity_vector = zeros(N_order * N_order, 1);
    NSP_num = 0;
    Neighbour_patch_vector = zeros(PATCH_ORDER * PATCH_ORDER * 3, N_order * N_order);      % 列向量是基底
    N_width = int32((N_order - 1) / 2);
    for new_image_row = image_row - N_width : image_row + N_width 
        for new_image_column = image_column - N_width : image_column + N_width 
            row_down = new_image_row - patch_width; row_up = new_image_row + patch_width;
            column_left = new_image_column - patch_width; column_right = new_image_column + patch_width;
            compare_patch = inpainting_image(row_down: row_up, column_left: column_right, 1:3);
             [cmp_black_area_row, ~] = FIND_BLACK(compare_patch, 1);        % row, column_lst
             if numel(cmp_black_area_row) == 0           % No black area, can be an atom
                 NSP_num  = NSP_num + 1;
                 Orijincompare_patch = reshape(compare_patch, PATCH_ORDER * PATCH_ORDER * 3, 1);
                 compare_patch = Orijincompare_patch;
                 compare_patch(patch_black_area_index) = 0;
                 patch_distance = Calculate_distance(patch, compare_patch);
                 Neighbour_patch_vector(:, NSP_num) =  Orijincompare_patch;
                 similarity =  exp(-patch_distance / 25);
                 similarity_vector(NSP_num) = similarity; 
             end
        end
    end
    similarity_vector(NSP_num + 1:end) = [];
    normalize_coefficient = sum(similarity_vector);
    similarity_vector = similarity_vector / normalize_coefficient; % Normalization
    Neighbour_patch_vector(:, NSP_num + 1:end) = [];
    structure_sparsity = norm(similarity_vector, 2) * sqrt(NSP_num) / N_order;
    patch_priority = structure_sparsity * patch_confidence;
else
    patch_priority = 0;
    similarity_vector = [];
    Neighbour_patch_vector = [];
end
end


% Function： to inpaint image based on the patch choose in Patch_select func
% Call format: Patch_inpaint(priority_row, priority_column, PATCH_ORDER)
% Input paratemter: Priority_row_sub, Priority_column_sub, Patch_order
% Output parameter: None, change global parameter inpainting_image
function Patch_inpaint(priority_row, priority_column, PATCH_ORDER)
global inpainting_image obj_similarity_vector obj_Neighbour_patch_vector
patch_width = (PATCH_ORDER - 1) / 2;
priority_patch = inpainting_image(priority_row - patch_width: priority_row + patch_width,...
    priority_column - patch_width: priority_column + patch_width, 1:3);
[priority_patch_black_area_index, ~] = FIND_BLACK(priority_patch, 0);   % index
priority_patch = reshape(priority_patch, PATCH_ORDER * PATCH_ORDER * 3, 1);
[~, index_lst] = sort(obj_similarity_vector, 'descend'); % similarity 越高越接近
N = 10;
if numel(index_lst) > N
    index_lst = index_lst(1:N);
else
    N = numel(index_lst);
end
fprintf('beforeN = %d \n', N);
% 非黑色区域
P_ = setdiff(1:PATCH_ORDER* PATCH_ORDER * 3, priority_patch_black_area_index);
P_plus = [];
for num = 1:size(obj_Neighbour_patch_vector, 2) 
    P_plus = [P_plus, P_ + (num-1) * (PATCH_ORDER * PATCH_ORDER * 3)];
end
Known_toblack_obj_Neibour_vector = obj_Neighbour_patch_vector;
Known_toblack_obj_Neibour_vector(P_plus) = 0 ; 
after_sum_obj_Neighbour_patch_vector = Known_toblack_obj_Neibour_vector * obj_similarity_vector;  % 求和相加
priority_patch(priority_patch_black_area_index) = 0;
Objective_function = [priority_patch; after_sum_obj_Neighbour_patch_vector * 0.5];      % pusaiT(Column vector)
Sm = [];
constraint_min_episilon = 100000000;
X = [];
while 1
    new_min_episilon = constraint_min_episilon;
    for i = 1:N
        sort_index = index_lst(i); % min distance
        similar_patch = obj_Neighbour_patch_vector(:, sort_index);
        Try_Sm = [Sm, similar_patch];
        P_row_pusaiP = similar_patch;  P_row_pusaiP(priority_patch_black_area_index) = 0;
        P_pusaiP = similar_patch;  P_pusaiP(P_) = 0; 
        D_pusaiP = [P_row_pusaiP; P_pusaiP];         % Approximate function(Column vector)
        TryX = [X, D_pusaiP] ; 
        repeat_obj_matrix = repmat(Objective_function, 1, size(TryX, 2));
        inv_Gram_matrix = pinv((repeat_obj_matrix - TryX)' * (repeat_obj_matrix - TryX));
        objective_vector = sum(inv_Gram_matrix, 2) / sum(sum(inv_Gram_matrix));
        Projection_orijin_matrix = Try_Sm * objective_vector;   % Column vector
        Projection_matrix = Projection_orijin_matrix; Projection_matrix(P_) = 0;
        consistency_episilon = 0.25 * Calculate_distance(Projection_matrix, after_sum_obj_Neighbour_patch_vector);
        Projection_matrix = Projection_orijin_matrix; Projection_matrix(priority_patch_black_area_index) = 0;
        in_episilon = Calculate_distance(Projection_matrix, priority_patch);
        episilon = max(in_episilon, consistency_episilon);
        if episilon < constraint_min_episilon
            constraint_min_episilon = episilon;
            min_index = sort_index;
            min_X = TryX;
        end
    end
    if N == 1 || new_min_episilon == constraint_min_episilon
        break
    elseif new_min_episilon ~= constraint_min_episilon
         Sm = [Sm, obj_Neighbour_patch_vector(:, min_index)];
         index_lst(index_lst == min_index) = [];
         X = min_X;
         N = N - 1;
         min_patch = Projection_orijin_matrix;
    end
end
min_patch(P_) = priority_patch(P_);
size(Projection_orijin_matrix)
% size(priority_patch)
size(min_patch)
result_Projection_matrix = reshape(min_patch, PATCH_ORDER, PATCH_ORDER, 3);
inpainting_image(priority_row - patch_width: priority_row + patch_width,...
    priority_column - patch_width: priority_column + patch_width, 1:3) = result_Projection_matrix;
disp('finish inpainting');
end
