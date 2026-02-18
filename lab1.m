% Lawrence Cai, Prof. Loh, ENGR 006, Spring 2026, Lab 1.
% Collaborators: Kinley Zangmo, Yinuo Gao.
function lab1()
    clc;
    % All the constant values: 
    % the radius of pulley; the weight of the center mass;
    % the positions of the 3 counterweights.
    r = 0.75;
    wt = 3.02;
    w2 = [-95.75, -68, 116.61417];
    w3 = [-95.25, 64, 116.69291];
    w4 = [94.5, -12.25, 116.73228];
    
    % The weights of the three counterweights in order of trial 1 and 2.
    weights_1 = [1.634, 1.616, 2.9025] + 0.055;
    weights_2 = [2.218, 1.616, 3.541] + 0.055;
    
    % The position of the center mass in all four trials: two before 
    % changing the masses of the counterweights, and two after.
    c1_1 = [0.25 -9 59.75];
    c1_2 = [4.5 -9.75 60.75];
    c2_1 = [16.5 -14.75 73.75];
    c2_2 = [7 -14.75 71.75];
    
    % The function that corrects the Z coordinate based on the original Z
    % and the original ground distance.
    % params: ground distance and the original Z.
    % return: the corrected Z.
    function new_Z = correctZ(gd, z)
        new_gd = gd - r; % Calculate the new distance by substracting r.
        alpha = atand(z/new_gd); % Calculate the angle of the smaller triangle.
        hyp = norm([new_gd, z]); % Calculate the hypotenuse based on the new distance and the original z.
        beta = asind(r/hyp); % Calculate the angle between the large and the small triangles.
        theta = alpha + beta; % Summing the angles gets the angle of the larger triangle.
        new_Z = gd*tand(theta); % Calculate the opposite of the larger triangle, which is the corrected Z coord.
    end
    
    % The function that calculate the unit vector that goes from c to w.
    % params: the uncorrected coordinates of the pulley, the coordinates of the center mass.
    % return: the unit vector lambda.
    function lam = findLam(w, c)
        x = w(1)-c(1); % Find the x and y component of the distance vector from the center mass to the weight.
        y = w(2)-c(2);
        gd = norm([x,y]); % Find the ground distance, which is the norm of x and y. 
        new_w = [w(1), w(2), correctZ(gd, w(3))] % Correct Z using the correctZ function to find the corrected Z.
        T = new_w-c; % Find the corrected tension vector.
        lam = T/norm(T); % Calculate the lambda by dividing the norm.
    end
    
    % The function that returns the magnitude of the three tension forces.
    % params: the position of the center mass.
    % return: the vector of the three tension forces in order.
    function tensions = findTensions(c)
        lam2 = findLam(w2, c); % Find the three corrected tension unit vectors.
        lam3 = findLam(w3, c);
        lam4 = findLam(w4, c);
        matrix = [lam2;lam3;lam4]'; % Assemble the matrix, and solve for the system using the known constant weight. 
        tensions = matrix\[0;0;wt];
        matrix = [lam2;lam3;lam4]% or we can use: tensions = inv(matrix)*[0;0;wt], as discussed on class.
    end
    
    % The function that returns the percent error as a vector.
    % params: the matrix of tensions and weights, where each row represent
    % the tensions and weights in each trial.
    function error = findPercentErrors(t, w)
        error = zeros(1, 4); % Since there are four trials, we set the size to be 4.
        for i = 1:4 % A for loop to iterate through each trial.
            
            % Each entry of the error is calculated by summing the
            % corresponding row of tensions and weights, find their
            % absolute value, and then divide by the ideal value which is
            % the sum of the tensions in that row/trial.
            sum_t = sum(t(i,:));
            sum_w = sum(w(i,:));
            error(i) = abs(sum_t - sum_w)/sum_w;
        end
    end
    
    % We assemble the tensions matrix and the weights matrix.
    tensions = [findTensions(c1_1)'; findTensions(c1_2)'; 
        findTensions(c2_1)'; findTensions(c2_2)'];
    weights = [weights_1; weights_1; weights_2; weights_2];

    % Find the percent_errors using the function and display it.
    percent_errors = findPercentErrors(tensions, weights);
    percent_error_1 = (percent_errors(1) + percent_errors(2))/2
    percent_error_2 = (percent_errors(3) + percent_errors(4))/2
    
end