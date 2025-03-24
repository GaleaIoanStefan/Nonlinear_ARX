clc; clear; close all;

% Loading the input data and storing it into easier-to-work-with variables
load("iddata-06.mat");
u_id = id.u;
y_id = id.y;
u_val = val.u;
y_val = val.y;

N_id = length(y_id);
N_val = length(y_val);

Ts = 0.01;
mmax = 3;
nmax = 5;

for m = 1:mmax
    for na = 1:nmax
        for nb = 1:nmax

            n = na+nb;
            nk = 1;

            % Reccursively calling the function recc to obtain all the possible
            % combinations of powers up until degree m
            combinations = recc(m,n);

            % Keeping only the valid combinations of powers, with the sum of the
            % indexes smaller than m
            k2 = 1;
            valid_comb = [];
            for k = 1:length(combinations)
                if sum(combinations(k,:)) <= m
                    valid_comb(k2,:) = combinations(k,:);
                    k2 = k2 + 1;
                end
            end

            % Prediction
            % Computing delayed versions of y and u and store them in vector d
            d_id = d_calc(u_id, y_id, na, nb);
            d_val = d_calc(u_val, y_val, na, nb);
            
            % Identification Data
            phi_id = phi_calc(N_id, n, d_id, valid_comb);
            theta = phi_id\y_id;
            ypred_id = phi_id * theta;

            % Validation Data
            phi_val = phi_calc(N_val, n, d_val, valid_comb);
            ypred_val = phi_val * theta;


            % Simulation
            % Identification Data
            ysim_id = zeros(N_id,1);
            ysim_id(1) = phi_id(1,:) * theta;
            for i = 2:N_id
                d_iter = d_sim(u_id, ysim_id, na, nb, i);
                phi_iter = phi_sim(n, d_iter, valid_comb);
                ysim_id(i) = phi_iter * theta;
            end

            % Validation Data
            ysim_val = zeros(N_val,1);
            ysim_val(1) = phi_val(1,:) * theta;
            for i = 2:N_val
                d_iter = d_sim(u_val, ysim_val, na, nb, i);
                phi_iter = phi_sim(n, d_iter, valid_comb);
                ysim_val(i) = phi_iter * theta;
            end
            
            % Computing the MSE for prediction and simulation
            mse_prediction_id(na,nb) = mean((y_id - ypred_id).^2);
            mse_simulation_id(na,nb) = mean((y_id - ysim_id).^2);
            mse_prediction_val(na,nb) = mean((y_val - ypred_val).^2);
            mse_simulation_val(na,nb) = mean((y_val - ysim_val).^2);

            % Forming matrices with all the values of the MSE, depending on
            % the parameters m, na, nb
            MSE_pred_id(5*(m-1)+na,nb) = mse_prediction_id(na,nb);
            MSE_sim_id(5*(m-1)+na,nb) = mse_simulation_id(na,nb);
            MSE_pred_val(5*(m-1)+na,nb) = mse_prediction_val(na,nb);
            MSE_sim_val(5*(m-1)+na,nb) = mse_simulation_val(na,nb);

            % Storing all the approximated output matrices into larger matrices so
            % that only the best approximator to be plotted afterwards
            start_index_id = (m-1)*na*nb*N_id + (na-1)*nb*N_id + (nb-1)*N_id + 1;
            end_index_id = start_index_id + N_id - 1;
            start_index_val = (m-1)*na*nb*N_val + (na-1)*nb*N_val + (nb-1)*N_val + 1;
            end_index_val = start_index_val + N_val - 1;
            
            id_pred(start_index_id:end_index_id, 1) = ypred_id;
            id_sim(start_index_id:end_index_id, 1) = ysim_id;

            val_pred(start_index_val:end_index_val, 1) = ypred_val;
            val_sim(start_index_val:end_index_val, 1) = ysim_val;
        end
    end
end


% Displaying Minimum MSE for each scenario in the command line
[min_pred_id, min_index_pred_id] = min(MSE_pred_id, [], "all");
[min_sim_id, min_index_sim_id] = min(MSE_sim_id, [], "all");
[min_pred_val, min_index_pred_val] = min(MSE_pred_val, [], "all");
[min_sim_val, min_index_sim_val] = min(MSE_sim_val, [], "all");


% Finding the minimum indexes
[m_pred_id, na_pred_id, nb_pred_id] = findParameters(min_index_pred_id, mmax, nmax);
[m_sim_id, na_sim_id, nb_sim_id] = findParameters(min_index_sim_id, mmax, nmax);
[m_pred_val, na_pred_val, nb_pred_val] = findParameters(min_index_pred_val, mmax, nmax);
[m_sim_val, na_sim_val, nb_sim_val] = findParameters(min_index_sim_val, mmax, nmax);


fprintf("\nMinimum MSE for Prediction on Identification Data\nMSE = %.6f \nm = %d, na = %d, nb = %d", min_pred_id, m_pred_id, na_pred_id, nb_pred_id);
fprintf("\n\nMinimum MSE for Simulation on Identification Data \nMSE = %.6f \nm = %d, na = %d, nb = %d", min_sim_id, m_sim_id, na_sim_id, nb_sim_id);

fprintf("\n\nMinimum MSE for Prediction on Validation Data \nMSE = %.6f \nm = %d, na = %d, nb = %d", min_pred_val, m_pred_val, na_pred_val, nb_pred_val);
fprintf("\n\nMinimum MSE for Simulation on Validation Data \nMSE = %.6f \nm = %d, na = %d, nb = %d\n", min_sim_val, m_sim_val, na_sim_val, nb_sim_val);


% 3D Plots
figure("Name","Prediction for Identification Set")
plots3D(MSE_pred_id, min_pred_id, m_pred_id, na_pred_id, nb_pred_id, nmax);
figure("Name","Simulation for Identification Set")
plots3D(MSE_sim_id, min_sim_id, m_sim_id, na_sim_id, nb_sim_id, nmax);
figure("Name","Prediction for Validation Set")
plots3D(MSE_pred_val, min_pred_val, m_pred_val, na_pred_val, nb_pred_val, nmax);
figure("Name","Simulation for Validation Set")
plots3D(MSE_sim_val, min_sim_val, m_sim_val, na_sim_val, nb_sim_val, nmax);

% Finding best predictions and simulations
ypred_id_best = extractBestValues(id_pred, m_pred_id, na_pred_id, nb_pred_id, N_id);
ysim_id_best = extractBestValues(id_sim, m_sim_id, na_sim_id, nb_sim_id, N_id);
ypred_val_best = extractBestValues(val_pred, m_pred_val, na_pred_val, nb_pred_val, N_val);
ysim_val_best = extractBestValues(val_sim, m_sim_val, na_sim_val, nb_sim_val, N_val);


% Comparing the best approximations
% Identification Data - Prediction and Simulation
figure
identification_data = iddata(y_id, u_id, Ts);
id_prediction_arx = iddata(ypred_id_best, u_id, Ts);
id_simulation_arx = iddata(ysim_id_best, u_id, Ts);
id_arx_model = arx(identification_data, [na_sim_id nb_sim_id nk]);
compare(identification_data, id_prediction_arx, id_simulation_arx, id_arx_model);
legend("Identification data", "id_{prediction}_{arx}", "id_{simulation}_{arx}", "id_{arx}_{model}")
title("Comparison for the Identification Data")

% Validation Data - Prediction and Simulation
figure
validation_data = iddata(y_val(n-1:end), u_val(n-1:end), Ts);
val_prediction_arx = iddata(ypred_val_best(n-1:end), u_val(n-1:end), Ts);
val_simulation_arx = iddata(ysim_val_best(n-1:end), u_val(n-1:end), Ts);
val_arx_model = arx(identification_data, [na_pred_val nb_pred_val nk]);
compare(validation_data, val_prediction_arx, val_simulation_arx, val_arx_model);
title("Comparison for the Validation Data")


% Functions

% Computing the phi matrix, formed by all the delays at all the possible
% combinations of powers, having the sum smaller than m
function phi = phi_calc(N, n, d, comb)
    for k = 1 : N
        for le = 1 : length(comb)
            prod = 1;
            for j = 1 : n
                prod = prod * d(k, j) ^ comb(le, j);
            end
            phi(k, le) = prod;
        end
    end
end


% Similar to the phi matrix, used for the calculus of the phi matrix for
% simulation nonlinear arx. The key difference is that only one line of phi
% is computed via a function call 
function phi = phi_sim(n, d, comb)
    for l = 1 : length(comb)
        prod = 1;
        for j = 1 : n
            prod = prod * d(j) ^ comb(l, j);
        end
        phi(l) = prod;
    end
end


% Computing the d matrix, formed by the delays of the output and input for
% all the given datapoints in the identification or validation set
function d = d_calc(u, y, na, nb)
    N = length(y);
    for i = 1:N
        for j = 1:na
            if i > j
                dy(i,j) = y(i-j);
            else
                dy(i,j) = 0;
            end
        end
        for j = 1:nb
            if i > j
                du(i,j) = u(i-j);
            else
                du(i,j) = 0;
            end
        end
    end
    d = [dy du];
end

% Similar to the d matrix, used for the calculus of the simulation delay for
% nonlinear arx. The key difference is that only one line of d is computed
% via a function call 
function d = d_sim(u, y, na, nb, i)
    for j = 1:na
        if i > j
            dy(j) = y(i-j);
        else
            dy(j) = 0;
        end
    end
    for j = 1:nb
        if i > j
            du(j) = u(i-j);
        else
            du(j) = 0;
        end
    end
    d = [dy du];
end

% Reccursive function that computes all the possible combinations of numbers
% up until a polynomial degree m, for a given number of terms n, which is
% equal to na + nb
function comb = recc(m, n)
    if n == 1
        comb = (0:m)';
    else
        comb2 = recc(m, n-1);
        comb = [];
        for k = 0:m
            new_comb = [k*ones(length(comb2), 1), comb2];
            comb = [comb; new_comb];
        end
    end
end

% Function to create 3D plots with the variation of the MSE with the
% parameters m, na and nb
function [] = plots3D(MSE, min, m, na, nb, nmax)
    [Y, X] = meshgrid(1:nmax, 1:nmax);
    for p = 1:3
        subplot(1,3,p)
        surf(X, Y, MSE((p-1)*nmax+1 : (p-1)*nmax+nmax,:));
        title("for m = " + p);
        xlabel("na"); ylabel("nb"); zlabel("MSE"); hold on
        if m == p
            plot3(na, nb, min, 'pr','MarkerSize', 15)
        end
    end
end

% Function to find the positions of the minimum values of the MSE from the
% 3D-array containing all the computed MSEs 
function [m, na, nb] = findParameters(min_index, mmax, nmax)
    m = floor(mod(min_index-1, mmax*nmax) / nmax) + 1;
    na = mod(min_index-1, nmax) + 1;
    nb = floor((min_index-1) / (mmax*nmax)) + 1;
end


% Function to extract the subset containing the best values from the data matrices for the
% specified parameters (m, na, nb, N)
function best_values = extractBestValues(ymatrix, m, na, nb, N)
    start_index = (m-1)*na*nb*N + (na-1)*nb*N + (nb-1)*N + 1;
    end_index = start_index + N - 1;
    best_values = ymatrix(start_index:end_index, 1);
end
