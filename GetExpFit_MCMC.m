function [output] = GetExpFit_MCMC(params0,step,x,y,s,n_steps)
%%%MCMC This function runs a Markov Chain Monte Carlo for growth data
%%%(absorbance over time) that has been pre-smoothed by a rolling median.
%%%   INPUT ARGUMENTS:
%%%   1. params0 = array; starting location in parameter space
%%%      [a,k,leftbound,rightbound]
%%%   2. step = array; step size for each param, same length as param0
%%%   3. x = array; x-values
%%%   4. y = array; y-values
%%%   5. s = array; error on y-values
%%%   6. n_steps = val; number of steps to take
%%%
%%%
%%%   OUTPUT ARGUMENTS:
%%%   output: a struct containing the following fields
%%%   1. outputput.params: the most optimized model parameters
%%%   2. output.loc = array; location in parameter space
%%%   3. output.chi2 = val; chi2 value of model
%%%   4. output.acceptVal = array; accept (1) or reject (0) proposed step

tolerence = 1e-4; % stop the MCMC loop if chi2 begins to converge

% pre-allocate arrays for storage
output.acceptVal = NaN(n_steps, 1);
output.loc = NaN(n_steps, length(params0));
output.chi2 = NaN(n_steps, 1);

% Define the chi2 function using matrix algebra
calculateChi2 = @(model_out, data, error) dot((model_out - data) , (model_out - data))/((sum(error))^2);



for  n =1:n_steps+1
    
    % 1. Evaluate the y-values of the model, given the current guess of the
    % parameter values, within the range of the bounds
    model = evalModel(params0, x, y, s);
    
    % 2. Get the chi2 value of the model y-values given current params
    chi2_old = calculateChi2(model.model, model.y_sub , model.s_sub);
    
    % 3. Randomly pick a param to perturb & return new param
    new_param = updateParams(params0, step);
    
    % 4. Caculate chi2 for new parameter + old (unmodified other) parameter,
    % and then accept or reject
    new_model = evalModel(new_param, x, y, s); % evaluate model with new parameters
    
    chi2_new = calculateChi2(new_model.model, new_model.y_sub , new_model.s_sub); % calculate new chi2
    
    % 5. Accept / reject
    [params0, output.acceptVal(n), output.chi2(n)]  = evalStep(chi2_new, chi2_old, new_param, params0);
    output.loc(n,:) = params0; % store the parameter at this step
    
    % 6.  stop this loop if there is no improvement in accuracy
    if (abs(chi2_old - chi2_new) < tolerence) && (chi2_old ~= chi2_new)
        fprintf('stopping loop at %i iterations \n', n)
        break
    end
    
    
end

% cleaning up outputs

output.loc = output.loc(1:n, :);
output.chi2 = output.chi2(1:n);
output.acceptVal = output.acceptVal(1:n);
output.params = params0;
end

function out = evalModel(params, x, y, s)
%%% evaluate the growth curve model
%%% Inputs:
%%% params: (5x1) array with model parameters to be optomized
%%% x: (1xn) array containing the domain of the model
%%% y: (1xn) array containing the corresponding data to x
%%% s: (1xn) array containing error on y-values 
%%%
%%% Outputs:
%%% out: a struct containing
%%% a. out.x_sub: an (1xn) array that is subset of domain
%%% b. out.y_sub: an (1xn) array containing the corresponding data for x_sub
%%% c. out.s_sub: an (1xn) array containing the corresponding data for
%%% s_sub
%%% d. out.model: an (1xn) array containing points evaluated at x_sub

% defining parameters from state vector
a = params(1);
loss_frequency = params(2);
left_bound = params(3);
right_bound = params(4);
offset = params(5);

left_index = round(left_bound*length(x)); %make percentage into integer number
right_index = round(right_bound*length(x)); %make percentage into integer number
x_sub = x(max(left_index, 1): min(right_index, length(x))); %subset of x-array using indices
y_sub = y(max(left_index, 1): min(right_index, length(y))); %subset of y-array using indices
s_sub = s(max(left_index, 1): min(right_index, length(s))); %subset of y-array using indices
model = a*exp(loss_frequency*x_sub)+offset; %Calculate y-values for sub-x-array

% put results into a struct
out.x_sub = x_sub;
out.y_sub = y_sub;
out.s_sub = s_sub;
out.model = model;

end
% end of evalModel function

function new_params = updateParams(old_params, step)
%%% makes new parameter vector by perturbing one parameter with gaussian noise
%%% inputs:
%%% old_params: (1x5) array with model parameters
%%% steps: (1x5) array containing the standard deviation of the perturbation size for each parameter
%%% outputs:
%%% new_parameters: (1x5) array containing new model parameters

% random integer w/in length of param array
r = randi([1 length(old_params)],1);
% calculate value of new parameter
new_val = old_params(r)+normrnd(0,step(r)); %orginal param + std dev


% define new param
new_params = old_params;
% change perturbed param with new value of that param
new_params(r) = new_val;
end
% end of makeNewParams function

function [params_out, accept, chi2_out] = evalStep(new_chi2, old_chi2, new_params, old_params)
%%% evaluates whether to keep or reject the step in the MCMC
%%% inputs:
%%% 1. new_chi2: scalar with the chi2 evaluated at evalModel(new_params)
%%% 2. old_chi2: scalar with the chi2 evaluated at evalModel(old_params)
%%% 3. new_params: (1x5) array containing new model parameters
%%% 4. (1x5) array containing old model parameters

% Case 1: accept the step
if new_chi2 < old_chi2 %&& ij_new < ij_old && kl_new < kl_old
    accept = 1; %accept
    %            loc = new_param;
    chi2_out = new_chi2;
    params_out = new_params; % update location
    % ij_old = ij_new
    % kl_old = kl_new
    
else %Case 2: need to calc. ratio of probabilities
    eval = exp((old_chi2 - new_chi2)/2);
    
    % Case 2A: then compare those numbers
    if eval > rand(1)
        accept = 1; %accept
        %                loc = new_param;
        chi2_out = new_chi2;
        params_out = new_params; %update location
        
        % Case 2B: reject step
    else
        accept = 0; %reject
        %                loc = params0; %stay in old place
        chi2_out = old_chi2;
        params_out = old_params;
    end
end
end
% end of function evalStep


