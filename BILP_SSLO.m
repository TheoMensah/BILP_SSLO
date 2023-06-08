function [sol,xFlag,OFV,Tm,EBV,Grade,stp,mm,mp] = BILP_SSLO(EBV,Grade,Cutoff)

%% Block Model and Parameter File Input:
% filein = -1;
% while filein < 0        % repeat until valid file name entered
%     [inputfile, inputpath, ~] = uigetfile('*.xlsx',...
%         'Select input file');  % get GUI file selector
%     filein = fopen([inputpath inputfile],'r');    % read only
%     if filein < 0, disp('! File not found or file sharing problem !'), end
% end
% fclose(filein);
%% Read Data
ts = tic;
% Para = xlsread(inputfile,'Cost_data', 'C6:C30');
% F = readmatrix(inputfile,'Sheet','Block_Model','Range','Q5:AJ44');%'data5.xlsx');
%F = readmatrix(inputfile,'Sheet','Block_Model','Range','N5:N679');%'data5.xlsx');
%F = readmatrix(inputfile,'Sheet','Block_Model','Range','N5:N1759');%'data5.xlsx');
%% Optimization Parameter Initialization:
 

%%Optimization Parameter Initialization:
%    Price = Para(2)/31.1035;        % Selling Price of the mineral
%      Cref = Para(4);                % Cost of selling + Refinery Cost 
%       Rec = Para(10);               % Recovery Matrix of the ore
%      Cmin = Para(6);                % Mining Cost per tonne of material
%      Cpro = Para(8);                % Processing cost per tonne of ore
%    BlkTon = Para(1)*10;             % Block Tonnage            
%    Cutoff = 4.5%3.5/31.1035; %Para(12)/31.10348;      % Stope Cutoff grade
    %Grade = reshape(F,27,65);       % Block Grades
%     Grade = F;                     % Block Grades
    % Grade = F(1:18,10:22); % inserted to make a smaller problem for testing
    % Grade = reshape(Grade,[],1);
%     [r,c] = size(Grade);            % Dimensions
 %% Economic Block Value calculation:   Equation 1
% EBV = (((Price-Cref)*Grade*Rec)-(Cmin+Cpro))*2000;
%EBV = EBV;
%% Geometrical and Geotechnical Parameter Initialization:
%     n = r;               % Number of blocks in the x direction in Model
%     J = c;                % Number of blocks in the y direction in Model
[n, J] = size(EBV);
% alpha = 2; %Para(14);        % Minimum Mining Width
% 
% gamma = 2; %Para(18);%5;     % Minimum Pillar Width
% 
%     S = 20; %Para(16);        % Number of stopes
%  beta1 = 3;
%  beta2 = 3;

alpha1 = 4; %Para(14);          % Minimum Mining Width in the x direction
alpha2 = 4; %Para(14);          % Minimum Mining Width in the z direction
beta1 = 4;                      % Maximum Mining Width in the x direction
beta2 = 4;                      % Maximum Mining Width in the z direction
gamma1 = 3; %Para(18);%5;       % Minimum Pillar Width in the x direction
gamma2 = 3; %Para(18);%5;       % Minimum Pillar Width in the z direction

    S = 20; %Para(16);           % Number of stopes
%% Mathematical Model:

% Decision Varibales Definition:
  prob = optimproblem('ObjectiveSense','maximize');
     y = optimvar('y',n,J,S,'Type','integer','LowerBound',0,'UpperBound',1);
     z1 = optimvar('z1',n,J,S,'Type','integer','LowerBound',0,'UpperBound',1);
     w1 = optimvar('w1',n,J,'Type','integer','LowerBound',0,'UpperBound',1);

     z2 = optimvar('z2',n,J,S,'Type','integer','LowerBound',0,'UpperBound',1);
     w2 = optimvar('w2',n,J,'Type','integer','LowerBound',0,'UpperBound',1);

% z1 = optimvar('z1',I,J,K,'Type','integer','LowerBound',0,'UpperBound',1);
% w1 = optimvar('w1',I,J,'Type','integer','LowerBound',0,'UpperBound',1);
% 
% z2 = optimvar('z2',I,J,K,'Type','integer','LowerBound',0,'UpperBound',1);
% w2 = optimvar('w2',I,J,'Type','integer','LowerBound',0,'UpperBound',1);

%%
 % Objective function : Equation 2 
prob.Objective = 0;
for k = 1:S
    prob.Objective = prob.Objective + sum(sum(y(:,:,k).*EBV));
end

%% Variable by variable constraints
prob.Constraints.cons1 = optimconstr(n,J,S);
prob.Constraints.cons2 = optimconstr(n,J,S);
prob.Constraints.cons3 = optimconstr(n,J,S);
prob.Constraints.cons4 = optimconstr(n,J,S);
prob.Constraints.cons5 = optimconstr(n,J,S);
prob.Constraints.cons6 = optimconstr(n,J,S);
prob.Constraints.eq6h = optimconstr(n,J,S);
prob.Constraints.eq6v = optimconstr(n,J,S);

for i=1:n
    for j=1:J
        for k = 1:S
            if i > 2
                prob.Constraints.cons1(i,j,k) = y(i,j,k) - y(i-1,j,k) <= z1(i,j,k); 
            end
            
            if j > 2
                prob.Constraints.cons2(i,j,k) = y(i,j,k) - y(i,j-1,k) <= z2(i,j,k);
            end

            if i > 2
                prob.Constraints.cons3(i,j,k) = 1 - y(i-1,j,k) >= z1(i,j,k); 
            end
            
            if j > 2
                prob.Constraints.cons4(i,j,k) = 1 - y(i,j-1,k) >= z2(i,j,k);
            end
            prob.Constraints.cons5(i,j,k) = sum(z1(max(i-alpha1 + 1,1):i,j,k)) <= y(i,j,k); % Equation 5
            prob.Constraints.cons6(i,j,k) = sum(z2(i,max(j-alpha2 + 1,1):j,k)) <= y(i,j,k); % Equation 6
            prob.Constraints.eq6h(i,j,k) = z1(i,j,k) <= beta1 - sum(y(min(n,i+1):min(n,i+beta1),j,k));
            prob.Constraints.eq6v(i,j,k) = z2(i,j,k) <= beta2 - sum(y(i,min(j+1,J):min(J,j+beta2),k));
        end
    end
end



%%
%Constraint:---- Mining Block Control: Equation 8
prob.Constraints.cons7  = sum(y,3) <= 1; %Ensure only 1 block per stope
% % prob.Constraints.cons200  = sum(z1,3) <= 1;
% % prob.Constraints.cons201  = sum(z2,3) <= 1;


%% Stope by stope constraints
prob.Constraints.cons8 = optimconstr(S,1); % Initialize Cut-off grade constraints
%prob.Constraints.cons100h = optimconstr(n,S);
prob.Constraints.cons100v = optimconstr(n,J,S);
for k = 1:S

 
    prob.Constraints.cons8(k) = sum(sum(Grade.*y(:,:,k))) - Cutoff*sum(sum(y(:,:,k))) >= 0; % Cut-off grade constraints (Equation 9)
    prob.Constraints.cons100v(:,:,k) = beta1 + beta2 >= sum(sum(z1(:,:,k),1)) + sum(sum(z2(:,:,k),2));
    %prob.Constraints.cons100h(:,k) = beta2 >= sum(sum(z2(:,:,k),2));

    %gamma >= beta1 >=9
%     prob.Constraints.cons100v(:,k) = beta1 >= sum(sum(z1(:,:,k),1)) - alpha;
%     prob.Constraints.cons100h(:,k) = beta2 >= sum(sum(z2(:,:,k),2)) - alpha;


end 
%% Block by block constraints
prob.Constraints.cons9  = optimconstr(n,J); % Initialize Equation 10
prob.Constraints.cons10 = optimconstr(n,J); % Initialize Equation 11
prob.Constraints.cons11 = optimconstr(n,J); % Initialize Equation 12
prob.Constraints.cons12 = optimconstr(n,J); % Initialize Equation 13
for i=1:n
    for j=1:J
        prob.Constraints.cons9(i,j)  = sum(w1(max(i-gamma1+1,1):i,j)) <=1- sum(y(i,j,:)); % Equation 10
        prob.Constraints.cons10(i,j) = sum(w2(i,max(j-gamma2+1,1):j)) <=1- sum(y(i,j,:)); % Equation 11
        if i > 1; prob.Constraints.cons11(i,j) = sum(y(i,j,:)) - sum(y(i - 1,j,:)) == sum(z1(i,j,:)) - w1(i,j); end % Equation 12
        if j > 1; prob.Constraints.cons12(i,j) = sum(y(i,j,:)) - sum(y(i,j - 1,:)) == sum(z2(i,j,:)) - w2(i,j); end % Equation 13
    end
end 


%% Solve problem
options             = optimoptions('intlinprog','MaxTime',9999999999999999999999);
[sol,OFV,xFlag] = solve(prob,'Options',options);

Tm = toc(ts);  %End Time
stp = S;
mm = alpha1;
mp = gamma1;
end 