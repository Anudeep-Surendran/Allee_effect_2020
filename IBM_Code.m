%% Script for simulating the individual based model(IBM) in the article 
%% Surendran et al. (2020), Population dynamics with spatial structure and Allee effects. 
%
%% Author: Anudeep Surendran
%          anudeep.surendran@hdr.qut.edu.au
%          School of Mathematical Sciences
%          Queensland University of Technology
%          Brisbane, Australia
%
%% Last update: 01 May 2020
%
%% Instructions:
%
% 1- This script outputs a Matlab data file ('IBM.mat') consisting of the data corresponds to the 
%    density of individuals as a function of time and pair correlation function (PCFs) 
%    computed at t=100 as a function of separation distance from IBM.
%
% 2- Change the model parameters below to simulate various cases considered in the article and any other new senarios.

%%.................Model parameters........................................ 
% Initial population size
Initial_size=400;
% Neighbour dependent interaction strengths
gammap=0.009;gammad=0.009;

% Spatial extent of interactions
sigmap=4.0;sigmad=4.0;

% Dispersal range
sigmadi=4.0;

% Mean and standard deviation movement distance
mus=0.4;sigmas=0.1;

% Intrinsic event rates
p=0.2;d=0.4;m=0.1;

%Domain Length
L=20.0;

% Parameters and variables introduced for simulation of the IBM and calculation of the PCF
N_realizations=200; % Number of realizations
t_final=100; % Final time
xi_final=8.0;xi_initial=0.08;xi_increment=0.04;dxi=0.2; % Parameters for binning the distance($|\xi|$) to calculate PCFs 
steps=5000000; % Number of steps used in the Gillespie algorithm

%.................Preallocating arrays.....................................

t_ar=0:.2:t_final; % Time points at which density of consumers and resources are recorded for averaging.
xi_ar=xi_initial:xi_increment:xi_final;% For binning the distance
Z_1_ar_final=zeros(size(t_ar));% Array to store density of individuals 
Z_1_ar=zeros(size(t_ar));% Array to store density of individuals (in each realizations)
C_ar_final=zeros(size(xi_ar));% Array to store PCF values.

%.....................The actual code begins from here.....................

% Loop for creating ensamble of realizations
for rlzn=1:N_realizations
           
    % Initializing the population and calculating initial rates 
    t=0; 
    num=Initial_size;
    x_initial=(rand(1,num)*L)+(-L/2.0);y_initial=(rand(1,num)*L)+(-L/2.0); % Individuals are randomely placed on LxL domain
    x=x_initial;y=y_initial; % x and y coordinates of the locations of individuals.
    
    % Pre allocating arrays
    distance_ar=zeros(num); % Distance matrix storing the separation distance of all individuals from each other
    M_ar=m*ones(1,num);P_ar=zeros(1,num);D_ar=zeros(1,num); % Array for storing the Movement, proliferation and death rates of each individuals    
    D_square_ar=zeros(1,num);
    count=1;Z_1_ar(count)=num/(L^2); % Computing initial density of consumers and resources    
    for i=1:num % In these loops, we calculate the neighbour dependent event rates  as a result of the interaction of individuals
        for j=1:num
            xi_x=periodic(x(j)-x(i),L);xi_y=periodic(y(j)-y(i),L);
            xi_square=xi_x^2+xi_y^2;  
            distance_ar(i,j)=sqrt(xi_square); % Distance between individuals            
            if (j ~= i)  % Neglecting the self pairs                
                D_ar(i)=D_ar(i)+w(xi_square,gammad,sigmad);
                P_ar(i)=P_ar(i)+w(xi_square,gammap,sigmap);
            end
        end
        D_square_ar(i)=d+(D_ar(i))^2;P_ar(i)=p+P_ar(i);
    end
%................Gillespie simulation starts from here.....................

    rn=rand(1,steps); % uniform random numbers to use in Gillespie algorithm
    for ns=1:steps    
        % Calculating net event rates and time increment        
		rate_P = sum(P_ar);              % Total proliferation rate
		rate_D = sum(D_square_ar);       % Total death rate
        rate_M = sum(M_ar);              % Total movement rate
		net_rate = rate_P+rate_D+rate_M; % Total events rate		
		dt=exprnd(1/net_rate);  % Time step for the event (exponentially distributed) 
        
        if (rn(ns) <= (rate_P/net_rate)) % If the event is proliferation
            cum_sum=cumsum(P_ar);% Cumulative sum of proliferation rates
            R=rand*rate_P; % Uniform random number on interval [0,rate_P]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event
            dispersal=mvnrnd(0.0,sigmadi,2);% New location of individual from bivariate normal distribution
            new_loc_x=periodic(x(choose)+dispersal(1),L);new_loc_y=periodic(y(choose)+dispersal(2),L);
            x=[x(1:choose) new_loc_x x(choose+1:end)];y=[y(1:choose) new_loc_y y(choose+1:end)];            
            M_ar=[M_ar(1:choose) m M_ar(choose+1:end)];
            P_ar=[P_ar(1:choose) 0.0 P_ar(choose+1:end)];
            D_ar=[D_ar(1:choose) 0.0 D_ar(choose+1:end)];
            D_square_ar=[D_square_ar(1:choose) 0.0 D_square_ar(choose+1:end)];
            row=zeros(1,num+1);coloumn=zeros(num+1,1);
            for i=1:num+1
                xi_x=periodic(x(i)-x(choose+1),L);xi_y=periodic(y(i)-y(choose+1),L);
                xi_square=xi_x^2+xi_y^2;
                row(i)=sqrt(xi_square);coloumn(i)=sqrt(xi_square);% To record distance from daughter individual                
                if (i ~= choose+1)
                    D_ar(choose+1)=D_ar(choose+1)+w(xi_square,gammad,sigmad);
                    P_ar(choose+1)=P_ar(choose+1)+w(xi_square,gammap,sigmap);
                    D_ar(i)=D_ar(i)+w(xi_square,gammad,sigmad);  
                    P_ar(i)=P_ar(i)+w(xi_square,gammap,sigmap); 
                    D_square_ar(i)=d+(D_ar(i))^2;                   
                end
            end
            D_square_ar(choose+1)=d+(D_ar(choose+1))^2;P_ar(choose+1)=P_ar(choose+1)+p;
            dummy=zeros(num+1);% Temporary 
            dummy(1:choose,1:choose)=distance_ar(1:choose,1:choose);
            dummy(1:choose,choose+2:end)=distance_ar(1:choose,choose+1:end);
            dummy(choose+2:end,1:choose)=distance_ar(choose+1:end,1:choose);
            dummy(choose+2:end,choose+2:end)=distance_ar(choose+1:end,choose+1:end);
            dummy(1:end,choose+1)=coloumn;dummy(choose+1,1:end)=row;
            distance_ar=dummy;
            num=num+1;% Updating the number of individuals after proliferation event
            
        elseif (rn(ns) > (rate_P/net_rate) && rn(ns) <= ((rate_P+rate_M)/net_rate)) % If the event is movement           
            cum_sum=cumsum(M_ar);% Cumulative sum of movement rate
            R=rand*rate_M; % Uniform random number on interval [0,rate_M]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event
            old_loc_x=x(choose);old_loc_y=y(choose); % Previous location of the individual chosen for movement event
            magnitude_xi=-1; % Computing the movement distance
            while (magnitude_xi < mus-(4.0*sigmas) || magnitude_xi > (mus+(4.0*sigmas))) % Truncation 
                magnitude_xi=normrnd(mus,sigmas);
            end
            angle_xi=2.0*pi*rand; % Uniformly distributed direction angle
            new_loc_x=x(choose)+(magnitude_xi*cos(angle_xi)); % Updating the new location after the movement
            new_loc_y=y(choose)+(magnitude_xi*sin(angle_xi)); %                      ""
            x(choose)=periodic(new_loc_x,L);y(choose)=periodic(new_loc_y,L);%        ""            
            M_ar(choose)=m;P_ar(choose)=0.0;D_ar(choose)=0.0; % Updating the event rates after the movement event
            for i=1:num % Updating the event rates after the movement event           
                xi_x=periodic(x(i)-x(choose),L);xi_y=periodic(y(i)-y(choose),L);
                xi_x_old=periodic(x(i)-old_loc_x,L);xi_y_old=periodic(y(i)-old_loc_y,L);
                xi_square=xi_x^2+xi_y^2;xi_square_old=xi_x_old^2+xi_y_old^2;
                distance_ar(choose,i)=sqrt(xi_square); % Distance between individuals
                distance_ar(i,choose)=sqrt(xi_square);                
                if (i ~= choose)
                    D_ar(choose)=D_ar(choose)+w(xi_square,gammad,sigmad);
                    D_ar(i)=D_ar(i)+w(xi_square,gammad,sigmad)-w(xi_square_old,gammad,sigmad);
                    P_ar(choose)=P_ar(choose)+w(xi_square,gammap,sigmap);
                    P_ar(i)=P_ar(i)+w(xi_square,gammap,sigmap)-w(xi_square_old,gammap,sigmap);
                    D_square_ar(i)=d+(D_ar(i))^2;
                end
            end
            P_ar(choose)=P_ar(choose)+p;D_square_ar(choose)=d+(D_ar(choose))^2;
            
        else % If the event is death
            cum_sum=cumsum(D_square_ar); % Cumulative sum of death rates
            R=rand*rate_D; % Uniform random number on interval [0,rate_D]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event            
            for i=1:num
                xi_x=periodic(x(choose)-x(i),L);xi_y=periodic(y(choose)-y(i),L);
                xi_square=xi_x^2+xi_y^2;
                if (i ~= choose)
                    D_ar(i)=D_ar(i)-w(xi_square,gammad,sigmad);
                    D_square_ar(i)=d+(D_ar(i))^2;
                    P_ar(i)=P_ar(i)-w(xi_square,gammap,sigmap);
                end
            end
            x=[x(1:choose-1) x(choose+1:end)];y=[y(1:choose-1) y(choose+1:end)];
            M_ar=[M_ar(1:choose-1) M_ar(choose+1:end)];P_ar=[P_ar(1:choose-1) P_ar(choose+1:end)];
            D_ar=[D_ar(1:choose-1) D_ar(choose+1:end)];D_square_ar=[D_square_ar(1:choose-1) D_square_ar(choose+1:end)];
            distance_ar=[distance_ar(1:choose-1,1:choose-1) distance_ar(1:choose-1,choose+1:end);
            distance_ar(choose+1:end,1:choose-1) distance_ar(choose+1:end,choose+1:end)];
            num=num-1; % Updating the number of individuals after death event
        end
        
        t=t+dt; % Updating the time
        if (num==0)
            Z_1_ar(count:end)=0.0;
            break
        end
        if (t >= t_ar(count))
            Z_1_ar(count)=num/(L^2); % Storing the desnsity of individuals at time t
            if (t >= t_final)% Computing the PCF, C(|\xi|, t), at t=t_final
                if (num>1)
                    C_ar=zeros(size(xi_ar));% Pre-allocation
                    for ab=1:length(xi_ar)
                        xi=xi_ar(ab);
                        count_C=0;
                        for i=1:num
                            for j=1:num
                                if (i ~= j)% Counting distances excluding the self distance (when i=j).
                                    xival=distance_ar(j,i);
                                    if (xival >= (xi-(dxi/2.0)) && xival <= (xi+(dxi/2.0))) % Binning the distances
                                        count_C=count_C+1.0;
                                    end
                                end
                            end
                        end
                        C_ar(ab)=count_C/((num*(num-1.0))*(2.0*pi*xi*dxi)/(L^2));% Normalizing the count of distance to plot PCFs                         
                    end
                    C_ar_final=C_ar_final+C_ar;
                end
                break
            end
            count=count+1;
        end
    end
    Z_1_ar_final=Z_1_ar_final+Z_1_ar;
end
% Averaged results from IBM
Z_1_ar_final=Z_1_ar_final/N_realizations;
C_ar_final=C_ar_final/N_realizations;
%...........................IBM simulation ends here.............................

% Writing the data into 'IBM.mat' file
z1_data=zeros(length(t_ar),2);z1_data(:,1)=t_ar;z1_data(:,2)=Z_1_ar_final;
pcf_data=zeros(length(xi_ar),2);pcf_data(:,1)=xi_ar;pcf_data(:,2)=C_ar_final;
save('IBM.mat','z1_data','pcf_data');


%.......................Main code ends here................................

 function a = periodic(x,L)
% Function to impliment periodic boundary condition.
% Input   :  x  =  x or y cordinate of an individual's location
%		     L	=  spatial domain length
%
% Output  :  impliment periodic boundary condition 
    if x<(-L/2.0)	
        a=x+L;
    elseif x>=(L/2.0)	
		a=x-L;
    else	
        a=x;
    end
end

function a = w(dsquare,gamma,sigma) 
% Function to calculate interaction kernel
% Input   :  dsquare =  distance square
%		     gamma   =  interaction strength
%            sigma   =  spatial extent of interaction
%
% Output  :  value of interaction kernel
a=gamma*exp(-(dsquare/(2.0*(sigma^2))));
end                 