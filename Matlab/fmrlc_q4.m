%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q4: Implementation of FMRLC
%Author: Gansheng TAN, aegean0045@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specifications:
% This script is based on the question 4 of the Intelligent control courses
% assignment. Please refer to question description for further information 
% about the system transfer function and refernce signal.
% The refreshing rate of FLC is as same as the input sampling rate.
% tstop is set to 15 second for a reason. If you are trying to modify this
% parameter, please change the enlarging scale of membership function for
% u(input to the target system) as well. The enlarging scale in this
% example is 5.
clear;
close all;
clc;
%% parameters
k_p=1.5;
zeta_p=0.707;
omega_p=1.2;
k_r=1.5;
a_r=1;
M = tf(k_r,[1 a_r]);
G= tf(k_p,[1 2*zeta_p*omega_p omega_p^2]);

% FLC design whose inputs are e and edot, output is u
% the universe of disclose is [-1,1]
num_mfs_e = 3;
num_mfs_edot = 3;
num_mfs_u = 3;
% gains
g_e = 0.8;
g_edot = 0.8;
g_u = 0.8;
% mfs: triangular mfs are used in this question
hwidth_mfs_e = 2/num_mfs_e; % half of the triangle base width
hwidth_mfs_edot = 2/num_mfs_edot;
hwidth_mfs_u = 2/num_mfs_u; % base width of the output mfs
center_mfs_e = [-1 0 1];
center_mfs_edot= [-1 0 1];
% rules
% e\edot -1   0    1
% -1      1  0.5   0
% 0      0.5  0   -0.5
% 1       0  -0.5  -1
rules= [ -1   -0.5   0;
        -0.5    0   0.5;
          0    0.5   1];
 
% fuzzy inverse model (Guessing) that takes ye and yedot as inputs and
% outputs p
g_ye =2;
g_yedot =1;
g_p = 0.4;
num_mfs_ye = 3;
num_mfs_yedot = 3;
num_mfs_p = 3;
hwidth_mfs_ye = 2/num_mfs_ye;
hwidth_mfs_yedot = 2/num_mfs_yedot;
hwidth_mfs_p = 2/num_mfs_p;
center_mfs_ye = [-1 0 1];
center_mfs_yedot = [-1 0 1];
inverseFLCrules = [-1  -0.5  0;
                  -0.5  0   0.5;
                    0  0.5   1];
%% knowledge-based modifier
memory_step = 1; %number of times that the modifier looks back

%% simulation
% initialization
begin_simulation=input('Press 1 to confirm the simulation?  \n (type 0 to exit) ');

if begin_simulation==1
t = 0;
index = 1;
tstop = 15; %in sec
int_step = 0.1; %integration step
sampling_time = 0.1; % the inverse of the sampling frequency
counter_int = 10; % begins at 0.1 and counts the times of the integration
                   % it determine the refreshing rate of FLC
e_before =0;
r_ini =0;
y=0; %initialize y
u=0; %initialize system input u
p=0;
ye=0;
yedot=0;
t_vector = linspace(0,tstop+sampling_time,tstop/sampling_time+2);
r = sin(0.6*t_vector);
%% reference model
[ym t_lsim] = lsim(M,r,t_vector);
% run the loop
while t<=tstop
    
%% FLC behavior
% manipulation the begining values
e(index)=r(index)-y(index);
edot(index) = (e(index)-e_before)/sampling_time;
e_before=e(index);
num_activated_mfs_e = 0;
num_activated_mfs_edot = 0;
e_times_ge=e*g_e;
edot_times_gedot = edot*g_edot;
if e_times_ge(index)<=center_mfs_e(1)   % extreme value
    activation_mfs_e=[1 0 0];
    lastID_mfs_e = 1;
    num_activated_mfs_e = 1+num_activated_mfs_e;
elseif e_times_ge(index)>=center_mfs_e(num_mfs_e)
    activation_mfs_e=[0 0 1];
    lastID_mfs_e = num_mfs_e;
    num_activated_mfs_e = 1+num_activated_mfs_e;
else
    for i =1:num_mfs_e
        if e_times_ge(index)<=center_mfs_e(i)
            activation_mfs_e(i)=max(0,1+(e_times_ge(index)-center_mfs_e(i))/hwidth_mfs_e);
            if activation_mfs_e(i)~=0
                num_activated_mfs_e = 1+num_activated_mfs_e;
                lastID_mfs_e=i;
            end
        else
            activation_mfs_e(i)=max(0,1+(center_mfs_e(i)-e_times_ge(index))/hwidth_mfs_e);
            if activation_mfs_e(i)~=0
                num_activated_mfs_e = 1+num_activated_mfs_e;
                lastID_mfs_e=i;
            end
        end
    end
end

% similar operation on edot, I could have done better if I used Python
if edot_times_gedot(index)<=center_mfs_edot(1)   % extreme value
    activation_mfs_edot=[1 0 0];
    lastID_mfs_edot = 1;
    num_activated_mfs_edot = 1+num_activated_mfs_edot;
elseif edot_times_gedot(index)>=center_mfs_edot(num_mfs_edot)
    activation_mfs_edot=[0 0 1];
    lastID_mfs_edot = num_mfs_edot;
    num_activated_mfs_edot = 1+num_activated_mfs_edot;
else
    for i =1:num_mfs_edot
        if edot_times_gedot(index)<=center_mfs_edot(i)
            activation_mfs_edot(i)=max(0,1+(edot_times_gedot(index)-center_mfs_edot(i))/hwidth_mfs_edot);
            if activation_mfs_edot(i)~=0
                num_activated_mfs_edot = 1+num_activated_mfs_edot;
                lastID_mfs_edot=i;
            end
        else
            activation_mfs_edot(i)=max(0,1+(center_mfs_edot(i)-edot_times_gedot(index))/hwidth_mfs_edot);
            if activation_mfs_edot(i)~=0
                num_activated_mfs_edot = 1+num_activated_mfs_edot;
                lastID_mfs_edot=i;
            end
        end
    end
end

u_mfs_implied = cell(num_activated_mfs_e*num_activated_mfs_edot,1);
u_universe_of_disclose = linspace(-1,1);
u_mfs_implied_count = 1;
% implication and defuzzification
for i = lastID_mfs_e-num_activated_mfs_e+1:lastID_mfs_e
    for j = lastID_mfs_edot-num_activated_mfs_edot+1:lastID_mfs_edot
        J = min([activation_mfs_e(i) activation_mfs_edot(j)]);% T0 conditional function, we take min as T-norm
        trapmf1 = rules(i,j)-hwidth_mfs_u;
        trapmf2 = trapmf1+J*(rules(i,j)-trapmf1);
        trapmf3 = 2*rules(i,j)-trapmf2;
        trapmf4 = rules(i,j)+hwidth_mfs_u;
        u_mfs_implied{u_mfs_implied_count} = trapmf(5*u_universe_of_disclose,[trapmf1 trapmf2 ...
            trapmf3 trapmf4]);
        u_mfs_implied_count=u_mfs_implied_count+1;
    end
end
if num_activated_mfs_e*num_activated_mfs_edot==1 %take care of single rule activation
    u_max_mf_implied = cell2mat(u_mfs_implied);
else
    u_max_mf_implied = max(cell2mat(u_mfs_implied));
end
u_defuzzified = defuzz(5*u_universe_of_disclose,u_max_mf_implied,'centroid');
% 5 is a parameter that enlarge the universe of disclose since p would
% modify the center of membership function for u
%% system response
u(index+1) = u_defuzzified*g_u;
t_vector_u = linspace (0,(index+1)*sampling_time,index+1);
[y t_lsim] = lsim(G,u,t_vector_u);
%% fuzzy inverse model - similar to fuzzy model
ye(index+1)=ym(index+1)-y(index+1);
yedot(index+1)=(ye(index+1)-ye(index))/sampling_time;
num_activated_mfs_ye = 0;
num_activated_mfs_yedot = 0;
ye_times_gye = ye*g_ye;
yedot_times_gyedot = yedot*g_yedot;
% input 1:ye
if ye_times_gye(index+1)<=center_mfs_ye(1)
    activation_mfs_ye = [1 0 0];
    lastID_mfs_ye = 1;
    num_activated_mfs_ye = 1+num_activated_mfs_ye;
elseif ye_times_gye(index+1)>=center_mfs_ye(num_mfs_ye)
    activation_mfs_ye = [0 0 1];
    lastID_mfs_ye = num_mfs_ye;
    num_activated_mfs_ye = 1+num_activated_mfs_ye;
else
    for i = 1:num_mfs_ye
        if ye_times_gye(index+1)<=center_mfs_ye(i)
            activation_mfs_ye(i) = max([0 1+(ye_times_gye(index+1)-...
                center_mfs_ye(i))/hwidth_mfs_ye]);
            if activation_mfs_ye(i)~=0
                num_activated_mfs_ye = 1+num_activated_mfs_ye;
                lastID_mfs_ye = i;
            end
        else
            activation_mfs_ye(i) = max([0 1+(-ye_times_gye(index+1)+...
                center_mfs_ye(i))/hwidth_mfs_ye]);
            if activation_mfs_ye(i)~=0
                num_activated_mfs_ye = 1+num_activated_mfs_ye;
                lastID_mfs_ye = i;
            end
        end
    end    
end
% input 2: yedot
if yedot_times_gyedot(index+1)<=center_mfs_yedot(1)
    activation_mfs_yedot = [1 0 0];
    lastID_mfs_yedot = 1;
    num_activated_mfs_yedot = 1+num_activated_mfs_yedot;
elseif yedot_times_gyedot(index+1)>=center_mfs_yedot(num_mfs_yedot)
    activation_mfs_yedot = [0 0 1];
    lastID_mfs_yedot = num_mfs_yedot;
    num_activated_mfs_yedot = 1+num_activated_mfs_yedot;
else
    for i = 1:num_mfs_yedot
        if yedot_times_gyedot(index+1)<=center_mfs_yedot(i)
            activation_mfs_yedot(i) = max([0 1+(yedot_times_gyedot(index+1)-...
                center_mfs_yedot(i))/hwidth_mfs_yedot]);
            if activation_mfs_yedot(i)~=0
                num_activated_mfs_yedot = 1+num_activated_mfs_yedot;
                lastID_mfs_yedot = i;
            end
        else
            activation_mfs_yedot(i) = max([0 1+(-yedot_times_gyedot(index+1)+...
                center_mfs_yedot(i))/hwidth_mfs_yedot]);
            if activation_mfs_yedot(i)~=0
                num_activated_mfs_yedot = 1+num_activated_mfs_yedot;
                lastID_mfs_yedot = i;
            end
        end
    end    
end
% defuzzification for inverse FLC
p_mfs_implied = cell(num_activated_mfs_ye*num_activated_mfs_yedot,1);
p_universe_of_disclose = linspace(-1,1);
p_mfs_implied_count = 1;
% implication and defuzzification
for i = lastID_mfs_ye-num_activated_mfs_ye+1:lastID_mfs_ye
    for j = lastID_mfs_yedot-num_activated_mfs_yedot+1:lastID_mfs_yedot
        J = min([activation_mfs_ye(i) activation_mfs_yedot(j)]);% T0 conditional function, we take min as T-norm
        trapmf1 = inverseFLCrules(i,j)-hwidth_mfs_p;
        trapmf2 = trapmf1+J*(inverseFLCrules(i,j)-trapmf1);
        trapmf3 = 2*inverseFLCrules(i,j)-trapmf2;
        trapmf4 = inverseFLCrules(i,j)+hwidth_mfs_p;
        p_mfs_implied{p_mfs_implied_count} = trapmf(p_universe_of_disclose,[trapmf1 trapmf2 ...
            trapmf3 trapmf4]);
        p_mfs_implied_count=p_mfs_implied_count+1;
    end
end
if num_activated_mfs_ye*num_activated_mfs_yedot==1 %take care of single rule activation
    p_max_mf_implied = cell2mat(p_mfs_implied);
else
    p_max_mf_implied = max(cell2mat(p_mfs_implied));
end
p_defuzzified = defuzz(p_universe_of_disclose,p_max_mf_implied,'centroid');

p(index+1)=p_defuzzified*g_p;

%% Knowledge-based modifier
for i =lastID_mfs_e-num_activated_mfs_e+1:lastID_mfs_e
    for j = lastID_mfs_edot-num_activated_mfs_edot+1:lastID_mfs_edot
        rules (i,j) = rules(i,j)+p(index+1);
    end
end
t=t+sampling_time;
index=index+1;
end
end

%% test function
%system performance
figure(1)
clf
time = linspace(0,t,index);
plot(time,r(1:length(y)),'k-',time,ym(1:length(y)),'k--',time,y,'k-.')
legend('r','ym','y')
title('Performance of reference model and target model')
xlabel('time in second')

% FLC behavior
figure(1)
clf
time = linspace(0,t,index);
plot(time(1:length(y)-1),e(1:length(y)-1),'k-',time(1:length(y)-1),edot(1:length(y)-1),...
    'k--',time(1:length(y)-1),u(1:length(y)-1),'k-.')
legend('e','edot', 'u')
title('FLC behavior')
xlabel('time in second')

% inverse FLC behavior
figure(1)
clf
time = linspace(0,t,index);
plot(time(1:length(y)),ye(1:length(y)),'k-',time(1:length(y)),yedot(1:length(y)),...
    'k--',time(1:length(y)),p(1:length(y)),'k-.')
legend('ye','yedot', 'p')
title('inverse FLC behavior')
xlabel('time in second')

%% Extensions
% for further information, please contact Gansheng TAN