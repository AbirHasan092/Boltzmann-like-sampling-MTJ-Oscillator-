%% SOT_Elliptic_Ensemble_Ic_always_on.m
% Ensemble stochastic LLG with:
%   - signed Oersted-ring current sweep Ic
%   - direct hard-axis spin-current pulse Is
%   - Ic is constant throughout each trajectory, including Is initialization
%   - Gaussian device and direction-offset parameter sets
%   - initial-condition sensitivity testing
%
% Spin-current equation:
% (1+a^2) dm/dt = -gamma0[m x H + a*m x(m x H)]
%                 -(1/qNs)m x(m x Is) +(a/qNs)m x Is
%
% Debashis/Ostwal direct-spin-current stochastic macrospin model.
% Run SOT_10percent_Elliptic_Parameter_set_generator.m first.

clear; clc; close all;

%% User controls
script_dir = fileparts(mfilename('fullpath'));
parameter_file = fullfile(script_dir, ...
    'INI_parameter_sets_elliptic.mat');
number_of_sweeps = 50;
representative_sweep = 1;
Ic_min = -25e-3;                  % [A], Debashis ring-current range
Ic_max =  25e-3;                  % [A]
Nic_points = 45;
Ic_trace_target = 10e-3;         % [A], representative trajectory only
Nens = 200;                       % Debashis simulation averaging count
base_seed = 120;

% Common initial state for device-parameter variation.
parameter_theta0_deg = 0;
parameter_phi0_deg   = 0;

% Initial-condition sensitivity test.
initial_theta_deg = [0];
initial_phi_deg   = [0];

% Hard-axis spin-current pulse strength. Is is the spin current appearing
% in the torque equation, not the charge current through the Ta Hall bar.
% The default is expressed relative to q*Ns*gamma0*Hk. A measured Hall-bar
% charge current must be converted to Is before it is used as an override.
Is_pulse_scale = 1.2;
Is_pulse_override = [0.18];           % [A]; [] uses the scale above

%% Constants and numerical controls
kB = 1.380649e-23;
mu0 = 4*pi*1e-7;
qe = 1.602176634e-19;
muB = 9.2740100783e-24;           % [A m^2]
dt = 1e-12;
Tsim = 50e-9;
Nstep = round(Tsim/dt);
time = (0:Nstep)*dt;
mz_switch_threshold = 0;

% Oersted-ring current remains the probability-sweep variable.
Ic_values = linspace(Ic_min,Ic_max,Nic_points);
Nic = numel(Ic_values);
Ic_plot_limits_mA = [min(Ic_values) max(Ic_values)]*1e3;
Ic_plot_ticks_mA = linspace(Ic_plot_limits_mA(1),Ic_plot_limits_mA(2),5);
[~,kTrace] = min(abs(Ic_values-Ic_trace_target));
Ic_trace_actual = Ic_values(kTrace);

%% Load parameter sets
if ~isfile(parameter_file)
    error(['Missing %s. Run generate_device_parameter_sets.m first, ' ...
        'with both scripts in the same folder.'],parameter_file);
end
S = load(parameter_file,'parameter_sets','baseline','variation_flags');
parameter_sets = S.parameter_sets;
baseline = S.baseline;
if ~isfield(baseline,'gamma') || ~isscalar(baseline.gamma) || ...
        ~isfinite(baseline.gamma) || baseline.gamma<=0
    error('The generated baseline must contain a positive gamma.');
end
gamma = baseline.gamma;            % [rad/(s T)]
gamma0 = mu0*gamma;                % [m/(A s)]
if ~isfield(S,'variation_flags')
    error(['The parameter file does not contain variation_flags. ' ...
        'Run the updated generate_device_parameter_sets.m first.']);
end
variation_flags = S.variation_flags;
required_flags = {'vary_Is_direction_with_device_variation', ...
    'vary_Ic_direction_with_device_variation','only_vary_Is_direction'};
for flag_index=1:numel(required_flags)
    if ~isfield(variation_flags,required_flags{flag_index})
        error('variation_flags is missing %s.',required_flags{flag_index});
    end
end
if number_of_sweeps>numel(parameter_sets)
    error('Requested %d sweeps, but the file contains only %d.', ...
        number_of_sweeps,numel(parameter_sets));
end
if representative_sweep<1 || representative_sweep>number_of_sweeps
    error('representative_sweep must be between 1 and number_of_sweeps.');
end

plot_parameter_boxplots(parameter_sets(1:number_of_sweeps),baseline);

Hk_baseline = baseline.Hk_Oe*1000/(4*pi);
if isempty(Is_pulse_override)
    Is_pulse = Is_pulse_scale*qe*baseline.Ns*gamma0*Hk_baseline;
else
    Is_pulse = Is_pulse_override;
end
if ~isscalar(Is_pulse) || ~isfinite(Is_pulse)
    error('The hard-axis Is pulse must be a finite scalar in amperes.');
end
fprintf('Baseline volume   = %.6e m^3\n',baseline.V);
fprintf('Baseline Ns       = %.6e\n',baseline.Ns);
fprintf('Baseline Hk       = %.1f mT = %.1f Oe\n', ...
    baseline.Hk_mT,baseline.Hk_Oe);
fprintf('Hard-axis Is pulse= %.6e A\n',Is_pulse);

Pplus_all = zeros(number_of_sweeps,Nic);
HOe_all = zeros(number_of_sweeps,Nic);
trace = struct('m',[],'Is',[],'term_precession',[], ...
    'term_field_damping',[],'term_Is_damping_like',[], ...
    'term_Is_alpha',[],'stored',false);

%% Device-parameter sweeps
for isweep=1:number_of_sweeps
    p=parameter_sets(isweep);
    fprintf('\n========== Sweep %d of %d ==========\n',isweep,number_of_sweeps);
    store_this_sweep=isweep==representative_sweep;
    [Pplus_all(isweep,:),HOe_all(isweep,:),candidate_trace] = ...
        run_probability_curve(p,Ic_values,Is_pulse,Nens, ...
        base_seed+100000*(isweep-1),dt,Nstep,kTrace,store_this_sweep, ...
        parameter_theta0_deg,parameter_phi0_deg,kB,mu0,qe,muB,gamma0, ...
        mz_switch_threshold);
    if store_this_sweep
        trace=candidate_trace;
    end
end

%% Figure 6: P(+z) versus Oersted-ring current Ic
figProbability=figure('Color','w', ...
    'Name','Figure 6 - Switching probability','NumberTitle','off');
clf(figProbability); hold on;
colors=turbo(number_of_sweeps);
h=gobjects(number_of_sweeps,1);
for isweep=1:number_of_sweeps
    c=colors(isweep,:);
    h(isweep)=plot(Ic_values*1e3,Pplus_all(isweep,:),'o-', ...
        'LineWidth',2.5,'MarkerSize',8,'Color',c, ...
        'MarkerFaceColor',c,'MarkerEdgeColor',c);
end
labels=arrayfun(@(s) sprintf('Sweep %d',s),1:number_of_sweeps, ...
    'UniformOutput',false);
legend(h,labels,'Location','best','FontSize',13);
xlabel('Oersted ring current I_c [mA]','FontSize',22);
ylabel('P(+z)','FontSize',22);
set(gca,'FontSize',22,'LineWidth',2.5,'TickDir','in');
xlim(Ic_plot_limits_mA); xticks(Ic_plot_ticks_mA); ylim([-0.05 1.05]);
box on; grid off;
% The representative trajectory point is intentionally not marked.

%% Initial-condition sensitivity at the selected parameter set
if numel(initial_theta_deg)~=numel(initial_phi_deg)
    error('initial_theta_deg and initial_phi_deg must have equal lengths.');
end
Nsensitivity=numel(initial_theta_deg);
Pplus_initial=zeros(Nsensitivity,Nic);
p=parameter_sets(representative_sweep);
for q=1:Nsensitivity
    fprintf('\nInitial-condition test %d/%d: theta=%.1f, phi=%.1f deg\n', ...
        q,Nsensitivity,initial_theta_deg(q),initial_phi_deg(q));
    Pplus_initial(q,:)=run_probability_curve(p,Ic_values,Is_pulse,Nens, ...
        base_seed+700000+100000*(q-1),dt,Nstep,kTrace,false, ...
        initial_theta_deg(q),initial_phi_deg(q),kB,mu0,qe,muB,gamma0, ...
        mz_switch_threshold);
end

figSensitivity=figure('Color','w', ...
    'Name','Figure 8 - Initial-condition sensitivity','NumberTitle','off');
clf(figSensitivity); hold on;
sens_colors=turbo(Nsensitivity);
hs=gobjects(Nsensitivity,1);
sens_labels=cell(Nsensitivity,1);
for q=1:Nsensitivity
    c=sens_colors(q,:);
    hs(q)=plot(Ic_values*1e3,Pplus_initial(q,:),'o-', ...
        'LineWidth',2.2,'MarkerSize',7,'Color',c, ...
        'MarkerFaceColor',c,'MarkerEdgeColor',c);
    sens_labels{q}=sprintf('\\theta_0=%.1f^\\circ, \\phi_0=%.1f^\\circ', ...
        initial_theta_deg(q),initial_phi_deg(q));
end
legend(hs,sens_labels,'Location','best','FontSize',13);
xlabel('Oersted ring current I_c [mA]','FontSize',22);
ylabel('P(+z)','FontSize',22);
set(gca,'FontSize',22,'LineWidth',2.5,'TickDir','in');
xlim(Ic_plot_limits_mA); xticks(Ic_plot_ticks_mA); ylim([-0.05 1.05]);
box on; grid off;
title(sprintf('Initial-condition sensitivity, parameter Sweep %d', ...
    representative_sweep));

%% Representative trajectory figures
mxTrace=trace.m(1,:); myTrace=trace.m(2,:); mzTrace=trace.m(3,:);
thetaTrace=acos(max(min(mzTrace,1),-1));
phiTrace=unwrap(atan2(myTrace,mxTrace));
p=parameter_sets(representative_sweep);

figCartesian=figure('Color','w','Name','Figure 1 - Cartesian trajectory', ...
    'NumberTitle','off'); clf(figCartesian);
plot(time*1e9,mxTrace,'LineWidth',1.1); hold on;
plot(time*1e9,myTrace,'LineWidth',1.1);
plot(time*1e9,mzTrace,'LineWidth',1.4);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
xlabel('Time [ns]'); ylabel('m_i');
legend('m_x','m_y','m_z','Location','best');
ylim([-1.05 1.05]); xlim([0 Tsim*1e9]); grid on;
title(sprintf('Trajectory, Sweep %d, I_c=%.3f mA', ...
    representative_sweep,Ic_trace_actual*1e3));

figAngles=figure('Color','w','Name','Figure 2 - Magnetization angles', ...
    'NumberTitle','off'); clf(figAngles);
subplot(2,1,1); plot(time*1e9,thetaTrace*180/pi,'LineWidth',1.3);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
ylabel('\theta [deg]'); grid on; title('Polar angle');
subplot(2,1,2); plot(time*1e9,phiTrace*180/pi,'LineWidth',1.3);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
xlabel('Time [ns]'); ylabel('\phi [deg]'); grid on; title('Azimuthal angle');

figPulse=figure('Color','w','Name','Figure 3 - Hard-axis spin current', ...
    'NumberTitle','off'); clf(figPulse);
plot(time*1e9,trace.Is(1,:),'LineWidth',1.4); hold on;
plot(time*1e9,trace.Is(2,:),'LineWidth',1.4);
plot(time*1e9,trace.Is(3,:),'LineWidth',1.4);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
xlabel('Time [ns]'); ylabel('I_s components [A]');
legend('I_{s,x}','I_{s,y}','I_{s,z}','Location','best');
grid on; title('Hard-axis spin-current initialization pulse');

figCombined=figure('Color','w','Name','Figure 4 - Magnetization and Is', ...
    'NumberTitle','off'); clf(figCombined);
subplot(2,1,1); plot(time*1e9,mxTrace,'LineWidth',1); hold on;
plot(time*1e9,myTrace,'LineWidth',1); plot(time*1e9,mzTrace,'LineWidth',1.3);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
ylabel('m_i'); legend('m_x','m_y','m_z','Location','best'); grid on;
subplot(2,1,2); plot(time*1e9,vecnorm(trace.Is,2,1),'LineWidth',1.4);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
xlabel('Time [ns]'); ylabel('|I_s| [A]'); grid on;

figSphere=figure('Color','w','Name','Figure 5 - Unit-sphere trajectory', ...
    'NumberTitle','off'); clf(figSphere);
[Xs,Ys,Zs]=sphere(50);
surf(Xs,Ys,Zs,'FaceAlpha',0.06,'EdgeAlpha',0.06); hold on;
plot3(mxTrace,myTrace,mzTrace,'LineWidth',1.3);
plot3(mxTrace(1),myTrace(1),mzTrace(1),'o','MarkerSize',8,'LineWidth',1.5);
plot3(mxTrace(end),myTrace(end),mzTrace(end),'s','MarkerSize',8,'LineWidth',1.5);
xlabel('m_x'); ylabel('m_y'); zlabel('m_z'); axis equal;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); grid on;
title(sprintf('Trajectory, Sweep %d, I_c=%.3f mA', ...
    representative_sweep,Ic_trace_actual*1e3));

%% Representative trajectory: four explicit LLG terms
% Each stored term is a vector contribution to dm/dt for ensemble member 1.
% Plot its Euclidean magnitude so all four contributions can be compared
% without cancellation between x, y, and z components.
termPrecessionMag = vecnorm(trace.term_precession,2,1);
termFieldDampingMag = vecnorm(trace.term_field_damping,2,1);
termIsDampingLikeMag = vecnorm(trace.term_Is_damping_like,2,1);
termIsAlphaMag = vecnorm(trace.term_Is_alpha,2,1);

figLLGTerms=figure('Color','w','Name', ...
    'Figure 9 - Four LLG terms','NumberTitle','off');
clf(figLLGTerms);
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(time*1e9,termPrecessionMag,'LineWidth',1.2);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
ylabel('|T_{prec}| [s^{-1}]'); grid on;
title('Field-precession term');

nexttile;
plot(time*1e9,termFieldDampingMag,'LineWidth',1.2);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
ylabel('|T_{damp}| [s^{-1}]'); grid on;
title('Gilbert field-damping term');

nexttile;
plot(time*1e9,termIsDampingLikeMag,'LineWidth',1.2);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
ylabel('|T_{Is,DL}| [s^{-1}]'); grid on;
title('Damping-like spin-current term');

nexttile;
plot(time*1e9,termIsAlphaMag,'LineWidth',1.2);
xline(p.tIs_on*1e9,'--'); xline(p.tIs_off*1e9,'--');
xlabel('Time [ns]'); ylabel('|T_{Is,\alpha}| [s^{-1}]'); grid on;
title('\alpha spin-current correction term');

save(fullfile(script_dir,'INI_zonly_parameter_variation_results_0.18A.mat'), ...
    'Pplus_all','Pplus_initial','Ic_values','HOe_all','Is_pulse', ...
    'parameter_sets','representative_sweep','initial_theta_deg', ...
    'initial_phi_deg','trace','variation_flags','time','dt', ...
    'Ic_trace_actual','Is_pulse_scale','Is_pulse_override');

%% Local simulation functions
function [Pplus,HOe_values,trace]=run_probability_curve(p,Ic_values,Is_pulse, ...
    Nens,seed,dt,Nstep,kTrace,store_trace,theta0_deg,phi0_deg, ...
    kB,mu0,qe,muB,gamma0,threshold)

Hk_eff=p.Hk_Oe*1000/(4*pi);     % [A/m], target effective PMA field
u=[0;0;1];
V=pi*(p.diameter_x/2)*(p.diameter_y/2)*p.tF;
Ns=p.Ms*V/muB;
sigma_H=sqrt(2*p.alpha*kB*p.T/(gamma0*mu0*p.Ms*V*dt));
[Nx,Ny,Nz]=ellipsoid_demag_factors( ...
    p.diameter_x,p.diameter_y,p.tF);
Ndemag=[Nx;Ny;Nz];
Ku=0.5*mu0*p.Ms*(Hk_eff+p.Ms*(Nz-Nx));

Is_direction=p.Is_direction_ideal(:)+p.Is_offset(:);
if norm(Is_direction)<1e-12
    error('Misaligned spin-current direction has zero magnitude.');
end
Is_direction=Is_direction/norm(Is_direction);

Oe_direction=[0;0;1]+p.Oe_offset(:);
if norm(Oe_direction)<1e-12
    error('Misaligned Oersted-field direction has zero magnitude.');
end
Oe_direction=Oe_direction/norm(Oe_direction);

theta0=deg2rad(theta0_deg); phi0=deg2rad(phi0_deg);
m0=[sin(theta0)*cos(phi0);sin(theta0)*sin(phi0);cos(theta0)];
m0=m0/norm(m0);

Nic=numel(Ic_values); Pplus=zeros(1,Nic); HOe_values=zeros(1,Nic);
trace.m=zeros(3,Nstep+1); trace.Is=zeros(3,Nstep+1);
trace.term_precession=zeros(3,Nstep+1);
trace.term_field_damping=zeros(3,Nstep+1);
trace.term_Is_damping_like=zeros(3,Nstep+1);
trace.term_Is_alpha=zeros(3,Nstep+1);
trace.stored=false;

for kIc=1:Nic
    Ic=Ic_values(kIc);
    H_Oe=p.Nturns*Ic*p.R_ring^2/ ...
        (2*(p.R_ring^2+p.z_offset^2)^(3/2));
    HOe_values(kIc)=H_Oe;
    Happ_constant=H_Oe*Oe_direction;
    rng(seed+kIc,'twister');
    m=repmat(m0,1,Nens);
    do_trace=store_trace && kIc==kTrace;
    if do_trace, trace.m(:,1)=m(:,1); end

    for n=1:Nstep
        tn=(n-1)*dt; tnp1=n*dt;
        Is_n=Is_pulse*double(tn>=p.tIs_on && tn<=p.tIs_off)*Is_direction;
        Is_np1=Is_pulse*double(tnp1>=p.tIs_on && tnp1<=p.tIs_off)*Is_direction;

        % Apply the swept Ic field throughout the trajectory: before,
        % during, and after hard-axis initialization by the Is pulse.
        % Both stochastic Heun stages use the same applied field.
        Happ_n=Happ_constant;
        Happ_np1=Happ_constant;

        Hth=sigma_H*randn(3,Nens);
        if do_trace
            [k1,Tprec_n,Tdamp_n,TIsDL_n,TIsAlpha_n]= ...
                llg_spin_current_rhs(m,Hth,p.Ms,Ku,u,Happ_n,Ndemag, ...
                p.alpha,gamma0,mu0,qe,Ns,Is_n);
        else
            k1=llg_spin_current_rhs(m,Hth,p.Ms,Ku,u,Happ_n,Ndemag, ...
                p.alpha,gamma0,mu0,qe,Ns,Is_n);
        end
        mPredict=normalize_columns(m+dt*k1);
        if do_trace
            [k2,Tprec_np1,Tdamp_np1,TIsDL_np1,TIsAlpha_np1]= ...
                llg_spin_current_rhs(mPredict,Hth,p.Ms,Ku,u,Happ_np1,Ndemag, ...
                p.alpha,gamma0,mu0,qe,Ns,Is_np1);
        else
            k2=llg_spin_current_rhs(mPredict,Hth,p.Ms,Ku,u,Happ_np1,Ndemag, ...
                p.alpha,gamma0,mu0,qe,Ns,Is_np1);
        end
        m=normalize_columns(m+0.5*dt*(k1+k2));
        if do_trace
            trace.m(:,n+1)=m(:,1);
            trace.Is(:,n+1)=Is_np1;
            % Store predictor terms at t_n. The corrector terms from the
            % last step provide the final sample at t_{Nstep+1}.
            trace.term_precession(:,n)=Tprec_n(:,1);
            trace.term_field_damping(:,n)=Tdamp_n(:,1);
            trace.term_Is_damping_like(:,n)=TIsDL_n(:,1);
            trace.term_Is_alpha(:,n)=TIsAlpha_n(:,1);
            if n==Nstep
                trace.term_precession(:,n+1)=Tprec_np1(:,1);
                trace.term_field_damping(:,n+1)=Tdamp_np1(:,1);
                trace.term_Is_damping_like(:,n+1)=TIsDL_np1(:,1);
                trace.term_Is_alpha(:,n+1)=TIsAlpha_np1(:,1);
            end
        end
    end
    Pplus(kIc)=mean(m(3,:)>threshold);
    fprintf(['Ic=%+7.3f mA | constant mu0 H_Oe=%+8.3f mT | ' ...
        'P(+z)=%.3f\n'], ...
        Ic*1e3,mu0*H_Oe*1e3,Pplus(kIc));
    if do_trace, trace.stored=true; end
end
if store_trace && ~trace.stored
    error('Representative trajectory was not recorded.');
end
end

function [dmdt,Tprecession,TfieldDamping,TIsDampingLike,TIsAlpha]= ...
    llg_spin_current_rhs(m,Hth,Ms,Ku,u,Happ,Ndemag, ...
    alpha,gamma0,mu0,qe,Ns,Is)
mdotu=sum(m.*u,1);
Hani=(2*Ku/(mu0*Ms))*u*mdotu;
Hdemag=-Ms*(Ndemag.*m);
Heff=Hani+Hdemag+Happ+Hth;
mxH=cross(m,Heff,1);
mxmxH=cross(m,mxH,1);
Is_matrix=repmat(Is,1,size(m,2));
mxIs=cross(m,Is_matrix,1);
mxmxIs=cross(m,mxIs,1);

% Four explicit vector contributions to dm/dt.
Tprecession = -gamma0*mxH/(1+alpha^2);
TfieldDamping = -gamma0*alpha*mxmxH/(1+alpha^2);
TIsDampingLike = -(1/(qe*Ns))*mxmxIs/(1+alpha^2);
TIsAlpha = (alpha/(qe*Ns))*mxIs/(1+alpha^2);
dmdt = Tprecession+TfieldDamping+TIsDampingLike+TIsAlpha;
end

function m=normalize_columns(m)
n=sqrt(sum(m.^2,1)); n(n<eps)=1; m=m./n;
end

function [Nx,Ny,Nz] = ellipsoid_demag_factors( ...
    diameter_x,diameter_y,thickness)
% Recompute shape factors from the same sampled geometry used for V.
semi_axes = [diameter_x diameter_y thickness]/2;
if any(~isfinite(semi_axes)) || any(semi_axes<=0)
    error('All ellipsoid dimensions must be positive and finite.');
end
scaled_axes = semi_axes/min(semi_axes);
axis_squared = scaled_axes.^2;
prefactor = prod(scaled_axes)/2;
x_min = log(min(axis_squared))-35;
x_max = log(max(axis_squared))+35;
N = zeros(1,3);
for axis_index=1:3
    integrand = @(x) demag_log_integrand( ...
        x,axis_squared,axis_index,prefactor);
    N(axis_index) = integral(integrand,x_min,x_max, ...
        'RelTol',1e-9,'AbsTol',1e-12);
end
N=max(N,0); N=N/sum(N);
Nx=N(1); Ny=N(2); Nz=N(3);
end

function value = demag_log_integrand(x,axis_squared,axis_index,prefactor)
s=exp(x);
delta=sqrt((s+axis_squared(1)).*(s+axis_squared(2)).* ...
    (s+axis_squared(3)));
value=prefactor.*s./((s+axis_squared(axis_index)).*delta);
end

function plot_parameter_boxplots(parameter_sets,baseline)
% Only independently sampled primitive values are shown. V, Ns, and
% Nx/Ny/Nz are derived, while pulse timing, Nturns, and z_offset are fixed.
names={'Ms','alpha','T','diameter x','diameter y','tF','Hk','R ring'};
fields={'Ms','alpha','T','diameter_x','diameter_y','tF', ...
    'Hk_Oe','R_ring'};
ns=numel(parameter_sets); X=zeros(ns,numel(fields));
for j=1:numel(fields)
    values=arrayfun(@(p) p.(fields{j}),parameter_sets).';
    X(:,j)=100*(values-baseline.(fields{j}))/abs(baseline.(fields{j}));
end
O=zeros(ns,6);
for s=1:ns
    O(s,:)=[parameter_sets(s).Is_offset parameter_sets(s).Oe_offset];
end
figBox=figure('Color','w','Name','Device parameter variation','NumberTitle','off');
clf(figBox); tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; boxplot(X,'Labels',names,'Whisker',1.5);
ylabel('Variation from mean [%]');
title('Parameter variation (box = 25th to 75th percentile)');
set(gca,'FontSize',10,'LineWidth',1.2,'TickLabelInterpreter','none');
xtickangle(35); box on; grid off;
nexttile; boxplot(O,'Labels',{'I_s offset x','I_s offset y','I_s offset z', ...
    'Oe offset x','Oe offset y','Oe offset z'},'Whisker',1.5);
ylabel('Direction offset'); title('Direction-offset variation');
set(gca,'FontSize',11,'LineWidth',1.2,'TickLabelInterpreter','none');
box on; grid off;
end
