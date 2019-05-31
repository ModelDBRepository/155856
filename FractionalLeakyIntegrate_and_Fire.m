%% Fractional Leaky Integrate-and-Fire Model
% Model files for the paper
% Wondimu Teka,  Toma M.Marinov, Fidel Santamaria (2014) Neuronal Spike Timing Adaptation Described with a Fractional Leaky Integrate-and-Fire Model. 
% PLoS Comput Biol 10(3): e1003526. doi:10.1371/journal.pcbi.1003526

clear
 
rseed=2342;

Ncells=1; %number of cells,  do not change this value, this code is not for network, it is only for one cell. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   parameters  values that can be varied. 

dt=0.1; %ms  time step. This value is used for  all the data. 
t=0:dt:1000; %ms    mostly is  up to  5000 ms (5sec)  or  10000 ms (10 sec) is used. 
                


 % alpha= ??; %  the fractional exponent (order)  it is  varied  below using a for loop

 
Rm=40; %Mohms
taum=20*ones([1 Ncells]); %ms membrane time constant
Cm=0.5; %nF  % this is not used since taum and Rm are used and note taum = Cm*Rm

refrac=8*ones([1 Ncells]); %ms absolute refarctory period 

v0= -70*ones([1 Ncells]);  % mV  initial value 
% initial voltage was varied: used values are v0= -70, v0= -75 v0= -58, v0= -45, 

vrest= -70;     % mV  the resting  potential which is used for reset too. 
vth= vrest + 20*ones([1 Ncells]); %  is  -50  mV,   threshold value of membrane voltage
 vpeak = vrest + 100*ones([1 Ncells]); % is  30 mV,   peak value of he voltage at spike


Iinjamplitude= 3;  % nA  injected current amplitude, see Iinj  for the whole injected current 

%  for the subthershold  use Iinjamplitude= 0.3; 


Namp=0; %noise amplitude, to add white Gaussian noise with this  standard devation Namp, which is noise amplitude

 %  This value should be changed only for figure 11, 
% Namp=1  %  For only figure 11,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  used for most of the simulations
Iinj=Iinjamplitude*ones(length(t),Ncells)+ Namp*(randn(length(t),Ncells));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the sinusoidal simulation with  non constant freqency, use the following  ZAP current
%  f=0.1*(2./(1+exp(-t./1500))-1).^3;  % frequency  in cycle  per ms. since time is in ms; the freqency increases from 0  Hz to  100 Hz
%  Iinj=0.3*sin((2*pi).*f.*t)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %  for figure 10 use the following current instead of the above
% for the sinusoidal simulation with constant freqency, use the following
%f = 0.003  %  this is in  cycle per ms, it is the same is     3 Hz 
%Iinj= 2 +2*sin((0:dt:t(end)).*2*pi*f)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names for parameters to pass them to the  main  function: runNetworkderivative

NetProp.Ncells=Ncells;
NetProp.Rm=Rm; %Mohms
NetProp.TauM=taum;
NetProp.Cm=Cm; %nF  % this is not used since taum and Rm are used.

NetProp.Refrac=refrac;
NetProp.vTh=vth;
NetProp.v0=v0;
NetProp.vrest=vrest;
NetProp.vpeak=vpeak;

NetProp.Noise=Namp;

NetProp.dt=dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Main integrator
%  Integrating  output for different alpha values using  the fractional derivative

b=1;

 for alpha=1:-0.2:0.2
out(b)=runNetworkderivative(NetProp,Iinj,t,rseed,alpha); % main output
                % Note:This has the following output out.v=v  for voltage; out.sp=sp for spike ; out.Memo2=Memo2 for memery; out.t=t for time;

alphavalue(b)=alpha

  ISI{b}=diff(out(b).t(logical(out(b).sp(1:end-1))));  %  calculates  ISI in msec
  firingrate{b} = 1000./ISI{b};  %in  Hz, calculates firing rate
  totalspike(b) = sum(out(b).sp); % for toltal number of spikes for  total time
  b= b+1
 end
 
  cv=['kbgmr'];
  c=1;
 
  for b=1:length(out)  % for  alpha =1.0
  subplot(5,1,c)

      trange=1:1000/dt;
 
  plot(out(b).t(trange), out(b).v(trange), cv(b), 'Linewidth',1.5)  %


set(gca, 'FontSize', 14, 'FontName', 'Helvetica')

 set(gca, 'LooseInset', get(gca,'TightInset'))
 
 
ylabel('V (mV)', 'fontsize',14, 'FontName', 'Helvetica');
  legendinfo2= sprintf('\\alpha= %0.1f\n', alphavalue(b));
  c=c+1;
   axis([-1  1000  -75   35])
   
   legend(legendinfo2,'Location','SouthWest', 'fontsize',11, 'FontName', 'Helvetica')
   clear legendinfo2;
end

  xlabel('Time (ms)', 'fontsize',14, 'FontName', 'Helvetica');

  
  
  
  
