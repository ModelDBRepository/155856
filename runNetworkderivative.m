
%% This is a  funcion for the  Fractional Leaky Integrate-and-Fire Model. It integrates the fractional derivative and  the  voltage v at each time t.

function out=runNetworkderivative(NetProp,Iinj,t,rseed,alpha)
rand('seed',rseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells=NetProp.Ncells;
taum=NetProp.TauM;
refrac=NetProp.Refrac/NetProp.dt;
vth=NetProp.vTh;
vpeak=NetProp.vpeak;
v0=NetProp.v0;
vrest=NetProp.vrest;
Namp=NetProp.Noise;
rm=NetProp.Rm;
dt=NetProp.dt;


v=vrest.*ones(length(t),Ncells);
v(1,:) = v0(1,:);
vrest=vrest.*ones(length(t),Ncells);

sp=zeros(length(t),Ncells);

isstillrefrac=zeros(1,Ncells);
rfcounter=zeros(1,Ncells);


Inhib = zeros(length(t),Ncells);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The weight for the  voltage memory trace of the fractional drivative for
% calculated here for the total time t for faster simulation
NN=length(t);
nn=1:NN-1;
WCoet=(NN+1-nn).^(1-alpha)-(NN-nn).^(1-alpha);


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for a=1:length(t)-1
     
       Ihere=Iinj(a,:); 
    


 if a<2   %for the first  step, it is a normal neuron
        vdummy =dt*(vrest(a,:)-v(a,:)+rm*Ihere)./taum + v(a,:);
 else
        
    
        %%%%% Fractional derivative starts  here 
        d2dM=v(1:a,:); % to call all past values of voltage  for fractioanl integration

          WCoe=WCoet(end-a+2:end);  % The weight for the  voltage memory trace of the fractional drivative  at each  tiime t 
        TeDi=d2dM(2:a,:)-d2dM(1:a-1,:); % Delta V (using all past values of V)  of  the  voltage memory trace of the fractional drivative  at each  tiime t 
        
         fraccalcu=WCoe*TeDi-d2dM(a,:);  %  The fraction derivative 
         
         kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term
        %%%%% Fractional derivative ends   here 
        
        
        
         % used for  the fractioanl  neuron, value of the voltage if it is  on   subthreshold.
        vdummy = kr*(vrest(a,:)-v(a,:)+rm*Ihere)./taum - fraccalcu; 


        Memo2(a,:)=fraccalcu+v(a,:); % to save memory over time
        
        
        if Memo2(a,:) > 0 &&  Memo2(a-1,:) < 0
            Inhib(a-1,:) = Memo2(a-1,:);   % to save the value of mmemory which has a positive feedback
        end

 end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    test=((vdummy>vth)+(isstillrefrac==1))>0 ;%v larger than vth or inside refrac period
    sp(a+1,:)=test.*(rfcounter==0); % a spike if outside refrac period and test is possitive
    
    v(a+1,:)=vrest(a,:).*(isstillrefrac==1) + (~test).*vdummy + sp(a+1,:)*vpeak; % save  V value
    isstillrefrac=sp(a+1,:)+test.*isstillrefrac.*(rfcounter<refrac); %initiallize refrac flag
    rfcounter=rfcounter+test.*(isstillrefrac==1); %you enter or continue rfcounter
    rfcounter=(~(test.*(rfcounter==refrac))).*rfcounter; %exit rfcounter
    isstillrefrac=isstillrefrac.*(rfcounter>0); %resetting the refrac flag
end
out.v=v;
out.sp=sp;
out.Memo2=Memo2;
out.t=t;
out.Inhib = Inhib;
