clear
clc

fprintf("FCC model for the atomic nucleus\n");

answer=inputdlg({'Insert the number of protons Z of the nucleus you want to analyze:'});
numProtons=str2double(answer(1,1));

answer=inputdlg({'Insert the number of neutrons N of the nucleus you want to analyze:'});
numNeutrons=str2double(answer(1,1));

%Standard FCC build-up%			 
fm=2.026/(2*sqrt(2));
kP=1;
kN=1;

for n=0:7                   %No. of shells = 8 for calculations on nuclei up to 240 protons and 240 neutrons 
    for l=n:-1:0            %L-subshells, higher Ls first: L = N,..,2,1,0  
        j=(l+1/2);          %J is tot. ang. mom. value, J = |L+spin| 
        for m=-j:j          %M = -J,.,-3/2,-1/2,1/2,3/2,.,J
            for i=-1:2:1    %isospin: neutrons = -1, protons = 1.
                x=2*abs(m)*power(-1,m-1/2); 
                y=abs(2*j+1-2*abs(m))*power(-1,n-m+1/2+i/2-abs(n-j+1));
                z=2*abs(n-j+1)*power(-1,abs(n-j+1)-i/2);
                if i==1
                    nucleon(kP+240,1)=x;
                    nucleon(kP+240,2)=y;
                    nucleon(kP+240,3)=z;
                    if mod(x+1,4)==0
                        nucleon_s(kP+240,1)=-1/2;  
                    else
                        nucleon_s(kP+240,1)=1/2;
                    end
                    nucleon_j(kP+240,1)=(abs(x)+abs(y)-1)/2;
                    nucleon_n(kP+240,1)=(abs(x)+abs(y)+abs(z)-3)/2;
                    nucleon_m(kP+240,1)=0.5*abs(x)*power(-1,0.5*(x-1));
                    kP=kP+1;
                else
                    nucleon(kN,1)=x;
                    nucleon(kN,2)=y;
                    nucleon(kN,3)=z;	
                    if mod(x+1,4)==0
                        nucleon_s(kN,1)=-1/2;  
                    else
                        nucleon_s(kN,1)=1/2;
                    end
                    nucleon_j(kN,1)=(abs(x)+abs(y)-1)/2;
                    nucleon_n(kN,1)=(abs(x)+abs(y)+abs(z)-3)/2;
                    nucleon_m(kN,1)=0.5*abs(x)*power(-1,0.5*(x-1));
                    kN=kN+1;
                end
            end
        end
    end
end

for i=1:240
    LofN(i,1)=i;        %Make changes in the List of Nucleons LofN[] if you want to change the build-up sequence
    LofN(i+240,1)=i+240;
end

%Calculation of Total Coulomb repulsion term
fccQ=0;   
if numProtons>=2 
    for j=1:numProtons-1
        for k=j+1:numProtons                
            fccQ=fccQ+1.439965/sqrt((fm*nucleon(LofN(j+240),1)-fm*nucleon(LofN(k+240),1))^2+(fm*nucleon(LofN(j+240),2)-fm*nucleon(LofN(k+240),2))^2+(fm*nucleon(LofN(j+240),3)-fm*nucleon(LofN(k+240),3))^2);
        end
    end
end

%Calculation of nuclear spin
pspin=0;
nspin=0;
if numProtons>0
    for k=1:numProtons
        pspin=sign(nucleon_s(LofN(k+240)))*nucleon_j(LofN(k+240))+pspin;
    end
end
if numNeutrons>0
    for k=1:numNeutrons
        nspin=sign(nucleon_s(LofN(k)))*nucleon_j(LofN(k))+nspin;
    end
end
fccSpin=abs(pspin+nspin);

%Claculation of number of bonds used for BE calculation considering the number of neighbors per nucleon (1-12)
c1=1.4;c2=2.8;c3=7.2516;c4=5.6;c5=7.2;c6=8.4;c7=9.8;c8=10.4;c9=12.6;c10=14;c11=15.4;c12=16.8;
NumberNucleons=zeros(13,2);
for i=1:numProtons
    first_neighbors=0;
    for j=1:numProtons
        dist=sqrt((nucleon(i+240,1)-nucleon(j+240,1))^2+(nucleon(i+240,2)-nucleon(j+240,2))^2+(nucleon(i+240,3)-nucleon(j+240,3))^2);
        if dist==2*sqrt(2)
           first_neighbors=first_neighbors+1;
        end
    end
    for j=1:numNeutrons
         dist=sqrt((nucleon(i+240,1)-nucleon(j,1))^2+(nucleon(i+240,2)-nucleon(j,2))^2+(nucleon(i+240,3)-nucleon(j,3))^2);
        if dist==2*sqrt(2)
           first_neighbors=first_neighbors+1;
        end
    end
    NumberNucleons(first_neighbors+1,2)=NumberNucleons(first_neighbors+1,2)+1;
end
for i=1:numNeutrons
    first_neighbors=0;
    for j=1:numProtons
        dist=sqrt((nucleon(i,1)-nucleon(j+240,1))^2+(nucleon(i,2)-nucleon(j+240,2))^2+(nucleon(i,3)-nucleon(j+240,3))^2);
        if dist==2*sqrt(2)
           first_neighbors=first_neighbors+1;
        end
    end
    for j=1:numNeutrons
         dist=sqrt((nucleon(i,1)-nucleon(j,1))^2+(nucleon(i,2)-nucleon(j,2))^2+(nucleon(i,3)-nucleon(j,3))^2);
        if dist==2*sqrt(2)
           first_neighbors=first_neighbors+1;
        end
    end
    NumberNucleons(first_neighbors+1,1)=NumberNucleons(first_neighbors+1,1)+1;
end
fccBE=(NumberNucleons(2,1)+NumberNucleons(2,2))*c1+(NumberNucleons(3,1)+NumberNucleons(3,2))*c2+(NumberNucleons(4,1)+NumberNucleons(4,2))*c3+(NumberNucleons(5,1)+NumberNucleons(5,2))*c4+(NumberNucleons(6,1)+NumberNucleons(6,2))*c5+(NumberNucleons(7,1)+NumberNucleons(7,2))*c6+(NumberNucleons(8,1)+NumberNucleons(8,2))*c7+(NumberNucleons(9,1)+NumberNucleons(9,2))*c8+(NumberNucleons(10,1)+NumberNucleons(10,2))*c9+(NumberNucleons(11,1)+NumberNucleons(11,2))*c10+(NumberNucleons(12,1)+NumberNucleons(12,2))*c11+(NumberNucleons(13,1)+NumberNucleons(13,2))*c12-fccQ;

%Draw the default FCC structure, with protons as yellow spheres and neutrons as blue spheres
figure
ball_size=1000;
bond_width=1;
for i=1:numProtons
    scatter3(fm*nucleon(i+240,1),fm*nucleon(i+240,2),fm*nucleon(i+240,3),ball_size,'yellow','filled')
    hold on
end
for i=1:numNeutrons
    scatter3(fm*nucleon(i,1),fm*nucleon(i,2),fm*nucleon(i,3),ball_size,'blue','filled')
    hold on
end

%Number of first, second and third neighbor bonds and BE
pp_bonds=zeros(3,1);
nn_bonds=zeros(3,1);
pn_bonds=zeros(3,1);
for i=1:numProtons-1
   for j=i+1:numProtons
       dist=sqrt((nucleon(i+240,1)-nucleon(j+240,1))^2+(nucleon(i+240,2)-nucleon(j+240,2))^2+(nucleon(i+240,3)-nucleon(j+240,3))^2);
       if dist==2*sqrt(2)
           pp_bonds(1,1)=pp_bonds(1,1)+1;
           plot3([fm*nucleon(i+240,1),fm*nucleon(j+240,1)],[fm*nucleon(i+240,2),fm*nucleon(j+240,2)],[fm*nucleon(i+240,3),fm*nucleon(j+240,3)],'k','LineWidth',bond_width)
           hold on
       elseif dist==4
           pp_bonds(2,1)=pp_bonds(2,1)+1;
       elseif dist==2*sqrt(6)
           pp_bonds(3,1)=pp_bonds(3,1)+1;
       end
   end
end
for i=1:numNeutrons-1
   for j=i+1:numNeutrons
       dist=sqrt((nucleon(i,1)-nucleon(j,1))^2+(nucleon(i,2)-nucleon(j,2))^2+(nucleon(i,3)-nucleon(j,3))^2);
       if dist==2*sqrt(2)
           nn_bonds(1,1)=nn_bonds(1,1)+1;
           plot3([fm*nucleon(i,1),fm*nucleon(j,1)],[fm*nucleon(i,2),fm*nucleon(j,2)],[fm*nucleon(i,3),fm*nucleon(j,3)],'k','LineWidth',bond_width)
           hold on
       elseif dist==4
           nn_bonds(2,1)=nn_bonds(2,1)+1;
       elseif dist==2*sqrt(6)
           nn_bonds(3,1)=nn_bonds(3,1)+1;
       end
   end
end
for i=1:numProtons
   for j=1:numNeutrons
       dist=sqrt((nucleon(i+240,1)-nucleon(j,1))^2+(nucleon(i+240,2)-nucleon(j,2))^2+(nucleon(i+240,3)-nucleon(j,3))^2);
       if dist==2*sqrt(2)
           pn_bonds(1,1)=pn_bonds(1,1)+1;
           plot3([fm*nucleon(i+240,1),fm*nucleon(j,1)],[fm*nucleon(i+240,2),fm*nucleon(j,2)],[fm*nucleon(i+240,3),fm*nucleon(j,3)],'k','LineWidth',bond_width)
           hold on
       elseif dist==4
           pn_bonds(2,1)=pn_bonds(2,1)+1;
       elseif dist==2*sqrt(6)
           pn_bonds(3,1)=pn_bonds(3,1)+1;
       end
   end
end
bonds=pp_bonds+nn_bonds+pn_bonds;
beta=(fccBE+fccQ)/bonds(1,1);
bind_en=beta*bonds(1,1);
title({['FCC lattice model for Z = ',num2str(numProtons),' and N = ',num2str(numNeutrons),' (default build-up)'] ['Total binding energy = ',num2str(fccBE),' MeV'] ['Total Coulomb repulsion energy = ',num2str(fccQ),' MeV']})
axis equal
axis square
xlabel('X [fm]')
ylabel('Y [fm]')
zlabel('Z [fm]')

%Randomization of external nucleons
closure_shells=zeros(8,1);
for i=1:8
   closure_shells(i)=i*(i+1);
   if i>1
      closure_shells(i)=closure_shells(i)+closure_shells(i-1);
   end
end
n_lowerN=-1;
n_lowerP=-1;
for i=1:8
   if numProtons>=closure_shells(i)
      n_lowerP=i-1; 
   end
   if numNeutrons>=closure_shells(i)
      n_lowerN=i-1; 
   end
end

if n_lowerP==-1
   excessP0=numProtons;
   fprintf("In the default build-up the protons are not closing any shell and the excess is %d\n",excessP0);
   n_core_fixedP=-1;
   n_maximum_permutationP=1;
   excessP=excessP0;
else
   excessP0=numProtons-closure_shells(n_lowerP+1);
   fprintf("In the default build-up the protons are closing shell %d and the excess is %d\n",n_lowerP,excessP0);
   
   answer=inputdlg({'Type the n-value of the shell to be kept as fixed interior core for protons:'});
   n_core_fixedP=str2double(answer(1,1));
   answer=inputdlg({'Type the n-value of the maximum shell for the randomization of external protons:'});
   n_maximum_permutationP=str2double(answer(1,1));
   
   excessP=numProtons-closure_shells(n_core_fixedP+1);
end

if n_lowerN==-1
   excessN0=numNeutrons;
   fprintf("In the default build-up the neutrons are not closing any shell and the excess is %d\n",excessN0);
   n_core_fixedN=-1;
   n_maximum_permutationN=1;
   excessN=excessN0;
else
   excessN0=numNeutrons-closure_shells(n_lowerN+1);
   fprintf("In the default build-up the neutrons are closing shell %d and the excess is %d\n",n_lowerN,excessN0);
   
   answer=inputdlg({'Type the n-value of the shell to be kept as fixed interior core for neutrons:'});
   n_core_fixedN=str2double(answer(1,1));
   answer=inputdlg({'Type the n-value of the maximum shell for the randomization of external neutrons:'});
   n_maximum_permutationN=str2double(answer(1,1));
   
   excessN=numNeutrons-closure_shells(n_core_fixedN+1);
end

coreP=numProtons-excessP;
coreN=numNeutrons-excessN;
maximumP=closure_shells(n_maximum_permutationN+1);
maximumN=closure_shells(n_maximum_permutationN+1);
occupable_positionsP=maximumP-coreP;
occupable_positionsN=maximumN-coreN;
number_different_structuresP=excessP*(occupable_positionsP-excessP);
number_different_structuresN=excessN*(occupable_positionsN-excessN);
build_upP=zeros(240,number_different_structuresP);
build_upN=zeros(240,number_different_structuresN);
for i=1:240
    build_upP(i,:)=i;
    build_upN(i,:)=i;
end
for s=1:number_different_structuresP
    change1P_index=coreP+1+fix((s-1)/(occupable_positionsP-excessP));
    change2P_index=coreP+excessP+1+mod(s-1,(occupable_positionsP-excessP));
    build_upP([change1P_index change2P_index],s)=build_upP([change2P_index change1P_index],s);
end
for s=1:number_different_structuresN
    change1N_index=coreN+1+fix((s-1)/(occupable_positionsN-excessN));
    change2N_index=coreN+excessN+1+mod(s-1,(occupable_positionsN-excessN));
    build_upN([change1N_index change2N_index],s)=build_upN([change2N_index change1N_index],s);
end
build_up=zeros(480,number_different_structuresP*number_different_structuresN);
for s=1:number_different_structuresP
    for r=1:number_different_structuresN
       for i=1:240
          build_up(i,(s-1)*number_different_structuresN+r)=build_upN(i,r);
          build_up(i+240,(s-1)*number_different_structuresN+r)=build_upP(i,s)+240;
       end
    end
end
feas_struct=1;
feasible_buildup(:,1)=LofN(:);
for s=1:number_different_structuresP*number_different_structuresN
    
    %Check if the s-th random build-up leads to same spin as the default
    feasible_spin=0;
    pspin_rand=0;
    nspin_rand=0;
    if numProtons>0
        for k=1:numProtons
            pspin_rand=sign(nucleon_s(build_up(k+240,s)))*nucleon_j(build_up(k+240,s))+pspin_rand;
        end
    end
    if numNeutrons>0
        for k=1:numNeutrons
            nspin_rand=sign(nucleon_s(build_up(k,s)))*nucleon_j(build_up(k,s))+nspin_rand;
        end
    end
    fccSpin_rand=abs(pspin_rand+nspin_rand);
    if fccSpin_rand==fccSpin
       feasible_spin=1;
    end
    
    %Check if the s-th random build-up leads to disconnected nucleons
    feasible_bonds=0;
    NumberNucleons=zeros(13,2);
    for i=1:numProtons
        first_neighbors=0;
        for j=1:numProtons
            dist=sqrt((nucleon(build_up(i+240,s),1)-nucleon(build_up(j+240,s),1))^2+(nucleon(build_up(i+240,s),2)-nucleon(build_up(j+240,s),2))^2+(nucleon(build_up(i+240,s),3)-nucleon(build_up(j+240,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        for j=1:numNeutrons
             dist=sqrt((nucleon(build_up(i+240,s),1)-nucleon(build_up(j,s),1))^2+(nucleon(build_up(i+240,s),2)-nucleon(build_up(j,s),2))^2+(nucleon(build_up(i+240,s),3)-nucleon(build_up(j,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        NumberNucleons(first_neighbors+1,2)=NumberNucleons(first_neighbors+1,2)+1;
    end
    for i=1:numNeutrons
        first_neighbors=0;
        for j=1:numProtons
            dist=sqrt((nucleon(build_up(i,s),1)-nucleon(build_up(j+240,s),1))^2+(nucleon(build_up(i,s),2)-nucleon(build_up(j+240,s),2))^2+(nucleon(build_up(i,s),3)-nucleon(build_up(j+240,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        for j=1:numNeutrons
             dist=sqrt((nucleon(build_up(i,s),1)-nucleon(build_up(j,s),1))^2+(nucleon(build_up(i,s),2)-nucleon(build_up(j,s),2))^2+(nucleon(build_up(i,s),3)-nucleon(build_up(j,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        NumberNucleons(first_neighbors+1,1)=NumberNucleons(first_neighbors+1,1)+1;
    end
    if NumberNucleons(1,1)==0 && NumberNucleons(1,2)==0
       feasible_bonds=1;
    end
    
    %Accept or reject the s-th build-up
    if feasible_spin==1 && feasible_bonds==1
       feasible_buildup(:,feas_struct+1)=build_up(:,s);
       feas_struct=feas_struct+1;
    end

    progress = (s / (number_different_structuresP*number_different_structuresN)) * 100;
    fprintf('Computing all feasible FCC structures - Fesible structures so far: %d, Progress: %.2f%%\n', feas_struct, progress);
end

feas_coulomb=zeros(1,feas_struct);
feas_BE=zeros(1,feas_struct);
for s=1:feas_struct   
    if numProtons>=2 
        for j=1:numProtons-1
            for k=j+1:numProtons                
                feas_coulomb(1,s)=feas_coulomb(1,s)+1.439965/sqrt((fm*nucleon(feasible_buildup(j+240,s),1)-fm*nucleon(feasible_buildup(k+240,s),1))^2+(fm*nucleon(feasible_buildup(j+240,s),2)-fm*nucleon(feasible_buildup(k+240,s),2))^2+(fm*nucleon(feasible_buildup(j+240,s),3)-fm*nucleon(feasible_buildup(k+240,s),3))^2);
            end
        end
    end
    NumberNucleons=zeros(13,2);
    for i=1:numProtons
        first_neighbors=0;
        for j=1:numProtons
            dist=sqrt((nucleon(feasible_buildup(i+240,s),1)-nucleon(feasible_buildup(j+240,s),1))^2+(nucleon(feasible_buildup(i+240,s),2)-nucleon(feasible_buildup(j+240,s),2))^2+(nucleon(feasible_buildup(i+240,s),3)-nucleon(feasible_buildup(j+240,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        for j=1:numNeutrons
            dist=sqrt((nucleon(feasible_buildup(i+240,s),1)-nucleon(feasible_buildup(j,s),1))^2+(nucleon(feasible_buildup(i+240,s),2)-nucleon(feasible_buildup(j,s),2))^2+(nucleon(feasible_buildup(i+240,s),3)-nucleon(feasible_buildup(j,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        NumberNucleons(first_neighbors+1,2)=NumberNucleons(first_neighbors+1,2)+1;
    end
    for i=1:numNeutrons
        first_neighbors=0;
        for j=1:numProtons
            dist=sqrt((nucleon(feasible_buildup(i,s),1)-nucleon(feasible_buildup(j+240,s),1))^2+(nucleon(feasible_buildup(i,s),2)-nucleon(feasible_buildup(j+240,s),2))^2+(nucleon(feasible_buildup(i,s),3)-nucleon(feasible_buildup(j+240,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        for j=1:numNeutrons
             dist=sqrt((nucleon(feasible_buildup(i,s),1)-nucleon(feasible_buildup(j,s),1))^2+(nucleon(feasible_buildup(i,s),2)-nucleon(feasible_buildup(j,s),2))^2+(nucleon(feasible_buildup(i,s),3)-nucleon(feasible_buildup(j,s),3))^2);
            if dist==2*sqrt(2)
               first_neighbors=first_neighbors+1;
            end
        end
        NumberNucleons(first_neighbors+1,1)=NumberNucleons(first_neighbors+1,1)+1;
    end
    feas_BE(1,s)=(NumberNucleons(2,1)+NumberNucleons(2,2))*c1+(NumberNucleons(3,1)+NumberNucleons(3,2))*c2+(NumberNucleons(4,1)+NumberNucleons(4,2))*c3+(NumberNucleons(5,1)+NumberNucleons(5,2))*c4+(NumberNucleons(6,1)+NumberNucleons(6,2))*c5+(NumberNucleons(7,1)+NumberNucleons(7,2))*c6+(NumberNucleons(8,1)+NumberNucleons(8,2))*c7+(NumberNucleons(9,1)+NumberNucleons(9,2))*c8+(NumberNucleons(10,1)+NumberNucleons(10,2))*c9+(NumberNucleons(11,1)+NumberNucleons(11,2))*c10+(NumberNucleons(12,1)+NumberNucleons(12,2))*c11+(NumberNucleons(13,1)+NumberNucleons(13,2))*c12;
end
feas_net_BE=feas_BE-feas_coulomb;

figure
histogram(feas_net_BE)
title(['Distribution of binding energies for ',num2str(feas_struct),' different FCC structures'])
xlabel('Net binding energy [MeV]')
ylabel('Occurrence [-]')
[max_BE,opt_struct]=max(feas_net_BE);

%Draw structure with maximum BE
figure
for i=1:numProtons
    scatter3(fm*nucleon(feasible_buildup(i+240,opt_struct),1),fm*nucleon(feasible_buildup(i+240,opt_struct),2),fm*nucleon(feasible_buildup(i+240,opt_struct),3),ball_size,'yellow','filled')
    hold on
end
for i=1:numNeutrons
    scatter3(fm*nucleon(feasible_buildup(i,opt_struct),1),fm*nucleon(feasible_buildup(i,opt_struct),2),fm*nucleon(feasible_buildup(i,opt_struct),3),ball_size,'blue','filled')
    hold on
end
for i=1:numProtons-1
   for j=i+1:numProtons
       dist=sqrt((nucleon(feasible_buildup(i+240,opt_struct),1)-nucleon(feasible_buildup(j+240,opt_struct),1))^2+(nucleon(feasible_buildup(i+240,opt_struct),2)-nucleon(feasible_buildup(j+240,opt_struct),2))^2+(nucleon(feasible_buildup(i+240,opt_struct),3)-nucleon(feasible_buildup(j+240,opt_struct),3))^2);
       if dist==2*sqrt(2)
           plot3([fm*nucleon(feasible_buildup(i+240,opt_struct),1),fm*nucleon(feasible_buildup(j+240,opt_struct),1)],[fm*nucleon(feasible_buildup(i+240,opt_struct),2),fm*nucleon(feasible_buildup(j+240,opt_struct),2)],[fm*nucleon(feasible_buildup(i+240,opt_struct),3),fm*nucleon(feasible_buildup(j+240,opt_struct),3)],'k','LineWidth',bond_width)
           hold on
       end
   end
end
for i=1:numNeutrons-1
   for j=i+1:numNeutrons
       dist=sqrt((nucleon(feasible_buildup(i,opt_struct),1)-nucleon(feasible_buildup(j,opt_struct),1))^2+(nucleon(feasible_buildup(i,opt_struct),2)-nucleon(feasible_buildup(j,opt_struct),2))^2+(nucleon(feasible_buildup(i,opt_struct),3)-nucleon(feasible_buildup(j,opt_struct),3))^2);
       if dist==2*sqrt(2)
           plot3([fm*nucleon(feasible_buildup(i,opt_struct),1),fm*nucleon(feasible_buildup(j,opt_struct),1)],[fm*nucleon(feasible_buildup(i,opt_struct),2),fm*nucleon(feasible_buildup(j,opt_struct),2)],[fm*nucleon(feasible_buildup(i,opt_struct),3),fm*nucleon(feasible_buildup(j,opt_struct),3)],'k','LineWidth',bond_width)
           hold on
       end
   end
end
for i=1:numProtons
   for j=1:numNeutrons
       dist=sqrt((nucleon(feasible_buildup(i+240,opt_struct),1)-nucleon(feasible_buildup(j,opt_struct),1))^2+(nucleon(feasible_buildup(i+240,opt_struct),2)-nucleon(feasible_buildup(j,opt_struct),2))^2+(nucleon(feasible_buildup(i+240,opt_struct),3)-nucleon(feasible_buildup(j,opt_struct),3))^2);
       if dist==2*sqrt(2)
           plot3([fm*nucleon(feasible_buildup(i+240,opt_struct),1),fm*nucleon(feasible_buildup(j,opt_struct),1)],[fm*nucleon(feasible_buildup(i+240,opt_struct),2),fm*nucleon(feasible_buildup(j,opt_struct),2)],[fm*nucleon(feasible_buildup(i+240,opt_struct),3),fm*nucleon(feasible_buildup(j,opt_struct),3)],'k','LineWidth',bond_width)
           hold on
       end
   end
end
title(['Optimal FCC lattice (highest BE) model for Z = ',num2str(numProtons),' and N = ',num2str(numNeutrons),' - BE = ',num2str(max_BE),' MeV'])
axis equal
axis square
xlabel('X [fm]')
ylabel('Y [fm]')
zlabel('Z [fm]')

answer=inputdlg({'Do you want to analyze the fission of the nucleus? (Y/N):'});
if strcmp(answer,'Y')

    %Fission for all FCC structures and collection of fission fragments
    n_fission_planes=21;
    beta_fission=beta;
    fragments_Z_N_A_coul_BE=zeros(n_fission_planes*2*feas_struct,5);
    
    %For very small nuclei you might need to reduce the number of crystallographic plans, or remove a-posteriori cases where two fragments have not been generated
    coeff_fiss_plane=[1 0 0 0;1 0 0 2;1 0 0 -2;0 1 0 0;0 1 0 2;0 1 0 -2;0 0 1 0;0 0 1 2;0 0 1 -2;1 1 1 0;1 1 1 2;1 1 1 -2;-1 1 1 0;-1 1 1 2;-1 1 1 -2;1 -1 1 0;1 -1 1 2;1 -1 1 -2;1 1 -1 0;1 1 -1 2;1 1 -1 -2];
    
    for s=1:feas_struct
        for f=1:n_fission_planes
            broken_bonds=0;
            closer_f1=0;
            closer_f2=0;
            
            pf1=0; %Number of protons in fragment 1
            pf2=0; %Number of protons in fragment 2
            nf1=0; %Number of neutrons in fragment 1
            nf2=0; %Number of neutrons in fragment 2
            
            Pf1_sacr=zeros(kP,3);
            Pf2_sacr=zeros(kP,3);
            Closer_f1_sacr=zeros(kP+kN,3);
            Closer_f2_sacr=zeros(kP+kN,3);
            
            for i=1:numProtons
               dist=abs(coeff_fiss_plane(f,1)*nucleon(feasible_buildup(i+240,s),1)+coeff_fiss_plane(f,2)*nucleon(feasible_buildup(i+240,s),2)+coeff_fiss_plane(f,3)*nucleon(feasible_buildup(i+240,s),3)-coeff_fiss_plane(f,4))/sqrt(coeff_fiss_plane(f,1)^2+coeff_fiss_plane(f,2)^2+coeff_fiss_plane(f,3)^2);
               if coeff_fiss_plane(f,1)*nucleon(feasible_buildup(i+240,s),1)+coeff_fiss_plane(f,2)*nucleon(feasible_buildup(i+240,s),2)+coeff_fiss_plane(f,3)*nucleon(feasible_buildup(i+240,s),3)<coeff_fiss_plane(f,4)
                  Pf1_sacr(pf1+1,:)=nucleon(feasible_buildup(i+240,s),:);
                  pf1=pf1+1;
                  if dist<2
                      Closer_f1_sacr(closer_f1+1,:)=nucleon(feasible_buildup(i+240,s),:);
                      closer_f1=closer_f1+1;
                  end
               else
                  Pf2_sacr(pf2+1,:)=nucleon(feasible_buildup(i+240,s),:);
                  pf2=pf2+1;
                  if dist<2
                      Closer_f2_sacr(closer_f2+1,:)=nucleon(feasible_buildup(i+240,s),:);
                      closer_f2=closer_f2+1;
                  end
               end
            end
            Pf1=Pf1_sacr(1:pf1,:);
            Pf2=Pf2_sacr(1:pf2,:);
            
            for i=1:numNeutrons
               dist=abs(coeff_fiss_plane(f,1)*nucleon(feasible_buildup(i,s),1)+coeff_fiss_plane(f,2)*nucleon(feasible_buildup(i,s),2)+coeff_fiss_plane(f,3)*nucleon(feasible_buildup(i,s),3)-coeff_fiss_plane(f,4))/sqrt(coeff_fiss_plane(f,1)^2+coeff_fiss_plane(f,2)^2+coeff_fiss_plane(f,3)^2);
               if coeff_fiss_plane(f,1)*nucleon(feasible_buildup(i,s),1)+coeff_fiss_plane(f,2)*nucleon(feasible_buildup(i,s),2)+coeff_fiss_plane(f,3)*nucleon(feasible_buildup(i,s),3)<coeff_fiss_plane(f,4)
                  nf1=nf1+1;
                  if dist<2
                      Closer_f1_sacr(closer_f1+1,:)=nucleon(feasible_buildup(i,s),:);
                      closer_f1=closer_f1+1;
                  end
               else
                  nf2=nf2+1;
                  if dist<2
                      Closer_f2_sacr(closer_f2+1,:)=nucleon(feasible_buildup(i,s),:);
                      closer_f2=closer_f2+1;
                  end
               end
            end
            
            Closer_f1=Closer_f1_sacr(1:closer_f1,:);
            Closer_f2=Closer_f2_sacr(1:closer_f2,:);

            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,1)=pf1;
            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,1)=pf2;
            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,2)=nf1;
            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,2)=nf2;
            
            if pf1>0 && pf2>0
                for i=1:pf1
                    for j=1:pf2
                        dist=fm*sqrt((Pf1(i,1)-Pf2(j,1))^2+(Pf1(i,2)-Pf2(j,2))^2+(Pf1(i,3)-Pf2(j,3))^2);
                        fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,4)=fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,4)+1.439965/dist;
                        fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,4)=fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,4)+1.439965/dist;
                    end
                end
            end
            
            if closer_f1>0 && closer_f2>0
                for i=1:closer_f1
                    for j=1:closer_f2
                        dist=sqrt((Closer_f1(i,1)-Closer_f2(j,1))^2+(Closer_f1(i,2)-Closer_f2(j,2))^2+(Closer_f1(i,3)-Closer_f2(j,3))^2);
                        if dist==2*sqrt(2)
                            broken_bonds=broken_bonds+1;
                        end
                    end
                end
            end
            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,5)=beta_fission*broken_bonds-fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f-1,4);
            fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,5)=beta_fission*broken_bonds-fragments_Z_N_A_coul_BE((s-1)*n_fission_planes*2+2*f,4);
            
            clear Pf1
            clear Pf2
            clear Closer_f1
            clear Closer_f2
        end
        progress = (s / feas_struct) * 100;
        fprintf('Performing fission of the nucleus (21 planes x %d FCC models) - Progress: %.2f%%\n', feas_struct, progress);
    end
    fragments_Z_N_A_coul_BE(:,3)=fragments_Z_N_A_coul_BE(:,1)+fragments_Z_N_A_coul_BE(:,2);
    
    figure
    histogram(fragments_Z_N_A_coul_BE(:,5))
    title(['Distribution of energy values (MeV) to induce fission along the fission planes for the ',num2str(feas_struct),' different FCC structures'])
    xlabel('Energy value [MeV]')
    ylabel('Occurrence [-]')
    
    answer=inputdlg({'Type the minimum energy value (in MeV) you are interested in to analyze the fission statistics:'});
    energy_threshold_min=str2double(answer);
    answer=inputdlg({'Type the maximum energy value (in MeV) you are interested in to analyze the fission statistics:'});
    energy_threshold_max=str2double(answer);
    
    clear F
    f_statistics=0;
    F_sacr=zeros(n_fission_planes*2*feas_struct,5);
    for i=1:n_fission_planes*2*feas_struct
       if fragments_Z_N_A_coul_BE(i,5)<=energy_threshold_max && fragments_Z_N_A_coul_BE(i,5)>=energy_threshold_min
          F_sacr(f_statistics+1,:)=fragments_Z_N_A_coul_BE(i,:);
          f_statistics=f_statistics+1;
          progress = (i / (n_fission_planes*2*feas_struct)) * 100;
          fprintf('Computing energy-weighted fission results - Progress: %.2f%%\n', progress);
       end
    end
    F=F_sacr(1:f_statistics,:);
    
    figure
    subplot(2,1,1)
    histogram(F(:,1))
    xlabel('Z')
    ylabel('Occurrence [-]')
    title(['Distribution of fragments atomic number Z between ',num2str(energy_threshold_min),' MeV and ',num2str(energy_threshold_max),' MeV'])
    subplot(2,1,2)
    histogram(F(:,3))
    xlabel('A')
    ylabel('Occurrence [-]')
    title(['Distribution of fragments mass number A between ',num2str(energy_threshold_min),' MeV and ',num2str(energy_threshold_max),' MeV'])
    
    prob=zeros(n_fission_planes*2*feas_struct,1);
    for i=1:size(F,1)
        prob(i)=1/F(i,5);
    end
    prob=prob/sum(prob);

    %If some energies are negative (or 0) use something like this instead:
    %prob=zeros(n_fission_planes*2*feas_struct,1);
    %minimum_energy=min(F(:,5));
    %for i=1:size(F,1)
    %    prob(i)=1/(F(i,5)+minimum_energy+0.001);
    %end
    %prob=prob/sum(prob);

    Z_fragm=zeros(numProtons+1,1);
    A_fragm=zeros(numProtons+numNeutrons+1,1);
    for i=1:size(F,1)
        Z_fragm(F(i,1)+1,1)=Z_fragm(F(i,1)+1,1)+prob(i);
        A_fragm(F(i,3)+1,1)=A_fragm(F(i,3)+1,1)+prob(i);
    end
    
    figure
    subplot(2,1,1)
    plot(linspace(0,numProtons,numProtons+1),Z_fragm)
    xlabel('Z')
    ylabel('Probability (energy)-weighted occurrence [-]')
    title(['Distribution of fragments atomic number Z between ',num2str(energy_threshold_min),' MeV and ',num2str(energy_threshold_max),' MeV'])
    subplot(2,1,2)
    plot(linspace(0,numProtons+numNeutrons,numProtons+numNeutrons+1),A_fragm)
    xlabel('A')
    ylabel('Probability (energy)-weighted occurrence [-]')
    title(['Distribution of fragments mass number A between ',num2str(energy_threshold_min),' MeV and ',num2str(energy_threshold_max),' MeV']) 

    answer=inputdlg({'Do you want to take also into account the fragment decay for the statistics (Y/N)?'});
    if strcmp(answer,'Y')
        Tt=[60;3600;86400;2629800;31557600;31557600000;31557600000000;31557600000000000];%half-life threshold in seconds
        stability_table=xlsread('Nuclear_stability.xlsx','Foglio1','B3:CP158');
        half_lives=xlsread('Nuclear_stability.xlsx','Foglio2','B3:CP158');
        Z_distribution=zeros(numProtons+1,size(Tt,1));
        A_distribution=zeros(numProtons+numNeutrons+1,size(Tt,1));
        for time_scales=1:size(Tt,1)
            F_modified=F;
            final_fragments=zeros(size(F_modified,1),4);
            n_fragm=0;
            for i=1:size(F_modified,1)
               if stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==0
                   n_fragm=n_fragm+1;
                   final_fragments(n_fragm,1:2)=F_modified(i,1:2);
                   final_fragments(n_fragm,3)=final_fragments(n_fragm,1)+final_fragments(n_fragm,2);
                   final_fragments(n_fragm,4)=F_modified(i,5);
               elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)>0
                   if (1-exp(-log(2)*Tt(time_scales)/half_lives(F_modified(i,2)+1,F_modified(i,1)+1)))<rand()
                       n_fragm=n_fragm+1;
                       final_fragments(n_fragm,1:2)=F_modified(i,1:2);
                       final_fragments(n_fragm,3)=final_fragments(n_fragm,1)+final_fragments(n_fragm,2);
                       final_fragments(n_fragm,4)=F_modified(i,5);
                   else
                       while stability_table(F_modified(i,2)+1,F_modified(i,1)+1)>0 && (1-exp(-log(2)*Tt(time_scales)/half_lives(F_modified(i,2)+1,F_modified(i,1)+1)))>=rand()
                            if stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==1
                                F_modified(i,1)=F_modified(i,1)-1;
                                F_modified(i,2)=F_modified(i,2)+1;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==2
                                F_modified(i,1)=F_modified(i,1)+1;
                                F_modified(i,2)=F_modified(i,2)-1;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==3
                                F_modified(i,2)=F_modified(i,2)-1;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==4
                                F_modified(i,2)=F_modified(i,2)-2;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==5
                                F_modified(i,1)=F_modified(i,1)-1;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==6
                                F_modified(i,1)=F_modified(i,1)-2;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==7
                                F_modified(i,1)=F_modified(i,1)-1;
                                F_modified(i,2)=F_modified(i,2)+1;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==8
                                F_modified(i,1)=F_modified(i,1)-3;
                            elseif stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==9
                                F_modified(i,1)=F_modified(i,1)-2;
                                F_modified(i,2)=F_modified(i,2)-2;                
                            end
                       end
                       if stability_table(F_modified(i,2)+1,F_modified(i,1)+1)==0 || (1-exp(-log(2)*Tt(time_scales)/half_lives(F_modified(i,2)+1,F_modified(i,1)+1)))<rand()
                           n_fragm=n_fragm+1;
                           final_fragments(n_fragm,1:2)=F_modified(i,1:2);
                           final_fragments(n_fragm,3)=final_fragments(n_fragm,1)+final_fragments(n_fragm,2);
                           final_fragments(n_fragm,4)=F_modified(i,5);
                       end
                   end
               end
            end
            cut=0;
            for i=1:size(final_fragments,1)
               if final_fragments(i,1)==0 && final_fragments(i,2)==0 && final_fragments(i,3)==0 && final_fragments(i,4)==0
                   if cut==0
                       cut=i;
                   end
               end
            end
            final_fragments=final_fragments(1:(cut-1),1:4);
            
            prob=zeros(size(final_fragments,1),1);
            for i=1:size(final_fragments,1)
                prob(i)=1/final_fragments(i,4);
            end
            prob=prob/sum(prob);

            %If some energies are negative (or 0) use something like this instead:
            %prob=zeros(size(final_fragments,1),1);
            %minimum_energy=min(final_fragments(:,4));
            %for i=1:size(F,1)
            %    prob(i)=1/(final_fragments(i,4)+minimum_energy+0.001);
            %end
            %prob=prob/sum(prob);

            Z_fragm=zeros(numProtons+1,1);
            A_fragm=zeros(numProtons+numNeutrons+1,1);
            for i=1:size(final_fragments,1)
                Z_fragm(final_fragments(i,1)+1,1)=Z_fragm(final_fragments(i,1)+1,1)+prob(i);
                A_fragm(final_fragments(i,3)+1,1)=A_fragm(final_fragments(i,3)+1,1)+prob(i);
            end
            Z_distribution(:,time_scales)=Z_fragm;
            A_distribution(:,time_scales)=A_fragm;
            clear final_fragments
            
            progress = (time_scales / size(Tt,1)) * 100;
            fprintf('Computing decay results for Tt = %d seconds - Progress: %.2f%%\n', Tt(time_scales), progress);
        
        end
        
        figure
        subplot(1,2,1)
        for time_scales=1:size(Tt,1)
            plot(linspace(0,numProtons,numProtons+1),Z_distribution(:,time_scales),'LineWidth',2)
            hold on
        end
        xlabel('Z')
        ylabel('Probability (energy)-weighted occurrence [-]')
        subplot(1,2,2)
        for time_scales=1:size(Tt,1)
            plot(linspace(0,numProtons+numNeutrons,numProtons+numNeutrons+1),A_distribution(:,time_scales),'LineWidth',2)
            hold on
        end
        xlabel('A')
        legend('T_{char} = 1 minute','T_{char} = 1 hour','T_{char} = 1 day','T_{char} = 1 month','T_{char} = 1 year','T_{char} = 1000 years','T_{char} = 1 million years','T_{char} = 1 billion years')
    end
end