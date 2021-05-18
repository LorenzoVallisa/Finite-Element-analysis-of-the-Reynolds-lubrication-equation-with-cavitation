%%%% PART 1
%%%% The penalty operator problem  applied to bearings lubrication
%%
clc
clear all
close all

% Variabile attiva esegue i plot
debug_plot = 0;
incr_debug=0;
cavitat_area_debug=0;
somma = debug_plot + incr_debug + cavitat_area_debug;
if somma >1
    msg = ' One debug test a time allowed ';
    error(msg);
    break;
end

Test = 'Test1';
eta=0.95;
toll = 1e-9;

for jj = 3:6
    
    nRef = jj;
    
    % Soluzione Iniziale
    
    [solutions,femregion,Matrices,Dati,A,M]=C_main2D(Test,nRef,'non_sym',0,[]);
    p1 = solutions.uh;
    
    if debug_plot
        figure;
        trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p1));
        xlabel('x [-]'); ylabel('y [-]'); zlabel('P [-]');   title([' Pressure distribution initial solution '] );
        colorbar
    end
    
    % Intorduzione obstacle
    
    my_coord = femregion.coord;
    y = my_coord(:,2);
    x = my_coord(:,1);
    obs_f = Dati.obstacle;
    obs_ev = eval(obs_f);
    
    
    if debug_plot
        figure;
        trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(obs_ev));
        xlabel('x [-]'); ylabel('y [-]'); zlabel('P [-]');  title([' Obstacle '] );
        colorbar
    end
    
    
    if strcmpi(Test,'Test2')
        
        exact = Dati.exact_sol;
        exact_solution = eval(exact);
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),exact_solution);
            xlabel('x [-]'); ylabel('y [-]'); zlabel('P[-]');   title([' Exact Solution (Test2) '] );
            colorbar
        end
        
    end
    
    
    i = 1;
    abres =[];
    abres(i) = 1;
    
    while (abres(i) >= toll)
        
        
        xx = find(p1>obs_ev);
        
        % Attivare solo in caso si voglia vedere evoluzione area di
        % cavitazione per un determinato numero di refinement (disattivare prima però altri plot)
        if cavitat_area_debug
            yy = find(p1 < -0.000001);
            val = femregion.coord(yy,:);
            subplot(2,5,i),plot(val(:,1),val(:,2),'.r','LineWidth',3),hold on,grid on,grid minor,ylim([0 6.28]),xlim([0 6.28]),xlabel('x [-]'),ylabel('y [-]')
            drawnow
        end
        
        
        [solutions,femregion,Matrices,Dati,A,M]=C_main2D(Test,nRef,'non_sym',1,xx);
        p2 = solutions.uh;
        
        % Relaxation
        
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p2));
            title([' Pressure distribution before relaxation '] ); xlabel('x-axis'); ylabel('y-axis'); zlabel('P[-]');
            drawnow;
        end
        
        
        p3 = p1*(1-eta) + (eta)*p2;
        
        if incr_debug
            errp = abs(p3 - p1);
            subplot(2,3,ii)
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(errp))
            colorbar
            zlim([0 0.28])
            xlabel('x[-]');
            ylabel('y[-]');
        end

        p1=p3;
        
        
        res = A* p1 - Matrices.b;
        
        
        abres = [abres norm(res,2)];
        
        
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p1));
            xlabel('x [-]'); ylabel('y [-]'); zlabel('P[-]'); title([' Pressure distribution lifted '] );
            colorbar
            drawnow;
        end
        
        i = i+1;
        
        
        
        
    end
    
    
    %     if strcmpi(Test,'Test2')
    %
    %         %Soluzione esatta disponibile solo per Test2
    %         save([ 'p'  num2str(jj) '3.mat'],'p1')
    %
    %     end
    
end



%%
%%% PART 2
%% Error and Convergence Analysis (NEL CASO IN CUI SI VOGLIA VEDERE L'ERRORE PER UN EPSILON FISSATO)->(SOLO Test2)


% NOTA_1: Una volta deciso un epsilon, bisogna ri-runnare la simulazione qui
% di sopra per ogni valore di nRef dato quell'epsilon 

%NOTA_2: Per analisi di convergenza usando metodi numerici vedi cartella
% Numerical Error
clc
close all
clear all

addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_3')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_4')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_5')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_6')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Errors')




uh3 = load('p33.mat');
u_ex3 = load('p_exact_3.mat');



uh4 = load('p43.mat');
u_ex4 = load('p_exact_4.mat');

uh5 = load('p53.mat');
u_ex5 = load('p_exact_5.mat');


uh6 = load('p63.mat');
u_ex6 = load('p_exact_6.mat');


ref3 = struct('uh',uh3.p1,'u_ex',u_ex3.p1);

ref4 = struct('uh',uh4.p1,'u_ex',u_ex4.p1);

ref5= struct('uh',uh5.p1,'u_ex',u_ex5.p1);

ref6 = struct('uh',uh6.p1,'u_ex',u_ex6.p1);


sol_struct_vect = [ref3 ref4 ref5 ref6];

final_ref = 6;


[table_error,rates]= test_convergence_1(sol_struct_vect,final_ref)

%%
%%% PART 3
%% Error and Convergence Analysis (MODIFIED TO TAKE INTO ACCOUNT OTHER EPSILON)->(SOLO Test2)

% NOTA: I file .mat inclusi nella cartella sono simulazioni già eseguite
% con epsilon diversi (h^2 è quello senza underscore). Se si vuole
% ri-runnare queste simulazioni bisogna prima cambiare il nome con cui
% vengono salvati i file (Line 130).  

clc
close all
clear all

addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_3')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_4')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_5')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Ref_6')
addpath('D:\POLIMI\Specialistica\2o semestre\Numerical analysis for partial differential equation\Progetto\Prova con CG\LAST_VERSION CODE\CG_FEM\Errors')




uh3_h2 = load('p33.mat');
u_ex3 = load('p_exact_3.mat');
uh3_h1 = load('p33_h1.mat');
uh3_h3 = load('p33_h3.mat');


uh4_h2 = load('p43.mat');
uh4_h1 = load('p43_h1.mat');
uh4_h3 = load('p43_h3.mat');

u_ex4 = load('p_exact_4.mat');

uh5_h2 = load('p53.mat');
uh5_h1 = load('p53_h1.mat');
uh5_h3 = load('p53_h3.mat');

u_ex5 = load('p_exact_5.mat');


uh6_h2 = load('p63.mat');
uh6_h1 = load('p63_h1.mat');
uh6_h3 = load('p63_h3.mat');

u_ex6 = load('p_exact_6.mat');

figure
for kkk = 1:3
    switch kkk
        case 1
            a = uh3_h1.p1;
            b = uh4_h1.p1;
            c = uh5_h1.p1;
            d = uh6_h1.p1;
        case 2
            a = uh3_h2.p1;
            b = uh4_h2.p1;
            c = uh5_h2.p1;
            d = uh6_h2.p1;
        case 3
            a = uh3_h3.p1;
            b = uh4_h3.p1;
            c = uh5_h3.p1;
            d = uh6_h3.p1;
    end
    
    
    ref3 = struct('uh',a,'u_ex',u_ex3.p1);
    
    ref4 = struct('uh',b,'u_ex',u_ex4.p1);
    
    ref5= struct('uh',c,'u_ex',u_ex5.p1);
    
    ref6 = struct('uh',d,'u_ex',u_ex6.p1);
    
    
    sol_struct_vect = [ref3 ref4 ref5 ref6];
    
    final_ref = 6;
    
    
    [hh,E_H1]= test_convergence_2(sol_struct_vect,final_ref)
    
    
    
    grid on
    if kkk ==1
        loglog(hh,E_H1,'-dm');
    elseif kkk == 2
        loglog(hh,E_H1,'-ob');
    elseif kkk == 3
        loglog(hh,E_H1,'-*k');
    end
    grid on
    hold on
    
end
grid on
loglog(hh,hh.^(1),'-sr');
hold on
xlabel('h'); ylabel('Error');
legend('||u-u_h||_{1,\Omega}  \epsilon ~ h^1',...
    '||u-u_h||_{1,\Omega}  \epsilon ~ h^2',...
    '||u-u_h||_{1,\Omega}  \epsilon ~ h^3',...
    'h',...
    'Location','BestOutside')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% EXTRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% PART 4
%% LAGRANGIAN AUGMENTED

% Metodo utile nel caso in cui si voglia forzare condizioni al bordo non
% lineari nel problema (i.e. studio di un differente fenomeno o di un
% fenomeno correlato alla cavitazione): riducendo infatti la distanza tra
% gli autovalori nel Penalty Problem aumentandone di conseguenza la
% stabilità del metodo (Nocedal : Numerical Optimization (P.511))

clc
clear all
close all

% Variabile attiva esegue i plot
debug_plot = 0;

Test = 'Test1';
eta=0.95;
toll = 1e-9;


% Initial Guess Lagrangian Multiplier (According to the guess convergence
% result may improve (detail not yet investaigated)

u = 0;


for jj = 3:6
    
    nRef = jj;
    
    % Soluzione Iniziale
    
    [solutions,femregion,Matrices,Dati,A,M]=C_main2D(Test,nRef,'non_sym',0,[],u);
    p1 = solutions.uh;
    u = zeros(size(p1,1),1);
    
    if debug_plot
        figure;
        trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p1));
        xlabel('x [-]'); ylabel('y [-]'); zlabel('P [-]');   title([' Pressure distribution initial solution '] );
        colorbar
    end
    
    % Intorduzione obstacle
    
    my_coord = femregion.coord;
    y = my_coord(:,2);
    x = my_coord(:,1);
    obs_f = Dati.obstacle;
    obs_ev = eval(obs_f);
    
    
    if debug_plot
        figure;
        trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(obs_ev));
        xlabel('x [-]'); ylabel('y [-]'); zlabel('P [-]');  title([' Obstacle '] );
        colorbar
    end
    
    
    if strcmpi(Test,'Test2')
        
        exact = Dati.exact_sol;
        exact_solution = eval(exact);
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),exact_solution);
            xlabel('x [-]'); ylabel('y [-]'); zlabel('P[-]');   title([' Exact Solution (Test2) '] );
            colorbar
        end
        
    end
    
    
    i = 1;
    abres =[];
    abres(i) = 1;
    
    while (abres(i) >= toll)
        
        
        xx = find(p1>obs_ev);
        
        % Attivare solo in caso si voglia vedere evoluzione area di
        % cavitazione per un determinato numero di refinement (disattivare prima però altri plot)
        
        %         yy = find(p1 < -0.000001);
        %         val = femregion.coord(yy,:);
        %         subplot(2,5,i),plot(val(:,1),val(:,2),'.r','LineWidth',3),hold on,grid on,grid minor,ylim([0 6.28]),xlim([0 6.28]),xlabel('x [-]'),ylabel('y [-]')
        %         drawnow
        
        u = M*u - M*(p1-obs_ev)/femregion.epsilon;
        [solutions,femregion,Matrices,Dati,A,M]=C_main2D(Test,nRef,'non_sym',1,xx,u);
        p2 = solutions.uh;
        
        % Relaxation 
        
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p2));
            title([' Pressure distribution before relaxation '] ); xlabel('x-axis'); ylabel('y-axis'); zlabel('P[-]');
            drawnow;
        end
        
        
        p3 = p1*(1-eta) + (eta)*p2;
        p1=p3;
        
        
        res = A* p1 - Matrices.b;
        
        
        abres = [abres norm(res,2)];
        
        
        if debug_plot
            figure;
            trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(p1));
            xlabel('x [-]'); ylabel('y [-]'); zlabel('P[-]'); title([' Pressure distribution lifted '] );
            colorbar
            drawnow;
        end
        
        i = i+1;
        
        
        
        
    end
    
    
    %     if strcmpi(Test,'Test2')
    %
    %         %Soluzione esatta disponibile solo per Test2
    %         save([ 'p'  num2str(jj) '3.mat'],'p1')
    %
    %     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EXTRA %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KKT
clear all
close all
clc


Test = 'Test1';
eta=0.95;
toll = 1e-9;
nRef = 3;

[cavitation_vector,solutions,femregion,Matrices,Dati,A,M]=C_main2D(Test,nRef,'non_sym',0,[]);
p1 = solutions.uh;
C = 0.001;
ss = (length(solutions.uh));
s = C*rand(ss,1);
S = diag((s));
lambda = C*rand(ss,1);
L = diag((lambda));

Big_A = [A -ones(ss,ss) zeros(ss,ss);
    -ones(ss,ss) zeros(ss,ss) ones(ss,ss);
    zeros(ss,ss) S L];
Big_b = - [A*p1 - Matrices.b - lambda; s - (p1 - cavitation_vector); L*S*ones(ss,1)];

cond(Big_A)

x = Big_A\Big_b;




