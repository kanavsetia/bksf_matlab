% =========================================================================
%         Jordan-Wigner and Bravyi Kitaev Simulation for fermions
%
% Author: 	Kanav Setia*
% Date:     2016/11/19
% Version: 	0.1
%
% *Whitfield Group, Department of Physics and Astronomy,
%  Dartmouth College
% =========================================================================
% HELP SECTION
%
% This is the main file for the Code. The first part resets the workspace.
% Then the Flags defined are used to determine which algorithm needs to
% called, or if comparison of both the algorithm is required. Plots
% generated at the end are also determined by these flags. Please
% initialize these flags appropriately as per your requirements. After
% these flags, all this file is doing is calling the functions and plotting
% the data. The code is well commented and should be easy to read.
%
% G matrix --> Matrix that defines the intraction between the modes.
%
% Inputs : You need to manualy enter the G matrices in qubinit.m file which
% is the initialization file and later on we will also add unit tests in
% this file. It is also required that the flags be initialized
% appropirately.
%
% Outputs: Data_JW contains the data structure returned by jordan wigner
% algorithm and Data_BK contains the data structure returned by
% bravyi-kitaev algorithm. To see the the organization of data in
% data structures, see help sections of algorithms.
% See also jorwig, bk, bk_general, Ajk2, hydrogen, qubinit.
% =========================================================================

% Reset the workspace
clc;				% clear command window
clear all; 			% clear all objects
close all;     		% close all figures

Flag_JW=true;       %Flag for running jordan-wigner algorithm
Flag_BK=true;       %Flag for running Bravyi-Kitaev (nearest neighbor
%                   interaction model) algorithm
Flag_BK_general=false;%Flag for running Bravyi-Kitaev (general) algorithm
Flag_3D_Animation=false; %Flag for generating 3D animation
Flag_GIF=false;     %Flag for generating gif file of 3D animation
Flag_JW_Plots=false; %Flag for generating Jordan Wigner results
Flag_BK_Plots=false; %Flag for generating Bravyi Kitaev results
Flag_Comp_plots=true; %Flag for comparing the results of two algorithms
Flag_2D_animation=false; %Flag for genrating animation in comparison plots
Flag_2D_gif=false;
Flag_hydrogen=false; %Run Script file for comparing BK2 hydrogen model
%with jordan wigner model

% =========================================================================
% Initialization for both the algorithms. Here we initialize the number of
% qubit and the G matrix
% =========================================================================
qubinit;
% =========================================================================
% Calling algorithms to get the data
% =========================================================================
if Flag_JW
    Data_JW=jorwig(tot_qub,G);
end

if Flag_BK
    Data_BK=bk(tot_qub,G); % bk is the Bravyi-kitaev algorithm where only
    % nearest neighbor hopping terms have been
    % considered.
end
if Flag_BK_general
    Data_BK=bk_general(tot_qub,G); %bk_general is the general algorithm
    %                                  where the hopping terms are determined
    %                                  by the G matrix
end
% =========================================================================
% Plots from Jordan Wigner
% =========================================================================
if Flag_JW_Plots
    JW_x=linspace(1,tot_qub,tot_qub);
    JW_y=Data_JW(1).expectation.probability;
    h=figure();
    h.Position=[1 1 2500 1000];
    title('Jordan-Wigner Algorithm')
    hold on
    subplot(2,2,1)
    plot(Data_JW(1).expectation.Time,Data_JW(1).expectation.tau,'b',...
        Data_JW(1).expectation.Time,Data_JW(1).expectation.nu,'--',...
        Data_JW(1).expectation.Time,Data_JW(1).expectation.tau+...
        Data_JW(1).expectation.nu,'r','linewidth',2)
    xlabel('Time')
    ylabel('Energy')
    title('JW-Energy variations')
    legend('Kinetic Energy','Potential Energy')
    subplot(2,2,2)
    hold on;
    legend_cell={};
    for d2=1:length(JW_y(1,:))
        plot(Data_JW(1).expectation.Time,JW_y(:,d2),'linewidth',2);
        legend_cell{d2}=strcat('particle',num2str(d2));
    end
    legend(legend_cell)
    xlabel('Time')
    ylabel('Probability')
    title('Time evolution of probabilities in different modes')
    subplot(2,2,3)
    plot(Data_JW(1).expectation.Time,Data_JW(1).expectation.n)
    xlabel('Time')
    ylabel('Number Operator')
    title('Total Particles')
end
% =========================================================================
% Plots from Bravyi Kitaev
% =========================================================================
if Flag_BK_Plots
    BK_x=linspace(1,tot_qub,tot_qub);
    BK_y=Data_BK(1).expectation.probability;
    h1=figure(2);
    h1.Position=[1 1 2500 1000];
    hold on
    title('Bravyi-Kitaev Algorithm')
    subplot(2,2,1)
    plot(Data_BK(1).expectation.Time,Data_BK(1).expectation.tau,'b',...
        Data_BK(1).expectation.Time,Data_BK(1).expectation.nu,'--',...
        Data_BK(1).expectation.Time,Data_BK(1).expectation.tau+...
        Data_BK(1).expectation.nu,'r','linewidth',2)
    xlabel('Time')
    ylabel('Energy')
    title('BK-Energy variations')
    legend('Kinetic Energy','Potential Energy','Total Energy')
    subplot(2,2,2)
    hold on;
    legend_cell={};
    for d2=1:length(BK_y(1,:))
        plot(Data_BK(1).expectation.Time,BK_y(:,d2),'linewidth',2);
        legend_cell{d2}=strcat('Mode ',num2str(d2));
    end
    legend(legend_cell)
    xlabel('Time')
    ylabel('Probability')
    title('Time evolution of probabilities in different modes')
    subplot(2,2,3)
    plot(Data_BK(1).expectation.Time,Data_BK(1).expectation.n)
    xlabel('Time')
    ylabel('Number Operator')
    title('Total Particles')
end
% =========================================================================
% Comparison
% =========================================================================
if Flag_Comp_plots
    h3=figure(3);
    h3.Position=[1 1 2500 1000];
    JW_x=linspace(1,tot_qub,tot_qub);
    JW_y=Data_JW(1).expectation.probability;
    subplot(2,2,1)
    % jwfi=stem(JW_x,JW_y(1,:),50,'r','MarkerFaceColor','r');
    jwfi=stem(JW_y(1,:),'filled');
    jwgca=gca;
    xlabel('Mode','FontSize',14.5,'FontWeight','bold')
    ylabel('Probability','FontSize',14.5,'FontWeight','bold')
    title('Time evolution of probabilities in different modes in JW',...
        'FontSize',14.5,'FontWeight','bold');
    axis([0 tot_qub+1 0 1.2]);
    jwgca.XTick=linspace(0,tot_qub+1, tot_qub+2);
    grid on;
    a=50;
    BK_x=linspace(1,tot_qub,tot_qub);
    BK_y=Data_BK(1).expectation.probability;
    subplot(2,2,2)
    % bkfi=stem(BK_x,BK_y(1,:),50,'b','MarkerFaceColor','b');
    bkfi=stem(JW_y(1,:),'r','filled');
    bkgca=gca;
    xlabel('Mode','FontSize',14.5,'FontWeight','bold')
    ylabel('Probability','FontSize',14.5,'FontWeight','bold')
    title('Time evolution of probabilities in different modes in BK'...
        ,'FontSize',14.5,'FontWeight','bold')
    axis([0 tot_qub+1 0 1.2]);
    bkgca.XTick=linspace(0,tot_qub+1, tot_qub+2);
    grid on;
    a=50;
%     subplot(2,2,3)
%     xlabel('Time','FontSize',14.5,'FontWeight','bold')
%     ylabel('Probability','FontSize',14.5,'FontWeight','bold')
%     title('Time evolution of probabilities in different modes in JW'...
%         ,'FontSize',14.5,'FontWeight','bold')
%     hold on;
%     legend_cell={};
%     for d2=1:length(JW_y(1,:))
%         plot(Data_JW(1).expectation.Time,JW_y(:,d2),'linewidth',4);
%         legend_cell{d2}=strcat('Mode ',num2str(d2));
%     end
%     legend(legend_cell)
    h4=figure(2);
    h4.Position=[1 1 2500 1000];
    xlabel('Time','FontSize',30,'FontWeight','bold')
    ylabel('Probability','FontSize',30,'FontWeight','bold')
    title('Time evolution of probabilities in different modes in JW'...
        ,'FontSize',30,'FontWeight','bold')
    hold on;
    legend_cell={};
    for d2=1:length(JW_y(1,:))
        plot(Data_JW(1).expectation.Time,JW_y(:,d2),'linewidth',8);
        legend_cell{d2}=strcat('Mode ',num2str(d2));
    end
    legend(legend_cell)
    
    h5=figure(4);
    h5.Position=[1 1 2500 1000];
    xlabel('Time','FontSize',30,'FontWeight','bold')
    ylabel('Probability','FontSize',30,'FontWeight','bold')
    title('Time evolution of probabilities in different modes in BK'...
        ,'FontSize',30,'FontWeight','bold')
    hold on;
    legend_cell={};
    for d2=1:length(BK_y(1,:))
        plot(Data_BK(1).expectation.Time,BK_y(:,d2),'linewidth',8);
        legend_cell{d2}=strcat('Mode ',num2str(d2));
    end
    legend(legend_cell)
%     
%     subplot(2,2,4)
%     xlabel('Time','FontSize',14.5,'FontWeight','bold')
%     ylabel('Probability','FontSize',14.5,'FontWeight','bold')
%     title('Time evolution of probabilities in different modes in BK'...
%         ,'FontSize',14.5,'FontWeight','bold')
%     hold on;
%     legend_cell={};
%     for d2=1:length(BK_y(1,:))
%         plot(Data_BK(1).expectation.Time,BK_y(:,d2),'linewidth',4);
%         legend_cell{d2}=strcat('Mode ',num2str(d2));
%     end
%     legend(legend_cell)
    filename = 'bkgif2d.gif';
    if (Flag_2D_gif)
        frame = getframe(h3);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif',...
            'Loopcount',inf,'DelayTime',0.00005);
    end
    
    if(Flag_2D_animation)
        for d1=1:length(JW_y(:,1))
            %     jwfi.XDataSource='JW_x';
            jwfi.YDataSource='JW_y(d1,:)';
            %     bkfi.XDataSource='BK_x';
            bkfi.YDataSource='BK_y(d1,:)';
            refreshdata(jwfi);
            refreshdata(bkfi);
            drawnow;
            if (Flag_2D_gif)
                if mod(d1,2)==0
                    frame = getframe(h3);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    imwrite(imind,cm,filename,'gif',...
                        'WriteMode','append','DelayTime',0.00005);
                end
            end
            %             pause(0.1)
        end
    end
end
if (Flag_3D_Animation)
    p=Data_BK(1).expectation.probability;
    np=size(p);
    r=linspace(0,1,40);
    X=[];
    Y=[];
    Z=[];
    X1=[];
    Y1=[];
    Z1=[];
    C=[];
    h1=figure(4);
    h1.Position=[1 1 2500 1000];
    filename = 'bkgif.gif';
    for d1=1:np(1)
        if d1==1
            for d2=1:np(2)
                [val ind]=min(abs(p(d1,d2)/2-r));
                for d3=1:ind
                    points = round(r(d3)*70);
                    if points ~=0
                        [X1,Y1,Z1] = sphere(points);
                        X=[X;(r(d3).*X1(:))+2*d2];
                        Y=[Y;r(d3).*Y1(:)];
                        Z=[Z;r(d3).*Z1(:)];
                        C=[C;d3.*ones(length(X1(:)),1)];
                    end
                end
            end
            h=scatter3(X,Y,Z,30,C,'.');
            s=gca;
            %             s.Color='black';
            view(-10,20)
            axis([-1 15 -1 2 -1 2])
            xlabel('Particle','FontSize',14.5,'FontWeight','bold')
            ylabel('Y','FontSize',14.5,'FontWeight','bold')
            zlabel('Z','FontSize',14.5,'FontWeight','bold')
            title('Evolution of probabilities in different modes for Bravyi-Kitaev Algorithm'...
                ,'FontSize',14.5,'FontWeight','bold');
            annotation('TextBox',[0.5  0.7 0.2 0.1],'String',[char(169)...
                ,'Kanav Setia, Whitfield Group, Dartmouth College'],...
                'FitBoxToText','on');
            s.XTick=[];
            s.YTick=[];
            s.ZTick=[];
            camlight('headlight')
            if (Flag_GIF)
                frame = getframe(h1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif',...
                    'Loopcount',inf,'DelayTime',0.00005);
            end
            
        else
            X=[];
            Y=[];
            Z=[];
            X1=[];
            Y1=[];
            Z1=[];
            C=[];
            for d2=1:np(2)
                [val ind]=min(abs(p(d1,d2)/2-r));
                for d3=1:ind
                    points = round(r(d3)*70);
                    if points ~=0
                        [X1,Y1,Z1] = sphere(points);
                        X=[X;r(d3).*X1(:)+2*d2];
                        Y=[Y;r(d3).*Y1(:)];
                        Z=[Z;r(d3).*Z1(:)];
                        C=[C;d3.*ones(length(X1(:)),1)];
                    else
                    end
                end
            end
            h.XDataSource='X';
            h.YDataSource='Y';
            h.ZDataSource='Z';
            h.CDataSource='C';
            refreshdata(h);
            drawnow;
            if (Flag_GIF)
                if mod(d1,2)==0
                    frame = getframe(h1);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    imwrite(imind,cm,filename,'gif',...
                        'WriteMode','append','DelayTime',0.00005);
                end
            end
        end
    end
end
if Flag_hydrogen
    prob=hydrogen(G);
    JW_y=Data_JW(1).expectation.probability;
    probnorm=norm(JW_y(:,4)-prob(:,4))+norm(JW_y(:,3)-prob(:,3))+...
        norm(JW_y(:,2)-prob(:,2))+norm(JW_y(:,1)-prob(:,1))
    %Comparison of Hydrogen and Jorwig
    nf=figure(4);
    nf.Position=[1 1 2500 1000];
    JW_y=Data_JW(1).expectation.probability;
    title('Time evolution of probabilities of finding particles in different modes'...
        ,'FontSize',12);
    figure('Position',[230 250 670 410])
    %     subplot(2,2,1)
    plot(Data_JW(1).expectation.Time,JW_y(:,1),'linewidth',6);
    hold on;
    plot(Data_JW(1).expectation.Time,prob(:,1), 'linewidth',6);
    legend('Mode 1 JW','Mode 1 BK')
    xlabel('Time','FontSize',14.5)
    ylabel('Probability','FontSize',14.5)
    title('Time evolution of probability of finding a particle in mode 1'...
        ,'FontSize',12);
    %     subplot(2,2,2)
    figure('Position',[230 250 670 410])
    plot(Data_JW(1).expectation.Time,JW_y(:,2),'linewidth',6);
    hold on;
    plot(Data_JW(1).expectation.Time,prob(:,2), 'linewidth',6);
    legend('Mode 2 JW','Mode 2 BK')
    xlabel('Time','FontSize',14.5)
    ylabel('Probability','FontSize',14.5)
    title('Time evolution of probability of finding a particle in mode 2'...
        ,'FontSize',12);
    %     subplot(2,2,3)
    figure('Position',[230 250 670 410])
    plot(Data_JW(1).expectation.Time,JW_y(:,3),'linewidth',6);
    hold on;
    plot(Data_JW(1).expectation.Time,prob(:,3), 'linewidth',6);
    legend('Mode 3 JW','Mode 3 BK')
    xlabel('Time','FontSize',14.5)
    ylabel('Probability','FontSize',14.5)
    title('Time evolution of probability of finding a particle in mode 3'...
        ,'FontSize',12);
    figure('Position',[230 250 670 410])
    %     subplot(2,2,4)
    plot(Data_JW(1).expectation.Time,JW_y(:,4),'linewidth',6);
    hold on;
    plot(Data_JW(1).expectation.Time,prob(:,4), 'linewidth',6);
    legend('Mode 4 JW','Mode 4 BK')
    xlabel('Time','FontSize',14.5)
    ylabel('Probability','FontSize',14.5)
    title('Time evolution of probability of finding a particle in mode 4'...
        ,'FontSize',12);
end