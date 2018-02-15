function plot_benchmark_signal(problem)
% plot_benchmark_signal.m
%
% Plot the results from test_all.sh
%
% Input: problem, where problem is a string of either
%
% problem20040
% problem260100
% problem5010
%
% If one runs
%
% make benchmark_small P=str NF=?? M=?? N=??
%
% with some other parameters, then call
%
% plot_benchmark(str)
%
% to see the results.
%
% Gives plot for
% 
% Optimal objective (to check we obtain comparable results)
% Timings
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, November 7, 2012

addpath ../data_process

problem_name = sprintf('../data/%s',problem);

filename_matlab = sprintf('../data/sol_%s_Matlab.mat',problem);

filename_pdip_slp_primal = sprintf('../data/sol_%s_primal',problem);
filename_pdip_slp_dual = sprintf('../data/sol_%s_dual',problem);;

filename_pdip_slp_primal_sd = sprintf('../data/sol_%s_primal_sd',problem);
filename_pdip_slp_dual_sd = sprintf('../data/sol_%s_dual_sd',problem);;

filename_cvxgen = sprintf('../data/sol_%s_cvxgen',problem);
filename_cvxgen_blas = sprintf('../data/sol_%s_cvxgen_blas',problem);

if(exist( sprintf('%s_%d',filename_cvxgen,1),'file'))
	CVXGEN = true;
else
	CVXGEN = false;
end

if(exist( sprintf('%s_%d',filename_cvxgen_blas,1),'file'))
	CVXGEN_BLAS = true;
else
	CVXGEN_BLAS = false;
end


load(filename_matlab);

% --- Process the data from pdip_slp C++ implementation ---

f = @(alpha,x,X,gamma) norm(x-X*alpha,1)+gamma*norm(alpha,1);

% For results
timing_cpp_pdip_slp_primal = zeros(no_of_frames,1);
timing_cpp_pdip_slp_dual = zeros(no_of_frames,1);
timing_cpp_pdip_slp_primal_sd = zeros(no_of_frames,1);
timing_cpp_pdip_slp_dual_sd = zeros(no_of_frames,1);

result_cpp_pdip_slp_primal = zeros(n,no_of_frames);
result_cpp_pdip_slp_dual = zeros(n,no_of_frames);
result_cpp_pdip_slp_primal_sd = zeros(n,no_of_frames);
result_cpp_pdip_slp_dual_sd = zeros(n,no_of_frames);

obj_cpp_pdip_slp_primal = zeros(no_of_frames,1);
obj_cpp_pdip_slp_dual = zeros(no_of_frames,1);
obj_cpp_pdip_slp_primal_sd = zeros(no_of_frames,1);
obj_cpp_pdip_slp_dual_sd = zeros(no_of_frames,1);

for p = 1:no_of_frames
	[k status time alpha] = ...
			readSolution(sprintf('%s_%d',filename_pdip_slp_primal,p));
	if status == 0
		fprintf('Problem number %d was not solved correctly by pdip_slp_primal\n',p);
	end
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
	timing_cpp_pdip_slp_primal(p) = time;
	result_cpp_pdip_slp_primal(:,p) = alpha;
	obj_cpp_pdip_slp_primal(p) = f(alpha,x,X,gamma);
	
end

for p = 1:no_of_frames
	[k status time alpha] = ...
			readSolution(sprintf('%s_%d',filename_pdip_slp_dual,p));
	if status == 0
		fprintf('Problem number %d was not solved correctly by pdip_slp_dual\n',p);
	end
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
	timing_cpp_pdip_slp_dual(p) = time;
	result_cpp_pdip_slp_dual(:,p) = alpha;
	obj_cpp_pdip_slp_dual(p) = f(alpha,x,X,gamma);
end

for p = 1:no_of_frames
	[k status time alpha] = ...
			readSolution(sprintf('%s_%d',filename_pdip_slp_primal_sd,p));
	if status == 0
		fprintf('Problem number %d was not solved correctly by pdip_slp_primal_sd\n',p);
	end
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
	timing_cpp_pdip_slp_primal_sd(p) = time;
	result_cpp_pdip_slp_primal_sd(:,p) = alpha;
	obj_cpp_pdip_slp_primal_sd(p) = f(alpha,x,X,gamma);
	
end

for p = 1:no_of_frames
	[k status time alpha] = ...
			readSolution(sprintf('%s_%d',filename_pdip_slp_dual_sd,p));
	if status == 0
		fprintf('Problem number %d was not solved correctly by pdip_slp_dual_sd\n',p);
	end
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
	timing_cpp_pdip_slp_dual_sd(p) = time;
	result_cpp_pdip_slp_dual_sd(:,p) = alpha;
	obj_cpp_pdip_slp_dual_sd(p) = f(alpha,x,X,gamma);
end

if(CVXGEN)
	timing_cvxgen = zeros(no_of_frames,1);

	result_cvxgen = zeros(n,no_of_frames);

	obj_cvxgen = zeros(no_of_frames,1);

	for p = 1:no_of_frames
		[k status time alpha] = ...
				readSolution(sprintf('%s_%d',filename_cvxgen,p));
		if status == 0
			fprintf('Problem number %d was not solved correctly by cvxgen\n',p);
		end
		[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
		timing_cvxgen(p) = time;
		result_cvxgen(:,p) = alpha;
		obj_cvxgen(p) = f(alpha,x,X,gamma);
	end
end

if(CVXGEN_BLAS)
	timing_cvxgen_blas = zeros(no_of_frames,1);

	result_cvxgen_blas = zeros(n,no_of_frames);

	obj_cvxgen_blas = zeros(no_of_frames,1);

	for p = 1:no_of_frames
		[k status time alpha] = ...
				readSolution(sprintf('%s_%d',filename_cvxgen_blas,p));
		if status == 0
			fprintf('Problem number %d was not solved correctly by cvxgen\n',p);
		end
		[X x gamma] = readProblem(sprintf('%s_%d',problem_name,p));
	
		timing_cvxgen_blas(p) = time;
		result_cvxgen_blas(:,p) = alpha;
		obj_cvxgen_blas(p) = f(alpha,x,X,gamma);
	end
end

	
% --- Plot the solutions ---------------------------------


mark = {'ko','kx','k*','k+','ks','kd','kv','k^','kh','kp'};

marker_size = 12;
line_width = 2;

no_of_points = 40;
points = round(linspace(1,no_of_frames,no_of_points));


figure(2)
clf
hMosek = plot(obj_Mosek,mark{1},'markersize',marker_size,'linewidt',line_width);
hold on
hCVX = plot(obj_CVX_SeDuMi,mark{2},'markersize',marker_size,'linewidt',line_width);
hMprimal = plot(obj_pdip_slp_primal,mark{3},'markersize',marker_size,'linewidt',line_width);
hMdual = plot(obj_pdip_slp_dual,mark{4},'markersize',marker_size,'linewidt',line_width);
hCprimal = plot(obj_cpp_pdip_slp_primal,mark{5},'markersize',marker_size,'linewidt',line_width);
hCdual = plot(obj_cpp_pdip_slp_dual,mark{6},'markersize',marker_size,'linewidt',line_width);

if(CVXGEN)
	hcvxgen = plot(obj_cvxgen,mark{7},'markersize',marker_size,'linewidt',line_width);
end

if(CVXGEN_BLAS)
	hcvxgen_blas = plot(obj_cvxgen_blas,mark{8},'markersize',marker_size,'linewidt',line_width);
end

hCprimalsd = plot(obj_cpp_pdip_slp_primal_sd,mark{9},'markersize',marker_size,'linewidt',line_width);
hCdualsd = plot(obj_cpp_pdip_slp_dual_sd,mark{10},'markersize',marker_size,'linewidt',line_width);

xlabel('frame')
ylabel('fs')

l = [hMosek hCVX hMprimal hMdual hCprimal hCdual hCprimalsd hCdualsd];
L =	{'Mosek','CVX+SeDuMi','Mprimal','Mdual','Cprimal','Cdual','Cprimal(sd)','Cdual(sd)'};

if(CVXGEN)
	l =[l hcvxgen];
	L = {L{:} 'CVXGEN'};
end

if(CVXGEN_BLAS)
	l =[l hcvxgen_blas];
	L = {L{:} 'CVXGEN_BLAS'};
end

legend(l,L);

figure(1)
clf
semilogy(timing_Mosek,'-k','markersize',marker_size,'linewidt',line_width);
hold on
semilogy(points,timing_Mosek(points),mark{1},'markersize',marker_size,'linewidt',line_width);

semilogy(timing_CVX_SeDuMi,'-k','markersize',marker_size,'linewidt',line_width);
semilogy(points,timing_CVX_SeDuMi(points),mark{2},'markersize',marker_size,'linewidt',line_width);

semilogy(timing_pdip_slp_primal,'-k','markersize',marker_size,'linewidt',line_width);
semilogy(points,timing_pdip_slp_primal(points),mark{3},'markersize',marker_size,'linewidt',line_width);

semilogy(timing_pdip_slp_dual,-k,'markersize',marker_size,'linewidt',line_width);
semilogy(points,timing_pdip_slp_dual(points),mark{4},'markersize',marker_size,'linewidt',line_width);

semilogy(timing_cpp_pdip_slp_primal,'-k','markersize',marker_size,'linewidt',line_width);
semilogy(points,timing_cpp_pdip_slp_primal(points),mark{5},'markersize',marker_size,'linewidt',line_width);

semilogy(timing_cpp_pdip_slp_dual,-k,'markersize',marker_size,'linewidt',line_width);
semilogy(timing_cpp_pdip_slp_dual(points),mark{6},'markersize',marker_size,'linewidt',line_width);

if( CVXGEN )
	semilogy(timing_cvxgen,'-k','markersize',marker_size,'linewidt',line_width);
	hcvxgen = semilogy(points,timing_cvxgen(points),mark{7},'markersize',marker_size,'linewidt',line_width);
end

if( CVXGEN_BLAS )
	semilogy(timing_cvxgen_blas,'-k','markersize',marker_size,'linewidt',line_width);
	hcvxgen_blas = semilogy(points,timing_cvxgen_blas(points),mark{8},'markersize',marker_size,'linewidt',line_width);
end


xlabel('frame')
ylabel('time[s]')

leg = legend(l,L);

%set(leg,'FontSize',20)


% ---Print average timings out for possible transfer to table----
fprintf('--------------------------------------------------------\n');
fprintf('                %s\n',problem)
fprintf('---------------Avg (min/max)----------------------------\n');
fprintf('CVX+SeDumi    %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_CVX_SeDuMi),1000*min(timing_CVX_SeDuMi),1000*max(timing_CVX_SeDuMi));
fprintf('Mosek         %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_Mosek),1000*min(timing_Mosek),1000*max(timing_Mosek));
fprintf('Mprimal       %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_pdip_slp_primal),1000*min(timing_pdip_slp_primal),1000*max(timing_pdip_slp_primal));
fprintf('Mdual         %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_pdip_slp_dual),1000*min(timing_pdip_slp_dual),1000*max(timing_pdip_slp_dual));
fprintf('Cprimal       %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_cpp_pdip_slp_primal), 1000*min(timing_cpp_pdip_slp_primal),1000*max(timing_cpp_pdip_slp_primal));
fprintf('Cdual         %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_cpp_pdip_slp_dual),1000*min(timing_cpp_pdip_slp_dual),1000*max(timing_cpp_pdip_slp_dual));
fprintf('Cprimal(sd)   %6.2f (%6.2f/%6.2f) [ms]\n',...
1000*mean(timing_cpp_pdip_slp_primal_sd),1000*min(timing_cpp_pdip_slp_primal_sd), 1000*max(timing_cpp_pdip_slp_primal_sd));
fprintf('Cdual(sd)     %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_cpp_pdip_slp_dual_sd),1000*min(timing_cpp_pdip_slp_dual_sd),1000*max(timing_cpp_pdip_slp_dual_sd));

if( CVXGEN )
fprintf('CVXGEN        %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_cvxgen),1000*min(timing_cvxgen),1000*max(timing_cvxgen));
end

if( CVXGEN_BLAS )
fprintf('CVXGEN_BLAS   %6.2f (%6.2f/%6.2f) [ms]\n',...
		1000*mean(timing_cvxgen_blas),1000*min(timing_cvxgen_blas),1000*max(timing_cvxgen_blas));
end


% --- minimum object calculations ---
if( CVXGEN )
	if (CVXGEN_BLAS)
		F = zeros(10,no_of_frames);
	else
		F = zeros(9,no_of_frames);
	end
else
	F = zeros(8,no_of_frames);
end


F(1,:) = obj_Mosek;
F(2,:) = obj_CVX_SeDuMi;
F(3,:) = obj_pdip_slp_primal;
F(4,:) = obj_pdip_slp_dual;
F(5,:) = obj_cpp_pdip_slp_primal;
F(6,:) = obj_cpp_pdip_slp_dual;
F(7,:) = obj_cpp_pdip_slp_primal_sd;
F(8,:) = obj_cpp_pdip_slp_dual_sd;

if( CVXGEN )
	F(9,:) = obj_cvxgen;
end

if( CVXGEN_BLAS )
	F(10,:) = obj_cvxgen_blas;
end

fprintf('Max relative deviation in the objective %2.4e\n',max((max(F)-min(F))./min(F)));
