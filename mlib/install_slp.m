function succ = install_slp(varargin)
% succ = install_slp(varargin)
%
% Installs the sparse linear prediction software using the default 
% BLAS and LAPACK libraries provided with Matlab. This can generate 
% some problems.  If so, please report this to the author.
%
% If you would like to link to another BLAS/LAPACK library, then
% provide the link line as e.g. install_slp('-lblas -llapack')
%
% AUTHOR:
%   Post-doc Tobias L. Jensen, Aalborg University, Denmark.
%     E-mail: tlj@es.aau.dk
%
% VERSION HISTORY:
%   0.1 [22-MAR-2013] Supports both windows and linux(unix) 
%                  based installation
%
% COPYRIGHT:
%   2012 Tobias L. Jensen
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE?2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
%   implied. See the License for the specific language governing
%   permissions and limitations under the License.


ext = mexext;
if strcmp(ext(end-1:end),'64')
	is64 = true;
	arg = '-largeArrayDims';
else
	arg = '';
	is64 = false;
end

CFLAGS = sprintf('-I../src %s', arg);

%-------------------------------------------------------------------------------
% BLAS and Lapack option
%-------------------------------------------------------------------------------

% This section is (almost) the same as done by T.A. Davies
% url http://www.cise.ufl.edu/research/sparse/SPQR/SPQR/MATLAB/spqr_make.m

% This is exceedingly ugly.  The MATLAB mex command needs to be told where to
% find the LAPACK and BLAS libraries, which is a real portability nightmare.
% The correct option is highly variable and depends on the MATLAB version.

if isempty(varargin)
	if ispc
	    if (verLessThan ('matlab', '6.5'))
	        % MATLAB 6.1 and earlier: use the version supplied in CHOLMOD
	        bl_line = '../../CHOLMOD/MATLAB/lcc_lib/libmwlapack.lib' ;
	    elseif (verLessThan ('matlab', '7.5'))
	        % use the built-in LAPACK lib (which includes the BLAS)
	        bl_line = 'libmwlapack.lib' ;
	    else
	        % need to also use the built-in BLAS lib 
	        bl_line = 'libmwlapack.lib libmwblas.lib' ;
	    end
	else
	    if (verLessThan ('matlab', '7.5'))
	        % MATLAB 7.5 and earlier, use the LAPACK lib (including the BLAS)
	        bl_line = '-lmwlapack' ;
	    else
	        % MATLAB 7.6 requires the -lmwblas option; earlier versions do not
	        bl_line = '-lmwlapack -lmwblas' ;
	    end
	end
	
	if (is64 && ~verLessThan ('matlab', '7.8'))
	    % versions 7.8 and later on 64-bit platforms use a 64-bit BLAS
			% fprintf ('with 64-bit BLAS\n') ;
	    bl_line = ['-DBLAS64 ' bl_line] ;
	end
else
	bl_line = varargin{1};
end

%MKL_STATIC=['-DMKL'...
%			' /opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_lp64.a'...
%      ' /opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_thread.a'...
%		  '	/opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_core.a'...
%			' /opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_lp64.a'...
%      ' /opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_thread.a'...
%		  '	/opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_core.a'...
%			' -L/opt/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64'...
%			' -liomp5 -lpthread -lm'];
% 
%MKL_DYNAMIC = ['-DMKL -I/opt/intel/composer_xe_2011_sp1.11.339/mkl/include/ '...
%			'-L/opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64 '...
%			'-L/opt/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64 '...
%			'-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm'];
%standard = '-DBLAS64 -lmwblas -lmwlapack';

cs_primal = sprintf(['mex %s -DMATLAB -DPRIMAL -DDOUBLE  slp_primal.cpp' ...
		' ../src/tools.cpp ../src/core_primal.cpp %s'], CFLAGS, bl_line);

cs_dual = sprintf(['mex %s -DMATLAB -DDUAL -DDOUBLE %s slp_dual.cpp' ...
		' ../src/tools.cpp ../src/core_dual.cpp %s'], CFLAGS, bl_line);

lc = true;
try
    if ispc
        system('copy slp.cpp slp_primal.cpp');
        system('copy slp.cpp slp_dual.cpp');
    else
        system('cp slp.cpp slp_primal.cpp');
        system('cp slp.cpp slp_dual.cpp');
    end

    fprintf('>> %s\n',cs_primal);eval(cs_primal);
		fprintf('>> %s\n',cs_dual);eval(cs_dual);
catch
    lc = false;
    fprintf('Compilation/linking of slp failed. Try follow the above instruction\n')
    fprintf('to correct the error.\n')
end

if lc
    try
        X = randn(10,5); x = randn(10,1); gamma = 1.0;
        [alpha slp_info] = slp(x,X,gamma);
    catch
        fprintf('There was an error during execution. This probably have something\n')
        fprintf('to do with the linking to the BLAS/LAPACK library and/or the\n')
        fprintf('PATH settings on this computer.\n')		
    end
    
    if exist('slp_info')
        if( slp_info.solved == 0)
            fprintf('The problem was not solved correctly. This probably have something\n')
            fprintf('to do with incompability of the used BLAS/LAPACK library and\n')
            fprintf('the data types used in slp. Try abother BLAS/LAPACK library.\n')
        end
    end
end
