%{
vis.vanGogh2 (manual) # conditions for the van Gogh 2 stimulus
-> vis.Condition
----
-> vis.vanGogh2Cache

pre_blank_period        :decimal(5,3)  #  (seconds)
duration                :decimal(5,3)  #  (seconds)

pattern_width           :smallint      #  pixel size of the resulting pattern
pattern_aspect          :float         #  the aspect ratio of the pattern
pattern_upscale         :tinyint       #  integer upscale factor of the pattern
pattern_aspect          :float         #  the aspect ratio of the pattern

Parameters for the Gabor filter. Filter patch is square.
gaborpatchsize          :smallint      #  size of the gabor filter = gaborpatchsize*pattern_width
gabor_wlscale           :float         #  wavelength scale: lambda = Ny/gabor_wlscale
gabor_envscale          : #  stddev of the Gabor's Gaussian envelope = Ny/cond.gabor_envscale
gabor_ell               : #  specifies the ellipticity of support of the Gabor

Parameters for the spatial Gaussian smoothing filter used to create orientation and
 contrast fields. Filter patch is square.
gaussfilt_scale         : #  size of filter = Ny*gaussfilt_scale
gaussfilt_istd          : #  inverse standard deviation 

Parameters for temporal filtering. We use a 2nd order IIR filter.
y(n) = 2(1-a)*y(n-1) - (1-a)^2*y(n-2) + a*x(n)
filt_noiseBW            : #  BW (in Hz) for filter for noise input 
filt_oriBW              : #  BW (in Hz) for filter for orientation field
filt_contBW             : #  BW (in Hz) for filter for external contrast field
filt_gammshape          : #  shape parameter of gamma distribution to generate inputs for external contrast field
filt_gammscale          : #  scale parameter of gamma distribution

%}

classdef vanGogh2 
    
    methods(Static)
        
        function test()
            fps = 60;
            cond.pre_blank_period   = 5.0;
            cond.noise_seed         = 100;
            cond.pattern_upscale    = 10;
            cond.pattern_width      = 32;
            cond.duration           = 60;
            cond.pattern_aspect     = 1.7;
            cond.gaborpatchsize     = 0.35; 
            cond.gabor_wlscale      = 4;
            cond.gabor_envscale     = 6;
            cond.gabor_ell          = 1;
            cond.gaussfilt_scale    = 2; 
            cond.gaussfilt_istd     = 3; 
            cond.filt_noiseBW       = 0.5;
            cond.filt_oriBW         = 0.5;
            cond.filt_contBW        = 0.5;
            cond.filt_gammshape     = 0.35;
            cond.filt_gammscale     = 2;
            
            % rng(cond.noise_seed)
            
            tic
            img = vanGogh2.make(cond, fps);
            toc
            
            v = VideoWriter('vanGogh2', 'MPEG-4');
            v.FrameRate = fps;
            v.Quality = 100;
            open(v), writeVideo(v, permute(img, [1 2 4 3])); close(v);

        end
        
        function [img, hash] = make(cond, fps)
            
            % video size
            Nframes = round((cond.duration + cond.pre_blank_period)*fps);
            Nx = cond.pattern_width;
            Ny = round(Nx/cond.pattern_aspect);
            
            %  Gabor filter size
            gy = round(cond.pattern_width*cond.gaborpatchsize);
            gx = round(cond.pattern_width*cond.gaborpatchsize);
            
            % Generating the meshgrid for the gabor filter
            ymin = -(gy-1)/2; ymax = (gy-1)/2;
            xmin = -(gx-1)/2; xmax = (gx-1)/2;
            [y,x] = meshgrid(linspace(ymin,ymax,gy),linspace(xmin,xmax,gx));
            
            % Gabor parameters
            lambda  = Ny/cond.gabor_wlscale;    % wavelength
            kappa   = 2*pi/lambda;              % wavenumber
            phi     = 0;                        % phase offset
            sigm    = Ny/cond.gabor_envscale;   % standard deviation of the Gaussian envelope
            gamma   = cond.gabor_ell;           % spatial aspect ratio of the Gabor 
            
            % Gaussian filter for smoothing
            fx          = Ny*cond.gaussfilt_scale;              % size of filter in pixels
            GaussFilt   = gausswin(fx,cond.gaussfilt_istd);     
            GaussFilt   = GaussFilt*GaussFilt'; 
            GaussFilt   = GaussFilt/sum(GaussFilt(:)); 
            
            % Initialize the orientation and contrast fields
            OMap    = zeros(Ny,Nx,Nframes);
            CMap    = zeros(Ny,Nx,Nframes);  
            CMapExt = zeros(Ny,Nx,Nframes);
            
            
            % Filter parameters
            % obtain the 2nd order filter parameter for the given filter BW
            wr_w = 2*pi*cond.filt_noiseBW/fps;
            a_w  = 2 - cos(wr_w) - sqrt((2 - cos(wr_w))^2 - 1);
            
            wr_o = 2*pi*cond.filt_oriBW/fps;
            a_o  = 2 - cos(wr_o) - sqrt((2 - cos(wr_o))^2 - 1);
            
            wr_c = 2*pi*cond.filt_contBW/fps;
            a_c  = 2 - cos(wr_c) - sqrt((2 - cos(wr_c))^2 - 1);
            
            a_w = 1 - a_w;
            a_o = 1 - a_o;
            a_c = 1 - a_c;
            
            % Initialize all the filter inputs
            Wold1 = randn(Ny + gy - 1, Nx + gx -1);  % white noise to be filtered
            Wold2 = randn(Ny + gy - 1, Nx + gx -1);  % white noise to be filtered
            
            OMapRawOld1 = randn(Ny,Nx) + 1i*randn(Ny,Nx);
            OMapRawOld2 = randn(Ny,Nx) + 1i*randn(Ny,Nx);
            
            k = cond.filt_gammshape/(gy*gx*2/a_c); % shape parameter for the gamma distribution 
            % k is proportional to 1/(no. of pixels being filtered in space and time)
            % 2/a_c term is because that is the rough time scale for temporal averaging with the 2nd order IIR filter
            CMapExtOld1 = gamrnd(k,cond.filt_gammscale,[Ny,Nx]);
            CMapExtOld2 = gamrnd(k,cond.filt_gammscale,[Ny,Nx]);
            
            % generate filter for upscaling
            f_up = cond.pattern_upscale;
            kernel_sigma = f_up;
            sz = f_up*[Ny,Nx];
            [fy,fx] = ndgrid(...
                (-floor(sz(1)/2):floor(sz(1)/2-0.5))*2*pi/sz(1), ...
                (-floor(sz(2)/2):floor(sz(2)/2-0.5))*2*pi/sz(2));

            fmask = exp(-(fy.^2 + fx.^2)*kernel_sigma.^2/2);
            fmask = ifftshift(fmask);
            
            % Output
            Y       = zeros(Ny,Nx,Nframes);
            YsVid   = zeros(Ny*f_up,Nx*f_up,Nframes);
            
            
            for tt = 1:Nframes
                
                % Temporal filtering of complex Gaussian noise
                OMapRawNew      = 2*(1-a_o)*OMapRawOld1 - (1-a_o)^2*OMapRawOld2 + a_o*(randn(Ny,Nx) + 1i*randn(Ny,Nx));
                % Spatial filtering to generate orientation field
                OMapFilt        = convn(OMapRawNew,GaussFilt,'same');
                OMap(:,:,tt)    = angle(OMapFilt);
                CMap(:,:,tt)    = abs(OMapFilt);
                
                % Temporal filtering of gamma distributed noise
                CMapExtNew      = 2*(1-a_c)*CMapExtOld1 - (1-a_c)^2*CMapExtOld2 + a_c*gamrnd(k,cond.filt_gammscale,[Ny,Nx]);
                % Spatial filtering to generate external contrast field
                CMapExtFilt     = convn(CMapExtNew,GaussFilt,'same');
                CMapExt(:,:,tt) = CMapExtFilt;
                
                % Temporal filtering of white noise to generate input to
                % the Gabor filters
                Wnew = 2*(1-a_w)*Wold1 - (1-a_w)^2*Wold2 + a_w*randn(Ny + gy - 1, Nx + gx -1);

                % Gabor filtering of temporally filtered noise input W
                for ii = 1:Ny
                    for jj = 1:Nx
                        % centered indices for the input 
                        lx = (gy-1)/2;
                        ly = (gx-1)/2;
                        ic = ii + lx; jc = jj + ly;
                        Wreq = Wnew(ic-lx:ic+ly,jc-lx:jc+ly);

                        theta = OMap(ii,jj,tt);                     % orientation
                        cont  = CMap(ii,jj,tt)*CMapExt(ii,jj,tt);   % contrast

                        yp = y*cos(theta) + x*sin(theta);
                        xp = -y*sin(theta) + x*cos(theta);
                        % Gabor defined at output pixel location (ii,jj) at time tt
                        F = cont*exp(-(yp.^2 + gamma^2*xp.^2)/(2*sigm^2)).*cos(kappa*yp + phi);  
                        Y(ii,jj,tt) = Wreq(:)'*F(:);
                    end
                end
                
                % Filter states updated
                OMapRawOld2 = OMapRawOld1;
                OMapRawOld1 = OMapRawNew;
                
                CMapExtOld2 = CMapExtOld1;
                CMapExtOld1 = CMapExtNew;
                
                Wold2       = Wold1;
                Wold1       = Wnew;
                
                % Upscaling
                YsVid(:,:,tt) = upscale(Y(:,:,tt), fmask, f_up);

            end
            
            Tstart  = 1 + fps*cond.pre_blank_period;
            YsVid   = YsVid(:,:,Tstart:end);
            
            K       = 0.0125; % hardcoded for now. LUT later
            img     = uint8(round(256*(YsVid/2/K + 0.5)));
            

            hash = 0;
        end
        
    end
end


function out = upscale(X, fmask, f_up)

    X = upsample(X', f_up, round(f_up/2))*f_up;
    X = upsample(X', f_up, round(f_up/2))*f_up;
    out = ifft2(fmask.*fft2(X));
       
end
                

