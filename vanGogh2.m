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
gaborpatchsize          :smallint      #  pixel size of the gabor filter
gabor_wlscale           :float         #  wavelength scale: lambda = Px/gabor_wlscale
gabor_envscale          : #  stddev of the Gabor's Gaussian envelope = Px/cond.gabor_envscale
gabor_ell               : #  specifies the ellipticity of support of the Gabor

Parameters for the spatial Gaussian smoothing filter used to create orientation and
 contrast fields. Filter patch is square.
gaussfilt_scale         : #  size of filter = Px*gaussfilt_scale
gaussfilt_istd          : #  inverse standard deviation 

Parameters for temporal filtering. We use a 2nd order IIR filter.
y(n) = 2(1-a)*y(n-1) - (1-a)^2*y(n-2) + a*x(n)
filt_noise              : #  filter for noise input
filt_orientation        : #  filter for orientation field
filt_contrast           : #  filter for external contrast field
filt_gammshape          : #  shape parameter of gamma distribution to generate inputs for external contrast field
filt_gammscale          : #  scale parameter of gamma distribution

%}

classdef vanGogh2 
    
    methods(Static)
        
        function test()
            fps = 60;
            cond.pre_blank_period   = 5.0;
            cond.noise_seed         = 100;
            cond.pattern_upscale    = 8;
            cond.pattern_width      = 32;
            cond.duration           = 30;
            cond.pattern_aspect     = 1.7;
            cond.gaborpatchsize     = 11; 
            cond.gabor_wlscale      = 4;
            cond.gabor_envscale     = 6;
            cond.gabor_ell          = 1;
            cond.gaussfilt_scale    = 2; 
            cond.gaussfilt_istd     = 3; 
            cond.filt_noise         = 0.06;
            cond.filt_orientation   = 0.06;
            cond.filt_contrast      = 0.06;
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
            Py = cond.pattern_width;
            Px = round(Py/cond.pattern_aspect);
            
            %  Gabor filter size
            gx = cond.gaborpatchsize;
            gy = cond.gaborpatchsize;
            
            % Generating the meshgrid for the gabor filter
            xmin = -(gx-1)/2; xmax = (gx-1)/2;
            ymin = -(gy-1)/2; ymax = (gy-1)/2;
            [x,y] = meshgrid(linspace(xmin,xmax,gx),linspace(ymin,ymax,gy));
            
            % Gabor parameters
            lambda  = Px/cond.gabor_wlscale;    % wavelength
            kappa   = 2*pi/lambda;              % wavenumber
            phi     = 0;                        % phase offset
            sigm    = Px/cond.gabor_envscale;   % standard deviation of the Gaussian envelope
            gamma   = cond.gabor_ell;           % spatial aspect ratio of the Gabor 
            
            % Gaussian filter for smoothing
            fx          = Px*cond.gaussfilt_scale;              % size of filter in pixels
            GaussFilt   = gausswin(fx,cond.gaussfilt_istd);     
            GaussFilt   = GaussFilt*GaussFilt'; 
            GaussFilt   = GaussFilt/sum(GaussFilt(:)); 
            
            % Initialize the orientation and contrast fields
            OMap    = zeros(Px,Py,Nframes);
            CMap    = zeros(Px,Py,Nframes);  
            CMapExt = zeros(Px,Py,Nframes);
            
            % Output
            Y = zeros(Px,Py,Nframes);
            
            % Filter parameters
            a_w = cond.filt_noise;
            a_o = cond.filt_orientation;
            a_c = cond.filt_contrast;
            
            % Initialize all the filter inputs
            Wold1 = randn(Px + gx - 1, Py + gy -1);  % white noise to be filtered
            Wold2 = randn(Px + gx - 1, Py + gy -1);  % white noise to be filtered
            
            OMapRawOld1 = randn(Px,Py) + 1i*randn(Px,Py);
            OMapRawOld2 = randn(Px,Py) + 1i*randn(Px,Py);
            
            k = cond.filt_gammshape/(gx*gy*2/a_c); % shape parameter for the gamma distribution 
            % k is proportional to 1/(no. of pixels being filtered in space and time)
            % 2/a_c term is because that is the rough time scale for temporal averaging with the 2nd order IIR filter
            CMapExtOld1 = gamrnd(k,cond.filt_gammscale,[Px,Py]);
            CMapExtOld2 = gamrnd(k,cond.filt_gammscale,[Px,Py]);
            
            
            for tt = 1:Nframes
                
                % Temporal filtering of complex Gaussian noise
                OMapRawNew      = 2*(1-a_o)*OMapRawOld1 - (1-a_o)^2*OMapRawOld2 + a_o*(randn(Px,Py) + 1i*randn(Px,Py));
                % Spatial filtering to generate orientation field
                OMapFilt        = convn(OMapRawNew,GaussFilt,'same');
                OMap(:,:,tt)    = angle(OMapFilt);
                CMap(:,:,tt)    = abs(OMapFilt);
                
                % Temporal filtering of gamma distributed noise
                CMapExtNew      = 2*(1-a_c)*CMapExtOld1 - (1-a_c)^2*CMapExtOld2 + a_c*gamrnd(k,cond.filt_gammscale,[Px,Py]);
                % Spatial filtering to generate external contrast field
                CMapExtFilt     = convn(CMapExtNew,GaussFilt,'same');
                CMapExt(:,:,tt) = CMapExtFilt;
                
                % Temporal filtering of white noise to generate input to
                % the Gabor filters
                Wnew = 2*(1-a_w)*Wold1 - (1-a_w)^2*Wold2 + a_w*randn(Px + gx - 1, Py + gy -1);

                % Gabor filtering of temporally filtered noise input W
                for ii = 1:Px
                    for jj = 1:Py
                        % centered indices for the input 
                        lx = (gx-1)/2;
                        ly = (gy-1)/2;
                        ic = ii + lx; jc = jj + ly;
                        Wreq = Wnew(ic-lx:ic+ly,jc-lx:jc+ly);

                        theta = OMap(ii,jj,tt);                     % orientation
                        cont  = CMap(ii,jj,tt)*CMapExt(ii,jj,tt);   % contrast

                        xp = x*cos(theta) + y*sin(theta);
                        yp = -x*sin(theta) + y*cos(theta);
                        % Gabor defined at output pixel location (ii,jj) at time tt
                        F = cont*exp(-(xp.^2 + gamma^2*yp.^2)/(2*sigm^2)).*cos(kappa*xp + phi);  
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

            end
            
            Tstart  = 1 + fps*cond.pre_blank_period;
            YsVid   = Y(:,:,Tstart:end);
            YsVid   = imresize(YsVid,cond.pattern_upscale,'lanczos3');
            % disp(std(YsVid(:)));
            % K       = 10*std(YsVid(:));
            K       = 12*0.00225; % hardcoded for now. LUT later
            img     = uint8(round(256*(YsVid/2/K + 0.5)));
             
            hash = 0;
        end
        
    end
end

