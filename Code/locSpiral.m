function [M W] = locSpiral(map,mask,pixel_size,minRes,maxRes,threshold,bw)
%locSpiral - method to local cryo-EM map enhancement
%
%   Method to improve the interpretability of cryo-EM maps affected by
%   inhomogeneous resolutions.
%
%   INPUT PARAMETERS:
%   map: cryo-EM map to be processed. This maps should be the average of
%   the half maps obtained after refinment or the final map without any
%   filtering postprocessing (masking, b-factor correction).
%
%   mask: Tight binary map that differenciates between the macromolecule (1)
%   and the background (0).
%
%   pixel_size: sampling rate of the cryo-EM map in Angstroms.
%
%   minRes: minimun resolution of the sweeping resolution range in Angstroms.
%   Usually a value of 25-35 Angstroms is good.
%
%   maxRes: maximum resolution of the sweeping resolution range. Provide here 
%   the obtained FSC global resolution value in Angstroms.
%
%   threshold: noise signficance in the comparision with noise. Excellent  
%   values are 0.9 or 0.95.
%
%   bw: bandwith of the bandpass filter applied in Fourier Space. For high
%   quality maps with global resolution (<4 Angstroms) and showing homogeneous 
%   distribution of resolutions a value of 3.5-5 is good. For maps presenting 
%   heterogeneous distributions of resolutions, or worse global resolution, 
%   values between 6-8 provide excellents results.
%
%   OUTPUT PARAMETERS:
%   M: Improved output map
%   W: Sum of all amplitude maps during the resolution sweeping that can be
%   consider as an quality map representing the amount of signal
%
%   If this method was useful to you, please cite the following paper:
%   S Kaur, J Gomez-Blanco, A Khalifa, S Adinarayanan, R Sanchez-Garcia,
%   D Wrapp, JS Mclellan, KH Bui, J Vargas, "Local computational methods to 
%   improve the interpretability and analysis of cryo-EM maps" BioRxiv, 
%   doi: https://doi.org/10.1101/2020.05.11.088013 (2020)
%
%   Javier Vargas, 01/06/20
%   jvargas@fis.ucm.es
%   Copyright 2020, Universidad Complutense de Madrid
%   $ Revision: 1.0.0.0
%   $ Date: 01/06/20
%
%Copyright 2020 Javier Vargas @UCM
%
%Redistribution and use in source and binary forms, with or without modification, 
%are permitted provided that the following conditions are met:
%
%1. Redistributions of source code must retain the above copyright notice, 
%this list of conditions and the following disclaimer.
%
%2. Redistributions in binary form must reproduce the above copyright notice, 
%this list of conditions and the following disclaimer in the documentation 
%and/or other materials provided with the distribution.
%
%3. Neither the name of the copyright holder nor the names of its contributors 
%may be used to endorse or promote products derived from this software without 
%specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
%USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

%The eco option is an aproximation for the spiral phase transform that
%works fine in all the cases tested and saves a lot of computational
%resources. If case you use a computer with limited RAM memory (<32Gb and
%box sizes smaller than 400x400x400px) use the eco option. If not, it is
%more appropiate to use eco = false.
eco = true;

M = map.*0;
W = map.*0;

[n n1 n2]=size(map);
maxFreq = floor(pixel_size*n/maxRes);
minFreq = floor(pixel_size*n/minRes);

% We create necessary staff for the bank of filters once:
[NR, NC, NZ]=size(map);
[u,v,w]=meshgrid(1:NC, 1:NR, 1:NZ);

%Temporal variables
u0=floor(mean(u(:))); v0=floor(mean(v(:))); w0=floor(mean(w(:)));
u=u-u0; v=v-v0; w=w-w0;

%Circular mask to be used in the quantile evaluation
maskC = sqrt(u.^2+v.^2+w.^2);
maskC = (maskC < u0);

clear u0 v0 w0;

[Theta Phi Radial] = cart2sph(u,v,w);
Hr=1-exp(-(u.^2+v.^2+w.^2)/(2*0.05^2)); %Gaussian DC filter with sigma=0.05

clear u v w Theta Phi;

%We obtain the Spiral Filter once:
if (eco)
    %Using the eco option and the approximation for the spiral phase
    %transform.
    Hs= spiralFilter(map);
    Hsx = nan;
    Hsy = nan;
    Hsz = nan;
else
    %Using the vectorized spiral:
    [Hsx Hsy Hsz]= spiralFilterVect(map);
    Hs = nan;
end

%We do one FFT of the map
C=fftn(map);
num = abs(maxFreq-minFreq);
freqs = linspace(minFreq,maxFreq,num);

parfor i=1:length(freqs)
    % Actual filter
    H=exp(-((Radial-freqs(i)).^2)/(2*bw^2)); %Annular filter with radius=R and sigma=S
    H = Hr.*H;
    H = H/sum(H(:));
    CH=C.*ifftshift(H);    
    ch=real(ifftn(CH));
    
    if (eco)
        %Using the eco option and the approximation for the spiral phase
        %transform.
        CH=CH.*ifftshift(Hs);
        s=abs(conj(ifftn(CH)));
    else
        %Using the vectorized spiral phase transform:
        CHx=((ifftn(CH.*ifftshift(Hsx))));
        CHx = CHx.*conj(CHx);
        CHy=((ifftn(CH.*ifftshift(Hsy))));
        CHy = CHy.*conj(CHy);
        CHz=((ifftn(CH.*ifftshift(Hsz))));
        CHz = CHz.*conj(CHz);
        s = sqrt(CHx+CHy+CHz);
    end
    
    %normalized igram
    cn=cos(atan2(s,ch));
    %modulation
    m=abs(ch+1i*s);
    
    q = quantile(m(((~mask).*maskC)>0.5),threshold);
    Cref = m./(m+q);

    m = double(m >= q);
    
    M = M+(cn.*m.*Cref);
    W = W+(m.*Cref);
    
end

M = M/maxFreq;
W = W/maxFreq;

M = M - mean(M(((~mask).*maskC)>0.5));
M = M ./ std(M(((~mask).*maskC)>0.5));

W = W - mean(W(((~mask).*maskC)>0.5));
W = W ./ std(W(((~mask).*maskC)>0.5));

%We impose posstivity
M(M<0) = 0;
end


function [Hu, Hv, Hw] = spiralFilterVect(c)

[NR, NC, NZ]=size(c);
[u,v,w]=meshgrid(1:NC, 1:NR, 1:NZ);

u0=floor(mean(u(:))); v0=floor(mean(v(:))); w0=floor(mean(w(:)));

u=u-u0;
v=v-v0;
w=w-w0;

Hu=(-1i*u)./sqrt((u.^2)+(v.^2)+(w.^2));
Hu(isnan(Hu))=0;

Hv=(-1i*v)./sqrt((u.^2)+(v.^2)+(w.^2));
Hv(isnan(Hv))=0;

Hw=(-1i*w)./sqrt((u.^2)+(v.^2)+(w.^2));
Hw(isnan(Hw))=0;

end

function H= spiralFilter(c)

[NR, NC, NZ]=size(c);
[u,v,w]=meshgrid(1:NC, 1:NR, 1:NZ);

u0=floor(mean(u(:))); v0=floor(mean(v(:))); w0=floor(mean(w(:)));

u=u-u0;
v=v-v0;
w=w-w0;

H=(u+1i*v+w)./abs(u+1i*v+w);
H(isnan(H))=0;

end