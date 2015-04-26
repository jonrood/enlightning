% ===========================================================================
% Copyright (C) 2015 Jon Rood.
%
% This file is part of Enlightning source code.
% 
% Enlightning source code is free software; you can redistribute it
% and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 3 of the License,
% or (at your option) any later version.
% 
% Enlightning source code is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
% ===========================================================================

function DrawSource(domainW,domainH,X,Y,alpha,srcAmp)

    %domainW,domainH are width and height of domain
    %X,Y are resolution of grid
    %alpha is the spread of the source
    %srcAmp is the amplitude of the source
    
    dx=domainW/X;
    dy=domainH/Y;
    H0=zeros(X,Y);
    
    fid = fopen('lightning.txt');
    S = fscanf(fid, '%f %f %f %f', [4 inf]);
    fclose(fid);
    
    for k=1:size(S,2)
        theta=-atan((S(4,k)-S(2,k))/(S(3,k)-S(1,k)));

        if((S(1,k)>S(3,k) && S(2,k)>S(4,k)) ||...
                (S(1,k)<=S(3,k) && S(2,k)>=S(4,k)))
            orientation=1;
        elseif((S(1,k)>=S(3,k) && S(2,k)<=S(4,k)) ||...
                (S(1,k)<S(3,k) && S(2,k)<S(4,k)))
            orientation=2;
        end

        for j=1:Y
            yWorld=dy*(j-0.5);
            for i=1:X
                xWorld=dx*(i-0.5);
                line0=tan(pi/2-theta)*(xWorld-S(1,k))+S(2,k);
                line1=tan(pi/2-theta)*(xWorld-S(3,k))+S(4,k);
                if((orientation==1 && yWorld>=line0)...
                        || (orientation==2 && yWorld<=line0))
                    fun=srcAmp*exp(-(log(2)/alpha^2)*...
                        ((xWorld-S(1,k))^2+(yWorld-S(2,k))^2));
                end
                if((orientation==1 && yWorld<=line1)...
                        || (orientation==2 && yWorld>=line1))
                    fun=srcAmp*exp(-(log(2)/alpha^2)*...
                        ((xWorld-S(3,k))^2+(yWorld-S(4,k))^2));
                end
                if((orientation==1 && yWorld<line0 && yWorld>line1)...
                        || (orientation==2 && yWorld>line0 &&...
                        yWorld<line1))
                    fun=srcAmp*exp(-(log(2)/alpha^2)*...
                        (sin(theta)*(xWorld-S(1,k))...
                        +cos(theta)*(yWorld-S(2,k)))^2);
                end
                H0(i,j)=H0(i,j)+fun*(1-(H0(i,j)/srcAmp));
            end
        end
    end
    %surf(H0);
    %shading interp;
    %axis([1 X 1 Y -0.001 srcAmp]);
    imagesc(rot90(H0));
    %axis off;
end