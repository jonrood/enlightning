% ===========================================================================
% Copyright (C) 2017 Jon Rood.
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

function LightningChannel(domainW,domainH,startX,startY,...
    maxIterations,branchLevels,meanAngle,write,overwrite,N,...
    branchyness,meanBranchLength,numAvgAngles)

    %domainW,domainH are width and height of domain
    %startX,startY are x,y start of lightning channel
    %maxIterations is maximum amount of main branch segment iterations
    %branchLevels is level of recursive branching
    %meanAngle is mean of random angles
    %write==1 write to file, write==0 don't write to file
    %overwrite==1 overwrite file, overwrite==0 append to file
    %N=20 from Roy and Ribner; bias towards vertical
    %branchyness is in [0,1], lower numbers mean less branching
    %meanBranchLength is average length of branches
    %numAvgAngles is number of previous angles to use in channel smoothing
    
    x=zeros(1,2);
    y=zeros(1,2);
    i=0;
    x(1)=startX;
    y(1)=startY;
    angles=zeros(numAvgAngles,1);
    
    if(overwrite==1)
        fid=fopen('lightning.txt', 'w');
        fclose(fid);
    end
    
    while((y(1)>0) && (i<maxIterations))
        theta=GetAngle(meanAngle,angles,N);
        
        length=GetLength();
        
        for n=numAvgAngles-1:-1:1
            angles(n+1,1)=angles(n,1);
        end
        
        angles(1,1)=theta;
        
        y(2)=y(1)-length*sind(90-theta);
        x(2)=x(1)+length*cosd(90-theta);
        
        if(y(2)<0)
            y(2)=0;
        end
        
        if(write==1)
            fid=fopen('lightning.txt', 'a');
            fprintf(fid, '%f %f %f %f\n',x(1),y(1),x(2),y(2));
            fclose(fid);
        end
        
        line(x, y);
        axis([0 domainW 0 domainH]);
        x(1)=x(2);
        y(1)=y(2);
        i=i+1;
        
        if(branchLevels>0)
            if(rand(1,1)<branchyness)
                branchMaxIterations=...
                    GetMaxIterations(meanBranchLength/2,branchLevels);
                LightningChannel(domainW,domainH,x(2),y(2),...
                    branchMaxIterations,branchLevels-1,meanAngle,...
                    write,0,N,branchyness/10,meanBranchLength/2,...
                    numAvgAngles);
            end
        end
    end
end

function iter=GetMaxIterations(meanBranchLength,branchLevels)
    iter=ceil(((meanBranchLength*branchLevels)/2)*(1+randn(1,1)));
end

function angle=GetAngle(meanAngle,angles,N)
    angle=-180:180;
    P=exp((-(angle-mean(angles)+mean(angles)/N).^2)./(meanAngle^2));
    P=P./trapz(angle,P(:));
    anglei=linspace(min(angle),max(angle),10000)';
    Pi=interp1(angle,P,anglei);
    cdf=cumtrapz(anglei,Pi);
    i=[true; not(diff(cdf)==0)];
    cdf=cdf(i);
    anglei=anglei(i);
    angle=interp1(cdf,anglei,rand(1,1));
end

function length=GetLength()
    length=3;
end
