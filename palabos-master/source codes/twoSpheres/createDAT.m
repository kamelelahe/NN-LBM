function createDAT(numFiles,path,baseInput,baseOutput)

tic
basename = [path '/' baseInput]; % base name of the bmp files
fid = fopen(baseOutput, 'w');    % open the output file to write in



%%%%%%%%%%%%%%%%%%%%% INLET SLICE %%%%%%%%%%%%%%%%%%%%%%%%
fname = [basename num2str(1, '%0.4i') '.bmp'];

fnameu = [basename num2str(2, '%0.4i') '.bmp'];

BB=imread(fname ,'BMP');
CC=imread(fnameu ,'BMP');
nx=size(BB,2)
ny=size(BB,1)
B=zeros(ny,nx);

wholeGeom=zeros(ny,nx,2);


wholeGeom(:,:,1)=BB;
wholeGeom(:,:,2)=CC;

indexMin=find(wholeGeom==0);
indexMax=find(wholeGeom>0);

wholeGeom(indexMin)=255;
wholeGeom(indexMax)=0;


%i
rA = circshift(wholeGeom,[0,1, 0]);
lA = circshift(wholeGeom,[0,-1, 0]);
fA = circshift(wholeGeom,[1,0, 0]);
bA = circshift(wholeGeom,[-1,0,0]);

rfA = circshift(wholeGeom,[1,1,0]);
rbA = circshift(wholeGeom,[-1,1,0]);
lfA = circshift(wholeGeom,[1,-1,0]);
lbA = circshift(wholeGeom,[-1,-1,0]);


%i+1
uA=circshift(wholeGeom,[0,0, 1]);

urA = circshift(wholeGeom,[0,1, 1]);
ulA = circshift(wholeGeom,[0,-1, 1]);
ufA = circshift(wholeGeom,[1,0, 1]);
ubA = circshift(wholeGeom,[-1,0,1]);

urfA = circshift(wholeGeom,[1,1,1]);
urbA = circshift(wholeGeom,[-1,1,1]);
ulfA = circshift(wholeGeom,[1,-1,1]);
ulbA = circshift(wholeGeom,[-1,-1,1]);

for i=1:nx
        for j=1:ny
            
            if (wholeGeom(j,i,1) == 255 && rA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && lA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && fA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && bA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,1) == 255 && rfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && rbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && lfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && lbA(j,i)==0)
                B(j,i) = 1;
            
                
            
            elseif (wholeGeom(j,i,1) == 255 && uA(j,i)==0)
                B(j,i) = 1;
                
            elseif (wholeGeom(j,i,1) == 255 && urA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && ulA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && ufA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && ubA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,1) == 255 && urfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && urbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && ulfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,1) == 255 && ulbA(j,i)==0)
                B(j,i) = 1;
                
                
            elseif (wholeGeom(j,i,1) == 0)
                B(j,i) = 0;
            else
                B(j,i)=2;

            end

        end
end

image(30*B)
axis equal
drawnow

'printing first slice'
fprintf(fid, '%i\n', B);


%%%%%%%%%%%%%%%%%%%%% INTERNAL SLICES %%%%%%%%%%%%%%%%%%%%%%%%

for ii=2:numFiles-1
    ii
fname = [basename num2str(ii, '%0.4i') '.bmp'];
fnamed =[basename num2str(ii-1, '%0.4i') '.bmp'];
fnameu = [basename num2str(ii+1, '%0.4i') '.bmp'];
AA=imread(fnamed ,'BMP');
BB=imread(fname ,'BMP');
CC=imread(fnameu ,'BMP');

wholeGeom=zeros(ny,nx,3);

wholeGeom(:,:,1)=AA;
wholeGeom(:,:,2)=BB;
wholeGeom(:,:,3)=CC;

indexMin=find(wholeGeom==0);
indexMax=find(wholeGeom>0);

wholeGeom(indexMin)=255;
wholeGeom(indexMax)=0;

%i
rA = circshift(wholeGeom,[0,1, 0]);
lA = circshift(wholeGeom,[0,-1, 0]);
fA = circshift(wholeGeom,[1,0, 0]);
bA = circshift(wholeGeom,[-1,0,0]);

rfA = circshift(wholeGeom,[1,1,0]);
rbA = circshift(wholeGeom,[-1,1,0]);
lfA = circshift(wholeGeom,[1,-1,0]);
lbA = circshift(wholeGeom,[-1,-1,0]);


%i-1
dA=circshift(wholeGeom,[0,0, -1]);

drA = circshift(wholeGeom,[0,1, -1]);
dlA = circshift(wholeGeom,[0,-1, -1]);
dfA = circshift(wholeGeom,[1,0, -1]);
dbA = circshift(wholeGeom,[-1,0,-1]);

drfA = circshift(wholeGeom,[1,1,-1]);
drbA = circshift(wholeGeom,[-1,1,-1]);
dlfA = circshift(wholeGeom,[1,-1,-1]);
dlbA = circshift(wholeGeom,[-1,-1,-1]);

%i+1
uA=circshift(wholeGeom,[0,0, 1]);

urA = circshift(wholeGeom,[0,1, 1]);
ulA = circshift(wholeGeom,[0,-1, 1]);
ufA = circshift(wholeGeom,[1,0, 1]);
ubA = circshift(wholeGeom,[-1,0,1]);

urfA = circshift(wholeGeom,[1,1,1]);
urbA = circshift(wholeGeom,[-1,1,1]);
ulfA = circshift(wholeGeom,[1,-1,1]);
ulbA = circshift(wholeGeom,[-1,-1,1]);

for i=1:nx
        for j=1:ny
            if (wholeGeom(j,i,2) == 255 && rA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && fA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && bA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,2) == 255 && rfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && rbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lbA(j,i)==0)
                B(j,i) = 1;
            
            
                
            elseif (wholeGeom(j,i,2) == 255 && dA(j,i)==0)
                B(j,i) = 1;
                
            elseif (wholeGeom(j,i,2) == 255 && drA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dbA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,2) == 255 && drfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && drbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlbA(j,i)==0)
                B(j,i) = 1;
                
            
            elseif (wholeGeom(j,i,2) == 255 && uA(j,i)==0)
                B(j,i) = 1;
                
            elseif (wholeGeom(j,i,2) == 255 && urA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && ulA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && ufA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && ubA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,2) == 255 && urfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && urbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && ulfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && ulbA(j,i)==0)
                B(j,i) = 1;
                
                
            elseif (wholeGeom(j,i,2) == 0)
                B(j,i) = 0;
            else
                B(j,i)=2;

            end
  
        end
end

image(30*B)
axis equal
drawnow

 fprintf(fid, '%i\n', B);
end






%%%%%%%%%%%%%%%%%%%%% OUTLET SLICES %%%%%%%%%%%%%%%%%%%%%%%%

fname = [basename num2str(numFiles, '%0.4i') '.bmp'];
fnamed =[basename num2str(numFiles-1, '%0.4i') '.bmp'];

AA=imread(fnamed ,'BMP');
BB=imread(fname ,'BMP');


wholeGeom=zeros(ny,nx,2);

wholeGeom(:,:,1)=AA;
wholeGeom(:,:,2)=BB;

indexMin=find(wholeGeom==0);
indexMax=find(wholeGeom>0);

wholeGeom(indexMin)=255;
wholeGeom(indexMax)=0;

%i
rA = circshift(wholeGeom,[0,1, 0]);
lA = circshift(wholeGeom,[0,-1, 0]);
fA = circshift(wholeGeom,[1,0, 0]);
bA = circshift(wholeGeom,[-1,0,0]);

rfA = circshift(wholeGeom,[1,1,0]);
rbA = circshift(wholeGeom,[-1,1,0]);
lfA = circshift(wholeGeom,[1,-1,0]);
lbA = circshift(wholeGeom,[-1,-1,0]);


%i-1
dA=circshift(wholeGeom,[0,0, -1]);

drA = circshift(wholeGeom,[0,1, -1]);
dlA = circshift(wholeGeom,[0,-1, -1]);
dfA = circshift(wholeGeom,[1,0, -1]);
dbA = circshift(wholeGeom,[-1,0,-1]);

drfA = circshift(wholeGeom,[1,1,-1]);
drbA = circshift(wholeGeom,[-1,1,-1]);
dlfA = circshift(wholeGeom,[1,-1,-1]);
dlbA = circshift(wholeGeom,[-1,-1,-1]);



for i=1:nx
        for j=1:ny
            if (wholeGeom(j,i,2) == 255 && rA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && fA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && bA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,2) == 255 && rfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && rbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && lbA(j,i)==0)
                B(j,i) = 1;
            
            
                
            elseif (wholeGeom(j,i,2) == 255 && dA(j,i)==0)
                B(j,i) = 1;
                
            elseif (wholeGeom(j,i,2) == 255 && drA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dbA(j,i)==0)
                B(j,i) = 1;
            
            elseif (wholeGeom(j,i,2) == 255 && drfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && drbA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlfA(j,i)==0)
                B(j,i) = 1;
            elseif (wholeGeom(j,i,2) == 255 && dlbA(j,i)==0)
                B(j,i) = 1;
                
           
                
            elseif (wholeGeom(j,i,2) == 0)
                B(j,i) = 0;
            else
                B(j,i)=2;

            end

        end
end

image(30*B)
axis equal
drawnow

'printing last slice'
fprintf(fid, '%i\n', B);


fclose(fid);

toc