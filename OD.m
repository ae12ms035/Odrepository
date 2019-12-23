% This is a function which calculates the Ocular Dominance Column Width and
% angle distribution of the orientation of columns with respect to V1-V2
% border.
%--------------------------------------------------------------------------
%                               Prerequisites
%--------------------------------------------------------------------------
% MATLAB 2013 (or higher versions) with Image Processing Toolbox
% Make sure you have the following supporting functions(available in the online repository) 
% in the directory(path) you are going to work. The image should also be in
% the same direcory.
% 1. flood (performs flood-filling)
% 2. pathfinder_general (identifies the local long-axis)
% 3. linorfit2 (required for line fitting along the long-axis of the ODC)
% 4. linorfitn (required for line fitting)
% 5. OrthoLine2 (calculates the orathogonal line from long axis of the ODC)
% 6. lengthfinder_newmin (calculates the width)
% 7. OD (the function performing the computation with the help of above functions)
%- -------------------------------------------------------------------------
%                               Function
%--------------------------------------------------------------------------
% function [MeanWidthBoth,STDWidthBoth,DI,W,w1,w2,O,o1,o2,H] = OD(file,CF);
% INPUT ARGUMENTS:
% 1. file : the address of the image file(OD map) to be analyzed
% Eg: If the image address is: C:\Users\Admin\Desktop\A\OD.png, then 
% >> file = 'C:\Users\Admin\Desktop\A\OD.png' 
% 2. CF is the pixel to micron conversion factor(a positve number) for the image. 
% Eg: >> CF = 20.94
% OUTPUT:
% 1.MeanWidthBoth - the mean width of the OD columns from both the eyes
% 2.STDWidthBoth - the standard deviation of the mean width of the OD columns from both the eyes
% 3.DI - the 'Diversity Index' of the OD column pattern
% 4.W - the array of all width measures from both the eyes
% 5.w1 - the array of all width measures from the white columns
% 6.w2 - the array of all width measures from the black columns
% 7.O - the array of all angle measures from both the eyes
% 8.o1 - the array of all angle measures from the white columns
% 9.o2 - the array of all angle measures from the black columns
% 10. H - a 1x6 array with some above information. Elements are: 
%MeanWidthBoth, STDWidthBoth, no. of width measures, mean orientation of OD 
%columns, no. of angle measures and Shannon's entropy of the angle distribution
%respectively
%--------------------------------------------------------------------------
%                               HOW TO USE
%--------------------------------------------------------------------------
% Enter the following in the command window with valid input arguments: 
% >> [MeanWidthBoth,STDWidthBoth,DI,W,w1,w2,O,o1,o2] = OD(file,CF);
%Example:
% >> [MeanWidthBoth,STDWidthBoth,DI,W,w1,w2,O,o1,o2] = OD('C:\Users\Admin\Desktop\A\OD.png',20.95,5,3);
% The program calculates numbers and writes to a folder named 'ImageName Analaysis Results'.
% where ImageName is the name of the image. For the above example it will
% be 'OD Analysis Results'. 
% It will have subfolders 'Eye1' and 'Eye2', with output files, data and graphs
% from eye 1(white columns) and eye 2(black columns) respectively. The
% third subfolder will be 'BothEyes', which has the combined data
% from both the eyes, and the width and angle histograms.
%
% Help the computer with the following steps, as 'OD' starts running. 
%--------------------------------------------------------------------------
%               STEPS WITH USER INPUT REQUIREMENTS
%--------------------------------------------------------------------------
% 0. Choose whether you want Excel outputs from the menu
% 1. Rotating the image for ROI selection: If no rotation is required,
% click 'NO' on the menu shown. If you would like to rotate thev image for 
% ROI selection, click on 'YES', in the menu. Then click on two points
% in the image which you think would be on a line with respect to which
% the image should be rotated. 
%
% 2. ROI Selection: If you want the entire image to be analysed, click
% 'NO', in the shown menu. If you want to study a specific region of the image,
% crop the rectangular region which you would like to analyze. Avoid
% regions with blood vessels, and irregular lighting for a better estimate. 
%
% 3. Edge-detection: You wil be shown three edge-dectected images along with 
% the original image for three sigma(the value of the standard deviation of the 
% gaussian window used for the purpose) values which was used for the process.
% You will be asked to identify the image with best-matching edges using the menu
% shown. Follow the instructions in the command window for choosing the right image. 
% 
% 4. Image Segmentation: You will be shown two figures: Original OD map and
% edge-map with edges(white) in a black bakcgruond. Click on the regions
% corresponding to the white regions in the OD map. This would flood fill
% the particular column. Do NOT click on the edges(or white) regions in the edge map
% If you did, you will be shown a Warning. Continue the process until no column
% (except very small ones) is left. 
%
%--------------------------------------------------------------------------

function [MeanWidthBoth,STDWidthBoth,DI,W,w1,w2,O,o1,o2,H] = OD(file,CF);
close all; clc;
disp('A pop-up MENU will appear on the upper-left corner of the screen. Use the MENU for selcting desired options!');
disp('Do you want an Excel output for this analysis?...Choose from the MENU...');
disp('Excel output is supported only on Windows Platform!');%Works only on Windows platform
xl = menu('Do you want an Excel output for this analysis?','YES','NO');
close all;
%--------------------------------------------------------------------------
%%%%%%%%%%USER INPUTS(Edit if required)
%--------------------------------------------------------------------------
nbins = 10;%No. of bins for the histograms
wlb = 100;%The lower bound for width calculation in microns
wub = 1000;%The upper bound for width calculation in microns
Window = 5;%Window size for computations(odd number > 1)
BrPtWindow = 3;%Window size to be excluded near branch points(odd number > 1)
%--------------------------------------------------------------------------
%End of USER INPUTS
%--------------------------------------------------------------------------
wsize = int16(floor(Window/2));
wsize2 = int16(floor(BrPtWindow/2));
[path,name,ext] = fileparts(file);
ImName = ([name ext]);
InputIm = imread(ImName);
%CF = input('Input the pixel to micron conversion factor: ');
%window = 5;

cd(path);
direc = [name ' Analysis Results'];
mkdir(direc);
cd(direc);
for n2 = 1:2
    mkdir([ 'Eye' num2str(n2)]);
end
cd(path);
LenUB = 500./CF;

%Rotation
%----------------------------------------------------------
%cd(direc);
currentimage  = InputIm;
imshow(currentimage);
disp('---------- STEP 1: Image Rotation for ROI selection (optional) ---------');
disp('Do you want to rotate this image for ROI selection?');
disp('Choose from the menu......');
rr = menu('Do you want to rotate this image for ROI selection?','YES','NO');

if rr == 1
disp('You need to specify two points which you think would be on a line w.r.to which rotation should take place')
disp('Choose two points on the image w.r.to which you would like to rotate the image....');
p = ginput(2);
clc;
angle1 = atan((p(2,2)-p(1,2))./(p(2,1)-p(1,1)))*180/pi;
close all;
RotImage = imrotate(currentimage,angle1,'bilinear','crop');
imshow(RotImage);
elseif rr == 2
    angle1 = 0;
    RotImage = currentimage;
end
clc;
disp('---------- STEP 2: V1-V2 Border Detection ---------');
disp('You need to identify two points lying on the V1-V2 border(on a line)');
disp('Now, Click on two points on the V1/V2 border....');
p2  = ginput(2);
Theta0 = -atan((p2(2,2)-p2(1,2))./(p2(2,1)-p2(1,1)));
if Theta0 < 0
        Theta0 = pi + Theta0;% 
else
        Theta0 = Theta0;%
end
close all;
clc;
imshow(RotImage);

%Cropping  
%------------------------------------------------------------
disp('---------- STEP 3: Cropping the ROI (Optional) ---------');
disp('Do you want to crop this image for ROI selection?');
disp('Choose from the menu......');
cr = menu('Do you want to crop this image for ROI selection?','YES','NO');
if cr == 1
disp('Keep the mouse pressed for selecting the ROI using the selection tool');
disp('Choose a rectangular region(ROI) for cropping....');
rect = getrect;
rect = round(rect);
close all;
clc;
MM = RotImage(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
% Display the subsetted image with appropriate axis ratio
figure; imshow(MM);
elseif cr == 2
  
    MM = RotImage(1:size(RotImage,1),1:size(RotImage,2));
    figure; imshow(MM);
end
        
% Write image to graphics file. 
imwrite(MM,[path,filesep direc filesep 'Eye' num2str(1) filesep 'SelectedROI.bmp']); 
imwrite(MM,[path, filesep direc filesep 'Eye' num2str(2) filesep 'SelectedROI.bmp']);
save([path,filesep direc filesep 'Eye' num2str(1) filesep 'CropXYtest.mat'],'MM','Theta0','angle1');
save([path,filesep direc filesep 'Eye' num2str(1) filesep 'V1V2Angle.mat'],'Theta0');
close all;
cd([path, filesep, direc , filesep 'Eye' num2str(1)]);
%Edge-Detection
%-------------------------------------------
im2filt  =imread('SelectedROI.bmp');
%     im2filt = imcomplement(im2filt);
%figure, imshow(im2filt);
h1 = fspecial('average');
ImFilt = imfilter(im2filt,h1,'symmetric');
im = ImFilt;
%figure, imshow(im);

for sigma = [2:0.25:20]
         ed1 = edge(im,'log',0,sigma);
         ed(1).f = ed1;
         ed2 = edge(im,'log',0,sigma+0.25);
          ed(2).f = ed2;
         ed3 = edge(im,'log',0,sigma+0.5);
           ed(3).f = ed3;
         edg1 = imfuse(im,ed1,'blend');
          edg(1).f = edg1;
          
         edg2 = imfuse(im,ed2,'blend');
          edg(2).f = edg2;
         edg3 = imfuse(im,ed3,'blend');
          edg(3).f = edg3;
%          imshow(ed);
        subplot(2,3,1), imshow(im2filt); title('Original');
        subplot(2,3,2), imshow(im2filt); title('Original');
        subplot(2,3,3), imshow(im2filt); title('Original');
        subplot(2,3,4), imshow(edg1); title(['Sigma = ', num2str(sigma)]);
        subplot(2,3,5), imshow(edg2); title(['Sigma = ', num2str(sigma+0.25)]);
        subplot(2,3,6), imshow(edg3); title(['Sigma = ', num2str(sigma+0.5)]);
        clc;
        disp('---------- STEP 4: Edge-detection ---------');
        disp('You are shown three edge-detected images(second row) of the OD images in the first row, with different "sigma"');
        disp('A best edge-detected image will have edges along the ODC borders with less undulations');  
        disp('_____________________________________________________________________________________');
        disp('Figure with best edges? Choose from the menu.....');
        disp('If none is best edge-detected image, click on "None, go to next iteration"');
        disp('TIP: Use "Maximize" option to see images and edges more clearly!');
        a = menu('Figure with best edges?', ['Sigma = ' num2str(sigma)], ['Sigma = ' num2str(sigma+0.25)], ['Sigma = ' num2str(sigma+0.5)], 'None, go to next iteration');
       
         if a ~= 4
%              ed = load(['ed' num2str(a) '.mat']);
%              edg = ['edg' num2str(a) '.mat'];
             imwrite(ed(a).f,[path filesep direc filesep 'Eye' num2str(1) filesep 'ODCBordersEye1.bmp']);
             imwrite(ed(a).f,[path, filesep direc filesep 'Eye' num2str(2) filesep 'ODCBordersEye2.bmp']);
             imwrite(edg(a).f,[path filesep direc filesep 'Eye' num2str(1) filesep 'ODCwithEdgeEye1.bmp']);
             imwrite(edg(a).f,[path, filesep direc filesep 'Eye' num2str(2) filesep 'ODCwithEdgeEye2.bmp']);
             save([path,filesep direc filesep 'Eye' num2str(1) filesep 'EdgeSigma.mat'],'sigma');
             break
         end

         if a  == 4
             close all;
             clc;
             sigma  = sigma + 0.25;
         end
         close all;
end
close all;
%Semi-automatic region denotation
%----------------------------------------------------------
ImEd = imread('ODCwithEdgeEye1.bmp');
ImEd= im2double(ImEd);
 
ImEdge = imread('ODCBordersEye1.bmp');
%     I = imfuse(ImFilt,ImEdge,'blend');
figure, imshow(ImEd);
    
CC = bwconncomp(ImEdge);
N = CC.NumObjects;
N = (ceil(N/2)+3);
ImEdge = ImEdge(2:end-1,2:end-1);
figure, imshow(ImEdge);
cd(path);
 clc;
%https://stackoverflow.com/questions/14238083/flood-fill-using-matlab
for i = 1:2*N
        disp('---------- STEP 5: Binary Conversion ---------');
       if i == 1
       disp('You are shown two figures for binary conversion using flood-filling:');
       disp('After you have flood-filled a fraction of the total columns, the program will ask you whether you have completed the process');
       disp('If asked choose the right answer from the menu, until you generate correct binary map.');
       disp('_________________________________________________________________________________________');
       disp('Figure(1): The edge-detected image & Figure(2): the edge-map');
       disp('TIP: You may want to move one figure to see the other one!');
       disp('_________________________________________________________________________________________');
       disp('Looking at Figure(1), do the following in Figure(2) ...');
       disp('Click on a point in Figure(2)which corresponds to white column of Figure(1)...');
       NN = ginput(1);
       NN = round(NN);
       elseif i == 2
       disp('You are now shown the flood-filled output image from the previous step');
       disp('OD columns show alternate bright and dark regions in the map');
       disp('We need to reconstruct the binary map using this property of ODCs');
       disp('So, click on the next possible region which should be filled with white');
       disp('_____________________________________________________________________________________');
       disp('TIP 1: No two adjacent column will be filled with white(or black)');
       disp('TIP 2: Use "Maximize" option to see the regions clearly....');
       disp('.......when flood-filling with smaller columns or columns close to the image boundary');
       NN = ginput(1);
       NN = round(NN);
       elseif i>2 
       disp('You are now shown the flood-filled output image from the previous step');
       disp('Click on the next possible region which should be filled with white');
       disp('_____________________________________________________________________________________');
       disp('TIP 1: No two adjacent column will be filled with white(or black)');
       disp('TIP 2: Use "Maximize" option to see the regions clearly....');
       disp('.......when flood-filling smaller columns or columns close to the image boundary');
       NN = ginput(1);
       NN = round(NN);  
       end
       if ImEdge(NN(1,2),NN(1,1)) == 1
           clc;
           disp('WARNING! Do NOT click on edges or filled white regions. Try again');
           disp('Choose a point within a region corresponding to white column...');
           NN = ginput(1);
           NN = round(NN);
       end
     
       img = flood(NN(1,1),NN(1,2),ImEdge,1);
       clc;
       close all;
       figure, imshow(img);
        ImEdge = img;
       if i > floor(0.45*N)
           disp('Did you complete? Press YES or NO in the menu');
           Complete = menu('Did you complete?', 'YES', 'NO');
           
             if Complete == 1
                 clc;
                 
                 
                break
             end
       end
   end
imwrite(ImEdge,[path filesep direc filesep 'Eye' num2str(1) filesep 'FilledODCEye1.bmp']);
imwrite(imcomplement(ImEdge),[path, filesep direc filesep 'Eye' num2str(2) filesep 'FilledODCEye2.bmp']);
close all;
cd([path, filesep direc]);

disp('Width measurement process about to start....');
%Core Computation
%-------------------------------------------------
G = zeros(2,7);

for n2 = 1:2
        cd([path, filesep direc, filesep 'Eye' num2str(n2)]);
        bw = imread(['FilledODCEye' num2str(n2) '.bmp']);%load the binary image
        TotObj = bwconncomp(bw);
        TotObjN = TotObj.NumObjects;
        ObjArea = regionprops(bw,'Area');%calculating individual object areas
        Area = [];
        for ar = 1:TotObjN
            Area(ar,:) = ObjArea(ar,1).Area;
        end
        MeanWArea = mean(Area);%mean area of a column
        ElimArea = 0.1*MeanWArea; %lower bound for area
        bw2 = bwareaopen(bw,round(ElimArea));%delete small objects
        imshowpair(bw,bw2,'montage');
        saveas(gcf,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'Binary OCD before (left) and after (right) small region elimination.bmp']);
        close all;
        imshow(bw2);
        imwrite(bw2,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'Binary OCD after small region elimination.bmp']);
        close all;
        cd ..
        
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%SCAN LENGTH FIXING
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
        UB = LenUB;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%SKELETONIZATION
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        bw = bw2;
        thick = bwmorph(imcomplement(bw),'thicken',inf);
        thick = imcomplement(thick);
        skel = bwmorph(thick,'skel',Inf);
        branch = bwmorph(skel,'branchpoints');
        [bry,brx] = find(branch);
       
        for br = 1:length(brx)
            if (bry(br)>wsize2 && brx(br)>wsize2)
            skel(bry(br)-wsize2:bry(br)+wsize2,brx(br)-wsize2:brx(br)+wsize2) = 0;
            end
        end
        figure, imshow(imfuse(bw,skel,'blend'));
        imwrite(imfuse(bw,skel,'blend'),[path, filesep direc filesep 'Eye' num2str(n2) filesep 'ODCwithSkeletonEye' num2str(n2) '.bmp']);
        close all;
        [yskel,xskel] = find(skel);  %SKELETON COORDINATES
        
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%OBJECT and boundary IDENTIFICATION
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        CC = bwconncomp(bw);
        ObjNum = CC.NumObjects;
        [B,L] = bwboundaries(bw,'holes');%IDENTIFYING THE BOUNDARY COORDINATES
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%OBJECT Labeling with colors - Region Denotation
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      
        L = bwlabel(bw);
        s = regionprops(L,'Centroid');
        figure, imshow(skel);
        hold on
        [B,L,n,A] = bwboundaries(bw,'holes');
        imshow(label2rgb(L, @jet, [.5 .5 .5]));

        hold on
        for k = 1:numel(s)
            c = s(k).Centroid;
            text(c(1), c(2), sprintf('%d' ,k ), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
        end           
        saveas(gcf,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'TaggedODColumnsEye' num2str(n2) '.png']);
        close all;
        cd ..
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%FINDING NEAREST NEIGHBOURS
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\   
        GlobWidth = [];
         Theta = [];
        ObWidth = zeros(length(xskel),ObjNum);
        point = 0;
        for k = 1:length(xskel)
             xsk = xskel(k);
             ysk = yskel(k);
             
             if (xsk>wsize && ysk>wsize && ysk <size(bw,1)-wsize && xsk < size(bw,2)-wsize)
                
                [xway,yway] = pathfinder_general(ysk,xsk,skel,wsize);
                if length(xway)>2 
                [xline,pline,m1] = OrthoLine2(xway,yway,UB);
                fig= figure(20);
%                 fig2 = figure(30);
               
                Theta(k,1) = ysk;
                Theta(k,2) = xsk;
                invtan = -m1;%-
                if invtan < 0
                    Theta(k,3) = pi + invtan-Theta0;% i-(atan(m1)*180/pi);%-ThetaV12(n1,2));
                     if Theta(k,3) < 0
                        Theta(k,3) = pi + Theta(k,3);
                    end
                else
                    Theta(k,3) = invtan-Theta0;%180-(atan(m1)*180/pi);%-ThetaV12(n1,2));
                     if Theta(k,3) < 0
                        Theta(k,3) = pi + Theta(k,3);
                    end
                end
                for ob = 1:ObjNum
                     enclosed_boundaries = find(A(:,ob));
                   
                    BB = B{ob,1};
                    if ~isempty(enclosed_boundaries)
                        
                        cat = vertcat(BB,B{enclosed_boundaries(:),1});
                        BB = cat;
                    end
%                     BB = B{ob,1};%BOUNDARY COORINATE OF OBJECT 'OB'
                    by = BB(:,1);            
                    bx = BB(:,2);
                    
                    [x1,y1,x2,y2,len] = lengthfinder_newmin(xline,pline,xsk,ysk,bx,by,bw,UB);%calculates the length
                    if (len ~=0 && y1 ~= 0 && y2 ~=0 && x1 ~= 0 && x2 ~=0)
                         figure(20), line([x1 x2], [-y1 -y2]);xlim([1 size(bw,2)]);ylim([-size(bw,1) 1]);hold on;
                         plot(xsk,-ysk,'r.');
                        GlobWidth(k,1) = ysk;
                        GlobWidth(k,2) = xsk;
                        GlobWidth(k,3) = y1;
                        GlobWidth(k,4) = x1;
                        GlobWidth(k,5) = y2;
                        GlobWidth(k,6) = x2;
                        GlobWidth(k,7) = len;
                        ObWidth(k,ob) = len;
                        point = point + 1;
                        disp('Figure(20): Skeleton points(red) and the measured widths(blue lines) of the ROI of...');
                        disp(['Image: ' ImName ' - Eye = ' num2str(n2) ', Column: ' num2str(ob) ', # of Points = ' num2str(point)])
                        
                    end
                end
                end
             end
        end
        clc;
        if n2 == 1
        disp('First eye measurements completed ....');
        elseif n2 == 2
        disp('Second eye measurements completed....')
        end
        saveas(fig,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'PlinesEye' num2str(n2) '.png']);
        save([path,filesep direc filesep 'Eye' num2str(n2) filesep 'ObWidth.mat'],'ObWidth');
        wid = nonzeros(GlobWidth(:,7))*CF;
        wid = wid(wid>wlb & wid <wub);
        Ang = Theta(any(Theta,2),:);
        z = Ang(:,3);
        MeanODW = mean(wid);
        StdODW = std(wid);
        summ = sum(wid);
        len = length(wid);
       
        [counts,centers] = hist(wid,nbins);
        
        figure, bar(centers,counts);legend(['ODC Width = ' num2str(MeanODW) '\mum \pm ' num2str(StdODW)]);xlabel('ODC Width (\mu m)'); ylabel('Frequency');
        saveas(gcf,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'WidthDistriHistogramEye' num2str(n2) '.png']);
        close all;
        
        [f,x] = ecdf(wid);
        figure; hist(wid);hold on;plot(x,f*max(counts),'g');legend(['ODC Width = ' num2str(MeanODW) '\mum \pm ' num2str(StdODW)]);xlabel('ODC Width (\mu m)'); ylabel('Frequency');
        saveas(gcf,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'WidthDistriCDFEye' num2str(n2) '.png']);
        close all;
        
%         MeanAngle = angle(sum(exp(i*z)))*180/pi;%
%         MeanAngleSq = angle(sum(exp(i*z.*2)))*180/pi;%
%         VarAngle = abs(sqrt(MeanAngleSq-(MeanAngle).^2));
        MeanAngle = atan2(sum(sin(z)),sum(cos(z)));
 %        MeanAngleSq = atan2(sum(sin(2*z)),sum(cos(2*z)));
%         VarAngle = abs(sqrt(MeanAngleSq-(MeanAngle).^2))*180/pi;
%      VarAngle = 1 - sqrt((sum(sin(z)))^2+(sum(cos(z)))^2)/length(z);
 %     StdAngle = sqrt(-2*log(1-VarAngle))*180/pi;
        
        
        [countsAng,centersAng] = hist(z*180/pi,[9:18:171]);
          countEn = nonzeros(countsAng);
        countEn = countEn/sum(countEn);
        entropy = -(log(countEn.'))*countEn;
        G(n2,1) = 1;
        G(n2,2) = MeanODW;
        G(n2,3) = StdODW;
        G(n2,4) = summ;
        G(n2,5) = len;
        G(n2,6) = MeanAngle*180/pi;
        G(n2,7) = entropy;
        G1 = [G(n2,2:3) G(n2,5:6) length(z) G(n2,7) (exp(G(n2,7))-1)/(nbins-1)]; 
        if xl == 1
            if ispc
                cd([path,filesep direc filesep 'Eye' num2str(n2)]);
                HeaderGlob = {'Y-Skel' 'X-Skel' 'Y-Border1' 'X-Border1' 'Y-Border2' 'X-Border2' 'Width(pixels)' 'Width(microns)'};
                GB = GlobWidth(any(GlobWidth,2),:);
                GBSort = [GB GB(:,7)*CF];
                OutputAll  = [HeaderGlob; num2cell(GBSort)];
                HeaderW = {'Width(microns)'};
                HeaderA = {'Angle(degrees)'};
                OutputSortedW = [HeaderW; num2cell(wid)];
                OutputSortedA = [HeaderA; num2cell(z*180/pi)];
                xlswrite(['CordinatesAndWidths_Eye' num2str(n2) '.xls'],OutputAll);
                xlswrite(['ODCWidths_Eye' num2str(n2) '.xls'],OutputSortedW);
                xlswrite(['ODCAngles_Eye' num2str(n2) '.xls'],OutputSortedA);
                HeaderE = {'MeanWidth(microns)' 'STD' '#of Width measures' 'Mean Angle(deg)' '#of angle measures' 'Shannon entropy' 'Div. Index'};
                OutputE  = [HeaderE; num2cell(G1)];
                xlswrite(['Summary_Eye_' num2str(n2) '.xls'],OutputE);
                cd ..
            end
        end
        if n2 == 1
            w1 = wid;
            o1 = z;
        end
        if n2 == 2
            w2 = wid;
            o2 = z;
        end
        figure, bar(centersAng,countsAng);legend(['ODC Angle Distribution - Mean = ' num2str(MeanAngle*180/pi) ' Entropy' num2str(entropy)]);xlabel('ODC Angle (Degrees)'); ylabel('Frequency');
        saveas(gcf,[path, filesep direc filesep 'Eye' num2str(n2) filesep 'AngleDistiHistogramEye' num2str(n2) '.png']);%xlim([-ThetaV12(n1,2)*180./pi 180+ThetaV12(n1,2)*180/pi]);
        save([path,filesep direc filesep 'Eye' num2str(n2) filesep 'GlobWidth.mat'],'GlobWidth','MeanODW','StdODW','wid','counts','centers','Theta','countsAng','centersAng','MeanAngle','entropy','G'); 
        close all;
        
end
mkdir([path, filesep direc filesep 'BothEyes']);
H = zeros(1,7);
W = [w1 ; w2];
O = [o1 ; o2];
MeanWidthBoth = mean(W);
STDWidthBoth = std(W);
[countsBoth,centersBoth] = hist(W,nbins);
[countsAngBoth,centersAngBoth] = hist(O*180/pi,[9:18:171]);
countEnBoth = nonzeros(countsAngBoth);
countEnBoth = countEnBoth/sum(countEnBoth);
entropyBoth = -(log(countEnBoth.'))*countEnBoth;

H(1,1) = mean(W);
H(1,2) = std(W);
H(1,3) = length(W);
H(1,4) = (atan2(sum(sin(O)),sum(cos(O))))*180/pi;
H(1,5) = length(O);
H(1,6) = entropyBoth;
DI = (exp(entropyBoth)-1)/(nbins-1);
H(1,7) = DI;
figure, bar(centersBoth,countsBoth);legend(['ODC Width = ' num2str(H(1,1)) '\mum \pm ' num2str(H(1,2))]);xlabel('ODC Width (\mu m)'); ylabel('Frequency');
saveas(gcf,[path, filesep direc filesep 'BothEyes' filesep 'WidthDistriHistogramBothEyes.png']);
close all;
     
figure, bar(centersAngBoth,countsAngBoth);legend(['ODC Angle Distribution - Mean = ' num2str(H(1,4)) ' Diversity Index' num2str(DI)]);xlabel('ODC Angle (Degrees)'); ylabel('Frequency');
saveas(gcf,[path, filesep direc filesep 'BothEyes' filesep 'AngleDistriHistogramBothEyes.png']); close all;
cd(path);
save([path,filesep direc filesep 'BothEyes' filesep 'CombinedWidth.mat'],'O','W','entropyBoth','H','DI'); 
if xl == 1
        cd([path,filesep direc filesep 'BothEyes']);
        if ispc
        
            HeaderBoth = {'MeanWidth(microns)' 'STD' '#of Width measures' 'Mean Angle(deg)' '#of angle measures' 'Shannon entropy' 'Div. Index'};
            OutputBoth  = [HeaderBoth; num2cell(H(any(H,2),:))];
            xlswrite(['Summary(BothEyes)_' name '.xls'],OutputBoth);
            HeaderWidth = {'Width(both, microns)'};
            OutputWidth  = [HeaderWidth; num2cell(W)];
            xlswrite(['ODCWidthMeasues(BothEyes)_' name '.xls'],OutputWidth);
            HeaderAngle = {'Angle(both, degrees)'};
            OutputAngle  = [HeaderAngle; num2cell(O*180/pi)];
            xlswrite(['ODCAngleMeasues(BothEyes)_' name '.xls'],OutputAngle);
        elseif ismac
            MeanWidthMicrons = H(1,1);
            STD = H(1,2);
            NoOfWidthMeasures = H(1,3);
            MeanAngleDeg = H(1,4);
            NoOfAngleMeasures = H(1,5);
            ShannonEntropy = H(1,6);
            DivIndex = H(1,7);
            Measure = {'Value'};
            T = table(Measure, MeanWidthMicrons,  STD, NoOfWidthMeasures, MeanAngleDeg, NoOfAngleMeasures,  ShannonEntropy, DivIndex);
            writetable(T, ['ODCAngleMeasues(BothEyes)_' name '.txt'] );
        end
end
disp('Information combined from both eyes. Calculations completed!')
cd(path);
end