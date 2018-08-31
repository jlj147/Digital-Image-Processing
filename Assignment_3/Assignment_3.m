%Jayce Jones
%EE 4323 Spring 2018
%Assignment_3

%This program creates a UI that allows filters to be applied to images and
%for the brightness and contrast of an image to be adjusted
%Select an image to edit by using the open image button
%Once the image is loaded a current image and preview image will be seen
%Image brightness can be adjusted by using the slider labled brightness
%Image contrast can be adjusted by using the slider labled contrast
%Different image enhancements can be chosen in the enhancement drop down
%After using a slider or picking an enhancement, changes made to the image
%can be viewed by first clicking 'Set Enhancement' and then by clicking
%'Preview Image'
%At this point the preview image will change to show the changes made to
%the current image, if one wishes to keep the changes the 'Apply Changes'
%button can be clicked
%Clicking the 'Apply Changes' button assigns the preview image to the current image
%Any enhancements set will now be applied to the current image and will appear in the
%preview image window
%To retrieve the original image at anytime the no filter option can be
%chosen from the enhancement drop down menu by following the steps as
%described above
%To save the current image the 'Save Image' button can be clicked

classdef Assignment_3 < handle
    
    %Defines properties for the image read in
    properties (GetAccess='private', SetAccess='private')
        original = []; %Stores the original image read in from the file
        current = []; %Stores the current image being manipulated
        preview = []; %Stores the temporary preview image after enhancements have been made
        choice = 0; %Variable for switch case choice
        brightness = .5; %Defines a brightness value for an image read in
        contrast = .5; %Defines a contrast value for an image read in
        enhance = 0; %Boolean value that changes when enahnce button is clicked
    end
   
methods
    %Function that creates the user interface
    function this=Assignment_3
        
        %Create figure window with the title Project #2 - Jayce Jones
        figure('Name','Final Project-Jayce Jones','NumberTitle','off');

        %Create a subplot to display the "current" image
        subplot(2,2,1);
        title('Current Image')
        
        %Create a subplot to dispaly the "preview" image
        subplot(2,2,2);
        title('Preview Image')
        
        %Create a pushbutton to select and open a new source image as the "current" image
        uicontrol('Style','pushbutton','String','Open Image','Position',[10,20,80,20],'Callback',...
        {@Assignment_3.OpenImage,this});
      
        %Create a pushbutton to save the "current" image as a file
        uicontrol('Style','pushbutton','String','Save Image','Position',[100,20,60,20],'Callback',...
        {@Assignment_3.SaveImage,this});
      
        %Create a pushbutton to apply whichever filter is currently applied to the "current" 
        %image to create the "preview" image
        uicontrol('Style','pushbutton','String','Preview Image','Position',[170,20,80,20],'Callback',...
        {@Assignment_3.PreviewImage,this});
      
        %Create a pushbutton to make the "preview" image into the new "current" image
        uicontrol('Style','pushbutton','String','Apply Changes','Position',[260,20,80,20],'Callback',...
        {@Assignment_3.ApplyChanges,this});
              
        %Create a pushbutton to apply the selected enhancement to the "current" image 
        uicontrol('Style','pushbutton','String','Set Enhancement','Position',[350,20,90,20],'Callback',...
        {@Assignment_3.SelectEnhancement,this});
        
        %Create a dropdown menu to display enhancements/filters
        uicontrol('Style','popup','String',...
        {'No Filter','Greyscale','Ideal Low Pass 3x3','Ideal Low Pass 5x5',' IdealLow Pass 7x7',' Ideal Low Pass 9x9',...
        'Global Histogram Equalization', 'Adaptive Histogram Equalization','High Boost Filter',...
        'Butterworth Low Pass','Gaussian Low Pass'...
        'Ideal High Pass','Butterworth High Pass','Gaussian High Pass',...
        'Ideal Bandpass Filter','Butterworth Bandpass Filter','Gaussian Bandpass Filter','Homomorphic Filter'...
        'Generate Mask','Morphological Erosion','Morphological Dilation','Morphological Opening','Morphological Closing'...
        'Morphological Boundary','Morphological Smoothing','Morphological Gradient','Object Recognition'},...
        'Position',[350,70,90,20],'Callback',...
        {@Assignment_3.SelectFilter,this});
    
        %Create a slider to change the brightness of the image
        uicontrol('Style','slider','Min',0,'Max',1,'Value',.5,'Position',[10,70,150,20],'Callback',...
                 {@Assignment_3.SelectBrightness,this});
        
        %Create a slider to change the contrast of the image     
        uicontrol('Style','slider','Min',0,'Max',1,'Value',.5,...
                 'Position',[180,70,150,20],'Callback',{@Assignment_3.SelectContrast,this});
        
        %Create a label for the brightness slider
        uicontrol('Style','text','Position',[10,90,150,20],'String','Brightness');
        
        %Create a label for the contrast slider        
        uicontrol('Style','Text','Position',[180,90,150,20],'String','Contrast');   
        
        %Create a lable for the enhancement drop down menu
        uicontrol('Style','Text','Position',[315,90,150,20],'String','Enhancement');
    end
end

methods (Static)
    
    %This function opens and image
    function this=OpenImage(~,~,this)
        
        %Define file types that can be loaded in
        file = uigetfile({...
        '*.jpg;*.jpeg;*.bmp;*.gif;*.tif;*.tiff;*.png;*.pbm;*.pgm;*.ppm;*.pnm;*.pcx;*.ras;*.xwd;*.cur;*.ico','All Image Types';...
        '*.jpg;*.jpeg','JPEG (*.jpg, *.jpeg)';...
        '*.bmp','BMP (*.bmp)';...
        '*.gif','GIF (*.gif)';...
        '*.tif;*.tiff','TIFF (*.tif, *.tiff)';...
        '*.png','PNG (*.png)';...
        '*.pbm;*.pgm;*.ppm;*.pnm','PNM (*.pbm, *.pgm, *.ppm, *.pnm)';...
        '*.pcx','PCX (*.pcx)';...
        '*.ras','RAS (*.ras)';...
        '*.xwd','XWD (*.xwd)';...
        '*.cur;*.ico','Cursors and Icons (*.cur, *.ico)';...
        '*.*','All files (*.*)'});  
    
        if (~file)
            return
        end
    
        image = imread(file); %Assigns the chosen picture to image
    
        %Assign the chosen image to the "current" array
        subplot(2,2,1)
        imshow(image);
        title('Current Image')
        axis on;
        this.current=image;
    
        %Assign the chosen image to the "preview" array
        subplot(2,2,2)
        imshow(image);
        title('Preview Image')
        axis on;
        this.preview=image;
        
        this.original=image;
    end
    
    %This function saves the image
    function this=SaveImage(~,~,this)
        
        %Define file types that an image can be saved as
        file=uiputfile ({...
        '*.jpg;*.jpeg;*.bmp;*.gif;*.tif;*.tiff;*.png;*.pbm;*.pgm;*.ppm;*.pnm;*.pcx;*.ras;*.xwd;*.cur;*.ico','All Image Types';...
        '*.jpg;*.jpeg','JPEG (*.jpg, *.jpeg)';...
        '*.bmp','BMP (*.bmp)';...
        '*.gif','GIF (*.gif)';...
        '*.tif;*.tiff','TIFF (*.tif, *.tiff)';...
        '*.png','PNG (*.png)';...
        '*.pbm;*.pgm;*.ppm;*.pnm','PNM (*.pbm, *.pgm, *.ppm, *.pnm)';...
        '*.pcx','PCX (*.pcx)';...
        '*.ras','RAS (*.ras)';...
        '*.xwd','XWD (*.xwd)';...
        '*.cur;*.ico','Cursors and Icons (*.cur, *.ico)';...
        '*.*','All files (*.*)'});
                     
        if (~file)
            return
        end
        
        %Write the curent image to the file
        imwrite(this.current,file);
        
    end
    
    %This functiom implements all enhancements/filters
    %The switch case corresponds to enhancement drop down menu selection
    %choice
    function this=PreviewImage(~,~,this)
        switch this.choice
            %No filter
            case 1
                DeleteSubplot(); %Delete any subplots that are not needed for selected enhancement
                this.preview=this.original; %Assign the preview image to be the original image
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
  
            %Convert image to greyscale
            case 2
                DeleteSubplot();
                
                %Try to convert an image to greyscale, if it already is
                %greyscale just assign it to the preview image
                try
                    this.preview=this.current;
                
                    %Enhancement
                    this.preview=rgb2gray(this.preview); %Convert an rgb image to grayscale
                    colormap(gray(256));                 %Apply the greyscale color map to the image
                
                    %Show enhancement in preview image
                    subplot(2,2,2);
                    imshow(this.preview);
                    title('Preview Image');
                    axis on;
                catch
                    this.preview=this.current;
                    
                    %Show enhancement in preview image
                    subplot(2,2,2);
                    imshow(this.preview);
                    title('Preview Image');
                    axis on;
                end
                   
            %Ideal Low Pass Filter 3x3
            case 3
                DeleteSubplot();
                this.preview=this.current;
                
                %Enhancement
                filter=ones(3,3)/9; %Create a normalized, 3-by-3, averaging filter
                this.preview=imfilter(this.preview,filter); %Apply the averaging filter to the grayscale image using imfilter
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Ideal Low Pass Filter 5x5
            case 4
                DeleteSubplot();
                this.preview=this.current;
                
                %Enhancement
                filter=ones(5,5)/25; %Create a normalized, 5-by-5, averaging filter
                this.preview=imfilter(this.preview,filter); %Apply the averaging filter to the grayscale image using imfilter
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Ideal Low Pass Filter 7x7
            case 5
                DeleteSubplot();
                this.preview=this.current;
                
                %Enhancement
                filter=ones(7,7)/49; %Create a normalized, 7-by-7, averaging filter
                this.preview=imfilter(this.preview,filter); %Apply the averaging filter to the grayscale image using imfilter
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Ideal Low Pass Filter 9x9
            case 6
                DeleteSubplot();
                this.preview=this.current;
                
                %Enhancement
                filter=ones(9,9)/81; %Create a normalized, 9-by-9, averaging filter
                this.preview=imfilter(this.preview,filter); %Apply the averaging filter to the grayscale image using imfilter   
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Global Histogram Equalization (Enhances contrast using histogram equalization)
            case 7
               DeleteSubplot();
               this.preview=this.current;
               
               %Enhancement               
               this.preview = histeq(this.preview); %Transforms the intensity of the image so that the histogram of the
                                                    %output intensity image bins approximately match the target histogram
               %Show enhancement in preview image
               subplot(2,2,2);
               imshow(this.preview);
               title('Preview Image');
               axis on;
               
            %Adaptive Histogram Equalization (Contrast-limited adaptive histogram equalization)
            case 8
               DeleteSubplot();
               
               %Try to convert image, if error is encountered than the
               %image is in color and needs to be converted, go into catch
               %and convert to grayscale then perform operation
               try
                   this.preview=this.current;
               
                    %Enhancement 
                    this.preview = adapthisteq(this.preview); %Enhances the contrast of a grayscale image by transforming
                                                              %the values using CLAHE
                    %Show enhancement in preview image
                    subplot(2,2,2);
                    imshow(this.preview);
                    title('Preview Image');
                    axis on;
               catch
                   this.current=rgb2gray(this.current);
                   this.preview=this.current;
                   
                    %Enhancement 
                    this.preview = adapthisteq(this.preview); %Enhances the contrast of a grayscale image by transforming
                                                              %the values using CLAHE
                    %Show enhancement in preview image
                    subplot(2,2,2);
                    imshow(this.preview);
                    title('Preview Image');
                    axis on;
               end
               
            %High Boost Filter (Enhance high frequency components while still keeping the low frequency components)
            case 9
                DeleteSubplot();
                this.preview=this.current;
                
                %Prompt user for input
                prompt='Enter a value for the boost scaling coefficient: ';
                k=input(prompt);
                   
                fgauss=fspecial('gaussian',[25,25],4.0);
                filtim=imfilter(this.preview,fgauss,'symmetric','conv');
                high_boost=this.preview+k*(this.preview-filtim);
                this.preview=high_boost;
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                           
            %Butterworth Low Pass Filter
            case 10
                DeleteSubplot();
                this.preview=this.current;
                
                h=size(this.preview,1);
                w=size(this.preview,2);
                
                %Prompt user for input
                prompt='Enter desired value for Butterworth order: ';
                n=input(prompt);
               
                prompt='Enter desired value for cutoff frequency: ';
                d=input(prompt);
                
                [x,y]=meshgrid(-floor(w/2):floor(w-1)/2,-floor(h/2):floor(h-1)/2);
                out=1./(1.+(d./(x.^2+y.^2).^0.5).^(2*n));
                
                this.preview=fftshift(fft2(this.preview));
                afhb=this.preview.*(1.0001-out);
                afhbi=ifft2(afhb);
                
                this.preview=afhbi;
                
                %Show enhancement in preview image
                subplot(2,2,2);
                fftshow(this.preview,'abs');
                title('Preview Image');
                axis on;
                
            
            %Gaussian Low Pass Filter
            case 11
                DeleteSubplot();
                this.preview=this.current;
                
                %Prompt user for input
                prompt='Enter a value for sigma: ';
                sigma=input(prompt); 
                
                %Use matlab function to apply a low pass gaussian filter to
                %the current image with a user specified sigma value
                this.preview=imgaussfilt(this.preview,sigma);
                
                %Apply to preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Ideal High Pass Filter    
            case 12
                DeleteSubplot();
                this.preview=this.current;
                filterhp=[-1 -1 -1;-1 8 -1;-1 -1 -1];
                this.preview=imfilter(this.preview,filterhp,'same');
                
                %Apply to preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Butterworth High Pass Filter
            case 13
                DeleteSubplot();
                this.preview=this.current;
                
                h=size(this.preview,1);
                w=size(this.preview,2);
                fftim=fftshift(fft2(double(this.preview)));
                [x,y]=meshgrid(-floor(w/2):floor(w/2)-1,-floor(h/2):floor(h/2)-1);
                
                %Prompt user for input
                prompt='Enter desired value for Butterworth order: ';
                n=input(prompt);
               
                prompt='Enter desired value for cutoff frequency: ';
                d=input(prompt);
                
                B=sqrt(2)-1;
                D=sqrt(x.^2+y.^2); %Define distance to center
                hhp=1./(1+B*((d./D).^(2*n)));
                out_spec_center=fftim.*hhp;
                out_spec=ifftshift(out_spec_center); %Uncenter the spectrum
                out=real(ifft2(out_spec)); %Take the inverse fft and get real components
                
                %Normalize and cast
                out=(out-min(out(:)))/(max(out(:))-min(out(:)));
                out=uint8(255*out); 
                this.preview=out;
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Gaussian High Pass Filter
            case 14
                DeleteSubplot();
                this.preview=this.current;
                
                %Prompt user for input
                prompt='Enter a value for sigma: ';
                sigma=input(prompt);
                
                %Subtract low pass filter from original image for high pass filter
                this.preview=this.preview-imfilter(this.preview,fspecial('gaussian',5,sigma));
     
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                             
            %Ideal Bandpass Filter
            case 15
                this.preview=this.current;
                imageI=this.preview;
                imageI=double(imageI);
                subplot(2,2,1)
                imshow(this.preview);
                title('Current Image');
                axis on;
                
                %Fourier Spectrum of image using fft
                imafft = fftshift(fft2(fftshift(imageI)));
                imafft2=fft2(imageI);
                imafft3=fftshift(imafft2);
                
                s=size(imageI);
                for i=1:s(1)
                    for j=1:s(2)
                        r=sqrt((i-1-s(1)/2)^2+(j-1-s(2)/2)^2);
                        if (r<30)
                            z(i,j)=0;
                        else
                            if (r>120)
                                z(i,j)=0;
                            else
                                z(i,j)=255;
                            end
                        end
                    end
                end
                
                %Apply Fourier transform to image and show enhancement in preview image
                subplot(2,2,2)
                imafft=imafft.*z/255;
                image_out=fftshift(ifft2(fftshift(imafft)));
                image_out=image_out-imageI;
                fftshow(image_out,'abs');
                title('Ideal Bandpass Filtered Image');
                
                subplot(2,2,3)
                fftshow(imafft3,'log');
                title('Fourier Spectrum of Image')
                
                subplot(2,2,4)
                fftshow(z,'log');
                title('Frequency Domain Filter');
                
            %Butterworth Bandpass Filter (Should be grayscale image)
            case 16
                this.preview=this.current;
                imageB=double(this.preview);
                [nx,ny]=size(imageB);
                u=imageB;
                imageB=uint8(u);
                fftu=fft2(u,2*nx-1,2*ny-1);
                fftu=fftshift(fftu);
                
                subplot(2,2,1)
                imshow(imageB,[]);
                title('Current Image');
                axis on;
                
                subplot(2,2,3)
                fftshow(fftu,'log')
                title('Fourier Spectrum of Image')
                
                %Initialize filters
                filter1=ones(2*nx-1,2*ny-1);
                filter2=ones(2*nx-1,2*ny-1);
                filter3=ones(2*nx-1,2*ny-1);
                
                %Prompt user for input
                prompt='Enter desired value for Butterworth order: ';
                n=input(prompt);
                
                prompt='Enter a value for high end cutoff frequency: ';
                h=input(prompt);
                
                prompt='Enter desired value for low end cutoff frequency: ';
                l=input(prompt);
                
                %Apply Butterworth Filter to the image
                for i=1:2*nx-1
                    for j=1:2*ny-1
                        dist=((i-(nx+1))^2+(j-(ny+1))^2)^.5;
                        filter1(i,j)=1/(1+(dist/h)^(2*n)); %Where 120 is high end cutoff
                        filter2(i,j)=1/(1+(dist/l)^(2*n));  %Where 30 is the low end cutoff
                        filter3(i,j)=1.0-filter2(i,j);
                        filter3(i,j)=filter1(i,j).*filter3(i,j);
                    end
                end
                
                %Update the image with passed frequencies
                filtered_imageB = fftu+filter3.*fftu;
                
                subplot(2,2,4)
                fftshow(filter3,'log')
                filtered_imageB=ifftshift(filtered_imageB);
                filtered_imageB=ifft2(filtered_imageB,2*nx-1,2*ny-1);
                filtered_imageB=real(filtered_imageB(1:nx,1:ny));
                filtered_imageB=uint8(filtered_imageB);
                title('Frequency Domain Filter');
                
                %Show enhancement in preview image
                subplot(2,2,2)
                imshow(filtered_imageB,[]);
                title('Buttersworth Bandpass Filtered Image');
                axis on;
                
                this.preview=filtered_imageB;
                
            %Gaussian Bandpass Filter
            case 17
                this.preview=this.current;
                imageG = double(this.preview);
                [nx,ny] = size(imageG);
                imageG = uint8(imageG);
                fftI = fft2(imageG,2*nx-1,2*ny-1);
                fftI = fftshift(fftI);

                subplot(2,2,1)
                imshow(imageG,[]);
                title('Original Image')
                axis on;

                subplot(2,2,3)
                fftshow(fftI,'log')
                title('Fourier Spectrum of Image')
                
                % Initialize filter
                filter1 = ones(2*nx-1,2*ny-1);
                filter2 = ones(2*nx-1,2*ny-1);
                filter3 = ones(2*nx-1,2*ny-1);
                
                %Prompt user for input
                prompt='Enter a value for high end cutoff frequency: ';
                h=input(prompt);
                
                prompt='Enter desired value for low end cutoff frequency: ';
                l=input(prompt);

                %Apply Gaussian Filter to the image
                for i = 1:2*nx-1
                    for j =1:2*ny-1
                        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
                        filter1(i,j) = exp(-dist^2/(2*h^2));
                        filter2(i,j) = exp(-dist^2/(2*l^2));
                        filter3(i,j) = 1.0 - filter2(i,j);
                        filter3(i,j) = filter1(i,j).*filter3(i,j);
                    end
                end
                
                % Update image with passed frequencies
                filtered_imageG = fftI + filter3.*fftI;

                subplot(2,2,4)
                fftshow(filter3,'log')
                title('Frequency Domain Filter')
                filtered_imageG = ifftshift(filtered_imageG);
                filtered_imageG = ifft2(filtered_imageG,2*nx-1,2*ny-1);
                filtered_imageG = real(filtered_imageG(1:nx,1:ny));
                filtered_imageG = uint8(filtered_imageG);

                %Show enhancement in preview image
                subplot(2,2,2)
                imshow(filtered_imageG,[])
                title('Gaussian Bandpass Filtered Image')
                axis on;
                
                this.preview=filtered_imageG;
               
            %Homomorphic Filter    
            case 18
                DeleteSubplot();
                this.preview=this.current;
                subplot(2,2,1);
                imshow(this.current);
                title('Current Image');
                axis on;
                
                %Apply 2D Fourier Transform on the image
                dim=this.preview;
                cim=double(dim);
                [~,~]=size(dim);
                cim=cim+1;
                lim=log(cim);
                fim=fft2(lim);
                
                %Prompt user for input
                prompt='Enter a value for sigma: ';
                sigma=input(prompt);
                
                prompt='Enter desired high end gamma value: ';
                hg=input(prompt);
                
                prompt='Enter desired low end gamma value: ';
                lg=input(prompt);
                
                dif=(hg-lg);
                [r,c]=size(fim);
                him=fim;
                
                for i=1:r
                    for j=1:c
                        p=-(((i-(r/2))^2+(j-(c/2))^2)/(2*(sigma^2)));
                        k=(1/2*3.14*(sigma^2));
                        term(i,j)=(1-k*exp(p));
                        
                    end
                end
                
                for i=1:r
                    for j=1:c
                        h(i,j)=(dif*term(i,j))+lg;
                    end
                end
                
                for i=1:r
                    for j=1:c
                        him(i,j)=fim(i,j)*h(i,j);
                    end
                end
                
                %Inverese 2D Fourier Transform & exponent of the result for the filtered image 
                ifim=ifft2(him);
                eim=exp(ifim);
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(uint8(eim));
                title('Homomorphic Filtered Image');
                axis on;
                this.preview=uint8(eim);
                
            %Generate Mask
            case 19
                DeleteSubplot();
                this.preview=this.current;
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                this.preview=mask;
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
             
            %Erosion
            case 20 
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Prompt user to select a structuring element
                prompt='Select a structuring element:\n 1 for Disk\n 2 for Square\n 3 for Rectnagle\n Enter Choice: ';
                x=input(prompt);
                
                %Apply erosion to the mask
                if x==1
                    r=7; %Create variable for adjustable structuring element
                    SE=strel('disk',r);
                    mask=imerode(mask,SE);
                    this.preview=mask;
                end
                
                if x==2
                    mn=[10,7];
                    SE=strel('rectangle',mn);
                    mask=imerode(mask,SE);
                    this.preview=mask;
                end
                
                if x==3
                    w=7;
                    SE=strel('square',w);
                    mask=imerode(mask,SE);
                    this.preview=mask;
                end
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Dilation
            case 21
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Prompt user to select a structuring element
                prompt='Select a structuring element:\n 1 for Disk\n 2 for Square\n 3 for Rectnagle\n Enter Choice: ';
                x=input(prompt);
                
                %Dilate the mask
                if x==1
                    r=7;
                    SE=strel('disk',r);
                    mask=imdilate(mask,SE);
                    this.preview=mask;
                end
                
                if x==2
                    mn=[10,7];
                    SE=strel('rectangle',mn);
                    mask=imdilate(mask,SE);
                    this.preview=mask;
                end
                
                if x==3
                    w=7;
                    SE=strel('square',w);
                    mask=imdilate(mask,SE);
                    this.preview=mask;
                end
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
            
            %Opening (Performs morphological opening (erosion followed by dilation))
            case 22
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Prompt user to select a structuring element
                prompt='Select a structuring element:\n 1 for Disk\n 2 for Square\n 3 for Rectnagle\n Enter Choice: ';
                x=input(prompt);
                
                %Open the mask
                if x==1
                    r=7;
                    SE=strel('disk',r);
                    mask=imerode(imdilate(mask,SE),SE);
                    this.preview=mask;
                end
                
                if x==2
                    mn=[10,7];
                    SE=strel('rectangle',mn);
                    mask=imerode(imdilate(mask,SE),SE);
                    this.preview=mask;
                end
                
                if x==3
                    w=7;
                    SE=strel('square',w);
                    mask=imerode(imdilate(mask,SE),SE);
                    this.preview=mask;
                end

                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
            
            %Closing (Performs morphological closing (dilation followed by erosion))
            case 23
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Prompt user to select a structuring element
                prompt='Select a structuring element:\n 1 for Disk\n 2 for Square\n 3 for Rectnagle\n Enter Choice: ';
                x=input(prompt);
                
                %Close the mask
                if x==1
                    r=7;
                    SE=strel('disk',r);
                    mask=imdilate(imerode(mask,SE),SE);
                    this.preview=mask;
                end
                
                if x==2
                    mn=[10,7];
                    SE=strel('rectangle',mn);
                    mask=imdilate(imerode(mask,SE),SE);
                    this.preview=mask;
                end
                
                if x==3
                    w=7;
                    SE=strel('square',w);
                    mask=imdilate(imerode(mask,SE),SE);
                    this.preview=mask;
                end
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
                
            %Boundary (Traces the exterior boundaries of objects, as well as boundaries of holes inside these objects)     
            case 24
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Trace boundaries of the mask
                [B,L]=bwboundaries(mask,'noholes');
                subplot(2,2,2);
                imshow(label2rgb(L, @jet, [.5 .5 .5]))
                this.preview = label2rgb(L, @jet, [.5 .5 .5]);
                title('Preview Image');
                axis on;
                hold on
                for k=1:length(B)
                    boundary=B{k};
                    plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
                end
                
            %Smoothing (Use morphological closing to join the circles in an image
            %together by filling in the gaps between them and by smoothing
            %their outer edges)
            case 25
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Prompt user to select a structuring element
                prompt='Select a structuring element:\n 1 for Disk\n 2 for Square\n 3 for Rectnagle\n Enter Choice: ';
                x=input(prompt);
                
                %Smooth the mask
                if x==1
                    r=7;
                    SE=strel('disk',r);
                    mask=imclose(gpuArray(this.preview),SE);
                    this.preview=mask;
                end
                
                if x==2
                    mn=[10,7];
                    SE=strel('rectangle',mn);
                    mask=imclose(gpuArray(this.preview),SE);
                    this.preview=mask;
                end
                
                if x==3
                    w=7;
                    SE=strel('square',w);
                    mask=imclose(gpuArray(this.preview),SE);
                    this.preview=mask;
                end
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
        
            %Gradient (The difference between the dilation and the erosion of a given image)
            case 26
                DeleteSubplot();
                this.preview=this.current;
                
                %Create mask of image
                image=this.preview;
                height=size(image,1);
                width=size(image,2);
                mask=ones(height,width);
                imagehsv=rgb2hsv(image);
                keyplane=imagehsv(:,:,1);
                keycolor=imagehsv(1,1);
                diffplane=abs(keyplane-keycolor);
                mask=double(diffplane<.04);
                
                %Take the gradient of the mask
                SE=strel(ones(3,3));
                basic_gradient=imdilate(mask,SE)-imerode(mask,SE);
                this.preview=basic_gradient;
                
                %Show enhancement in preview image
                subplot(2,2,2);
                imshow(this.preview);
                title('Preview Image');
                axis on;
              
            %Object Recognition (Creates multiple figures)    
            case 27
                DeleteSubplot();
                
                %Create red mask
                this.preview=this.current;
                [~,this.preview]=createMaskRed(this.preview);
                figure('Name','Final Project-Jayce Jones-Red Mask','NumberTitle','off');
                imshow(this.preview);
                
                %Create green mask
                this.preview=this.current;
                [~,this.preview]=createMaskGreen(this.preview);
                figure('Name','Final Project-Jayce Jones-Green Mask','NumberTitle','off');
                imshow(this.preview);
                
                %Create blue mask
                this.preview=this.current;
                [~,this.preview]=createMaskBlue(this.preview);
                figure('Name','Final Project-Jayce Jones-Blue Mask','NumberTitle','off');
                imshow(this.preview);
                
                %Create yellow mask
                this.preview=this.current;
                [~,this.preview]=createMaskYellow(this.preview);
                figure('Name','Final Project-Jayce Jones-Yellow Mask','NumberTitle','off');
                imshow(this.preview);
                
                %Create pink mask
                this.preview=this.current;
                [~,this.preview]=createMaskPink(this.preview);
                figure('Name','Final Project-Jayce Jones-Light Blue Mask','NumberTitle','off');
                imshow(this.preview);
                
                %Create light blue mask
                this.preview=this.current;
                [~,this.preview]=createMaskLightBlue(this.preview);
                figure('Name','Final Project-Jayce Jones-Light Blue Mask','NumberTitle','off');
                imshow(this.preview);
        end
            
        %Brightness (Function that allows for the image brightness top be adjusted)
        if (this.brightness>0.5)
            this.preview=imadjust(this.current,[0,1],[2*this.brightness-1 1]);
            subplot(2,2,2);
            imshow(this.preview);
            title('Preview Image');
            axis on;
        elseif(this.brightness<0.5)
            this.preview=imadjust(this.current,[0,1],[0 2*this.brightness]);
            subplot(2,2,2);
            imshow(this.preview);
            title('Preview Image');
            axis on;
        end
  
        %Contrast (Function that allows for the image contrast to be adjusted)
        if (this.contrast>0.5)
            this.preview=imadjust(this.current,[0 2*this.contrast-1],[0,1]);
            subplot(2,2,2);
            imshow(this.preview);
            title('Preview Image');
            axis on;
       elseif (this.contrast<0.5)
            this.preview=imadjust(this.current,[2*this.contrast 1],[0,1]);
            subplot(2,2,2);
            imshow(this.preview);
            title('Preview Image');
            axis on;
        end
    end
    
    %ApplyChanges callback function for the apply changes button
    function this=ApplyChanges(~,~,this)
        this.current=this.preview;
        subplot(2,2,1);
        imshow(this.current);
        title('Current Image');
        axis on;
    end
    
    %SelectFilter callback function for the select filter button
    function this=SelectFilter(obj,~,this)
        this.choice=get(obj,'Value');
    end
    
    %SelectBrightness callback function for the select brightness slider
    function this=SelectBrightness(obj,~,this)
        
        this.brightness=get(obj,'Value');
    end
    
    %SelectContrast callback function for the select contrast slider
    function this=SelectContrast(obj,~,this)
        this.contrast=get(obj,'Value');
    end
    
    %SelectEnhancement callback function for the select enhancement button
    function this=SelectEnhancement(obj,~,this)
        this.enhance=get(obj,'Value');
    end
    
end
end
