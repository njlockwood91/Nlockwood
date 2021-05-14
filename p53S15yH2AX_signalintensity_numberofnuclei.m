clear all 
close all 
addpath(genpath('/Users/lockwon/Documents/MATLAB/Patient Cells'))
cd('/Users/lockwon/Documents/MATLAB/Patient Cells/270421_patientcells/21488_yh2ax555/Patient/UT')
addpath(genpath('/Users/lockwon/Documents/MATLAB/imreadBF_2/'))

directory = cd;
d=dir(directory);

%LOADING IN THE DATA
experimentname=['270421_patientcells_quantyp53S15_pernuc_PatientUT.mat']; %EDIT THIS TO SOMETHING USEFUL 
imagenumber=10;  %MAKE SURE YOU EDIT THIS
diff=length(d)-imagenumber; %Difference in the directory and number of images. 

for i=1:imagenumber
      dapi_im(:,:,i)=imreadBF(d(i+diff).name,1,1,1);
      p53_im(:,:,i)=imreadBF(d(i+diff).name,1,1,3);
      S15_im(:,:,i)=imreadBF(d(i+diff).name,1,1,2);
end

for j=1:size(dapi_im,3)
    %run through each image one at a time

    dapi1=uint16(dapi_im(:,:,j));
    nuclei_bin=imfill(im2bw(dapi1,graythresh(dapi1)),'holes');
    
    %labelling the individual nuclei
    dapilabel = bwlabel(nuclei_bin,4);
    dapi_props=regionprops(dapilabel,'Area');
    dapiA=[dapi_props.Area];
    keepdapiIndices=find(dapiA >= 150);
    keepdapiCell=find(dapiA <= 10000);
    dapi_mask=ismember(dapilabel,keepdapiIndices);
    dapi_mask_cell1=ismember(dapilabel,keepdapiCell);
    toobig=(dapi_mask-dapi_mask_cell1)<1;
    dapi_mask_cell=dapi_mask;
    dapi_mask_cell(~toobig)=0;
    dapifiltered=nuclei_bin.*dapi_mask; 
    dapicellfiltered=nuclei_bin.*dapi_mask_cell;
    dapifiltered2=imclearborder(dapifiltered);
    dapicellfiltered2=imclearborder(dapicellfiltered);
    label = bwlabel(dapifiltered2,4);
    
    

    S15_foci_blur = medfilt2(S15_im(:,:,j),[2 2]);
    p53_blur = medfilt2(p53_im(:,:,j),[2 2]);

    %for each nuclei
    for i = 1:max(label(:))
        clear nuclei yh2ax_foci_seg yh2ax_foci_seg2 yh2ax_int nuc_size 
        nuclei=label==i;
        
        
        %gives S15 foci for just one nuclei
        S15_foci_seg = S15_foci_blur;
        S15_foci_seg(~nuclei) = 0; 
        S15_foci_seg2 = S15_foci_seg(:); 
        
        %sum of S15 intensity
        S15_int=sum(S15_foci_seg2);
        
        %gives p53 for just one nuclei
        p53_seg = p53_blur;
        p53_seg(~nuclei) = 0; 
        p53_seg2 = p53_seg(:); 
        
        %sum of p53 intensity
        p53_int=sum(p53_seg2);
 
        
        %per nucleus
        S15_int_all(j,i)=S15_int;
        p53_int_all(j,i)=p53_int;
        
 
        
        
    end%of each nuclei
    
      %nuclei number
        nucl_number(j)=max(label(:));
    
   %now working on each image. 
      
  
end

S15notzero=nonzeros(S15_int_all');

p53notzero=nonzeros(p53_int_all');

mean_S15pernuc_allimages=mean(S15notzero); 
sdS15_pernuc=std(S15notzero);
mean_p53pernuc_allimages=mean(p53notzero); 
sdp53_pernuc=std(p53notzero);

sprintf('The mean p53 signal per cell over all images is %f with a SD of %f',mean_p53pernuc_allimages, sdp53_pernuc)
sprintf('The mean percentage of s15 signal per cell over all images is %f with a SD of %f',mean_S15pernuc_allimages, sdS15_pernuc)


sprintf('The number of nuclei over all images is %f',length(S15notzero))

cd('/Users/lockwon/Documents/MATLAB/Patient Cells/270421_patientcells/21488_yh2ax555/Patient/UT_Results') %MAKE SURE PATHS RIGHT
save (experimentname, 'mean_S15pernuc_allimages', 'sdS15_pernuc','mean_p53pernuc_allimages', 'sdp53_pernuc');
    
