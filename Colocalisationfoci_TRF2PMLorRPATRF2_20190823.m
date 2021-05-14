
%OPEN IMAGE FILE - tiffs
clear all 
close all 
addpath(genpath('/Users/lockwon/Documents/MATLAB/APBquantification/'))
addpath(genpath('/Users/lockwon/Documents/MATLAB/imreadBF_2/'))
addpath(genpath('/Users/lockwon/Documents/MATLAB/'))
cd('/Users/lockwon/Documents/MATLAB/APBquantification/050521_U2OS_BLMkd_TRF2488_PML555/Chk1 ICRF')

directory = cd;
d=dir(directory);
minArea=1; maxArea=20; 

%LOADING IN THE DATA

experimentname=['05052021_quantTRF2PMLColocalisation_Chk1 ICRF_20190823script_98PML99TRF2perc.mat']; %EDIT THIS TO SOMETHING USEFUL 
imagenumber=10;  
diff=length(d)-imagenumber; 

for i=1:imagenumber
    clear dapi_stack dapi_MIP PML_stack PML_MIP TRF2_stack pTRF2_MIP metadata
    metadata=imreadBFmeta_editedforTS(d(i+diff).name); 
    zsize=metadata.zsize;
    
    dapi_stack(:,:,:)=imreadBF(d(i+diff).name,[1:zsize],1,1);
    dapi_MIP(:,:)=max(dapi_stack,[],3);
    dapi_im(:,:,i)=dapi_MIP;
    
    PML_stack(:,:,:)=imreadBF(d(i+diff).name,[1:zsize],1,3);
    PML_MIP(:,:)=max(PML_stack,[],3);
    PML_im(:,:,i)=PML_MIP;
    
    TRF2_stack(:,:,:)=imreadBF(d(i+diff).name,[1:zsize],1,2);
    TRF2_MIP(:,:)=max(TRF2_stack,[],3);
    TRF2_im(:,:,i)=TRF2_MIP;

     
end
 

for j=1:size(dapi_im,3)

    dapi1=uint16(dapi_im(:,:,j));
    nuclei_bin=imfill(im2bw(dapi1,graythresh(dapi1)),'holes');
  
    %labelling the individual nuclei
    dapilabel = bwlabel(nuclei_bin,4);
    dapi_props=regionprops(dapilabel,'Area');
    dapiA=[dapi_props.Area];
    keepdapiIndices=find(dapiA >= 1000);
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
    labelcell=bwlabel(dapicellfiltered2,4); 
    
    nuclei_number(:,j)=max(labelcell(:));

    TRF2_foci_blur = medfilt2(TRF2_im(:,:,j),[2 2]);
    PML_foci_blur = medfilt2(PML_im(:,:,j),[2 2]);

    
    %for each nuclei
    for i = 1:max(label(:))
        clear nuclei TRF_foci_seg TRF2_foci_seg2 PML_foci_seg PML_foci_seg2 TRF2_var pml_var TRF2_foci PML_foci foci_props
        nuclei=label==i;
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TRF2 focifigure(
       TRF2_foci_seg = TRF2_foci_blur;
        TRF2_foci_seg(~nuclei) = 0; %gives just one nuclei
        TRF2_foci_seg2 = TRF2_foci_seg; 
        if max(TRF2_foci_seg2(:)) > 0
            TRF2_foci_seg2(TRF2_foci_seg == 0) = [];
            TRF2_var=var(double(TRF2_foci_seg2(:)));
        TRF2_percentile=prctile(TRF2_foci_seg2,97.5);
        TRF2_focimedian=prctile(TRF2_foci_seg2,50);
        TRF2_percentilehigh=prctile(TRF2_foci_seg2,99.75); %%%%%%%
        else
            TRF2_foci_seg2=0;
            TRF2_var=0;
        end
        
        if TRF2_var < 5 %IF no foci
            TRF2_foci = 0;
        elseif max(TRF2_foci_seg2(:)) < 20 %if no foci just variable noise
            TRF2_foci =0;
        elseif isempty(TRF2_foci_seg2) >0 %if totally empty
            TRF2_foci=0;
        elseif TRF2_focimedian > 50; %if high background
            TRF2_foci=(TRF2_foci_seg>TRF2_percentilehigh);
        else
        TRF2_foci= (TRF2_foci_seg>TRF2_percentile);
        end
        
        %label and filter topo foci 
       TRF2_FOCI_labeled=bwlabel(TRF2_foci,4);
        TRF2_FOCI_props=regionprops(TRF2_FOCI_labeled,'Area'); 
        TRF2_A=[TRF2_FOCI_props.Area];
        TRF2_keepIndices=find(TRF2_A >= minArea & TRF2_A <= maxArea);
        TRF2_FOCI_mask=ismember(TRF2_FOCI_labeled,TRF2_keepIndices);
        TRF2_FOCIfiltered=TRF2_FOCI_labeled.*TRF2_FOCI_mask; %
        TRF2_numFOCI=bwconncomp(TRF2_FOCIfiltered);
        TRF2_foci_number(i)=TRF2_numFOCI.NumObjects;
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          
        %pml foci 
        PML_foci_seg = PML_foci_blur;
        PML_foci_seg(~nuclei) = 0; %gives just one nuclei
        PML_foci_seg2 = PML_foci_seg; 
        if max(PML_foci_seg2(:)) > 0
            PML_foci_seg2(PML_foci_seg == 0) = [];
            PML_var=var(double(PML_foci_seg2(:)));
        PML_percentile=prctile(PML_foci_seg2,99.5);
        PML_focimedian=prctile(PML_foci_seg2,50);
        PML_percentilehigh=prctile(PML_foci_seg2,98); %%%%%%%
        else
            PML_foci_seg2=0;
            PML_var=0;
        end
        
        if PML_var < 5 %IF no foci
            PML_foci = 0;
        elseif max(PML_foci_seg2(:)) < 20 %if no foci just variable noise
            PML_foci =0;
        elseif isempty(PML_foci_seg2) >0 %if totally empty
            PML_foci=0;
        elseif PML_focimedian > 50; %if high background
            PML_foci=(PML_foci_seg>PML_percentilehigh);
%        
        else
        PML_foci= (PML_foci_seg>PML_percentile);
        end
        
        %label and filter pml foci 
        PML_FOCI_labeled=bwlabel(PML_foci,4);
        PML_FOCI_props=regionprops(PML_FOCI_labeled,'Area'); 
        PML_A=[PML_FOCI_props.Area];
        PML_keepIndices=find(PML_A >= minArea & PML_A <= maxArea);
        PML_FOCI_mask=ismember(PML_FOCI_labeled,PML_keepIndices);
        PML_FOCIfiltered=PML_FOCI_labeled.*PML_FOCI_mask; %
        PML_numFOCI=bwconncomp(PML_FOCIfiltered);
        PML_foci_number(i)=PML_numFOCI.NumObjects;
      
          
          
%%%% coloc per nuclei

colocalised_foci=PML_FOCIfiltered.*TRF2_FOCIfiltered;
colocFOCI=bwconncomp(colocalised_foci);
colocFOCI_number(i)=colocFOCI.NumObjects; %number of colocalised foci per nuclei in the image


if PML_foci_number(i)>0
    perccoloc(i)= (colocFOCI_number(i)/PML_foci_number(i))*100;%percent of colocalisation as a funciton of PML foci per nuclei
else 
    perccoloc(i)= 0;
end
  
    end %of each nuclei
    
   %now working on each image. 
   mean_perccolocperimage(j)=mean(perccoloc(:)); 
    mean_numbercolocFOCIperimage(j) = mean(colocFOCI_number(:));
    mean_perccolocperimage(j)=mean(perccoloc(:)); %mean colocalised for each image
   mean_pmlnumberpernuclei(j)=mean(sum(PML_foci_number(:))/length(PML_foci_number)); %mean pml number per cell for each image.
   mean_trf2numberpernuclei(j)=mean(sum(TRF2_foci_number(:))/length(TRF2_foci_number)); %mean pml number per cell for each image.

   %    
   %mean pml foci per nuclear area
   mean_pml_nucleararea(j)=sum(PML_foci_number(:))/sum(dapiA);

   %mean trf2 foci per nuclear area
   mean_trf2_nucleararea(j)=sum(TRF2_foci_number(:))/sum(dapiA);


    
end
mean_numbercoloc_allimages=mean(mean_numbercolocFOCIperimage(:)); 
numberofcells=sum(nuclei_number);
mean_pmlnuclarea_allimages=mean(mean_pml_nucleararea(:));
mean_trf2nuclarea_allimages=mean(mean_trf2_nucleararea(:));
mean_pmlnucl_allimages=mean(mean_pmlnumberpernuclei(:));
mean_trf2nucl_allimages=mean(mean_trf2numberpernuclei(:));
mean_perccoloc_allimages=mean(mean_perccolocperimage(:));

sdperccoloc=std(mean_perccolocperimage(:));

sd_numberccoloc=std(mean_numbercolocFOCIperimage(:));
sdpmlarea=std(mean_pml_nucleararea(:));
sdtrf2area=std(mean_trf2_nucleararea(:));
sdpmlnuc=std(mean_pmlnumberpernuclei(:));
sdtrf2nuc=std(mean_trf2numberpernuclei(:));

sprintf('The mean number of colocalised foci per cell over all images is %f with a SD of %f',mean_numbercoloc_allimages, sd_numberccoloc)
sprintf('The mean percentage of colocalised foci per cell over all images is %f with a SD of %f',mean_perccoloc_allimages, sdperccoloc)
sprintf('The mean number of PML foci per nuclear area over all images is %f with a SD of %f',mean_pmlnuclarea_allimages, sdpmlarea)
sprintf('The mean number of TRF2 foci per nuclear area over all images is %f with a SD of %f',mean_trf2nuclarea_allimages, sdtrf2area)
sprintf('The mean number of PML foci per nucleus over all images is %f with a SD of %f',mean_pmlnucl_allimages, sdpmlnuc)
sprintf('The mean number of TRF2 foci per nuclear area over all images is %f with a SD of %f',mean_trf2nucl_allimages, sdtrf2nuc)
sprintf('The number of nuclei counted was %f',numberofcells)
cd('/Users/lockwon/Documents/MATLAB/APBquantification/050521_U2OS_BLMkd_TRF2488_PML555/BLM ICRF Results') %MAKE SURE PATHS RIGHT
save (experimentname, 'mean_numbercoloc_allimages', 'sd_numberccoloc','numberofcells','mean_pmlnuclarea_allimages','mean_trf2nuclarea_allimages','sdpmlarea','sdtrf2area');

