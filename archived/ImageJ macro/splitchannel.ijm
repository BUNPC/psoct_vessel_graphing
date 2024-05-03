//SplitChannelMacro 

        inputFolder= getDirectory("Choose a Directory"); 
        print(inputFolder); 
        images= getFileList(inputFolder); 

        setBatchMode(true); //batch mode 
        for (i=0; i<images.length; i++) { 
                if(endsWith(images[i],"_Ch1_MIP.tif")){
                	inputPath= inputFolder + images[i]; 
                	open (inputPath); 
                	imagesName=getTitle(); 
                	info=split(imagesName,'-');
                	index=split(info[3],'_');
                	a=index[1];
                	b=substring(a,7,10);
                	saveAs("tiff", inputFolder + b + "_ch1_mip.tif");	
                	close(); 	 
                }

        
} 
setBatchMode(false); 