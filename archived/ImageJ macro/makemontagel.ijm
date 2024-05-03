//SplitChannelMacro 

        inputFolder= getDirectory("Choose a Directory"); 
        images= getFileList(inputFolder); 


        for (y=20; y>0; y--) {
        	for (x=23; x>0; x--) {
	        	for (j=0; j<images.length; j++) {
	        		index = (x-1)*20 + y;
	        		filename = d2s(index,0);
	        		filename = filename + "_aip.tif";
	                if(images[j] == filename){
	                	inputPath= inputFolder + images[j]; 
	                	open(inputPath);
	                	//print(inputPath);
	                } 
	        	}
        	}    
		} 
		run("Images to Stack", "name=day0 title=[]");
		run("Make Montage...", "columns=23 rows=20 scale=0.25");
		//saveAs("Tiff", "C:/Users/jryang/Desktop/Data/0629TDE/day0/2x/50um/aip/vol1/Montage.tif");
		//close();
