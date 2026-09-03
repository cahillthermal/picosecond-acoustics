%initial

function[x,y,a,da,ndata,istrt,yoffset]=initial(filename,paramname)


%	/* ASSUME 
%		CRAP AT TOP(ASSUME 8 LINES OF JUNK)is gone
%		PUT DATA IN "DOGFILE" FOR USE IN PROGRAM*/
	k=0;
    maxstr=500;
	f2=fopen(filename,'r');
	%f1=fopen(dogfile,'w');
	%for i=1:8
    %    fgets(strng,maxstr,f2);
    %end
	
    data=fscanf(f2,'%f',[2,1125]);  %fix this to be max
    %data=flipud(data);   %stuff from labview comes out reversed
	
	fclose(f2);
	paramfile=fopen(paramname,'r');
    
	%f3=fopen("dogfile","r");
    
%    	/*INITIALIZE COUNTING PARAMETERS*/
	ymax=0.0;
	i=0;
	tcnt=0;
	count=0;
    tstart=300;
    tend=4200;
    yoffset=0;
    
    %/*USE TIME BEFORE PULSE HITS TO FIND OFFSET OF DATA SET-THEN READ
	%	IN DATA WHERE IT WILL BE COMPARED TO THE FIT FUNCTION*/
	j=0;
    for i=1:1125
        if (data(1,i)>-200 && data(1,i)<-10) 
			count=count+1;
			array(count)=data(2,i);
			yoffset=yoffset+data(2,i);
        end     
        
		if (data(1,i)>=0 && data(1,i)<=tend)
			j=j+1;
			if(data(1,i)>=tstart && data(1,i)<=tend) 
                tcnt=tcnt+1;
            end
			x(j)=data(1,i);
			y(j)=data(2,i);
        end
    end
	yoffset=yoffset/count;
    if (yoffset < 0)
        y=y-yoffset;
        yoffset=0;
    end
    
    a=fscanf(paramfile,'%f',[1,5]);
    da=a*.001;
    
    fclose(paramfile);
    
    ndata=tcnt;
    istrt=j-tcnt;
    
    
    
