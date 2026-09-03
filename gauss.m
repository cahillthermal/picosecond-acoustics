%gauss

n_steps=100;
     
      step=1.0/sqrt(a_n_steps);
      x=0;
      for n=1:n_steps
        y=rand;
        if (y > 0.5) then 
          x=x+step; 
        else if
          x=x-step;
        end
      end 
      
      end  