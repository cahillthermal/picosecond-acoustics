%gaussrand


function[x]=gaussrand
n_steps=100;
     
      step=1.0/sqrt(n_steps);
      x=0;
      for n=1:n_steps
        y=rand;
        if (y > 0.5)
          x=x+step; 
        else
          x=x-step;
        end
      end