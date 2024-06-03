function r = ImageNorm(s,minv,maxv)  


  r=(s-minv)/(maxv-minv);
  r(r<0)=0;
  r(r>1)=1;
end