function [u, v] = getgrid(x, y, xmin, xmax, xres, ymin, ymax, yres)

   %// Find how many grid points we are from the center
   u=round(x/xres);     
   v=round(y/yres);

   %// Add the center point of grid to each offset
   u=u+(-xmin/xres)+mod(1+(xmax-xmin)/xres,2);
   v=v+(-ymin/yres)+mod(1+(ymax-ymin)/yres,2);

end