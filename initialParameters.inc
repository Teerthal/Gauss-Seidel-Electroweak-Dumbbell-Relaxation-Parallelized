c Here we list the problem-specific parameters. These will be used
c in the initial conditions, in detecting features, and in main.f.
c latx, dx have been defined in parameters.inc and we use those here.
      integer im,jm
      real*8 xm,ym,zm,vxm,vym,str_rad,szoff
      parameter(im=0,jm=0)
c location of monopole (zm is specified as an argument when executing):
      parameter(xm=(float(im)+0.5)*dx,ym=(float(jm)+0.5)*dx)
      parameter(vxm=0.,vym=0.)
c Radius of the string profile
      parameter(str_rad=(1.)*dx)
      parameter(szoff=0.*dx)