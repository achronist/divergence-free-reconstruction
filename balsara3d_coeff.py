def balsara3d_coeff(Dx,Dy,Dz,
                    Bxp,DyBxp,DzBxp,Bxm,DyBxm,DzBxm,
                    Byp,DxByp,DzByp,Bym,DxBym,DzBym,
                    Bzp,DxBzp,DyBzp,Bzm,DxBzm,DyBzm):

    '''Divergence-free reconstruction in axis-aligned hexahedron given
    bilinear variation of normal magnetic field on faces
    
   Arguments
   Dx    size of cell in x direction
   Dy    size of cell in y direction
   Dz    size of cell in z direction
   Bxp,DyBxp,DzBxp    linear variation of Bx on face xp of cell
   Bxm,DyBxm,DzBxm    linear variation of Bx on face xm of cell
   Byp,DxByp,DzByp    linear variation of By on face yp of cell
   Bym,DxBym,DzBym    linear variation of By on face ym of cell
   Bzp,DxBzp,DyBzp    linear variation of Bz on face zp of cell
   Bzm,DxBzm,DyBzm    linear variation of Bz on face zm of cell

    Returns
    (a0,ax,ay,az,axx,axy,axz,b0,bx,by,bz,bxy,byy,byz,c0,cx,cy,cz,cxz,cyz,czz)
        Polynomial coefficients for reconstruction in cell
        Polynomial arguments to be offset so that
        -Dx/2<x<Dx/2, -Dy/2<y<Dy/2, -Dz/2<z<Dz/2

     Reference: Balsara2004
     Notation DyBxp is described in Balsara2004 eqn (3.2)
    '''
    
    # Determine divergence-free quadratic polynomial
    # matching linear variation of B on faces of quad

    ax = (Bxp - Bxm)/Dx
    ay = 0.5*(DyBxp/Dy + DyBxm/Dy)
    az = 0.5*(DzBxp/Dz + DzBxm/Dz)
    axy = (DyBxp/Dy - DyBxm/Dy)/Dx
    axz = (DzBxp/Dz - DzBxm/Dz)/Dx

    # (a,b,c,x,y,z) -> (b,c,a,y,z,x)
    by = (Byp - Bym)/Dy
    bz = 0.5*(DzByp/Dz + DzBym/Dz)
    bx = 0.5*(DxByp/Dx + DxBym/Dx)
    byz = (DzByp/Dz - DzBym/Dz)/Dy
    bxy = (DxByp/Dx - DxBym/Dx)/Dy

    # (a,b,c,x,y,z) -> (b,c,a,y,z,x)
    cz = (Bzp - Bzm)/Dz
    cx = 0.5*(DxBzp/Dx + DxBzm/Dx)
    cy = 0.5*(DyBzp/Dy + DyBzm/Dy)
    cxz = (DxBzp/Dx - DxBzm/Dx)/Dz
    cyz = (DyBzp/Dy - DyBzm/Dy)/Dz

    # Cross-talk
    axx = -0.5*(bxy + cxz)
    byy = -0.5*(cyz + axy)
    czz = -0.5*(axz + byz)

    # Offset
    a0 = 0.5*(Bxp + Bxm) - 0.25*axx*Dx**2
    b0 = 0.5*(Byp + Bym) - 0.25*byy*Dy**2
    c0 = 0.5*(Bzp + Bzm) - 0.25*czz*Dz**2

    return (a0,ax,ay,az,axx,axy,axz,
            b0,bx,by,bz,bxy,byy,byz,
            c0,cx,cy,cz,cxz,cyz,czz)
