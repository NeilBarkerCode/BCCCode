# See excel sheet for definitions

def BCCFEA01(X, Y, Z, R, OF, File):

    # Packages to be imported
    import numpy as np

    # BCCEquation packing equation
    BCCEquation = (4 / (np.sqrt(3)))*R*OF

    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA01-'
                    + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-' + str(OF) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA01-'
                                + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-' + str(OF) + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                               'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim spheres at cube corners
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates internal spheres for BCCEquation structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    SpaceClaimFile.close()

    # write model information to a text file
    TextFile.write('The Function used: BCCFEA01' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'Number of spheres in the X direction ' + str(X)
                   + '\n' + 'Number of spheres in the Y direction ' + str(Y) + '\n'
                   + 'Number of spheres in the Z direction ' + str(Z) + '\n' + 'The sphere radius '
                   + str(R) +'mm' + '\n' + 'The sphere overlapping fraction ' + str(OF) + '\n' + '\n')


    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = (X-1)*BCCEquation
    YDim = (Y-1)*BCCEquation
    ZDim = (Z-1)*BCCEquation

    TextFile.write('The dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCornerSpheres = (X - 1) * (Y - 1) * (Z - 1)
    TotalNumberOfSpheres = NumberOfCentreSpheres + NumberOfCornerSpheres
    NumberOfWindows = NumberOfCentreSpheres * 8

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) +'\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) +'\n')
    TextFile.write('The total number of windows - ' + str(NumberOfWindows) +'\n')

    # distance between centre of 2 spheres
    r = R
    D = (2*R) * OF
    # Window radius of the intersecting spheres
    WindowRadius = float((1 / (2 * D)) * np.sqrt((-D + r - R) * (-D - r + R) * (-D + r + R) * (D + r + R)))

    # finding the volume of the removed lens common to the intersecting spheres
    LensVolume = float((np.pi * (4 * R + D) * ((2 * R - D) ** 2)) / 12)
    TotalLensVolume = LensVolume*NumberOfWindows

    # total volume of the voids within the cube
    CubeVoidVolume = TotalNumberOfSpheres * (float((4 * np.pi * R ** 3) / 3)) - TotalLensVolume


    TextFile.write('The window radius is ' + str(WindowRadius) + 'mm' + '\n')
    TextFile.write('The single lens volume - ' + str(LensVolume) + 'mm^3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm^3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm^3' + '\n')

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1) * (Y - 1) * (Z - 1) * BCCEquation**3)
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm^3' + '\n')

    # Volume of porous material
    VolumeOfPorousMaterial = TotalCubeVolume-CubeVoidVolume
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) +'\n')

    TextFile.close()

def BCCFEA02(Xmm, Ymm, Zmm, R, P, File):

    # Packages to be imported
    import numpy as np

    # Massive equation that calculates the overlapping fraction from the sphere radius and the porosity argument.
    # See Maple file in the python function folder to see how this was completed.
    OF = abs(1/2/(np.pi*3**(1/2)+4*P)*(-2*np.pi*(3*3**(1/2)-(-3*(7*np.pi*3**(1/2)-36*P)/(np.pi*3**(1/2)
        +4*P))**(1 /2))*(np.pi*3**(1/2)+4*P)**2)**(1/3)+2*np.pi*3**(1/2)/(-2*np.pi*(3*3**(1/2)
        -(-3*(7*np.pi*3**(1/2)-36*P)/(np.pi*3**(1/2)+4*P))**(1/2))*(np.pi*3**(1/2)+4*P)**2)**(1/3))

    # BCCEquation packing equation
    BCCEquation = (4 / (np.sqrt(3))) * R * OF

    # converting the dimensions in mm to number of spheres
    X = int(Xmm / BCCEquation) + 2  # number of spheres in the X direction
    Y = int(Ymm / BCCEquation) + 2  # number of spheres in the Y direction
    Z = int(Zmm / BCCEquation) + 2  # number of spheres in the Z direction


    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA02-'
                    + str(Xmm) + '-' + str(Ymm) + '-' + str(Zmm) + '-' + str(R) + '-' + str(P) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA02-'
                          + str(Xmm) + '-' + str(Ymm) + '-' + str(Zmm) + '-' + str(R) + '-' + str(P) + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                         'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim spheres at cube corners
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates internal spheres for BCCEquation structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    SpaceClaimFile.close()

    # write model information to a text file
    TextFile.write('The Function used: BCCFEA02' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'X direction dimensions ' + str(Xmm)
                   + '\n' + 'Y direction dimensions ' + str(Ymm) + '\n'
                   + 'Z direction dimensions ' + str(Zmm) + '\n' + 'The sphere radius '
                   + str(R) + 'mm' + '\n' + 'The porosity is ' + str(P) + '\n' + '\n')

    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = (X - 1) * BCCEquation
    YDim = (Y - 1) * BCCEquation
    ZDim = (Z - 1) * BCCEquation

    TextFile.write('The actual dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCornerSpheres = (X - 1) * (Y - 1) * (Z - 1)
    TotalNumberOfSpheres = NumberOfCentreSpheres + NumberOfCornerSpheres
    NumberOfWindows = NumberOfCentreSpheres * 8

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) + '\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) + '\n')
    TextFile.write('The total number of windows - ' + str(NumberOfWindows) + '\n')

    # distance between centre of 2 spheres
    r = R
    D = (2 * R) * OF
    # Window radius of the intersecting spheres
    WindowRadius = float((1 / (2 * D)) * np.sqrt((-D + r - R) * (-D - r + R) * (-D + r + R) * (D + r + R)))

    # finding the volume of the removed lens common to the intersecting spheres
    LensVolume = float((np.pi * (4 * R + D) * ((2 * R - D) ** 2)) / 12)
    TotalLensVolume = LensVolume * (NumberOfWindows)

    # total volume of the voids within the cube
    CubeVoidVolume = TotalNumberOfSpheres * (float((4 * np.pi * R ** 3) / 3)) - TotalLensVolume

    TextFile.write('The window radius is ' + str(WindowRadius) + 'mm' + '\n')
    TextFile.write('The single lens volume - ' + str(LensVolume) + 'mm**3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm**3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm**3' + '\n')

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1) * (Y - 1) * (Z - 1) * BCCEquation ** 3)
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm**3' + '\n')

    # Volume of porous material
    VolumeOfPorousMaterial = TotalCubeVolume - CubeVoidVolume
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)
    TextFile.write('The BCC overlapping fraction is - ' + str(OF) + '\n')

    TextFile.close()


def BCCFEA03(X, Y, Z, R, Rc, File):

    # Packages to be imported
    import numpy as np

    # BCCEquation packing equation
    BCCEquation = (4 / (np.sqrt(3)))*R

    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA03-'
                    + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-'+ str(Rc) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA03-'
                                + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-' + str(Rc) + '-' + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                               'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim corner spheres
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates centre spheres for BCC structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + Rc) + '), MM(' + str(
                        xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # closes the file just written
    SpaceClaimFile.close()



    # write model information to a text file
    TextFile.write('The Function used: BCCFEA03' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'Number of spheres in the X direction ' + str(X)
                   + '\n' + 'Number of spheres in the Y direction ' + str(Y) + '\n'
                   + 'Number of spheres in the Z direction ' + str(Z) + '\n' + 'The corner sphere radius '
                   + str(R) + 'mm' + '\n' + 'The centre sphere radius ' + str(Rc) + '\n' + '\n')

    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = (X - 1) * BCCEquation
    YDim = (Y - 1) * BCCEquation
    ZDim = (Z - 1) * BCCEquation

    TextFile.write('The dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCornerSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCentreCornerWindows = NumberOfCornerSpheres * 8

    # The total number of centre centre windows changes with geometry as they are shared with connecting units.
    # In addition, the external windows need to be calculated separately as the lenz volume needs to be halfed

    if (2 * Rc) > BCCEquation:
        NumberOfCentreCentreWindows = (X * ((Y-1) * (Z-1)) + Y * ((X-1) * (Z-1)) + Z * ((X-1) * (Y-1)))
        ExternalNumberOfCentreCentreWindows = ((X-1)*(Y-1)*2)+((Y-1)*(Z-1)*2)+((Z-1)*(X-1)*2)
        InternalNumberOfCentreCentreWindows = NumberOfCentreCentreWindows - ExternalNumberOfCentreCentreWindows
    else:
        NumberOfCentreCentreWindows = 0
        ExternalNumberOfCentreCentreWindows = 0
        InternalNumberOfCentreCentreWindows = 0

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) + '\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) + '\n')
    TextFile.write('The total number of centre centre windows - ' + str(NumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of internal centre centre windows - '
                   + str(InternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of external centre centre windows - '
                   + str(ExternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of centre corner windows - ' + str(NumberOfCentreCornerWindows) + '\n')

    # calculation of the lens between centre centre sphere interaction

    # Distance between centre centre spheres
    DCentreCentre = BCCEquation
    DCentreCorner = 2*R

    # if the centre spheres do not overlap it returns a complex value. This If statemtent returns the Centre
    # Centre window radius as a zero if this is the case
    if (2*Rc)>BCCEquation:
         CentreCentreWindowRadius = float((1 / (2 * DCentreCentre)) * np.sqrt((-DCentreCentre + Rc - Rc)
                                    * (-DCentreCentre - Rc + Rc) * (-DCentreCentre + Rc + Rc)
                                    * (DCentreCentre + Rc + Rc)))
    else:
        CentreCentreWindowRadius = 0

    CentreCornerWindowRadius = float((1 / (2 * DCentreCorner)) * np.sqrt((-DCentreCorner + R - Rc)
                            * (-DCentreCorner - R + Rc) * (-DCentreCorner + R + Rc) * (DCentreCorner + R + Rc)))

    # finding the volume of the removed lens common to the centre centre spheres.
    # If the centre centre window radius is zero the lenz volume is zero.
    if CentreCentreWindowRadius == 0:
        CentreCentreLensVolume = 0
    else:
        CentreCentreLensVolume = float((np.pi * (4 * Rc + DCentreCentre) * ((2 * Rc - DCentreCentre) ** 2)) / 12)

    TotalCentreCentreLensVolume = float((CentreCentreLensVolume * InternalNumberOfCentreCentreWindows)
                                        + ((CentreCentreLensVolume * ExternalNumberOfCentreCentreWindows)/2))

    # finding the volume of the removed lens common to the centre corner spheres
    CentreCornerLensVolume = float(np.pi * ((Rc + R - DCentreCorner) ** 2) * (DCentreCorner ** 2 + 2 * DCentreCorner
                            * R- 3 * R ** 2 + 2 * DCentreCorner * Rc + 6 * R * Rc - 3 * Rc ** 2) / (12 * DCentreCorner))
    TotalCentreCornerLensVolume = CentreCornerLensVolume * NumberOfCentreCornerWindows

    # The total lens volume within the material.
    TotalLensVolume = TotalCentreCentreLensVolume + TotalCentreCornerLensVolume

    # Total volume of the voids within the cube
    CubeVoidVolume = (NumberOfCentreSpheres * (4 * np.pi * Rc ** 3) / 3) + (NumberOfCornerSpheres
                    * (4 * np.pi * R ** 3) / 3) - TotalLensVolume

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1) * (Y - 1) * (Z - 1) * BCCEquation ** 3)

    # Volume of porous material
    VolumeOfPorousMaterial = TotalCubeVolume - CubeVoidVolume

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)


    # Information for the text file


    TextFile.write('The centre centre window radius is ' + str(CentreCentreWindowRadius) + 'mm' + '\n')
    TextFile.write('The centre corner window radius is ' + str(CentreCornerWindowRadius) + 'mm' + '\n')
    TextFile.write('The single centre centre lens volume - ' + str(CentreCentreLensVolume) + 'mm^3' + '\n')
    TextFile.write('The single centre corner lens volume - ' + str(CentreCornerLensVolume) + 'mm^3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm^3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) + '\n')

    TextFile.close()

def BCCFEA04(Xmm, Ymm, Zmm, R, P, File):

    # Packages to be imported
    import numpy as np
    from sympy import symbols, pi, Eq, sqrt, solve

    U, W, V = symbols('U,W,V')

    U = P
    V = R

    eq = Eq(U, (-1.0 / 64.0 * pi * (-3.0 * W ** 4.0 + 12.0 * W ** 3.0 * V - 18.0 * W ** 2.0 * V ** 2.0 + V ** 4.0) * 3.0 ** (1.0 / 2.0) / V ** 4.0))
    X1 = solve(eq, W)

    Rc = float(X1[1])

    # BCCEquation packing equation
    BCCEquation = float((4 / (sqrt(3)))*R)

    # converting the dimensions in mm to number of spheres
    X = int(Xmm / BCCEquation) + 2  # number of spheres in the X direction
    Y = int(Ymm / BCCEquation) + 2  # number of spheres in the Y direction
    Z = int(Zmm / BCCEquation) + 2  # number of spheres in the Z direction


    # equation that calculates centre radius by using the corner radius and porosity


    # create text file with information in it

    TextFile = open(str(File) + 'BCCFEA04-'+ str(Xmm) + '-' + str(Ymm) + '-'
                    + str(Zmm) + '-' + str(R) + '-'+ str(P) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA04-'+ str(Xmm) + '-' + str(Ymm)
                          + '-' + str(Zmm) + '-' + str(R) + '-' + str(P) + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                               'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim corner spheres
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates centre spheres for BCC structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + Rc) + '), MM(' + str(
                        xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # closes the file just written
    SpaceClaimFile.close()


    # write model information to a text file
    TextFile.write('The Function used: BCCFEA04' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'X direction dimensions ' + str(Xmm) + 'mm'
                   + '\n' + 'Y direction dimensions ' + str(Ymm) + 'mm' + '\n'
                   + 'Z direction dimensions ' + str(Zmm) + 'mm' + '\n' + 'The corner sphere radius '
                   + str(R) + 'mm' + '\n' + 'The porosity is ' + str(P) + '\n' + '\n')

    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = float((X - 1.0) * BCCEquation)
    YDim = float((Y - 1.0) * BCCEquation)
    ZDim = float((Z - 1.0) * BCCEquation)

    TextFile.write('The actual dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = float((X - 1.0) * (Y - 1.0) * (Z - 1.0))
    NumberOfCornerSpheres = float((X - 1.0) * (Y - 1.0) * (Z - 1.0))
    NumberOfCentreCornerWindows = float(NumberOfCornerSpheres * 8.0)

    # The total number of centre centre windows changes with geometry as they are shared with connecting units.
    # In addition, the external windows need to be calculated separately as the lenz volume needs to be halfed


    NumberOfCentreCentreWindows = 0.0
    ExternalNumberOfCentreCentreWindows = 0.0
    InternalNumberOfCentreCentreWindows = 0.0

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) + '\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) + '\n')
    TextFile.write('The total number of centre centre windows - ' + str(NumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of internal centre centre windows - '
                   + str(InternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of external centre centre windows - '
                   + str(ExternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of centre corner windows - ' + str(NumberOfCentreCornerWindows) + '\n')

    # calculation of the lens between centre centre sphere interaction

    # Distance between centre centre spheres
    DCentreCentre = float(BCCEquation)
    DCentreCorner = float(2.0*R)

    # if the centre spheres do not overlap it returns a complex value. This If statemtent returns the Centre
    # Centre window radius as a zero if this is the case
    CentreCentreWindowRadius = 0

    CentreCornerWindowRadius = float((1.0 / (2.0 * DCentreCorner)) * sqrt((-DCentreCorner + R - Rc)
                                * (-DCentreCorner - R + Rc) * (-DCentreCorner + R + Rc) * (DCentreCorner + R + Rc)))

    # finding the volume of the removed lens common to the centre centre spheres.
    # If the centre centre window radius is zero the lenz volume is zero.

    CentreCentreLensVolume = 0.0


    TotalCentreCentreLensVolume = float((CentreCentreLensVolume * InternalNumberOfCentreCentreWindows) \
                                  + ((CentreCentreLensVolume * ExternalNumberOfCentreCentreWindows)/2))

    # finding the volume of the removed lens common to the centre corner spheres
    CentreCornerLensVolume = float(np.pi * ((Rc + R - DCentreCorner) ** 2.0) * (DCentreCorner ** 2.0 + 2.0 * DCentreCorner
                            * R- 3.0 * R ** 2.0 + 2.0 * DCentreCorner * Rc + 6.0 * R * Rc - 3.0 * Rc ** 2.0) / (12.0 * DCentreCorner))
    TotalCentreCornerLensVolume = float(CentreCornerLensVolume * NumberOfCentreCornerWindows)

    # The total lens volume within the material.
    TotalLensVolume = float(TotalCentreCentreLensVolume + TotalCentreCornerLensVolume)

    # Total volume of the voids within the cube
    CubeVoidVolume = float((NumberOfCentreSpheres * (4.0 * np.pi * Rc ** 3.0) / 3.0) + (NumberOfCornerSpheres
                    * (4.0 * np.pi * R ** 3.0) / 3.0) - TotalLensVolume)

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1.0) * (Y - 1.0) * (Z - 1.0) * BCCEquation ** 3.0)

    # Volume of porous material
    VolumeOfPorousMaterial = float(TotalCubeVolume - CubeVoidVolume)

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)


    # Information for the text file

    TextFile.write('The centre centre window radius is ' + str(CentreCentreWindowRadius) + 'mm' + '\n')
    TextFile.write('The centre corner window radius is ' + str(CentreCornerWindowRadius) + 'mm' + '\n')
    TextFile.write('The single centre centre lens volume - ' + str(CentreCentreLensVolume) + 'mm^3' + '\n')
    TextFile.write('The single centre corner lens volume - ' + str(CentreCornerLensVolume) + 'mm^3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm^3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) + '\n')

    TextFile.close()

def BCCFEA05(Xmm, Ymm, Zmm, R, File):

    # Packages to be imported
    import numpy as np

    # BCCEquation packing equation
    BCCEquation = (4 / (np.sqrt(3)))*R

    #Calculating the centre sphere radius so all windows (centrecentre and centre corner) are equal
    Rc = 1.27232156058011*R


    # converting the dimensions in mm to number of spheres
    X = int(Xmm / BCCEquation) + 2  # number of spheres in the X direction
    Y = int(Ymm / BCCEquation) + 2  # number of spheres in the Y direction
    Z = int(Zmm / BCCEquation) + 2  # number of spheres in the Z direction

    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA05-'
                    + str(Xmm) + '-' + str(Ymm) + '-' + str(Zmm) + '-' + str(R) +  '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA05-'
                                + str(Xmm) + '-' + str(Ymm) + '-' + str(Zmm) + '-' + str(R) + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                               'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim corner spheres
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates centre spheres for BCC structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + Rc) + '), MM(' + str(
                        xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # closes the file just written
    SpaceClaimFile.close()



    # write model information to a text file
    TextFile.write('The Function used: BCCFEA05' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'X direction dimensions ' + str(Xmm)
                   + '\n' + 'Y direction dimensions ' + str(Ymm) + '\n'
                   + 'Z direction dimensions ' + str(Zmm) + '\n' + 'The corner sphere radius '
                   + str(R) + 'mm' + '\n' + '\n' + '\n')

    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = (X - 1) * BCCEquation
    YDim = (Y - 1) * BCCEquation
    ZDim = (Z - 1) * BCCEquation

    TextFile.write('The actual dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCornerSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCentreCornerWindows = NumberOfCornerSpheres * 8

    # The total number of centre centre windows changes with geometry as they are shared with connecting units.
    # In addition, the external windows need to be calculated separately as the lenz volume needs to be halfed

    if (2 * Rc) > BCCEquation:
        NumberOfCentreCentreWindows = (X * ((Y-1) * (Z-1)) + Y * ((X-1) * (Z-1)) + Z * ((X-1) * (Y-1)))
        ExternalNumberOfCentreCentreWindows = ((X-1)*(Y-1)*2)+((Y-1)*(Z-1)*2)+((Z-1)*(X-1)*2)
        InternalNumberOfCentreCentreWindows = NumberOfCentreCentreWindows - ExternalNumberOfCentreCentreWindows
    else:
        NumberOfCentreCentreWindows = 0
        ExternalNumberOfCentreCentreWindows = 0
        InternalNumberOfCentreCentreWindows = 0

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) + '\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) + '\n')
    TextFile.write('The total number of centre centre windows - ' + str(NumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of internal centre centre windows - '
                   + str(InternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of external centre centre windows - '
                   + str(ExternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of centre corner windows - ' + str(NumberOfCentreCornerWindows) + '\n')

    # calculation of the lens between centre centre sphere interaction

    # Distance between centre centre spheres
    DCentreCentre = BCCEquation
    DCentreCorner = 2*R

    # if the centre spheres do not overlap it returns a complex value. This If statemtent returns the Centre
    # Centre window radius as a zero if this is the case
    if (2*Rc)>BCCEquation:
         CentreCentreWindowRadius = float((1 / (2 * DCentreCentre)) * np.sqrt((-DCentreCentre + Rc - Rc)
                                    * (-DCentreCentre - Rc + Rc) * (-DCentreCentre + Rc + Rc)
                                    * (DCentreCentre + Rc + Rc)))
    else:
        CentreCentreWindowRadius = 0

    CentreCornerWindowRadius = float((1 / (2 * DCentreCorner)) * np.sqrt((-DCentreCorner + R - Rc)
                            * (-DCentreCorner - R + Rc) * (-DCentreCorner + R + Rc) * (DCentreCorner + R + Rc)))

    # finding the volume of the removed lens common to the centre centre spheres.
    # If the centre centre window radius is zero the lenz volume is zero.
    if CentreCentreWindowRadius == 0:
        CentreCentreLensVolume = 0
    else:
        CentreCentreLensVolume = float((np.pi * (4 * Rc + DCentreCentre) * ((2 * Rc - DCentreCentre) ** 2)) / 12)

    TotalCentreCentreLensVolume = float((CentreCentreLensVolume * InternalNumberOfCentreCentreWindows)
                                        + ((CentreCentreLensVolume * ExternalNumberOfCentreCentreWindows)/2))

    # finding the volume of the removed lens common to the centre corner spheres
    CentreCornerLensVolume = float(np.pi * ((Rc + R - DCentreCorner) ** 2) * (DCentreCorner ** 2 + 2 * DCentreCorner
                            * R- 3 * R ** 2 + 2 * DCentreCorner * Rc + 6 * R * Rc - 3 * Rc ** 2) / (12 * DCentreCorner))
    TotalCentreCornerLensVolume = CentreCornerLensVolume * NumberOfCentreCornerWindows

    # The total lens volume within the material.
    TotalLensVolume = TotalCentreCentreLensVolume + TotalCentreCornerLensVolume

    # Total volume of the voids within the cube
    CubeVoidVolume = (NumberOfCentreSpheres * (4 * np.pi * Rc ** 3) / 3) + (NumberOfCornerSpheres
                    * (4 * np.pi * R ** 3) / 3) - TotalLensVolume

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1) * (Y - 1) * (Z - 1) * BCCEquation ** 3)

    # Volume of porous material
    VolumeOfPorousMaterial = TotalCubeVolume - CubeVoidVolume

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)


    # Information for the text file

    TextFile.write('The centre centre window radius is ' + str(CentreCentreWindowRadius) + 'mm' + '\n')
    TextFile.write('The centre corner window radius is ' + str(CentreCornerWindowRadius) + 'mm' + '\n')
    TextFile.write('The single centre centre lens volume - ' + str(CentreCentreLensVolume) + 'mm^3' + '\n')
    TextFile.write('The single centre corner lens volume - ' + str(CentreCornerLensVolume) + 'mm^3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm^3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) + '\n')

    TextFile.close()
def BCCFEA07(X, Y, Z, R, Rc1, Rc2, File):
    import numpy as np

    # BCCEquation packing equation

    BCCEquation = (4 / (np.sqrt(3))) * R

    # calculate the centre sphere difference between Rc1 and Rc2 along the Z direction
    if X > 2:
        RcInterval = (Rc2 - Rc1) / (X - 2)
    else:
        RcInterval = Rc1

    # populates an array of all the centre spheres radiusus
    xyz5 = ([])
    for i in range(X - 1):
        for ii in range((Z - 1) * (Y - 1)):
            xyz5 = np.append(xyz5, (np.array([Rc1 + i * RcInterval])))

    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA07-' + str(X) + '-' + str(Y) + '-' + str(Z)
                    + '-' + str(R) + '-' + str(Rc1) + '-' + str(Rc2) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA07-' + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-'
                          + str(Rc1) + '-' + str(Rc2) + '-' + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' + 'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write('result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str(
        (X - 1) * BCCEquation) + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation)
                         + ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim corner spheres
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1])
                                     + '), MM(' + str(xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM('
                                     + str(xyz1[1]) + '), MM(' + str(xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates centre spheres for BCC structures.

    # the count gives a count of all the centre spheres, This is used in collecting the sphere radius from the array
    # through the for loop. The count is reduced by one as the array has the initial element removed as it is a zero
    count = 0
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                count += 1
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM('
                                     + str(xyz1[1]) + '), MM(' + str(xyz1[2]) + ')), Point.Create(MM('
                                     + str(xyz1[0] + xyz5[count - 1]) + '), MM(' + str(xyz1[1]) + '), MM('
                                     + str(xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # closes the file just written
    SpaceClaimFile.close()

    ####################################################################################################################

    # text file information

    # write model information to a text file
    TextFile.write('The Function used: BCCFEA07' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'Number of spheres in the X direction ' + str(X) + '\n'
                   + 'Number of spheres in the Y direction ' + str(Y) + '\n' + 'Number of spheres in the Z direction '
                   + str(Z) + '\n' + 'The corner sphere radius ' + str(R) + 'mm' + '\n'
                   + 'The first set of centre sphere radiuses ' + str(Rc1) + 'mm' + '\n'
                   + 'The last set of centre sphere radiuses ' + str(Rc2) + 'mm' + '\n' + '\n')

    # calculating the cube dimensions
    XDim = (X - 1) * BCCEquation
    YDim = (Y - 1) * BCCEquation
    ZDim = (Z - 1) * BCCEquation

    TextFile.write('The dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = '
                   + str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # Calculating the centre sphere radiusus between Rc1 and Rc2

    CentreRadiuses = ([])
    for iii in range(X - 1):
        CentreRadiuses = np.append(CentreRadiuses, (np.array([Rc1 + iii * RcInterval])))

    # create for loop iiii that runs through the cnetrecentre array
    # iiii = [0]
    LensVolumeTotalArray = ([])
    CentreSphereVolumeArray = ([])
    CentreCornerWindowRadiusArray = ([])
    FBTBExternalRadiusArray = ([])
    LHSExternalLensVolumeTotalArray = ([])
    RHSExternalLensVolumeTotalArray = ([])

    # Calculating all the lens volumes

    for iiii in range(len(CentreRadiuses)):

        #  Lens volume calculation of the centre corner interaction.
        if R + CentreRadiuses[iiii] > 2 * R:
            TotalCentreCornerLensVolume = (np.pi * (CentreRadiuses[iiii] - R) ** 2 * (
                    -3 * CentreRadiuses[iiii] ** 2 + 10 * CentreRadiuses[iiii] * R + 5 * R ** 2) / (24 * R)) * (
                                                  Z - 1) * (Y - 1) * 8
            CentreCornerWindowRadius = (1.0 / (2.0 * 2.0 * R)) * np.sqrt(
                (-2.0 * R + R - CentreRadiuses[iiii]) * (-2.0 * R - R + CentreRadiuses[iiii]) * (
                        -2.0 * R + R + CentreRadiuses[iiii]) * (2.0 * R + R + CentreRadiuses[iiii]))
        else:
            TotalCentreCornerLensVolume = 0
            CentreCornerWindowRadius = 0
        LensVolumeTotalArray = np.append(LensVolumeTotalArray, TotalCentreCornerLensVolume)

        # Front back top bottom window (FBTB) external only
        if (2 * CentreRadiuses[iiii]) > BCCEquation:
            ExternalNumberOfFBTBWindows = (((Y - 1) + (Z - 1)) * 2)
            ExternalRadiusFBTBWindows = 0.5 * np.sqrt(
                BCCEquation ** 2 * (-BCCEquation ** 2 + 4 * CentreRadiuses[iiii] ** 2)) / BCCEquation
            TotalExternalFBTBLensVolume = (1.0 / 6.0) * np.pi * (-0.5 * BCCEquation + CentreRadiuses[iiii]) ** 2.0 * (
                    BCCEquation + 4.0 * CentreRadiuses[iiii]) * (2.0 * Y + 2.0 * Z - 4.0)
        else:
            ExternalNumberOfFBTBWindows = 0
            ExternalRadiusFBTBWindows = 0
            TotalExternalFBTBLensVolume = 0
        LensVolumeTotalArray = np.append(LensVolumeTotalArray, TotalExternalFBTBLensVolume)

        #  Front, back, top, bottom, window(FBTB), internal only
        if (2 * CentreRadiuses[iiii]) > BCCEquation:
            InternalNumberOfFBTBWindows = ((Y - 1) * (Z - 2)) + ((Z - 1) * (Y - 2))
            InternalRadiusFBTBWindows = 0.5 * np.sqrt(
                BCCEquation ** 2 * (-BCCEquation ** 2 + 4 * CentreRadiuses[iiii] ** 2)) / BCCEquation
            TotalInternalFBTBLensVolume = (1.0 / 3.0) * np.pi * (-0.5 * BCCEquation + CentreRadiuses[iiii]) ** 2.0 * (
                    BCCEquation + 4.0 * CentreRadiuses[iiii]) * ((Y - 2) * (Z - 1) + (Y - 1) * (Z - 2))
        else:
            InternalNumberOfFBTBWindows = 0
            InternalRadiusFBTBWindows = 0
            TotalInternalFBTBLensVolume = 0
        LensVolumeTotalArray = np.append(LensVolumeTotalArray, TotalInternalFBTBLensVolume)

        # Left Hand Side (LHS) window External only
        if (2 * CentreRadiuses[iiii]) > BCCEquation:
            ExternalNumberOfLHSWindows = (Z - 1) * (Y - 1)
            ExternalRadiusLHSWindows = 0.5 * np.sqrt(
                BCCEquation ** 2 * (-BCCEquation ** 2 + 4 * CentreRadiuses[iiii] ** 2)) / BCCEquation
            TotalExternalLHSLensVolume = np.pi * (BCCEquation - 2 * CentreRadiuses[iiii]) ** 2 * (
                    BCCEquation + 4 * CentreRadiuses[iiii]) * (Y - 1) * (Z - 1) / 24
        else:
            ExternalNumberOfLHSWindows = 0
            ExternalRadiusLHSWindows = 0
            TotalExternalLHSLensVolume = 0
        LHSExternalLensVolumeTotalArray = np.append(LHSExternalLensVolumeTotalArray, TotalExternalLHSLensVolume)

        #  Left Hand Side (LHS) windows internal only.
        if iiii + (iiii - 1) >= 0 and (CentreRadiuses[iiii] + CentreRadiuses[iiii - 1]) > BCCEquation:
            InternalNumberOfLHSWindows = (Z - 1) * (Y - 1)
            InternalRadiusLHSWindows = 0.5 * np.sqrt((-BCCEquation - CentreRadiuses[iiii] + CentreRadiuses[iiii - 1])
                                                     * (-BCCEquation + CentreRadiuses[iiii] - CentreRadiuses[iiii - 1])
                                                     * (-BCCEquation + CentreRadiuses[iiii] + CentreRadiuses[iiii - 1])
                                                     * (BCCEquation + CentreRadiuses[iiii] + CentreRadiuses[
                iiii - 1])) / BCCEquation
            TotalInternalLHSLensVolume = np.pi * (Y - 1) * (Z - 1) * (
                    -BCCEquation + CentreRadiuses[iiii] + CentreRadiuses[iiii - 1]) ** 2 * (
                                                 BCCEquation ** 2 + 2 * BCCEquation * CentreRadiuses[
                                             iiii] + 2 * BCCEquation * CentreRadiuses[iiii - 1] - 3 *
                                                 CentreRadiuses[iiii] ** 2 + 6 * CentreRadiuses[iiii] *
                                                 CentreRadiuses[iiii - 1] - 3 * CentreRadiuses[iiii - 1] ** 2) / (
                                                 12 * BCCEquation)

        else:
            InternalNumberOfLHSWindows = 0
            InternalRadiusLHSWindows = 0
            TotalInternalLHSLensVolume = 0
        LensVolumeTotalArray = np.append(LensVolumeTotalArray, TotalInternalLHSLensVolume)

        # Right Hand Side (RHS) windows External only CORRECt CHECKED WITH MATHCAD
        if (2 * CentreRadiuses[iiii]) > BCCEquation:
            ExternalNumberOfRHSWindows = (Z - 1) * (Y - 1)  # correct
            ExternalRadiusRHSWindows = 0.5 * np.sqrt(
                BCCEquation ** 2.0 * (-BCCEquation ** 2.0 + 4.0 * CentreRadiuses[iiii] ** 2.0)) / BCCEquation
            TotalExternalRHSLensVolume = np.pi * (BCCEquation - 2.0 * CentreRadiuses[iiii]) ** 2.0 * (
                    BCCEquation + 4.0 * CentreRadiuses[iiii]) * (Y - 1) * (Z - 1) / 24.0
        else:
            ExternalNumberOfRHSWindows = 0
            ExternalRadiusRHSWindows = 0
            TotalExternalRHSLensVolume = 0
        RHSExternalLensVolumeTotalArray = np.append(RHSExternalLensVolumeTotalArray, TotalExternalRHSLensVolume)

        # add all the lens volume elements

        # calculating the volume of the centre spheres
        CentreSphereVolumeArray = np.append(CentreSphereVolumeArray,
                                            ((Z - 1) * (Y - 1) * (CentreRadiuses[iiii]) ** 3 * (4.0 / 3.0) * np.pi))

    #  calculating the total lens volume to be removed from the spheres to obtain the void area
    TotalLensVolume = np.sum(LensVolumeTotalArray) + LHSExternalLensVolumeTotalArray[0] + \
                      RHSExternalLensVolumeTotalArray[(len(CentreRadiuses) - 1)]



    # calculate the volume of the centre spheres
    TotalCentreSphereVolume = np.sum(CentreSphereVolumeArray)

    # total volume of the corner spheres
    CornerSphereTotalVolume = (X - 1) * (Y - 1) * (Z - 1) * (4.0 / 3.0) * np.pi * R ** 3

    # Volume of cube without any pores in mm3
    TotalCubeVolume = (X - 1) * (Y - 1) * (Z - 1) * BCCEquation ** 3.0

    # calculating the volume of voids in the cube
    CubeVoidVolume = TotalCentreSphereVolume + CornerSphereTotalVolume - TotalLensVolume

    # calculate the porosity
    porosity = CubeVoidVolume / TotalCubeVolume

    # calculate the material to check against SpaceClaim
    TotalMaterial = TotalCubeVolume - CubeVoidVolume

    TextFile.write('Centre Sphere Radius along the X axis - ' + str(CentreRadiuses) + 'mm' + '\n')
    TextFile.write('Total solid cube volume - ' + str(TotalCubeVolume) + 'mm^3' + '\n')
    TextFile.write('Total void volume - ' + str(CubeVoidVolume) + 'mm^3' + '\n')
    TextFile.write('Porosity - ' + str(porosity) + '\n')

    SpaceClaimFile.close()
    TextFile.close()