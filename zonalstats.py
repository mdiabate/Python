
"""
Created on Thu Jan 05 14:08:26 2017
@author: Mouhamad
"""

import gdal, ogr, osr, numpy
import sys, os
from rasterstats import zonal_stats
import shp2kml
from math import ceil

gdal.UseExceptions()
osr.UseExceptions()

def grid(image, gridSize):

    #Open data
    raster = gdal.Open(image)
    proj = raster.GetProjection()

    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xmin = transform[0]
    ymax = transform[3]
    xmax = xmin + transform[1] * raster.RasterXSize
    ymin = ymax + transform[5] * raster.RasterYSize

    gridWidth = float(gridSize)
    gridHeight = float(gridSize)

    #Get Rows and columns
    rows = ceil((ymax - ymin)/gridHeight)
    cols = ceil((xmax - xmin)/gridWidth)

    #start grid shell
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax - gridHeight

    #create Output Files
    path = os.path.dirname(image)
    basename = os.path.splitext(os.path.basename(image))[0]
    fishnet = os.path.join(path, basename + '_ZonalStats.shp')

    if os.path.exists(fishnet):
        os.remove(fishnet)

    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    outDataSource = outDriver.CreateDataSource(fishnet)
    outLayer = outDataSource.CreateLayer('grid', geom_type = ogr.wkbPolygon)
    featureDefn = outLayer.GetLayerDefn()

    minfield = ogr.FieldDefn('MIN', ogr.OFTReal)
    outLayer.CreateField(minfield)

    maxfield = ogr.FieldDefn('MAX', ogr.OFTReal)
    outLayer.CreateField(maxfield)

    meanfield = ogr.FieldDefn('MEAN', ogr.OFTReal)
    outLayer.CreateField(meanfield)

    medianfield = ogr.FieldDefn('MEDIAN', ogr.OFTReal)
    outLayer.CreateField(medianfield)

    stdevfield = ogr.FieldDefn('STDEV', ogr.OFTReal)
    outLayer.CreateField(stdevfield)

    #create Grid Cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        #reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature = None

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    #create prj file for shapefile
    file = open(fishnet[:-4] + '.prj', 'w')
    file.write(proj)
    file.close()
    #save and close Datasources
    outDataSource = None
    return fishnet


#Function to fill in the constantarray with raster values
def blit(dest, src, loc):
    pos = [i if i >= 0 else None for i in loc]
    neg = [-i if i < 0 else None for i in loc]
    target = dest[[slice(i,None) for i in pos]]
    src = src[[slice(i, j) for i,j in zip(neg, target.shape)]]
    target[[slice(None, i) for i in src.shape]] = src
    return dest

def zstats(feat, fishnet, image):

    # Open data
    raster = gdal.Open(image)
    shp = ogr.Open(fishnet)
    lyr = shp.GetLayer()
    banddataraster = raster.GetRasterBand(1)
    ndv = banddataraster.GetNoDataValue()
    if ndv is None:
        ndv = -9999
    dataraster = numpy.array(banddataraster.ReadAsArray().astype(numpy.float32))

    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]

    #create masked array with extent of shapefile and values of raster
    x_min, x_max, y_min, y_max = lyr.GetExtent()
    cols = abs(int((x_max - x_min) / pixelHeight)+1)
    rows = abs(int((y_max - y_min) / pixelWidth)+1)
    array_9999 = numpy.full((rows, cols), ndv, dtype=numpy.float32)

    #fill in values from the raster
    arraytess = blit(array_9999, dataraster, (0,0))

    # Reproject vector geometry to same projection as raster
    sourceSR = lyr.GetSpatialRef()
    targetSR = osr.SpatialReference()
    targetSR.ImportFromWkt(raster.GetProjectionRef())
    coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)

    # Get extent of feat
    geom = feat.GetGeometryRef()
    geom.Transform(coordTrans)
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []; pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                    lon, lat, z = ring.GetPoint(p)
                    pointsX.append(lon)
                    pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []; pointsY = []
        for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)
    else:
        sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")

    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1
    xend = xoff + xcount
    yend = yoff + ycount

    #extract by feat extent
    zoneraster = arraytess[yoff:yend, xoff:xend]
    zoneraster[zoneraster==ndv]=numpy.nan

    stats = [float(numpy.nanmin(zoneraster)),
            float(numpy.nanmax(zoneraster)),
            float(numpy.nanmean(zoneraster)),
            float(numpy.nanmedian(zoneraster)),
            float(numpy.nanstd(zoneraster))]

    return stats

    #return zoneraster #stats

def main(image, gridSize=2):
    print ' ....Starting Zone Stats calculation...'
    fishnet = grid(image, gridSize)
    shp = ogr.Open(fishnet, update = True)
    lyr = shp.GetLayer()

    #stats = zonal_stats(fishnet, image)

    featList = range(lyr.GetFeatureCount())

    #Loop through each feature of the shapefile

    for FID in featList:

        feature = lyr.GetFeature(FID)

        #print feature

        stats = zstats(feature, fishnet, image)
        #print stats

        feature.SetField("MIN", stats[0])
        feature.SetField("MAX", stats[1] )
        feature.SetField("MEAN", stats[2])
        feature.SetField("MEDIAN", stats[3])
        feature.SetField("STDEV", stats[4])

        lyr.SetFeature(feature)
        #print  'COMPLETED'
        feature = None
    shp = None
    print 'Zonal Stats shapefile created '

    shp2kml.main(fishnet)
    print 'KML shapefile created '

if __name__ == "__main__":

    image = r"C:\Data\BGNIR\Corn\output\BGNIR_Corn_75lat_75fron_ENDVI.tif"
    #shp = None

    main(image, gridSize = 150)
