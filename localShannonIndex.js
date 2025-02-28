//This code is developed to run in the enviroment of Google Earth Engine code editor. 
// Define the geometry (in geographic coordinates) base in the polygon define to the study area where the analysis will be performed. 
var areaOfInterest = ee.Geometry.Polygon([
  [-1.529085, 46.87], [-2.229717, 47.50492] ....
]);
// Load the Dynamic World dataset for the year 2022. This dataset provides global land cover //classifications
var dataset = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterDate('2022-01-01', '2022-12-31') // Selects images from 2022
  .filterBounds(areaOfInterest); // Filters only images within the study area
// Extract the land cover classification band and compute the most common class over time
var landCover = dataset.select('label').mode().clip(areaOfInterest);
// Apply a mask to remove pixels with no data (values of 0 indicate missing data)
var landCoverMasked = landCover.updateMask(landCover.gt(0)).unmask(0); 
// Fills no-data pixels with 0
// Define a kernel (neighborhood window) for spatial analysis
var kernelSize = 7; // size of the moving window 
var kernel = ee.Kernel.square(kernelSize / 2 + 3, 'pixels', false); // Adds a 3-pixel buffer to the window
// Function to calculate the Shannon diversity index in a moving window
var calculateShannonOptimized = function(image) {
  // Compute the total number of pixels in the moving window
  var totalPixels = ee.Image(1).reduceNeighborhood({
    reducer: ee.Reducer.sum(), 
    kernel: kernel
  });
  // Compute the proportion of each land cover class within the moving window
  // This generates a set of images where each band represents the proportion of a land cover class
  var proportions = ee.Image.cat([
    image.eq(0).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_0'),
    image.eq(1).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_1'),
    image.eq(2).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_2'),
    image.eq(3).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_3'),
    image.eq(4).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_4'),
    image.eq(5).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_5'),
    image.eq(6).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_6'),
    image.eq(7).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_7'),
    image.eq(8).reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: ker-nel}).divide(totalPixels).rename('class_8')
  ]);
  // Compute the Shannon diversity index using the equation.
  var shannonIndex = proportions.expression(
    '-1 * (class_0 * log(class_0 + 1e-6) + ' +
    'class_1 * log(class_1 + 1e-6) + ' +
    'class_2 * log(class_2 + 1e-6) + ' +
    'class_3 * log(class_3 + 1e-6) + ' +
    'class_4 * log(class_4 + 1e-6) + ' +
    'class_5 * log(class_5 + 1e-6) + ' +
    'class_6 * log(class_6 + 1e-6) + ' +
    'class_7 * log(class_7 + 1e-6) + ' +
    'class_8 * log(class_8 + 1e-6))', {
      'class_0': proportions.select('class_0'),
      'class_1': proportions.select('class_1'),
      'class_2': proportions.select('class_2'),
      'class_3': proportions.select('class_3'),
      'class_4': proportions.select('class_4'),
      'class_5': proportions.select('class_5'),
      'class_6': proportions.select('class_6'),
      'class_7': proportions.select('class_7'),
      'class_8': proportions.select('class_8')
    });
  // Masking to ensure only valid areas are considered
  return shannonIndex.updateMask(totalPixels.gt(0)).max(0);
};
// Generate a grid over the study area to divide it into smaller analysis units
var scale = 50000; // Defines the size of each block (50 km)
var projection = ee.Projection('EPSG:3035'); // Uses the LAEA projection system
var grid = areaOfInterest.coveringGrid({
  proj: projection,
  scale: scale
});
// Convert the grid into a FeatureCollection
var gridFeatures = ee.FeatureCollection(
  grid.map(function(geometry) {
    return ee.Feature(geometry);
  })
);
// Display the grid on the map
Map.centerObject(areaOfInterest, 10);
Map.addLayer(gridFeatures, {}, 'Grid');
// Process each block separately and export results
gridFeatures.toList(gridFeatures.size()).evaluate(function(tiles) {
  tiles.forEach(function(tile, index) {
    var block = ee.Feature(tile).geometry();
    // Add a buffer to the block to minimize edge effects
    var bufferedBlock = block.buffer(scale / 10); // 10% buffer
    // Clip the land cover map to the block area
    var landCoverMaskedBlock = landCoverMasked.clip(bufferedBlock);

    // Compute the Shannon index for the block
    var shannonResultBlock = calculateShannonOptimized(landCoverMaskedBlock);
    // Clip the result to match the original block area
    var croppedResultBlock = shannonResultBlock.clip(block);
    // Export the computed Shannon index as a GeoTIFF file
    Export.image.toDrive({
      image: croppedResultBlock,
      description: 'Shannon_Index_Block_' + index,
      folder: 'EarthEngine_Exports',
      scale: 10, // Uses 10-meter resolution
      region: block,
      crs: 'EPSG:3035', 
      maxPixels: 1e12 
    });
  });
});
