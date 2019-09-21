/*
  @author Martin Jung - m.jung@sussex.ac.uk

Full extraction routine of EVI composite values
for all PREDICTS MLI buffered sites.
To be run in Google Earth Engine!

Desc: 
  Uses cFMASK filtered ( clouds (4) and shadows (2) out )
Landsat TOA estimates to calculate EVI
Reduces all estimates foreach given polygon.
*/
var l8_led = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
l4_led = ee.ImageCollection("LANDSAT/LT04/C01/T1_SR"),
l5_led = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
l7_led = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR");
// From Forum: Timeout tile 600 seconds per feature, 2 hours per aggregation and 10 minutes per tile.
// Coordinates and parameters 
// This points to a an asset table th euser has to upload. Alternatively a shapefile (replace line below)
var fullsites = ee.FeatureCollection('ft:1Oet2yGWvldNoVx8A6ZlCi75Ks_lrEG0dG9oD6k2j'); // All sites

// Filter full_sites radial to only those considered
fullsites = ee.FeatureCollection(fullsites).filterMetadata('Mx_lnr_','not_greater_than',3000);
var distinctSS = fullsites.distinct(['SS']);

//Map.addLayer(fullsites , {color: 'F40500'});

//////// Other Parameters //////////
var keep_properties = ['SS','SSBS']; // name of properties to keep from your feature such as id, class etc.
var export_geometry = false; // export the geometry in the csv
var scale = 30; // Resolution over which image collection should be reduced, 30m = native scale, 100 ~ avgML
var monthlyComp = false; // Make monthly composites
var startyear = 1982; //when to start
var endyear = 2014; // when to end
// max(), mean(), median(), min(), mode(), or(), product(), sum(), stdDev()
var reduce = ee.Reducer.mean();//.unweighted(); // Reducer for image collection
var what = "EVI_LEDAPS";
var folder = 'LandsatPREDICTS_circ';
var fm = "cfmask";
// Expressions
var f_evi = '2.5 * ((nir - red) / (nir + 2.4 * red + 1))'; // EVI2 formula
var f_ndmi = '((nir - swir) / (nir + swir) )'; // NDMI
// ################################################################### //
//                Function and main CODE starts here                   //
// ################################################################### //
//Function to mask excess EVI values defined as > 1.2 and < -1.0
var maskExcess = function(image) {
  var hi = image.lte(1.01);
  var lo = image.gte(-1.01);
  var masked = image.mask(hi.and(lo));
  return image.mask(masked);
};
// Function to remove clouds - expects the new SR data to have a cfmask layer
// 126122017 - Adapted to work with LT1
var maskclouds = function(scene) {
  // The "pixel_qa" band has various flags encoded in
  // different bits.  We extract some of them as individual mask bands.
  // Note: Cloud masking information is present in "pixel_qa" 
  // pixel_qa Bit 1: Clear pixel indicator.
  // pixel_qa Bit 2: Water indicator.
  // pixel_qa Bit 3: Cloud shadow indicator.
  // pixel_qa Bit 5: Cloud indicator.
  // pixel_qa Bits 6-7: Cloud confidence.
  // Fill = https://explorer.earthengine.google.com/#detail/LANDSAT%2FLE07%2FC01%2FT1_SR
    var clear = scene.select('pixel_qa').bitwiseAnd(2).neq(0);
    clear = scene.updateMask(clear);
    return(clear);
};
// Create TimeBand
function createTimeBand(image) {
  return image.addBands(image.metadata('system:time_start'));
}
// Remove bad edges around collections (- 1km)
var buffEdges = function(image){
  return image.clip(image.geometry().buffer(-1000));
};

var dateRanges = ee.List.sequence(0, ((endyear - startyear)*12)-1,1);
function makeMonthly(num) {
  // Move the start up by 32 days for each period, making sure to be in front of the start of the period
  var startDate = ee.Date(startyear+'-01-01');
  //var start = startDate.advance(ee.Number(num).multiply(31).subtract(1), 'day');
  // Move to the other end of the 32 day period
  //var end = start.advance(ee.Number(32),'day');
  var start = startDate.advance(ee.Number(num),'month');
  // advance by a month and then one day back to get full month
  var end = start.advance(ee.Number(1), 'month').advance(ee.Number(-1),'day');
  // Filter to the date range
  var filtered = Collection.filterDate(start, end);
  // Get the max of the selected period.
  var monthlyMeans = filtered.max();
  // Set the metadata date for this composite as the start of the month
  monthlyMeans = monthlyMeans.set('system:time_start',start.millis());
  return monthlyMeans;
}

// -------------------------------- //
// First Loop # BEST NOT TO PROCESS IT ALL AT ONCE. MAKE GROUPS
// To run, change increase run by one, execute and export!
// TBD: Convert into a python script for automatic task export
// RUN = Group Name
var runs = ee.List.sequence(1,distinctSS.size(),1); // Full list
var run = 0;// The run definining id. 
var batchsize = 20 ;// Go in batches of 20 (Limit of concurrent running tasks), thus 34+1 times total
var fails = ['11_227','11_233']; // The studies failed during export
  for (var f = run*batchsize ; f <= ((run+1)*batchsize)-1; f += 1){
    var ff = distinctSS.toList(1000).get(f); // Get Study SS
    var SS = ee.Feature(ff).get('SS');
    // Filter to sites
    var sites = fullsites.filter(ee.Filter.eq("SS",SS)).select(['SS','SSBS']);
    
    // ----------------------------------------- //
      // Filter the layers and set bounds
    // get the LC8 collection
    var L8 = l8_led
    //    .filterDate(year+'-01-01', year+'-12-31')
    .filterDate('1970-01-01', '2014-01-01') // Filter to up to latest sampling
    .filterBounds(sites)
    //    .filter(ee.Filter.lt('WRS_ROW', 122)) // Filter night times out -> http://tinyurl.com/hp6rap9
    .map(maskclouds)
    .map(createTimeBand);
    
    // Clip to feature collection geometry
    L8 = L8.map(function(i){return i.clipToCollection(sites);});    
    // get the LE7 collection
    var L7 = l7_led
    //    .filterDate(year+'-01-01', year+'-12-31')
    .filterDate('1970-01-01', '2014-01-01') // Filter to up to latest sampling
    .filterBounds(sites)
    //    .filter(ee.Filter.lt('WRS_ROW', 122)) // Filter night times out -> http://tinyurl.com/hp6rap9
    .map(maskclouds)
    .map(createTimeBand);
    
    // Clip to feature collection geometry
    L7 = L7.map(function(i){return i.clipToCollection(sites);});    
    
    // get the LE5 collection
    var L5 = l5_led
    //    .filterDate(year+'-01-01', year+'-12-31')
    .filterDate('1970-01-01', '2014-01-01') // Filter to up to latest sampling
    .filterBounds(sites)
    //    .filter(ee.Filter.lt('WRS_ROW', 122)) // Filter night times out -> http://tinyurl.com/hp6rap9
    .map(maskclouds)
    .map(createTimeBand);
    // Clip to feature collection geometry
    L5 = L5.map(function(i){return i.clipToCollection(sites);});    
    
    // get the LE5 collection
    var L4 = l4_led
    //    .filterDate(year+'-01-01', year+'-12-31')
    .filterDate('1970-01-01', '2014-01-01') // Filter to up to latest sampling
    .filterBounds(sites)
    //    .filter(ee.Filter.lt('WRS_ROW', 122)) // Filter night times out -> http://tinyurl.com/hp6rap9
    .map(maskclouds)
    .map(createTimeBand);
    // Clip to feature collection geometry
    L4 = L4.map(function(i){return i.clipToCollection(sites);});    
    
    
    // Landsat 8 EVI
    // SR data needs to be multiplied with 1000 before calculation !!!!
      var L8_evi = L8.map(
        function(image) {
          var evi = image.expression(
            f_evi,
            {
              red: image.select('B4').multiply(0.0001),    // 620-670nm, RED
              nir: image.select('B5').multiply(0.0001),    // 841-876nm, NIR
              swir: image.select('B6').multiply(0.0001),
              blue: image.select('B2').multiply(0.0001)    // 459-479nm, BLUE
            });
          // Rename that band to something appropriate
          var dimage = ee.Date(ee.Number(image.get('system:time_start'))).format();
          return evi.select([0], ['evi']).set({'datef': dimage,'landsat': 'L8','system:time_start': ee.Number(image.get('system:time_start'))});
          //    .addBands(image.metadata('system:time_start'));
        }
      );
    
    // calculate LE7 EVI
    var L7_evi = L7.map(
      function(image) {
        var evi = image.expression(
          f_evi,
          {
            red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
            nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
            swir: image.select('B5').multiply(0.0001),
            blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
          });
        // Rename that band to something appropriate
        var dimage = ee.Date(ee.Number(image.get('system:time_start'))).format();
        return evi.select([0], ['evi']).set({'datef': dimage,'landsat': 'L7','system:time_start': ee.Number(image.get('system:time_start'))});
        //    .addBands(image.metadata('system:time_start'));
      }
    );
    
    // calculate LE5 EVI
    var L5_evi = L5.map(
      function(image) {
        var evi = image.expression(
          f_evi,
          {
            red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
            nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
            swir: image.select('B5').multiply(0.0001),
            blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
          });
        // Rename that band to something appropriate
        var dimage = ee.Date(ee.Number(image.get('system:time_start'))).format();
        return evi.select([0], ['evi']).set({'datef': dimage,'landsat': 'L5','system:time_start': ee.Number(image.get('system:time_start'))});
        //    .addBands(image.metadata('system:time_start'));
      }
    );
    
    // Landsat 4
    var L4_evi = L4.map(
      function(image) {
        var evi = image.expression(
          f_evi,
          {
            red: image.select('B3').multiply(0.0001),    // 620-670nm, RED
            nir: image.select('B4').multiply(0.0001),    // 841-876nm, NIR
            swir: image.select('B5').multiply(0.0001),
            blue: image.select('B1').multiply(0.0001)    // 459-479nm, BLUE
          });
        // Rename that band to something appropriate
        var dimage = ee.Date(ee.Number(image.get('system:time_start'))).format();
        return evi.select([0], ['evi']).set({'datef': dimage,'landsat': 'L4','system:time_start': ee.Number(image.get('system:time_start'))});
        //    .addBands(image.metadata('system:time_start'));
      }
    );
    
    
    // Merge the EVI collections
    // this collection is sorted by time 
    var Collection = ee.ImageCollection(L8_evi.merge(L7_evi))
    .sort('system:time_start',true);
    Collection = ee.ImageCollection(Collection.merge(L5_evi))
    .sort('system:time_start',true);
    Collection = ee.ImageCollection(Collection.merge(L4_evi))
    .sort('system:time_start',true);
    
    //print(Collection);
    if(monthlyComp){
      // Make a list of months to be applied
      Collection = ee.ImageCollection(dateRanges.map(makeMonthly));
    }
    
    // Run function
    var extract = function(img,sites){
      var extracted = img.reduceRegions(sites, reduce, scale);
      return extracted;
    };
    
    var results = Collection.map(function(img){
      return extract(img,sites);
    }).flatten();
    
    // This will create a table with one row per feature and all dates in the image stack
    Export.table.toDrive({
      collection: results.select(["SSBS","mean","system:time_start","id"], null, export_geometry),
      folder: folder,
      description: run+'_'+f+'__'+what+'_MEAN_30_POLY',
      fileNamePrefix: run+'_'+f+'__'+what+'_MEAN_30_POLY',
      fileFormat: 'geoJSON'
    });
    
  }