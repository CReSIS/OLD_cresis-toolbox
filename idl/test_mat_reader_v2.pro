; Example batch script for reading Matlab -v7.3 CReSIS echogram files
;
; Modify the "fn" variable below to point to your CSARP echogram
; @test_mat_reader_v2

; ========================================================
; Setup and compile
.reset_session

; ========================================================
; Load data
fn = '/cresis/snfs1/dataproducts/public/data/snow/2017_Greenland_P3/CSARP_deconv/20170309_01/Data_20170309_01_100.mat'

file_id = H5F_OPEN(fn)

dataset_id_data = H5D_OPEN(file_id, '/Data')
data = H5D_Read(dataset_id_data)
H5D_CLOSE, dataset_id_data
dataset_id_lat = H5D_OPEN(file_id, '/Latitude')
lat = H5D_Read(dataset_id_lat)
H5D_CLOSE, dataset_id_lat
dataset_id_lon = H5D_OPEN(file_id, '/Longitude')
lon = H5D_Read(dataset_id_lon)
H5D_CLOSE, dataset_id_lon
dataset_id_elev = H5D_OPEN(file_id, '/Elevation')
elev = H5D_Read(dataset_id_elev)
H5D_CLOSE, dataset_id_elev
dataset_id_gps_time = H5D_OPEN(file_id, '/GPS_time')
gps_time = H5D_Read(dataset_id_gps_time)
H5D_CLOSE, dataset_id_gps_time

H5F_CLOSE, file_id

; ========================================================
; Plot Data
size_data = size(data)
device,decomposed=0
window,0,retain=2,xsize=size_data[2],ysize=size_data[1]
wset,0
color_table = 255-[[indgen(255)],[indgen(255)],[indgen(255)]]
tvlct,color_table
tvscl,transpose(10*alog10(data)),/order

; ========================================================
; Plot geographic information
window,1,retain=2,xsize=500,ysize=500
wset,1
tvlct,255-color_table
plot,lon,lat,xrange=[min(lon),max(lon)],yrange=[min(lat),max(lat)],xtitle='Longitude (deg)',ytitle='Latitude (deg)'

window,2,retain=2,xsize=500,ysize=300
wset,2
plot,elev,yrange=[min(elev),max(elev)],xtitle='Index',ytitle='Elevation (m)'

window,3,retain=2,xsize=500,ysize=300
wset,3
plot,gps_time,yrange=[min(gps_time),max(gps_time)],xtitle='Index',ytitle='GPS Time (sec)',title='Seconds since Jan 1, 1970 00:00:00'
