fibers = {'A': 'SCI',
          'B': 'FP'}
slice_names = {
    0 : '_SLICE1',
    1 : '_SLICE2'

}
data_format = "L2"#Can either be original or L2
slices= [0,1]
cam_names = {
    '1' : '_C1',
    '2' : '_C2'
}
cam_range = {
    '1' : [0, 90],
    '2' : [90, 170]
}
cam_num = 2
extnames = {
    'SCIDATA' : '_FLUX',
    'ERRDATA' : '_VAR',
    'WAVEDATA_VAC_BARY' : '_WAVE',
    'QUALDATA' : '_QUALDATA',
    'DLLDATA_VAC_BARY' : '_DISP',

}