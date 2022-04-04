function transport_num(c,disp,dx,dt){
    // This function approximates advection and dispersion numerically by a FVM
  
    // Move concentration by 1 cell
    for (let i = 0; i < c.length-1; i++) {
      c[i+1] = c[i]
    }
    // First cell gets inlet concentration
    c[0] = 1
    var Jd = []
    // Dispersive fluxes between cells
    for (let i = 0; i < c.length+1; i++) {
      if (i == 0) {
        Jd[i] = 0
      } else if (i == c.length) {
        Jd[i] = Jd[i-1]
      } else {
        Jd[i] = (c[i-1] - c[i])/dx*disp
      }
    }
  
    for (let i = 0; i < Jd.length-1; i++) {
      c[i] = c[i] + dt/dx * (Jd[i] - Jd[i+1])
    }
  
    return c
  }

// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var ymin= source1.data['ymin']
var ymax= source1.data['ymax']
var x2  = source2.data['x2']
var y2  = source2.data['y2']
var x3  = source3.data['xBTC']
var y3  = source3.data['yBTC']

var rg_CP     = rg_CP.active                      // [0]
var rg_SType  = rg_ST.active                      // [0]

// Values needed for all models
const col_len   = col_len_sl.value;                 // [m]
const rad       = col_rad_sl.value;                 // [m]
const reac_l    = Math.exp(reac_sl.value[0])/3600;  // [1/s]
const reac_h    = Math.exp(reac_sl.value[1])/3600;  // [1/s]
const disp_l    = Math.exp(disp_sl.value[0])/3600;  // [m2/s]
const disp_h    = Math.exp(disp_sl.value[1])/3600;  // [m2/s]
const Q         = flow_sl.value/1000/1000/3600;     // [m3/s]
const n         = poros_sl.value;                   // [-]
const t_inj     = pulse_inj_sl.value                // [s]
var xBTC        = x3[0];                            // [m]

// Derived entities
const A       = math.PI * rad**2;             // [m2]
const vel     = Q/A;                          // [m/s]
const sep_vel = vel / n                       // [m/s]
const reac    = (reac_l + reac_h)/2           // [1/s] 
const Dis     = (disp_l + disp_h)/2           // [m2/s]  
const PS      = col_len * A * n               // [m3]
const PV      = col_len/sep_vel               // [s] VEL oder SEP_VEL?
const c0      = 1;                            // [-] 

// Time span list
var tsp = []

// Discretize space (upper plot) and time (lower plot)
for (let j = 0; j < x.length; j++) {
  x[j] = -0.02*col_len + 1.02*col_len/x.length * j;
}
for (let j = 0; j < x2.length; j++) {
  tsp[j] = x2[j] * PV;
}

// x.length is number of spatial nodes
// x2.length is number of temporal nodes
// Idea: Work with independant data sources for numerical model and take closest value between 2 points and map it to the source


// This discretization follows the scheme of the GW transport code
// Requirement 1: Neumann-Number = 1/3
// Requirement 2: Courant-Number = 1
const dx    = 3*Dis/sep_vel
const dt    = 3*Dis/sep_vel**2
const t_end = PV * (x2.length-1)
  
// Lists for spatial and temporal span
var x_num = []
var t_num = []

for (let i = 0; ((i+1)*dx) < col_len; i++) {
    x_num[i] = (i+1)*dx
}
for (let i = 0; ((i+1)*dt) < t_end; i++) {
    t_num[i] = (i+1)*dt
}
console.log("Spatial Discretization is " + dx + " m\n"
              +"Temporal Discretization is " + dt + " s\n" 
              +"Seepage velocity is "+ sep_vel + "m/s\n"
              +"Dispersion Coefficient is "+ Dis + "m2/s\n"
              +"There are " + x_num.length + " spatial nodes\n"
              +"There are " + t_num.length + " temporal nodes\n"
              +"The Courant Number equals " + dt * sep_vel /dx
)

// Initialize list of concentrations
var c_array = Array(x_num.length).fill(0)
console.log(c_array)


for (let i = 0; i < t_num.length; i++) {
    // Set inlet concentration 
    if (rg_CP == 0) {
        var c_in = c0
    } else if (rg_CP == 1 && t_num[i] < t_inj) {
        var c_in = c0
    } else {
        var c_in = 0
    }

    // Transport by 1 cell
    //c_array = transport_num()
}
if (rg_SType == 0) { //Linear Sorption

}