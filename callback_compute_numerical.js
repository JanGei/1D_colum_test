function transport_num(c,c_in,disp,dx,dt){
    // This function approximates advection and dispersion numerically by a FVM
    
    // Move concentration by 1 cell
    for (let i = 1; i < c.length; i++) {
      c[c.length-i] = c[c.length-i-1]
    }
    // First cell gets inlet concentration
    c[0] = c_in
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

function total_conc(c,rg_SType,s,rho_s,poros,K_Fr,Fr_n) {
  // This function sums up the concentration in the aqueous and solid phase per cell
  var c_tot = Array(c.length).fill(0)
  if (rg_SType == 2) {
    for (let i = 0; i < c.length; i++) {
      c_tot[i] = (rho_s*(1-poros)*K_Fr*c[i]**(Fr_n-1)+ poros)*c[i]
    }
  } else {
    for (let i = 0; i < c.length; i++) {
      c_tot[i] = s[i]*rho_s*(1-poros) + poros*c[i]
    }
  }
  
  return c_tot
}

function zeros(dimensions) {
  // Credit: https://newbedev.com/javascript-initialize-2d-array-with-0-code-example
  var array = [];

  for (var i = 0; i < dimensions[0]; ++i) {
      array.push(dimensions.length == 1 ? 0 : zeros(dimensions.slice(1)));
  }

  return array;
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
const poros     = poros_sl.value;                   // [-]
const t_inj     = pulse_inj_sl.value                // [s]
var xBTC        = x3[0];                            // [m]
const rho_s     = rho_s_sl.value                    // [kg/m3]
const Kd        = Kd_sl.value                       // [m3/kg]
const K_ads     = Kads_sl.value                     // [mol/m3]
const s_max     = s_max_sl.value                    // [mol/kg]
const K_Fr      = K_Fr_sl.value                     // [mmol^(1-n)*L^(n) / kg]
const Fr_n      = Fr_n_sl.value                     // [-]

// Derived entities
const A       = math.PI * rad**2;             // [m2]
const vel     = Q/A;                          // [m/s]
const sep_vel = vel / poros                   // [m/s]
const reac    = (reac_l + reac_h)/2           // [1/s] 
const Dis     = (disp_l + disp_h)/2           // [m2/s]  
const PS      = col_len * A * poros           // [m3]
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

// Initialize lists
var c_array     = Array(x_num.length).fill(0)
var s_array     = Array(x_num.length).fill(0)
var c_tot_array = zeros([t_num.length,x_num.length]).fill(0)
var s_tot_array = zeros([t_num.length,x_num.length]).fill(0)


for (let i = 0; i < 6; i++) { // needs to be t_num.length
  // Set inlet concentration 
  if (rg_CP == 0) {
      var c_in = c0
  } else if (rg_CP == 1 && t_num[i] < t_inj) {
      var c_in = c0
  } else {
      var c_in = 0
  }

  // Transport
  //console.log(c_array)
  c_array = transport_num(c_array,c_in,Dis,dx,dt)
  //console.log(c_array)

  // Sorption
  if (rg_SType == 0) { // Linear Sorption
      var c_tot_lin = total_conc(c_array,rg_SType,s_array,rho_s,poros)

      for (let j = 0; j < c_array.length; j++) {
        c_array[j] = c_tot_lin[j]/(Kd*(1-poros)*rho_s+poros)
        s_array[j] = c_array[j]*Kd
      }

  } else if (rg_SType == 1) { // Langmuir Sorption
      var c_tot_lang = total_conc(c_array,rg_SType,s_array,rho_s,poros)

      for (let j = 0; j < c_array.length; j++) {
        var beta_lang  = (1-poros)*rho_s*s_max + poros*K_ads - c_tot_lang[j]
        var gamma_lang = -c_tot_lang[j]*K_ads

        c_array[j] = (-beta_lang + math.sqrt(beta_lang**2 - 4*poros*gamma_lang))/(2*poros)
        s_array[j] = s_max*c_array[j]/(K_ads+c_array[j]) //[mol/kg]
      }

  } else if (rg_SType == 2) { // Freundlich Sorption
      var c_tot_Fr = total_conc(c_array,rg_SType,s_array,rho_s,poros,K_Fr,Fr_n)
      //console.log(c_array)
      for (let j = 0; j < c_array.length; j++) {
        // Picard Iteration: Guess that everything is in the aqueous phase
        var c_old = c_tot_Fr[j]/poros
        // Loop until convergence criterion is met
        while (math.abs(c_old-c_array[j])>1e-9) {
          c_old = c_array[j]
          console.log(c_array[j])
          c_array[j] = c_tot_Fr[j]/(rho_s*(1-poros)*K_Fr*c_array[j]**(Fr_n-1)+poros)
          console.log(c_tot_Fr[j],rho_s,poros,K_Fr,Fr_n,poros)
          console.log(c_array[j])
        }
        s_array[j] = K_Fr * c_array[j]**Fr_n //[mmol/kg] --> different from other sorption types
      }
      //console.log(c_array)
  }

  // Store results
  for (let j = 0; j < c_array.length; j++) {
    c_tot_array[i,j] = c_array[j]
    s_tot_array[i,j] = s_array[j]
  }
  if (i==5){
    console.log(c_array)
  }
}

console.log(c_tot_array[5,0],c_tot_array[5,1],c_tot_array[5,2],c_tot_array[5,3],c_tot_array[5,4])
