function getc(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg) {
  var c = []
  var cmin = []
  var cmax = []
  var tau = []
  for (let i = 0; i < x.length; i++) { 
      if (x[i] <= 0) {
        c[i] = 1
        cmin[i] = 1
        cmax[i] = 1
      } else {
        var intlist = []
          for (let j = 0; j < Lcube1.length; j++) {
            var r_intermed = reac_l + (reac_h-reac_l)*Lcube1[j]
            var D_intermed = disp_l + (disp_h-disp_l)*Lcube2[j]
            var gam_intermed = get_gamma(r_intermed,D_intermed,vel)
            if (rg == 1  && t>t_inj) {
              // Pulse injection
              if (t<t_inj) {
                tau = 0
              } else {
                tau = t_inj
              }
              intlist[j] = 1/2 * ((1-math.erf((x[i]-vel*t) / math.sqrt(4*D_intermed*t)))-(1-math.erf((x[i]-vel*(t-tau))/math.sqrt(4*D_intermed*(t-tau)))));
            } else {
              // Cotinuous injection
              intlist[j] = 1/2 * Math.exp(x[i]*vel/(2*D_intermed))*(Math.exp((-x[i])*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]-vel*t*gam_intermed)/math.sqrt(4*D_intermed*t)))+math.exp(x[i]*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]+vel*t*gam_intermed)/math.sqrt(4*D_intermed*t))));
            }
          }
        c[i]    = math.mean(intlist)
        cmin[i] = math.min(intlist)
        cmax[i] = math.max(intlist)
    }
  }
  return [c, cmin, cmax]
}

// console.log() is accessable through F12


function getc_BTC(xBTC,vel,tsp,gam,t_inj,D) {
  const c = []
  var tau = []
  for (let i = 0; i < tsp.length; i++) {
      if (rg==1 && tsp[i] > t_inj){
        // Pulse injection |this produces negative values in the begining
        if (tsp[i]<t_inj) {
          tau = 0
        } else {
          tau = t_inj
        }
        c[i] = 1/2 * ((1-math.erf((xBTC-vel*tsp[i]) / math.sqrt(4*D*tsp[i])))-(1-math.erf((xBTC-vel*(tsp[i]-tau))/math.sqrt(4*D*(tsp[i]-tau)))));
      } else {
        // Continuous injection
        c[i] = 1/2 * Math.exp(xBTC*vel/(2*D))*(Math.exp((-xBTC)*vel*gam/(2*D))*(1-math.erf((xBTC-vel*tsp[i]*gam)/math.sqrt(4*D*tsp[i])))+math.exp(xBTC*vel*gam/(2*D))*(1-math.erf((xBTC+vel*tsp[i]*gam)/math.sqrt(4*D*tsp[i]))));
        }
  }
  return c
}

function get_gamma(reac,Dis,sep_vel) {
  var res = []
  res = Math.sqrt(1 + 4 * reac * Dis / sep_vel**2)
  return res
}

function transport_num(c,disp,dx,dt){
  // This function approximates advection and dispersion numerically by a FVM

  // Move concentration by 1 cell
  for (let i = 0; i < c.length-1; i++) {
    c[i+1] = c[i]
  }
  // First cell gets inlet concentration
  c[1] = 1
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
}

// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var ymin= source1.data['ymin']
var ymax= source1.data['ymax']
var x2  = source2.data['x2']
var y2  = source2.data['y2']

// Extracting slider values
const col_len   = col_len_sl.value;                 // [m]
const xBTC      = xBTC_sl.value;                    // [m]
const rad       = col_rad_sl.value;                 // [m]
const reac_l    = Math.exp(reac_sl.value[0])/3600;  // [1/s]
const reac_h    = Math.exp(reac_sl.value[1])/3600;  // [1/s]
const disp_l    = Math.exp(disp_sl.value[0])/3600;  // [m2/s]
const disp_h    = Math.exp(disp_sl.value[1])/3600;  // [m2/s]
const Q         = flow_sl.value/1000/1000/3600;     // [m3/s]
const n         = poros_sl.value;                   // [-]
const tPV       = Math.exp(pore_vol_sl.value);      // [-]
const t_inj     = pulse_inj_sl.value                // [s]
const rg        = rbgr.active                       // [false]

// Initializing empty lists
var c = []
var cmin = []
var cmax = []
var cBTX = []
var tsp = []

const c0      = 1;                            // [-]
const reac    = (reac_l + reac_h)/2           // [1/s]
const Dis     = (disp_l + disp_h)/2           // [m2/s]
const A       = math.PI * rad**2;             // [m2]
const vel     = Q/A;                          // [m/s]
const sep_vel = vel * n                       // [m/s]   
const PS      = col_len * A * n               // [m3]
const PV      = PS / (A*vel)                  // [s] VEL oder SEP_VEL?
const t       = tPV * PV                      // [s]

const gam     = Math.sqrt(1 + 4 * reac * Dis / sep_vel**2)  

for (let j = 0; j < x.length; j++) {
  x[j] = -0.02*col_len + 1.02*col_len/x.length * j;
}
for (let j = 0; j < x2.length; j++) {
  tsp[j] = x2[j] * PV;
}

console.log(tsp)

// This if statement has no meaning besides preventing a Type Error <-- why is that? It doesnt work without it
if (1<2){ 
  [c, cmin, cmax] = getc(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg)
  cBTX = getc_BTC(xBTC,sep_vel,tsp,gam,t_inj,Dis) 
}

    
// We have to loop through all indicies
for (let i = 0; i < c.length; i++) {
  y[i] = c[i]
  ymin[i] = cmin[i]
  ymax[i] = cmax[i]
}

for (let i = 0; i < x.length; i++) {
  y2[i] = cBTX[i]
}


// Update Sliders
pore_vol_sl.title = 'Pore Volume (1PV =' + (PV/3600).toFixed(2) +'h)';
xBTC_sl.end       = col_len_sl.value;
BTClocation.location = xBTC

if (rg==0) {
  pulse_inj_sl.visible = false
} else {
  pulse_inj_sl.visible = true
}

console.log(x2,y2)
console.log(x,y)

source1.change.emit();
source2.change.emit();