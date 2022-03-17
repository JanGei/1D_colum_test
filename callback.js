function getc_cont(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h) {
  var c = []
  var cmin = []
  var cmax = []
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

            intlist[j] = 1/2 * Math.exp(x[i]*vel/(2*D_intermed))*(Math.exp((-x[i])*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]-vel*t*gam_intermed)/math.sqrt(4*D_intermed*t)))+math.exp(x[i]*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]+vel*t*gam_intermed)/math.sqrt(4*D_intermed*t))));
          }
        c[i]    = math.mean(intlist)
        cmin[i] = math.min(intlist)
        cmax[i] = math.max(intlist)
    }
  }
  return [c, cmin, cmax]
}

function getc_pulse(x,vel,t,w,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h) {
  var c = []
  var cmin = []
  var cmax = []
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

            intlist[j] = 1/2 * (math.erf((x[i] + w/2 - vel*t)/(math.sqrt(4*D_intermed*t))) - math.erf((x[i]- w/2 - vel*t)/(math.sqrt(4*D_intermed*t))))*math.exp(-r_intermed*x[i]/vel);
          }
        c[i]    = math.mean(intlist)
        cmin[i] = math.min(intlist)
        cmax[i] = math.max(intlist)
    }
  }
  return [c, cmin, cmax]
}

// console.log() is accessable through F12


function getc_cont_BTC(xBTC,vel,tsp,gam,D) {
  const c = []
  for (let i = 0; i < tsp.length; i++) {
      c[i] = 1/2 * Math.exp(xBTC*vel/(2*D))*(Math.exp((-xBTC)*vel*gam/(2*D))*(1-math.erf((xBTC-vel*tsp[i]*gam)/math.sqrt(4*D*tsp[i])))+math.exp(xBTC*vel*gam/(2*D))*(1-math.erf((xBTC+vel*tsp[i]*gam)/math.sqrt(4*D*tsp[i]))));
  }
  return c
}

function getc_pulse_BTC(xBTC,vel,tsp,reac,w,D) {
  const c = []
  for (let i = 0; i < tsp.length; i++) {
      c[i] = 1/2 * (math.erf((xBTC + w/2 - vel*tsp[i])/(math.sqrt(4*D*tsp[i]))) - math.erf((xBTC- w/2 - vel*tsp[i])/(math.sqrt(4*D*tsp[i]))))*math.exp(-reac*xBTC/vel);
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


// Extracting slider values
const col_len   = col_len_sl.value;
const xBTC      = xBTC_sl.value;
const rad       = col_rad_sl.value;
const reac_l    = Math.exp(reac_sl.value[0])/3600;
const reac_h    = Math.exp(reac_sl.value[1])/3600;
const disp_l    = Math.exp(disp_sl.value[0])/3600;
const disp_h    = Math.exp(disp_sl.value[1])/3600;
const Q         = flow_sl.value;
const n         = poros_sl.value;
const tPV       = Math.exp(pore_vol_sl.value);
const t_w       = pulse_inj_sl.value
const rg        = rbgr.active

// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var ymin= source1.data['ymin']
var ymax= source1.data['ymax']
var x2  = source2.data['x2']
var y2  = source2.data['y2']

// Concentration lists for different parameter combinations
var c_ll= []
var c_lh= []
var c_hl= []
var c_hh= []
var c = []
var cmin = []
var cmax = []
var cBTX = []

// List for timespan (needed for BTC)
const tsp = []
const c0  = 1;

const reac= (reac_l + reac_h)/2
const Dis = (disp_l + disp_h)/2

const A       = math.PI * rad**2;            
const vel     = Q/1000/1000/3600/A;
const sep_vel = vel * n  
const PS      = col_len * A * n 
const PV      = PS / (A*vel) //normal velocity is needed, not seepage velocity tp obtain flux 
const t       = tPV * PV
const w       = sep_vel * t_w

const gam     = Math.sqrt(1 + 4 * reac * Dis / sep_vel**2)  
const gam_ll  = Math.sqrt(1 + 4 * reac_l * disp_l / sep_vel**2)  
const gam_lh  = Math.sqrt(1 + 4 * reac_l * disp_h / sep_vel**2)  
const gam_hl  = Math.sqrt(1 + 4 * reac_h * disp_l / sep_vel**2)  
const gam_hh  = Math.sqrt(1 + 4 * reac_h * disp_h / sep_vel**2)  


for (let j = 0; j < x.length; j++) {
  x[j] = -0.02*col_len + 1.02*col_len/x.length * j;
}
for (let j = 0; j < x2.length; j++) {
  tsp[j] = x2[j] * PV;
}


if (rg==0) {
  // Continuous Injection
  [c, cmin, cmax] = getc_cont(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h)
  cBTX = getc_cont_BTC(xBTC,sep_vel,tsp,gam,Dis)
  pulse_inj_sl.visible = false 
} else {
  // Pulse Injection
  [c, cmin, cmax] = getc_pulse(x,vel,t,w,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h)
  cBTX = getc_pulse_BTC(xBTC,sep_vel,tsp,reac,w,Dis)
  pulse_inj_sl.visible = true
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

source1.change.emit();
source2.change.emit();