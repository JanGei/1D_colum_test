function getc(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg_CP,Ret) {
  var c = []
  var cmin = []
  var cmax = []
  var cloQ = []
  var cupQ = []
  var vel_r = vel/Ret
  //var r_mean = (reac_l+reac_h)/2/Ret
  //var D_mean = (disp_l+disp_h)/2/Ret
  //var H_mean = 2*r_mean*D_mean/vel_r**2
  //var gam_mean = get_gamma(r_mean,D_mean,vel)
  for (let i = 0; i < x.length; i++) { 
      if (x[i] <= 0) {
        if (rg_CP == 1 && t > t_inj) {
          c[i] = 0
          cmin[i] = 0
          cmax[i] = 0
        } else {
          c[i] = 1
          cmin[i] = 1
          cmax[i] = 1
        }
        
      } else {
        var intlist = []
          for (let j = 0; j < Lcube1.length; j++) {
            var r_intermed = (reac_l + (reac_h-reac_l)*Lcube1[j])/Ret
            var D_intermed = (disp_l + (disp_h-disp_l)*Lcube2[j])/Ret
            var H_intermed = 2*r_intermed*D_intermed/vel_r**2
            //var gam_intermed = get_gamma(r_intermed,D_intermed,vel)
            if (rg_CP == 1) {
              // Pulse injection
              if (t<=t_inj) {
                //intlist[j] = 1/2 * (1-math.erf((x[i]-vel*t) / math.sqrt(4*D_intermed*t)))
                // eq 8 in Runkler 1996 (O'Loughlin and Bowmer)
                intlist[j] = 1/2 * ( math.exp(-r_intermed*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_intermed))/(2*math.sqrt(D_intermed*t)))) )
              } else {
                //intlist[j] = 1/2 * ((1-math.erf((x[i]-vel*t) / math.sqrt(4*D_intermed*t)))-(1-math.erf((x[i]-vel*(t-t_inj))/math.sqrt(4*D_intermed*(t-t_inj)))));
                // eq 10 in Runkler 1996 (O'Loughlin and Bowmer)
                intlist[j] = 1/2 * math.exp(-r_intermed*x[i]/vel_r) * ( (1-math.erf((x[i]-vel_r*t*(1+H_intermed))/(2*math.sqrt(D_intermed*t)))) - (1-math.erf((x[i]-vel_r*(t-t_inj)*(1+H_intermed))/(2*math.sqrt(D_intermed*(t-t_inj))))) )
              }
            } else {
              // Cotinuous injection
              //intlist[j] = 1/2 * Math.exp(x[i]*vel/(2*D_intermed))*(Math.exp((-x[i])*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]-vel*t*gam_intermed)/math.sqrt(4*D_intermed*t)))+math.exp(x[i]*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]+vel*t*gam_intermed)/math.sqrt(4*D_intermed*t))));
              // eq 8 in Runkler 1996 (O'Loughlin and Bowmer)
              intlist[j] = 1/2 * ( math.exp(-r_intermed*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_intermed))/(2*math.sqrt(D_intermed*t)))) )
            }
          }
        // Main line with mean values of dispersion and reaction
        //if (rg_CP==1){
          // Pulse injection
          //if (t<=t_inj) {
            //c[i] = 1/2 * (1-math.erf((x[i]-vel*t) / math.sqrt(4*D_mean*t)))
            // eq 8 in Runkler 1996 (O'Loughlin and Bowmer)
            //c[i] = 1/2 * ( math.exp(-r_mean*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_mean))/(2*math.sqrt(D_mean*t)))) )
          //} else {
            //c[i] = 1/2 * ((1-math.erf((x[i]-vel*t) / math.sqrt(4*D_mean*t)))-(1-math.erf((x[i]-vel*(t-t_inj))/math.sqrt(4*D_mean*(t-t_inj)))))
            // eq 10 in Runkler 1996 (O'Loughlin and Bowmer)
            //c[i] = 1/2 * math.exp(-r_mean*x[i]/vel_r) * ( (1-math.erf((x[i]-vel_r*t*(1+H_mean))/(2*math.sqrt(D_mean*t)))) - (1-math.erf((x[i]-vel_r*(t-t_inj)*(1+H_mean))/(2*math.sqrt(D_mean*(t-t_inj))))) )
          //}
        //} else {
          // Continuous injection
          // eq 8.66 in hydrogeology script (Ogata Banks)
          //c[i] = 1/2 * Math.exp(x[i]*vel/(2*D_mean))*(Math.exp((-x[i])*vel*gam_mean/(2*D_mean))*(1-math.erf((x[i]-vel*t*gam_mean)/math.sqrt(4*D_mean*t)))+math.exp(x[i]*vel*gam_mean/(2*D_mean))*(1-math.erf((x[i]+vel*t*gam_mean)/math.sqrt(4*D_mean*t))));
          // eq 8 in Runkler 1996 (O'Loughlin and Bowmer)
          //c[i] = 1/2 * ( math.exp(-r_mean*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_mean))/(2*math.sqrt(D_mean*t)))) )
        //}
        c[i]    = math.median(intlist)
        cmin[i] = math.min(intlist)
        cmax[i] = math.max(intlist)
        cloQ[i] = math.quantileSeq(intlist, 0.25)
        cupQ[i] = math.quantileSeq(intlist, 0.75)
    }
  }
  return [c, cmin, cmax, cloQ, cupQ]
}

// console.log() is accessable through F12


function getc_BTC(xBTC,vel,tsp,gam,t_inj,D_mean,r_mean,H_mean,Ret) {
  const c = []
  var vel_r = vel/Ret
  for (let i = 0; i < tsp.length; i++) {
      if (rg_CP==1){
        // Pulse injection 
        if (tsp[i]<=t_inj) {
          //c[i] = 1/2 * (1-math.erf((xBTC-vel_r*tsp[i]) / math.sqrt(4*D*tsp[i])))
          c[i] = 1/2 * ( math.exp(-r_mean*xBTC/vel_r) * (1-math.erf((xBTC - vel_r*tsp[i]*(1+ H_mean))/(2*math.sqrt(D_mean*tsp[i])))) )
        } else {
          //c[i] = 1/2 * ((1-math.erf((xBTC-vel_r*tsp[i]) / math.sqrt(4*D*tsp[i])))-(1-math.erf((xBTC-vel_r*(tsp[i]-t_inj))/math.sqrt(4*D*(tsp[i]-t_inj)))));
          c[i] = 1/2 * math.exp(-r_mean*xBTC/vel_r) * ( (1-math.erf((xBTC-vel_r*tsp[i]*(1+H_mean))/(2*math.sqrt(D_mean*tsp[i])))) - (1-math.erf((xBTC-vel_r*(tsp[i]-t_inj)*(1+H_mean))/(2*math.sqrt(D_mean*(tsp[i]-t_inj))))) )
        }
      } else {
        // Continuous injection
        // eq 8.66 in hydrogeology script (Ogata Banks)
        //c[i] = 1/2 * Math.exp(xBTC*vel_r/(2*D_mean))*(Math.exp((-xBTC)*vel_r*gam/(2*D_mean))*(1-math.erf((xBTC-vel_r*tsp[i]*gam)/math.sqrt(4*D_mean*tsp[i])))+math.exp(xBTC*vel_r*gam/(2*D_mean))*(1-math.erf((xBTC+vel_r*tsp[i]*gam)/math.sqrt(4*D_mean*tsp[i]))));
        // eq 8 in Runkler 1996 (O'Loughlin and Bowmer)
        c[i] = 1/2 * ( math.exp(-r_mean*xBTC/vel_r) * (1-math.erf((xBTC - vel_r*tsp[i]*(1+ H_mean))/(2*math.sqrt(D_mean*tsp[i])))) )
        }
  }
  return c
}

function get_gamma(r_mean,D_mean,sep_vel) {
  var res = []
  res = Math.sqrt(1 + 4 * r_mean * D_mean / sep_vel**2)
  return res
}

function expo(x, f) {
  return Number.parseFloat(x).toExponential(f);
  //https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Global_Objects/Number/toExponential
}

// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var ymin= source1.data['ymin']
var ymax= source1.data['ymax']
var yloQ= source1.data['yloQ']
var yupQ= source1.data['yupQ']
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
const Kd        = Kd_sl.value                       // [m3/kg]
const rho_s     = rho_s_sl.value                    // [kg/m3]

// Derived entities
const A       = math.PI * rad**2;             // [m2]
const vel     = Q/A;                          // [m/s]
const sep_vel = vel / n                       // [m/s]
const r_mean  = (reac_l+reac_h)/2             // [1/s]
const D_mean  = (disp_l+disp_h)/2             // [m2/s]
const H_mean  = 2*r_mean*D_mean/sep_vel**2    // [m/s]
const PS      = col_len * A * n               // [m3]
const PV      = col_len/sep_vel               // [s] VEL oder SEP_VEL?
const c0      = 1;                            // [-] 

// Retardation factor
if (rg_SType == 0) {
  var Ret = 1
} else {
  var Ret = 1 + (1-n)/n * Kd * rho_s
}

// Time span list
var tsp = []

// Discretize space (upper plot) and time (lower plot)
for (let j = 0; j < x.length; j++) {
  x[j] = -0.02*col_len + 1.02*col_len/x.length * j;
}
for (let j = 0; j < x2.length; j++) {
  tsp[j] = x2[j] * PV;
}

// Fix point draw tool to x-axis and limit its range on x-axis
y3[0] = 0
if (x3[0]<=0.001) {
  x3[0] = 0.01
} else if (x3[0] > col_len) {
  x3[0] = col_len
}

// Time 
const tPV       = Math.exp(pore_vol_sl.value);      // [-]
const t         = tPV * PV                          // [s]

const gam     = Math.sqrt(1 + 4 * r_mean * D_mean / sep_vel**2) 

// Initializing empty lists
var c = []
var cmin = []
var cmax = []
var cloQ = []
var cupQ = []
var cBTX = []

// This if statement has no meaning besides preventing a Type Error <-- why is that? It doesnt work without it
if (1<2){ 
  [c, cmin, cmax, cloQ, cupQ] = getc(x,sep_vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg_CP,Ret)
  cBTX = getc_BTC(xBTC,sep_vel,tsp,gam,t_inj,D_mean,r_mean,H_mean,Ret) 
}


// Update sources
for (let i = 0; i < c.length; i++) {
  y[i] = c[i]
  ymin[i] = cmin[i]
  ymax[i] = cmax[i]
  yloQ[i] = cloQ[i]
  yupQ[i] = cupQ[i]
}
for (let i = 0; i < x.length; i++) {
  y2[i] = cBTX[i]
}

// Update Sliders
pore_vol_sl.title = 'Pore Volume (1PV =' + (PV/3600).toFixed(2) +'h)';
BTCp.title.text   = 'Breakthrough Curve at x = ' + xBTC.toFixed(3) + ' m (Drag diamond in upper plot to change)'

// Displaying correct sliders
if (rg_CP==0) {
  pulse_inj_sl.visible = false
} else {
  pulse_inj_sl.visible = true
}
if (rg_SType == 1) {
  rho_s_sl.visible = true
  Kd_sl.visible = true
} else {
  rho_s_sl.visible = false
  Kd_sl.visible = false
}

// Displaying correct units
var r_format  = r_dict[r_us.value]
var D_format  = D_dict[D_us.value]
var fl_format = fl_dict[fl_us.value]
reac_sl.format = r_format
disp_sl.format = D_format
flow_sl.format = fl_format

reac_sl.format.change.emit();
disp_sl.format.change.emit();
flow_sl.format.change.emit();
source1.change.emit();
source2.change.emit();