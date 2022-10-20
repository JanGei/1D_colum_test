// Function: Compute concentration in column
function getc(x,vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg_CP,Ret,r_mean,D_mean,H_mean) {
  // Initializing empty concentration lists
  var c = []
  var cmin = []
  var cmax = []
  var cloQ = []
  var cupQ = []

  // Seepage velocity
  var vel_r = vel/Ret

  // Variables needed, if Ogata-Banks (1962) solution is used
  //var r_mean = (reac_l+reac_h)/2/Ret
  //var D_mean = (disp_l+disp_h)/2/Ret
  //var H_mean = 2*r_mean*D_mean/vel_r**2
  //var gam_mean = get_gamma(r_mean,D_mean,vel)

  for (let i = 0; i < x.length; i++) { 
      // "Behind" the column: concentration is either 1 or 0, depending on injction time
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

      // "In" the column
      } else {
        // Intermediate list to store all solutions of the latin hypercube
        var intlist = []
          for (let j = 0; j < Lcube1.length; j++) {
            // Computing values associated to the latin hypercube
            var r_intermed = (reac_l + (reac_h-reac_l)*Lcube1[j])/Ret
            var D_intermed = (disp_l + (disp_h-disp_l)*Lcube2[j])/Ret
            var H_intermed = 2*r_intermed*D_intermed/vel_r**2

            // Variable needed, if Ogata-Banks (1962) solution is used
            // var gam_intermed = get_gamma(r_intermed,D_intermed,vel)

            // Checking for pulse injection
            if (rg_CP == 1) {
              // Pulse injection
              if (t<=t_inj) {
                // Runkler (1996, Eq. 8, O'Loughlin and Bowmer)
                intlist[j] = 1/2 * (math.exp(-r_intermed*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_intermed))/(2*math.sqrt(D_intermed*t)))))
                // Ogata-Banks (1962) 
                // intlist[j] = 1/2 * (1-math.erf((x[i]-vel*t) / math.sqrt(4*D_intermed*t)))
              } else {
                // Runkler (1996, Eq. 10, O'Loughlin and Bowmer)
                intlist[j] = 1/2 * math.exp(-r_intermed*x[i]/vel_r) * ( (1-math.erf((x[i]-vel_r*t*(1+H_intermed))/(2*math.sqrt(D_intermed*t)))) - (1-math.erf((x[i]-vel_r*(t-t_inj)*(1+H_intermed))/(2*math.sqrt(D_intermed*(t-t_inj))))))
                // Ogata-Banks (1962) 
                //intlist[j] = 1/2 * ((1-math.erf((x[i]-vel*t) / math.sqrt(4*D_intermed*t)))-(1-math.erf((x[i]-vel*(t-t_inj))/math.sqrt(4*D_intermed*(t-t_inj)))));
              }
            } else {
              // Cotinuous injection
              // Runkler (1996, Eq. 8, O'Loughlin and Bowmer)
              intlist[j] = 1/2 * (math.exp(-r_intermed*x[i]/vel_r) * (1-math.erf((x[i] - vel_r*t*(1+ H_intermed))/(2*math.sqrt(D_intermed*t)))))
              // Ogata-Banks (1962) 
              //intlist[j] = 1/2 * Math.exp(x[i]*vel/(2*D_intermed))*(Math.exp((-x[i])*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]-vel*t*gam_intermed)/math.sqrt(4*D_intermed*t)))+math.exp(x[i]*vel*gam_intermed/(2*D_intermed))*(1-math.erf((x[i]+vel*t*gam_intermed)/math.sqrt(4*D_intermed*t))));
            }
          }
        
        // Obtain median, min, max, upper and lower quartile from intermediate list
        c[i]    = math.median(intlist)
        cmin[i] = math.min(intlist)
        cmax[i] = math.max(intlist)
        cloQ[i] = math.quantileSeq(intlist, 0.25)
        cupQ[i] = math.quantileSeq(intlist, 0.75)
    }
  }
  return [c, cmin, cmax, cloQ, cupQ]
}

// Function: Breakthrough Curve (BTC)
function getc_BTC(xBTC,vel,tsp,t_inj,D_mean,r_mean,H_mean,Ret) {
  // Initializing empty concentration list
  const c = []

  // Seepage velocity
  var vel_r = vel/Ret

  // Variables needed, if Ogata-Banks (1962) solution is used
  //const gam     = Math.sqrt(1 + 4 * r_mean * D_mean / sep_vel**2) 

  for (let i = 0; i < tsp.length; i++) {

      // Pulse injection 
      if (rg_CP==1){

        if (tsp[i]<=t_inj) {
          // Runkler (1996, Eq. 8, O'Loughlin and Bowmer)
          c[i] = 1/2 * ( math.exp(-r_mean*xBTC/vel_r) * (1-math.erf((xBTC - vel_r*tsp[i]*(1+ H_mean))/(2*math.sqrt(D_mean*tsp[i])))))
          // Ogata-Banks (1962) 
          //c[i] = 1/2 * (1-math.erf((xBTC-vel_r*tsp[i]) / math.sqrt(4*D*tsp[i])))
        } else {
          // Runkler (1996, Eq. 10, O'Loughlin and Bowmer)
          c[i] = 1/2 * math.exp(-r_mean*xBTC/vel_r) * ( (1-math.erf((xBTC-vel_r*tsp[i]*(1+H_mean))/(2*math.sqrt(D_mean*tsp[i])))) - (1-math.erf((xBTC-vel_r*(tsp[i]-t_inj)*(1+H_mean))/(2*math.sqrt(D_mean*(tsp[i]-t_inj))))) )
          // Ogata-Banks (1962) 
          //c[i] = 1/2 * ((1-math.erf((xBTC-vel_r*tsp[i]) / math.sqrt(4*D*tsp[i])))-(1-math.erf((xBTC-vel_r*(tsp[i]-t_inj))/math.sqrt(4*D*(tsp[i]-t_inj)))));
        }

      // Continuous injection
      } else {
        // Runkler (1996, Eq. 8, O'Loughlin and Bowmer)
        c[i] = 1/2 * ( math.exp(-r_mean*xBTC/vel_r) * (1-math.erf((xBTC - vel_r*tsp[i]*(1+ H_mean))/(2*math.sqrt(D_mean*tsp[i])))))
        // Ogata-Banks (1962) 
        //c[i] = 1/2 * Math.exp(xBTC*vel_r/(2*D_mean))*(Math.exp((-xBTC)*vel_r*gam/(2*D_mean))*(1-math.erf((xBTC-vel_r*tsp[i]*gam)/math.sqrt(4*D_mean*tsp[i])))+math.exp(xBTC*vel_r*gam/(2*D_mean))*(1-math.erf((xBTC+vel_r*tsp[i]*gam)/math.sqrt(4*D_mean*tsp[i]))));
        }
  }
  return c
}

// Function: Obtain gamma for Ogata-Banks solution
function get_gamma(r_mean,D_mean,sep_vel) {
  var res = []
  res = Math.sqrt(1 + 4 * r_mean * D_mean / sep_vel**2)
  return res
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


// Extracting values from callback
// For variable identification, refer to file 1D_column_test.py
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

// Derived values
const A       = math.PI * rad**2;             // [m2]   Area
const vel     = Q/A;                          // [m/s]  Velocity
const sep_vel = vel / n                       // [m/s]  Seepage velocity
const r_mean  = (reac_l+reac_h)/2             // [1/s]  1st order reaction constant
const D_mean  = (disp_l+disp_h)/2             // [m2/s] Dispersion coefficient
const H_mean  = 2*r_mean*D_mean/sep_vel**2    // [m/s]  Variable for Ogata-Banks
const r_med   = reac_l+(reac_h-reac_l)*math.median(Lcube1)   // [1/s]  Median 1st order reaction constant
const D_med   = disp_l+(disp_h-disp_l)*math.median(Lcube2)   // [m2/s] Median Dispersion coefficient
const H_med   = 2*r_med*D_med/sep_vel**2      // [m/s]  Median Variable for Ogata-Banks
const PS      = col_len * A * n               // [m3]   Pore space
const PV      = col_len/sep_vel               // [s]    Time equivalent for one pore volume
const c0      = 1;                            // [-]    Normed inlet concentration

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
  if (j == 0) {
    x[j] = - 0.005 * col_len 
  } else {
    x[j] = 0.005 * col_len * j;
}}

for (let j = 0; j < x2.length; j++) {
  tsp[j] = x2[j] * PV;
}

// Fix point draw tool to x-axis and limit range on x-axis
y3[0] = 0
if (x3[0]<=0.001) {
  x3[0] = 0.01
} else if (x3[0] > col_len) {
  x3[0] = col_len
}

// Compute time equivalent of selected pore volume 
const tPV       = Math.exp(pore_vol_sl.value);      // [-]
const t         = tPV * PV                          // [s]

// Initializing empty lists
var c = []
var cmin = []
var cmax = []
var cloQ = []
var cupQ = []
var cBTX = []

// This if statement has no meaning besides preventing a Type Error
if (1<2){ 
  [c, cmin, cmax, cloQ, cupQ] = getc(x,sep_vel,t,Lcube1,Lcube2,reac_l,reac_h,disp_l,disp_h,t_inj,rg_CP,Ret,r_mean,D_mean,H_mean)
  cBTX = getc_BTC(xBTC,sep_vel,tsp,t_inj,D_mean,r_mean,H_mean,Ret) 
}


// Updating sources
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

// Updating Sliders
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