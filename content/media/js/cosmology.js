/***************************************************
 * General Utilities
 **************************************************/


function print(s) {
    document.write(s + "<br/>");
}

function log10(x) {
    return Math.log(x) / Math.LN10;
}

function sinh(x) {
    return 0.5*(Math.exp(x) - Math.exp(-x));
}

// Emulate Numpy's useful linspace() function
function linspace(start, end, npts) {
    var vals = [];
    var diff = (end - start) / (npts - 1);
    for (var i = 0; i < npts; i++)
	vals.push(start + i*diff);
    return vals;
}

// Emulate Python's useful map() function
function map(func, vals) {
    var result = [];
    for (var i = 0; i < vals.length; i++)
	result.push(func(vals[i]));
    return result;
}

// Emulate Python's useful zip function
function zip(arrs) {
    var nrows, ncols, i, j;
    var row, result = [];
    nrows = arrs[0].length;
    ncols = arrs.length;
    for (i = 0; i < nrows; i++) {
	row = [];
	for (j = 0; j < ncols; j++)
	    row.push(arrs[j][i]);
	result.push(row);
    }
    return result;
}


// integrate between function values y(x) and x axis using trapezoidal rule.
// Must provide x and y values at each point.
function trapz(x, y) {
    var sum = 0.0;
    for (var i = 1; i < y.length; i++)
	sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);
    return sum * 0.5;
}


// set all the parameters when the page loads
function init_params() {
    // speed of light in km/s (exact)
    Ckms = 299792.458;
    seconds_per_Gyr = 3600 * 24 * 365.25 * 1e9;
    // 1 parsec in m (from wikipedia, which cites P. Kenneth Seidelmann,
    // Ed. (1992). Explanatory Supplement to the Astronomical
    // Almanac. Sausalito, CA: University Science Books. p. 716 and
    // s.v. parsec in Glossary.)
    PC = 3.08567782e16;

    // rough WMAP7 + BAO + H0 parameters
    H0 = 70;   // present day Hubble constant,km/s/Mpc
    Om = 0.27;   // matter density
    Ol = 1 - Om;   // dark energy density

    W = -1.;    // pressure / density for dark energy

    zmin = 0;
    zmax = 5;

    Ok = 0;                       // curvature energy density
    Dh = Ckms / H0;               // Hubble distance in Mpc

    n_integration_steps = 1000;

    Dcmin = Dc(zmin);
    Dcmax = Dc(zmax);

    Th = 1. / H0 * PC * 1.e3 / seconds_per_Gyr;   // Hubble time in Gyr
}

// display the current parameters in the form
function print_params() {
    document.getElementById("printOl").innerHTML = Ol.toFixed(2);
    document.getElementById("printDc1").innerHTML = Dcmin.toFixed(0);
    document.getElementById("printDc2").innerHTML = Dcmax.toFixed(0);
    document.parform.H0.value = H0;
    document.parform.Om.value = Om;
    document.parform.w.value = W;
    document.parform.zmin.value = zmin;
    document.parform.zmax.value = zmax;
}

// update the stored parameters using the values in the form.
function update_params() {
    H0 = parseFloat(document.parform.H0.value);
    if (H0 < 0) H0 = 0.;
    if (isNaN(H0)) H0 = 70.;
    Om = parseFloat(document.parform.Om.value);
    if (Om < 0) Om = 0;
    if (Om > 1) Om = 1;
    if (isNaN(Om)) Om = 0.27;
    zmin = parseFloat(document.parform.zmin.value);
    zmax = parseFloat(document.parform.zmax.value);
    if (isNaN(zmin)) zmin = 0.;
    if (isNaN(zmax)) zmax = 5.;
    if (zmin > zmax) {
	var temp = zmin;
	zmin = zmax;
	zmax = temp;
    }
    if (zmin < 0)  zmin = 0;
    if (zmax > 100) zmax = 100;
    Ol = 1 - Om;
    W = parseFloat(document.parform.w.value);
    if (isNaN(W)) W = -1.;
    Dh = Ckms / H0;               // Hubble distance in Mpc
    Th = 1. / H0 * PC * 1.e3 / seconds_per_Gyr;   // Hubble time in Gyr

    Dcmin = Dc(zmin);
    Dcmax = Dc(zmax);
}

function doplots() {
    update_params();
    print_params();
    var z, comdist;
    var zvals = linspace(zmin, zmax, 50);

    var p = zip([zvals, map(Tl, zvals)]);

    // Note the dollar sign here refers to the jquery object; plot is
    // a method added to the jquery object by flot.

    // jQuery("#plot_Tl") calls the jquery object, which returns the document
    // element with id "plot_Tl" (see the jquery documentation for more
    // details).
    jQuery.plot(jQuery("#plot_Tl"), [{data:p, shadowSize: 1}]);

    var p1 = zip([zvals, map(function (z) {return Dc(z)*Math.PI/180/60;},zvals)]);
    var p2 = [];
    for (var i = 0; i < p1.length; i++){
	z = zvals[i];
	comdist = p1[i][1];
	p2.push([z, comdist / (1 + z)]);
    }
    jQuery.plot(jQuery("#plot_dist"),[{label:"Comoving", data:p1, shadowSize: 1},
			    {label:"Proper", data:p2, shadowSize: 1}],
			    {legend: {position:"nw",backgroundOpacity:0}}
	  );

    p = zip([zvals, map(Hz, zvals)]);
    jQuery.plot(jQuery("#plot_H"), [{data:p, shadowSize: 1}]);

    p = zip([zvals, map(function (z) {return distmod(z)-21-2.5*log10(1+z);}, zvals)]);
    // dodgy heuristic
    if (zmin < 1e-8)
	p = p.slice(1);
    jQuery.plot(jQuery("#plot_Lstar"), [{data:p, shadowSize: 1}]);
}

// This is used to intercept a keypress event in the form fields,
// and plot if it's enter, otherwise ignore it.
function intercept_enter(event) {
  if (event && event.keyCode == 13) {
      doplots();
      return false;
  }
  else
      return true;
}


// Eqn 14 from Hogg, adding 'w' to Omega lambda.
function efunc(z) {
    var zp1 = 1. + z;
    return Math.sqrt(Om*Math.pow(zp1, 3) + Ok*zp1*zp1 +
		     Ol*Math.pow(zp1, 3. * (1. + W)));
}

function inv_efunc(z) {
    return 1. / efunc(z);
}

// Function for integration to find lookback time. Eqn 30 from Hogg.
function tfunc(z) {
    var zp1 = 1. + z;
    return 1. / (zp1 * efunc(z));
}

// Function for integration to find the absorption distance.
function xfunc(z) {
    var zp1 = 1. + z;
    return zp1*zp1 / efunc(z);
}

// Lookback time in Gyr
function Tl(z) {
    var zvals = linspace(0, z, n_integration_steps);
    return Th * trapz(zvals, map(tfunc, zvals));
}

// Comoving distance in Mpc
function Dc(z) {
    var zvals = linspace(0, z, n_integration_steps);
    return Dh * trapz(zvals, map(inv_efunc, zvals));
}

function Hz(z) {
    return H0 * efunc(z);
}

// Returns the transverse comoving distance in Mpc.  Dm*dtheta
// gives the transverse comoving distance between two objects
// separated by an angle dtheta radians at redshift z. Note Dm = Dc if
// Ok is zero (as in current lambda CDM model).
function Dm(z) {
    var comdist = Dc(z);
    if (Ok == 0)
        return comdist;
    var sqrtOk = Math.sqrt(Math.abs(Ok));
    if (Ok > 0)
        return comdist / sqrtOk * Math.sinh(sqrtOk * comdist / Dh);
    else
        return comdist / sqrtOk * Math.sin(sqrtOk * comdist / Dh);
}

// Returns the luminosity distance in Mpc.  (Relationship
// between bolometric flux and bolometric luminosity.)
function  Dl(z){
    return (1. + z) * Dm(z);
}

// Returns the distance modulus (apparent magnitude - absolute
// magnitude for an object at redshift z).
function distmod(z) {
	return 5. * log10(self.Dl(z) * 1e5);
}

// Return the absorption distance from redshift
// z2 to z1. Dimensionless (despite being called a distance)
function absorpdist(z) {
    var zvals =	linspace(0, z, n_integration_steps);
    return trapz(zvals, map(xfunc, zvals));
}
