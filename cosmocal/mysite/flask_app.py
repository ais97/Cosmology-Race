from flask import Response
from flask import Flask, render_template, request
import io
from astropy import constants as const
from scipy import integrate
from astropy import units as u
import numpy as np
import random
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from processing import age_of_universe, distance, temperature
app = Flask(__name__)
app.config["DEBUG"] = True

@app.route('/')
def student():
	return render_template('basic_layout.html')
	#return render_template('speech_recognition.html')
	#return render_template('index.html')

@app.route('/theory')
def theory():
	return render_template('theory.html')

@app.route('/memes')
def memes():
	return render_template('memes.html')


@app.route('/result',methods = ['POST', 'GET'])
def result():
	if request.method == 'POST':
		errors = ""
		H0 = 69.6
		z = 1.0
		omega_m0 = 0.286
		omega_lam0 = 0.713914
		omega_k0 = 0.0
		omega_rad0 = 8.4e-5
		try:
			H0 = float(request.form["H0"])
		except:
			errors += "{!r} input to Hubble Parameter is not a number.\n".format(request.form["H0"])
		try:
			z = float(request.form["z"])
			if(z<0):
				errors += "{!r} input to Redshift cannot be less than 0.\n".format(request.form["z"])
		except:
			errors += "{!r} input to Redshift is not a number.\n".format(request.form["z"])
		try:
			omega_k0 = float(request.form["omega_k0"])
		except:
			errors += "{!r} input to Curvature Parameter is not a number.\n".format(request.form["omega_k0"])
		try:
			omega_rad0 = float(request.form["omega_rad0"])
			if(omega_rad0<0):
				errors += "{!r} input to Radiation Density Parameter cannot be less than 0.\n".format(request.form["omega_rad0"])
		except:
			errors += "{!r} input to Radiation Density Parameter is not a number.\n".format(request.form["omega_rad0"])
		try:
			omega_m0 = float(request.form["omega_m0"])
			if(omega_m0<0):
				errors += "{!r} input to Matter Density Parameter cannot be less than 0.\n".format(request.form["omega_m0"])
		except:
			errors += "{!r} input to Matter Density Parameter is not a number.\n".format(request.form["omega_m0"])
		try:
			omega_lam0 = float(request.form["omega_lam0"])
		except:
			errors += "{!r} input to Dark Energy Density Parameter is not a number.\n".format(request.form["omega_lam0"])
		try:
			if(not np.allclose(omega_lam0+omega_m0+omega_rad0+omega_k0,1.)):
				errors += "Sum of density parameters not equal to 1. Hence, varying the Curvature parameter appropriately.\n"
				omega_k0 = 1 - (omega_lam0+omega_m0+omega_rad0)
		except:
			pass
		t,t0,t_lb = age_of_universe(z,H0,omega_m0, omega_rad0, omega_lam0, omega_k0)
		rc,dA,dL,dAs,Hl = distance(z,H0,omega_m0, omega_rad0, omega_lam0, omega_k0)
		T = temperature(z)
		#memes = ["age_of_universe.jpeg","avinash_amulya_universe.jpeg","cmb_fluctuations.jpeg","cosmology_expansion_of_universe.jpeg","cosmonaut.jpeg","curvature.png","dark_energy_not_constant.jpeg","different_distances_meme.jpeg","error_meme.jpeg","gravity_emergent_phenomenon.jpeg","have_you_ever_looked_up_at_the_sky.jpeg","homogeneous_isotropic.jpeg","ias_meme.jpeg","kbc.jpeg","multiverse.jpeg","peter_parker_meme.jpeg","random.jpeg","rohith_meme.png","so_what_can_i_do.jpeg","t&jcosmo_meme.gif"];memes = ["/static/"+meme for meme in memes];n = random.randint(0,len(memes));
		return render_template('result.html',
		errors=errors,
		H0 = H0,
		omega_m0=omega_m0,
		omega_rad0=omega_rad0,
		omega_lam0=omega_lam0,
		omega_k0=omega_k0,
		H0_str = str(H0),
		omega_m0_str=str(omega_m0),
		omega_rad0_str=str(omega_rad0),
		omega_lam0_str=str(omega_lam0),
		omega_k0_str=str(omega_k0),
		z=z,
		t0=t0,
		t=t,
		t_lb=t_lb,
		rc=rc,
		dA=dA,
		dL=dL,
		dAs=dAs,
		Hl=Hl,
		T=T,
		)

@app.route("/matplot-as-image-<string:H0_str>_<string:omega_m0_str>_<string:omega_rad0_str>_<string:omega_lam0_str>_<string:omega_k0_str>.png")
def plot_png(H0_str="69.6",omega_m0_str="0.286",omega_rad0_str="8.4e-5",omega_lam0_str="0.713914",omega_k0_str="0.0"):
	def plot_decoder(plot_code):
		new_tick_locations = np.array([t_arr[0]])
		new_ticks_labels = [z_str[0]]
		if(plot_code == "111"):
			if(np.log10(z_d)>0.75):
				new_tick_locations = np.append(new_tick_locations, t_d)
				if(z_d>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_d))
				else:
					new_ticks_labels.append("{}".format(np.round(z_d,3)))
			if(np.abs(np.log10(z_e)-np.log10(z_d))>1.0):
				new_tick_locations = np.append(new_tick_locations, t_e)
				if(z_e>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_e))
				else:
					new_ticks_labels.append("{}".format(np.round(z_e,3)))
			if(np.abs(np.log10(z_c)-np.log10(z_e))>1.0):
				new_tick_locations = np.append(new_tick_locations, t_c)
				if(z_c>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_c))
				else:
					new_ticks_labels.append("{}".format(np.round(z_c,3)))
		elif(plot_code=="011"):
			if(np.log10(z_d)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_d)
				if(z_d>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_d))
				else:
					new_ticks_labels.append("{}".format(np.round(z_d,3)))
			if(np.abs(np.log10(z_e)-np.log10(z_d))>0.5):
				new_tick_locations = np.append(new_tick_locations, t_e)
				if(z_e>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_e))
				else:
					new_ticks_labels.append("{}".format(np.round(z_e,3)))
		elif(plot_code=="110"):
			if(np.log10(z_e)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_e)
				if(z_e>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_e))
				else:
					new_ticks_labels.append("{}".format(np.round(z_e,3)))
			if(np.abs(np.log10(z_c)-np.log10(z_e))>0.5):
				new_tick_locations = np.append(new_tick_locations, t_c)
				if(z_c>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_c))
				else:
					new_ticks_labels.append("{}".format(np.round(z_c,3)))
		elif(plot_code=="101"):
			if(np.log10(z_d)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_d)
				if(z_d>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_d))
				else:
					new_ticks_labels.append("{}".format(np.round(z_d,3)))
			if(np.abs(np.log10(z_c)-np.log10(z_d))>0.5):
				new_tick_locations = np.append(new_tick_locations, t_c)
				if(z_c>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_c))
				else:
					new_ticks_labels.append("{}".format(np.round(z_c,3)))
		elif(plot_code=="001"):
			if(np.log10(z_d)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_d)
				if(z_d>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_d))
				else:
					new_ticks_labels.append("{}".format(np.round(z_d,3)))
		elif(plot_code=="010"):
			if(np.log10(z_e)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_e)
				if(z_e>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_e))
				else:
					new_ticks_labels.append("{}".format(np.round(z_e,3)))
		elif(plot_code=="100"):
			if(np.log10(z_c)>0.5):
				new_tick_locations = np.append(new_tick_locations, t_c)
				if(z_c>1000.0):
					new_ticks_labels.append("{:.2e}".format(z_c))
				else:
					new_ticks_labels.append("{}".format(np.round(z_c,3)))
		new_tick_locations = np.append(new_tick_locations, t_arr[-1])
		new_ticks_labels.append(z_str[-1])
		return new_tick_locations,new_ticks_labels

	H0 = float(H0_str)
	omega_m0 = float(omega_m0_str)
	omega_rad0 = float(omega_rad0_str)
	omega_lam0 = float(omega_lam0_str)
	omega_k0 = float(omega_k0_str)
	plt.style.use('dark_background')
	E = lambda z: omega_rad0*(1+z)**4 + omega_m0*(1+z)**3 + omega_lam0 + omega_k0*(1+z)**2

	def z_integrand(z):
		Z = -1.0/(H0*(1+z)*np.sqrt(E(z)))
		return Z

	def modified_z_integrand(z):
		Z = -s_to_Gyr/(H0*(1+z)*np.sqrt(E(z)))
		return Z

	pc_to_m = const.pc.value
	H0 = H0*(1e3/(1e6*pc_to_m))
	omega0 = {"rad": omega_rad0, "m": omega_m0, "lambda": omega_lam0}
	G = const.G.value
	rho_c0 = (3*(H0**2))/(8*np.pi*G)
	rho0 = dict(zip(list(omega0.keys()),[x*rho_c0 for x in list(omega0.values())]))
	s_to_yr = (1*u.second).to(u.yr).value
	s_to_Gyr = s_to_yr/1e9
	Gyr_to_s = 1.0/s_to_Gyr
	t0,dt0 = integrate.quad(z_integrand, 0, np.inf)
	t0 *= -(s_to_yr/1e9)
	dt0 *= -(s_to_yr/1e9)
	z_c = None
	z_d = None
	z_e = None
	z_arr = np.logspace(-2,5,1000)
	z_arr = np.append(z_arr,0)
	if(rho0["rad"] != 0):
		z_c = (rho0["m"]/rho0["rad"])-1.0
		if(z_c>=0):
			t_c,dt_c = integrate.quad(modified_z_integrand,0,z_c)
			t_c = t0 + t_c
			dt_c = dt0 + dt_c
			z_arr = np.append(z_arr,z_c)
		else:
			z_c = None
	if(rho0["m"] != 0):
		z_d = np.cbrt(rho0["lambda"]/rho0["m"])-1.0
		if(z_d >= 0):
			t_d,dt_d = integrate.quad(modified_z_integrand,0,z_d)
			t_d = t0 + t_d
			dt_d = dt0 + dt_d
			z_arr = np.append(z_arr,z_d)
		else:
			z_d = None
	if(rho0["rad"] != 0):
		z_e = (rho0["lambda"]/rho0["rad"])**(0.25)-1.0
		if(z_e >= 0):
			t_e,dt_e = integrate.quad(modified_z_integrand,0,z_e)
			t_e = t0 + t_e
			dt_e = dt0 + dt_e
			z_arr = np.append(z_arr,z_e)
		else:
			z_e = None
	z_arr = np.sort(z_arr)
	t_arr = np.array([integrate.quad(modified_z_integrand,0,z)[0] for z in z_arr]) + t0
	z_str = []
	for i in z_arr:
		if(i < 1.0):
			z_str.append("{}".format(np.round(i,2)))
		else:
			z_str.append("{:.2e}".format(i))

	f = {"rad": lambda z: rho0["rad"]*(1+z)**4,
	"m": lambda z: rho0["m"]*(1+z)**(3.0),
	"lambda": lambda z: rho0["lambda"],
	}

	rho = {}
	for i in rho0.keys():
		rho[i] = [f[i](z) for z in z_arr]

	plot_z_c = False
	plot_z_d = False
	plot_z_e = False
	ticks_code = ""
	if(z_c is not None): # radiation to matter
		idx = np.where(t_arr==t_c)[0][0]
		if(rho["rad"][idx+1] == max(rho["m"][idx+1],rho["rad"][idx+1],rho["lambda"][idx+1])):
			plot_z_c = True
			ticks_code += "1"
	else:
		ticks_code += "0"
	if(z_e is not None): # radiation to dark energy
		idx = np.where(t_arr==t_e)[0][0]
		if(rho["rad"][idx+1] == max(rho["m"][idx+1],rho["rad"][idx+1],rho["lambda"][idx+1])):
			plot_z_e = True
			ticks_code += "1"

	else:
		ticks_code += "0"
	if(z_d is not None): # matter to dark energy
		idx = np.where(t_arr==t_d)[0][0]
		if(rho["m"][idx+1] == max(rho["m"][idx+1],rho["rad"][idx+1],rho["lambda"][idx+1])):
			plot_z_d = True
			ticks_code += "1"
	else:
		ticks_code += "0"
	#new_tick_locations = np.array([t_arr[0],t_arr[-1]])
	#new_xtick_labels = [z_str[0],z_str[-1]]
	new_tick_locations, new_xtick_labels = plot_decoder(ticks_code)
	fs = 15
	clrs = ['orange','greenyellow','cyan']
	label = ["Radiation Density","Matter Density","Dark Energy Density"]
	fig = plt.figure(figsize=(12,8))
	ax1 = fig.add_subplot(111)
	for i,j in enumerate(rho.keys()):
		if(rho0[j]!=0.0):
			ax1.plot(np.log10(t_arr),np.log10(rho[j]),label=label[i],c=clrs[i],lw=3)
	if(plot_z_c):
		ax1.axvline(x=np.log10(t_c), label=r'Tranisition time, $t_{\gamma,m}$' + ' = {:.3e} Gyrs'.format(t_c),c='pink',ls='--')
	if(plot_z_d):
		ax1.axvline(x=np.log10(t_d), label=r'Tranisition time, $t_{m,\Lambda}$' + ' = {} Gyrs'.format(np.round(t_d,3)),c='yellow',ls='--')
	if(plot_z_e):
		ax1.axvline(x=np.log10(t_e), label=r'Tranisition time, $t_{\gamma,\Lambda}$' + ' = {} Gyrs'.format(np.round(t_e,3)),c='magenta',ls='--')
	ax1.axvline(x=np.log10(t_arr[0]), label=r'Present time, $t_0$ = {} Gyrs'.format(np.round(t_arr[0],3)),c='white',ls='--')
	plt.legend(loc="best",fontsize=fs)
	ax1.set_xlabel("Time since Big Bang, t (in Gyr)",fontsize=fs)
	ax1.set_ylabel(r"Density $\rho$ (in $kg/m^{3}$)",fontsize=fs)
	ax1.set_title("Evolution of Density with Cosmic Time",fontsize=fs+5)
	ax1.tick_params(axis='both', which='major', labelsize=fs)
	new_tick_labels = [r"$10^{%s}$" % str(int(lbl)) for lbl in ax1.get_xticks()]
	ax1.set_xticklabels(new_tick_labels)
	new_tick_labels = [r"$10^{%s}$" % str(int(lbl)) for lbl in ax1.get_yticks()]
	ax1.set_yticklabels(new_tick_labels)
	ax2 = ax1.twiny()
	ax1.grid(True)
	new_tick_locations = np.log10(new_tick_locations)
	ax2.set_xticks(new_tick_locations)
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xticklabels(new_xtick_labels)
	ax2.tick_params(axis='both', which='major', labelsize=fs)
	ax2.tick_params(axis='both', which='minor', labelsize=fs)
	ax2.set_xlabel("Redshift, z",fontsize=fs)
	output = io.BytesIO()
	FigureCanvasAgg(fig).print_png(output)
	return Response(output.getvalue(), mimetype="image/png")

if(__name__== '__main__'):
   app.run(debug = True)