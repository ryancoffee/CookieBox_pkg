#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel7-gcc48-opt/bin/python

import psana;
import numpy as np;
from numpy.fft import fft as FFT;
from numpy.fft import ifft as IFFT;
from numpy.fft import fftfreq as FREQS;
import matplotlib.pyplot as plt;
from cmath import rect;

nprect = np.vectorize(rect);

dirstr = 'check_data/'
skipshots = 20;
skipsteps = 1;
num = 0.0;

runstrs = ['172','164'];
expstrs = [str('amoi0314'),str('amoi0314')];

t0s = np.zeros((0,1),dtype = float);

for i in range(len(runstrs)):
	runstr = runstrs[i];
	expstr = expstrs[i];
	dsourcestr = 'exp=' + expstr + ':run=' + runstr;
	print("running:\t",dsourcestr)
	ds = psana.DataSource(dsourcestr);
	print(psana.DetNames('detectors'));
	det = psana.Detector('AmoETOF.0:Acqiris.0');
	for nevent,evt in enumerate(ds.events()):
		#acq_d = det.waveform(evt);
		acq_t = det.wftime(evt)
		#print(acq_t);
		t0s = np.row_stack((t0s,acq_t[0,0]));
		if nevent % 1000 ==0 :
			nbins = 100;
			hist,bins = np.histogram(t0s,bins=nbins);
			filename = dirstr + "-".join([expstr,runstr,'t0s_hist.dat']);
			headerstr = np.array_str(bins[:3]) + '...' + np.array_str(bins[-4:]);
			np.savetxt(filename,hist.T,fmt='%i',header = headerstr);
			filename = dirstr + "-".join([expstr,runstr,'t0s.dat']);
			headerstr = np.array_str(acq_t[0,:3]) + '...' + np.array_str(acq_t[0,-4:]);
			np.savetxt(filename,t0s,fmt='%3e',header = headerstr);
			print('updated at nevent = %i' % nevent);
		if nevent > 1e5:
			raise SystemExit;
	print('Done saving for run ',runstr);
print('Done.');
"""
	GDdet = psana.Detector('FEEGasDetEnergy'); #This I hope is the FEE Gas Detector readings
	EBdet = psana.Detector('EBeam'); #This I hope is the BLD data
	evr = psana.Detector('NoDetector.0:Evr.0');
	cd = psana.Detector('ControlData');

	y_init = 0;
	y_final = 0;
	R = np.zeros((0,nsamples),dtype=float);
	R_back = np.zeros((0,nsamples),dtype=float);
	F_abs = np.zeros((0,nsamples),dtype=float);
	F_arg = np.zeros((0,nsamples),dtype=float);
	D = np.zeros((0,nrolls+2),dtype=float);
	d_data = np.zeros(nrolls+2,dtype=float);
	P = np.zeros((0,1),dtype=float);
	sumsignal = np.zeros((1,nsamples),dtype=float);

	'''The third edition: take the average of each step and convert both axes to the right units. '''
		reference_img = np.zeros(nsamples);
		nrefshots = int(0);
		nsumevents = int(0);
		for nstep,step in enumerate(run.steps()):
			print("checking step ",nstep);
			if nstep%skipsteps==0:
				pvList = cd().pvControls();
				ft_abs = np.zeros((1,nsamples),dtype=float);
				ft_arg = np.zeros((1,nsamples),dtype=float);

				for pv in pvList:
					if y_init == 0:
						y_init = pv.value()
					y_final = pv.value()	
					print('Step', nstep, 'name/value',pv.name(),pv.value());
				for nevent,evt in enumerate(step.events()):
					if (printsample and nevent == 20):
						print('printing image');
						img = det.image(evt);
						filename=dirstr + expstr + '_r' + runstr + '_image200.dat';
						np.savetxt(filename,img,fmt='%.6e');
						printsample = False;

					ec = evr.eventCodes(evt)
					if 162 in ec: 
						if subref:
							img = det.image(evt);
							if (img is None):
								continue;
							if nrefshots==0:
								reference_img = np.sum(img[vwin[0]:vwin[1],:],axis=0)/num;
							else:
								reference_img *= (1.-ratio);
								reference_img += (ratio)*np.sum(img[vwin[0]:vwin[1],:],axis=0)/num;
							nrefshots += 1;
						continue;
					if nevent%skipshots == 0:
						lineout = np.zeros(nsamples,dtype=float);
						roller = np.zeros((nrolls,nsamples),dtype=float);
						lineoutFT = np.zeros(nsamples,dtype=complex);
						rollerFT = np.zeros((nrolls,nsamples),dtype=complex);
						dargs = np.zeros((nsamples-1,nrolls),dtype=float);
						darghist = np.zeros((nhist,nrolls),dtype=float);
	
						#print('Fetching event number',nevent);
    						img = det.image(evt);  #img is not areal image, it is a matrix
						ebResults = EBdet.get(evt);
						gdResults = GDdet.get(evt);
						#print(ebResults,gdResults);
 						if (img is None):
							continue;
						if (subref and 162 in ec):
							continue;
						try:
							lineout = (np.sum(img[vwin[0]:vwin[1],:],axis=0)/num) ;
							if subref:
								lineout = lineout - reference_img;
						except:
							continue
						if nsumevents == 0:
							sumsignal = lineout;
							nsumevents += 1;
						else:
							sumsignal *= float(nsumevents)/(nsumevents+1);
							sumsignal += lineout/(nsumevents+1);
							nsumevents + 1;
						R = np.row_stack((R,lineout));#it stacks from top to bottom

						for r in range(nrolls):
							#roller[r,:] = np.roll(np.copy(lineout) , (r*len(lineout))//nrolls);
							roller[r,:] = np.roll(np.copy(lineout), r*len(lineout)//nrolls);

						lineoutFT = FFT(np.roll(lineout,len(lineout)//2) );
						rollerFT = FFT(roller,axis=1);
						
						%Try dividing by the sum of all signals from the entire scan.
						%This might get rid of the etalon
						#ft_abs = np.copy(np.abs(lineoutFT));
						#ft_arg = np.copy(np.angle(lineoutFT));
						ft_abs = np.copy(np.abs(rollerFT[1,:]));
						ft_arg = np.copy(np.angle(rollerFT[1,:]));
						dargs = np.diff( np.unwrap( np.angle( rollerFT ) , axis = 1 ), axis =1 );

						#lineoutback = np.real( IFFT(nprect( w_weights, np.angle(rollerFT[0,:]) )) );
						#lineoutback = np.real( IFFT(nprect(w_weights,ft_arg)) );
						lineoutback = np.real( IFFT(nprect(ft_abs,ft_arg)) );
						#lineoutback = np.real( IFFT(rollerFT[1,:]) );
						R_back = np.row_stack((R_back,lineoutback));
						
						avg = np.average(dargs,axis=1,weights = ft_abs[1:]);
						#avg = rewrap(avg,-2.*np.pi/nrolls,0.,2.*np.pi/nrolls);
						d_data[0] = y_final*delayscales[i];
						d_data[1:-1] = avg[:];
						#d_data[-1] = chooseslope(avg,nrolls);
						d_data[-1] = timsChoice(avg,nrolls);
						#print(d_data)

						eb_data = (ebResults.ebeamL3Energy() , ebResults.ebeamCharge(), ebResults.ebeamEnergyBC1(), ebResults.ebeamEnergyBC2(), ebResults.ebeamLTU250(), ebResults.ebeamLTU450(), ebResults.ebeamLTUAngX(), ebResults.ebeamLTUAngY(), ebResults.ebeamLTUPosX(), ebResults.ebeamLTUPosY(), ebResults.ebeamUndAngX(), ebResults.ebeamUndAngY(), ebResults.ebeamUndPosX(), ebResults.ebeamUndPosY(), ebResults.ebeamPkCurrBC1(), ebResults.ebeamEnergyBC1(), ebResults.ebeamPkCurrBC2(), ebResults.ebeamEnergyBC2(), ebResults.ebeamDumpCharge());
						gd_data = ( gdResults.f_11_ENRC(), gdResults.f_12_ENRC(), gdResults.f_21_ENRC(), gdResults.f_22_ENRC(), gdResults.f_63_ENRC(), gdResults.f_64_ENRC() );
						p_data = np.concatenate((eb_data,gd_data));
						if len(P) < len(p_data):
							P = np.zeros((0,len(p_data)),dtype=float);
						D = np.row_stack((D,d_data));
						P = np.row_stack((P,p_data));

			F_abs = np.row_stack((F_abs,ft_abs));
			F_arg = np.row_stack((F_arg,ft_arg));

			if (nstep%10 == 0):
				filename="%s/%s_r%s_%i_%i_matrix.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
				np.savetxt(filename,R.T,fmt='%i');
				filename="%s/%s_r%s_%i_%i_matrix_back.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
				np.savetxt(filename,R_back.T,fmt='%i');
				filename="%s/%s_r%s_%i_%i_steps_ft_abs.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
				np.savetxt(filename,F_abs.T,fmt='%.6e');
				filename=dirstr + expstr + '_r' + runstr + '_delays.dat';
				np.savetxt(filename,D,fmt='%.6e');
				filename=dirstr + expstr + '_r' + runstr + '_sumsignal.dat';
				np.savetxt(filename,sumsignal,fmt='%.6e');
	#for plot
	y_dim = int(np.shape(R)[0]);
	x_dim = int(np.shape(R)[1]);
	lam = i2lam(np.arange(x_dim,dtype=float));

	print('\n\n\t\tWatch out!  Transposing ouput')
	if subref:
		runstr += "_refsub";
	filename="%s/%s_r%s_%i_%i_matrix.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
	np.savetxt(filename,R.T,fmt='%i');
	filename="%s/%s_r%s_%i_%i_matrix_back.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
	np.savetxt(filename,R_back.T,fmt='%i');
	filename="%s/%s_r%s_%i_%i_steps_ft_abs.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
	np.savetxt(filename,F_abs.T,fmt='%.6e');
	filename="%s/%s_r%s_%i_%i_steps_ft_arg.dat" % (dirstr,expstr,runstr,vwin[0],vwin[1]);
	np.savetxt(filename,F_arg.T,fmt='%.6e');
	filename=dirstr + expstr + '_r' + runstr + '_params.dat';
	np.savetxt(filename,P,fmt='%.6e');
	filename=dirstr + expstr + '_r' + runstr + '_delays.dat';
	np.savetxt(filename,D,fmt='%.6e');
	filename=dirstr + expstr + '_r' + runstr + '_wavelegths.dat';
	np.savetxt(filename,lam,fmt='%.6e');
	filename=dirstr + expstr + '_r' + runstr + '_sumsignal.dat';
	np.savetxt(filename,sumsignal,fmt='%.6e');
"""
