import psana
import psmon.config as psmonConfig # based on pyQTgraph
import psmon.publish as psmonPublish
import psmon.plots as psmonPlots

class Plotter(object):
	def __init__ (self):
		self.m_test = 0
		psmonPublish.init()
		self.eventcounter = 0
	def event (self,evt,env):
		self.eventcounter += 1
	def endrun (self,evt,env):
		multi = psmonPlots.MultiPlot(self.eventcounter, 'MULTI', ncols=2)
		self.m_test = 5
		for chan in range(4):
			keystring = "acq_" + str(chan)
			matrix = evt.get(psana.ndarray_float64_2,keystring)
			title = keystring
			subplot = psmonPlots.Image(self.eventcounter,title,matrix)
			multi.add(subplot)

		psmonPublish.send('MULTI',multi)

