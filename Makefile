all: run.py
	amuse run.py

clean:
	rm particles/*.hdf5
	rm trajectory.png
	rm trajectory_zoom.png
	rm energy.png
